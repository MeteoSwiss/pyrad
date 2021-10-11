#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_precipitation_comparison
================================================

This program compares radar data with a point measurement sensor.

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import atexit
import glob
from warnings import warn
import os

import numpy as np

from pyrad.io import read_ts_cum
from pyrad.graph import plot_scatter_comp
from pyrad.util import compute_1d_stats

print(__doc__)


def main():
    """
    """
    param_vec = ['RR_Z', 'RR_hydro']
    smn_station_vec = ['CIM', 'MAG', 'OTL']
    tstart = '20180401'
    tend = '20180430'

    np_radar_min = 6
    np_sensor_min = 6
    min_val = 0.2

    fbase = '/data/pyrad_products/mals_loc_dataquality/'
    img_ext = 'png'
    avg_time = 3600

    print("====== precipitation comparison started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== comparison finished: ")

    startdate = datetime.datetime.strptime(tstart, '%Y%m%d')
    enddate = datetime.datetime.strptime(tend, '%Y%m%d')

    ndays = (enddate - startdate).days + 1
    print('Number of days to process: '+str(ndays)+'\n\n')

    for param in param_vec:
        ts_vec = np.array([])
        val_radar = np.ma.array([])
        np_radar = np.array([])
        val_sensor = np.ma.array([])
        np_sensor = np.array([])

        for station in smn_station_vec:
            for day in range(ndays):
                current_date = startdate + datetime.timedelta(days=day)
                day_dir = current_date.strftime("%Y-%m-%d")
                daybase = current_date.strftime("%Y%m%d")

                fpath = (fbase+day_dir+'/rg'+station+'_'+param+'/RRcum' +
                         str(avg_time)+'s/')
                fname = glob.glob(
                    fpath+daybase+'_'+str(avg_time) +
                    's_acc_ts_comp_POINT_MEASUREMENT_*.csv')
                if not fname:
                    warn('No file found in '+fpath)
                    continue
                else:
                    (ts_aux, np_radar_aux, radar_value_aux, np_sensor_aux,
                     sensor_value_aux) = read_ts_cum(fname[0])
                    ts_vec = np.append(ts_vec, ts_aux)
                    val_radar = np.ma.append(val_radar, radar_value_aux)
                    np_radar = np.append(np_radar, np_radar_aux)
                    val_sensor = np.ma.append(val_sensor, sensor_value_aux)
                    np_sensor = np.append(np_sensor, np_sensor_aux)

        # filter out undesired data
        ind = np.where(np.logical_and(
            np.logical_and(
                np_radar >= np_radar_min, np_sensor >= np_sensor_min),
            np.logical_and(val_sensor >= min_val, val_radar >= min_val)))[0]

        val_sensor = val_sensor[ind]
        val_radar = val_radar[ind]

        # compute statistics
        stats = compute_1d_stats(val_sensor, val_radar)

        # create output image
        fpath = fbase+'RR/'
        if os.path.isdir(fpath):
            pass
        else:
            os.makedirs(fpath)

        figfname = [
            startdate.strftime('%Y%m%d')+'-'+enddate.strftime('%Y%m%d')+'_' +
            str(avg_time)+'s_acc_ts_comp_'+param+'.'+img_ext]

        for i in range(len(figfname)):
            figfname[i] = fpath+figfname[i]

        labelx = 'RG (mm)'
        labely = 'Radar (mm)'
        titl = (str(avg_time)+' s Acc. Comp. '+startdate.strftime('%Y%m%d') +
                '-'+enddate.strftime('%Y%m%d'))

        metadata = (
            'npoints: '+str(stats['npoints'])+'\n' +
            'NB: '+'{:.2f}'.format(float(stats['NB']))+'\n' +
            'corr: '+'{:.2f}'.format(float(stats['corr']))+'\n' +
            'RMS: '+'{:.2f}'.format(float(stats['RMS']))+'\n' +
            'Nash: '+'{:.2f}'.format(float(stats['Nash']))+'\n')

        plot_scatter_comp(
            val_sensor, val_radar, figfname, labelx=labelx,
            labely=labely, titl=titl, axis='equal', metadata=metadata,
            dpi=300)


def _print_end_msg(text):
    """
    prints end message

    Parameters
    ----------
    text : str
        the text to be printed

    Returns
    -------
    Nothing

    """
    print(text + datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))


# ---------------------------------------------------------
# Start main:
# ---------------------------------------------------------
if __name__ == "__main__":
    main()
