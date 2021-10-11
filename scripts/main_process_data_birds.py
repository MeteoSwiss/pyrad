#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
Pyrad: The MeteoSwiss Radar Processing framework
================================================

Welcome to Pyrad!

This program processes bird data

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import os
import glob
from warnings import warn

from pyrad.flow.flow_control import main as pyrad_main
from pyrad.io import get_fieldname_pyart
from pyrad.io import read_profile_ts
from pyrad.graph import get_field_name, _plot_time_range

from pyart.config import get_metadata

print(__doc__)


def main():
    """
    """
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    # positional arguments
    parser.add_argument(
        'proc_cfgfile', type=str, help='name of main configuration file')
    parser.add_argument(
        'starttime', type=str,
        help=('starting time of the data to be processed. ' +
              'Format ''YYYYMMDDhhmmss'''))
    parser.add_argument(
        'endtime', type=str,
        help='end time of the data to be processed. Format ''YYYYMMDDhhmmss''')

    # keyword arguments
    parser.add_argument(
        '--cfgpath', type=str,
        default=os.path.expanduser('~')+'/pyrad/config/processing/',
        help='configuration file path')

    parser.add_argument(
        '--storepath', type=str,
        default='/store/msrad/radar/pyrad_products/rad4alp_birds_PHA/',
        help='Base data storing path')

    parser.add_argument(
        '--hres', type=int, default=200, help='Height resolution [m]')

    args = parser.parse_args()

    print("====== PYRAD data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD data processing finished: ")

    print('config path: '+args.cfgpath)
    print('config file: '+args.proc_cfgfile)
    print('start time: '+args.starttime)
    print('end time: '+args.endtime)

    proc_starttime = datetime.datetime.strptime(
        args.starttime, '%Y%m%d%H%M%S')
    proc_endtime = datetime.datetime.strptime(
        args.endtime, '%Y%m%d%H%M%S')
    cfgfile_proc = args.cfgpath+args.proc_cfgfile

    pyrad_main(cfgfile_proc, starttime=proc_starttime, endtime=proc_endtime)

    # Plot time-height
    file_base = args.storepath
    hres = args.hres

    datatype_list = [
        'dBZc', 'eta_h', 'bird_density', 'WIND_SPEED', 'WIND_DIRECTION',
        'wind_vel_h_u', 'wind_vel_h_v', 'wind_vel_v']

    startdate = proc_starttime.replace(hour=0, minute=0, second=0, microsecond=0)
    enddate = proc_endtime.replace(hour=0, minute=0, second=0, microsecond=0)
    ndays = int((enddate-startdate).days)+1
    for datatype in datatype_list:
        flist = []
        for i in range(ndays):
            time_dir = (
                proc_starttime+datetime.timedelta(days=i)).strftime('%Y-%m-%d')

            filepath = (
                file_base+time_dir+'/VAD/PROFILE_WIND/' +
                '*_wind_profile_VAD_WIND_hres'+str(hres)+'.csv')
            labels = [
                'u_wind', 'std_u_wind', 'np_u_wind',
                'v_wind', 'std_v_wind', 'np_v_wind',
                'w_wind', 'std_w_wind', 'np_w_wind',
                'mag_h_wind', 'dir_h_wind']
            label_nr = 0
            if datatype == 'dBZc':
                filepath = (
                    file_base+time_dir+'/velFilter/PROFILE_dBZc/' +
                    '*_rhi_profile_*_dBZc_hres'+str(hres)+'.csv')
                labels = [
                    '50.0-percentile', '25.0-percentile', '75.0-percentile']

                # dBZ mean data
                # filepath = (
                #     file_base+time_dir+'/velFilter/PROFILE_dBZc_mean/' +
                #     '*_rhi_profile_*_dBZc_hres'+str(hres)+'.csv')
                # labels = [
                #     'Mean', 'Min', 'Max']

                # dBZ linear mean data
                # filepath = (
                #     file_base+time_dir+'/velFilter/PROFILE_dBZc_linear_mean/' +
                #     '*_rhi_profile_*_dBZc_hres'+str(hres)+'.csv')
                # labels = [
                #     'Mean', 'Min', 'Max']

                # dBZ before filtering with fitted velocity
                # filepath = (
                #     file_base+time_dir+'/echoFilter/PROFILE_dBZc/' +
                #     '*_rhi_profile_*_dBZc_hres'+str(hres)+'.csv')
                # labels = [
                #     '50.0-percentile', '25.0-percentile', '75.0-percentile']
                #
                # dBZ before filtering with fitted velocity. Linear mean
                # filepath = (
                #     file_base+time_dir+'/echoFilter/PROFILE_dBZc_linear_mean/' +
                #     '*_rhi_profile_*_dBZc_hres'+str(hres)+'.csv')
                # labels = [
                #     'Mean', 'Min', 'Max']
            elif datatype == 'eta_h':
                filepath = (
                    file_base+time_dir+'/vol_refl/PROFILE/' +
                    '*_rhi_profile_*_eta_h_hres'+str(hres)+'.csv')
                labels = [
                    '50.0-percentile', '25.0-percentile', '75.0-percentile']

                # mean data
                # filepath = (
                #     file_base+time_dir+'/vol_refl/PROFILE_mean/' +
                #     '*_rhi_profile_*_eta_h_hres'+str(hres)+'.csv')
                # labels = [
                #     'Mean', 'Min', 'Max']
            elif datatype == 'bird_density':
                filepath = (
                    file_base+time_dir+'/bird_density/PROFILE/' +
                    '*_rhi_profile_*_bird_density_hres'+str(hres)+'.csv')
                labels = [
                    '50.0-percentile', '25.0-percentile', '75.0-percentile']

                # mean data
                # filepath = (
                #     file_base+time_dir+'/bird_density/PROFILE_mean/' +
                #     '*_rhi_profile_*_bird_density_hres'+str(hres)+'.csv')
                # labels = [
                #     'Mean', 'Min', 'Max']
            elif datatype == 'WIND_SPEED':
                label_nr = 9
            elif datatype == 'WIND_DIRECTION':
                label_nr = 10
            elif datatype == 'wind_vel_h_u':
                label_nr = 0
            elif datatype == 'wind_vel_h_v':
                label_nr = 3
            elif datatype == 'wind_vel_v':
                label_nr = 6

            flist_aux = glob.glob(filepath)
            if not flist_aux:
                warn('No profile files found in '+filepath)
                continue
            flist.extend(flist_aux)

        if not flist:
            warn('No profile files found')
            continue
        flist.sort()

        field_name = get_fieldname_pyart(datatype)
        field_dict = get_metadata(field_name)
        titl = 'bird retrieval '+args.starttime+'\n'+get_field_name(
            field_dict, field_name)

        tbin_edges, hbin_edges, np_ma, data_ma, t_start = read_profile_ts(
            flist, labels, hres=hres, label_nr=label_nr)

        basepath_out = os.path.dirname(flist[0])
        fname = (
            basepath_out+'/'+args.starttime+'_TIME_HEIGHT_' +
            datatype+'_hres'+str(hres)+'.png')

        vmin = vmax = None
        _plot_time_range(
            tbin_edges, hbin_edges/1000., data_ma, field_name, [fname],
            titl=titl, figsize=[10, 8], vmin=vmin, vmax=vmax, dpi=72)

        print("----- plot to '%s'" % fname)

        # Plot number of points
        field_dict = get_metadata('number_of_samples')
        titl = 'bird retrieval '+args.starttime+'\n'+get_field_name(
            field_dict, 'number_of_samples')

        fname = (
            basepath_out+'/'+args.starttime+'_TIME_HEIGHT_' +
            datatype+'nsamples_hres'+str(hres)+'.png')

        vmin = vmax = None
        _plot_time_range(
            tbin_edges, hbin_edges/1000., np_ma, 'number_of_samples', [fname],
            titl=titl, figsize=[10, 8], vmin=vmin, vmax=vmax, dpi=72)

        print("----- plot to '%s'" % fname)


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
