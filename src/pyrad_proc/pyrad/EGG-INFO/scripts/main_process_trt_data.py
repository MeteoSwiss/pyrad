#!/home/daniel/anaconda3/bin/python
# -*- coding: utf-8 -*-

"""
================================================
main_process_trt_data.py
================================================

This program processes TRT data obtaining plots of time series of
parameters over the TRT cell trajectory

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import glob
import os
from shutil import copy
from warnings import warn

import numpy as np

from pyrad.io import read_trt_traj_data, write_trt_cell_scores
from pyrad.io import write_trt_cell_lightning
from pyrad.util import belongs_roi_indices
from pyrad.graph import plot_timeseries, plot_scatter_comp, plot_pos

print(__doc__)


def main():
    """
    """
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    # positional arguments
    parser.add_argument(
        'days', nargs='+', type=str,
        help='Dates to process. Format YYYY-MM-DD')

    # keyword arguments
    parser.add_argument(
        '--trtbase', type=str,
        default='/store/msrad/radar/trt/',
        help='name of folder containing the TRT cell data')

    parser.add_argument(
        '--lon', type=str,
        default='8.9000010,9.2000000,9.4999970,9.4999970,8.9000010',
        help=('longitude of the points defining the perimeter of the area ' +
              'of interest'))

    parser.add_argument(
        '--lat', type=str,
        default='47.0000030,47.0000030,47.0000030,47.5999930,47.5999930',
        help=('latitude of the points defining the perimeter of the area ' +
              'of interest'))

    args = parser.parse_args()

    print("====== TRT cell processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== comparison finished: ")

    time_dir_list = args.days
    lons = args.lon.split(',')
    lats = args.lat.split(',')

    if np.size(lons) != np.size(lats):
        warn(
            str(np.size(lons))+' longitudes but '+str(np.size(lats)) +
            ' latitudes. Their number must be equal')
        return

    lon_list = []
    lat_list = []
    for i, lon in enumerate(lons):
        lon_list.append(float(lon))
        lat_list.append(float(lats[i]))

    roi = {
        'lon': lon_list,
        'lat': lat_list
    }

    # List for collection of max data
    cell_ID_max_list = []
    nflashes_max_list = []
    area_flash_max_list = []
    flash_density_max_list = []
    time_flash_density_max_list = []
    flash_density_max_rank_list = []
    rank_max_list = []
    time_rank_max_list = []

    # List for collection of flashes data
    cell_ID_list = np.ma.asarray([], dtype=int)
    time_list = np.ma.asarray([], dtype=datetime.datetime)
    lon_list = np.ma.asarray([], dtype=float)
    lat_list = np.ma.asarray([], dtype=float)
    flash_density_list = np.ma.asarray([], dtype=float)
    rank_flash_density_list = np.ma.asarray([], dtype=float)
    area_list = np.ma.asarray([], dtype=float)
    nflash_list = np.ma.asarray([], dtype=int)

    for i, time_dir in enumerate(time_dir_list):
        data_input_path = args.trtbase+time_dir+'/TRTC_cell/'
        data_output_base = args.trtbase+time_dir+'/TRTC_cell_plots/'

        flist = glob.glob(data_input_path+'*.trt')

        for fname in flist:
            print('Reading TRT trajectory file '+fname)
            (traj_ID, yyyymmddHHMM, lon, lat, _, _, _, area, vel_x, vel_y,
             det, RANKr, CG_n, CG_p, CG, _, ET45, ET45m, ET15, ET15m, VIL,
             maxH, maxHm, POH, _, _, _, _) = read_trt_traj_data(fname)

            inds, is_roi = belongs_roi_indices(lat, lon, roi)

            if is_roi == 'None':
                continue
            elif is_roi == 'Some' and len(lat[inds]) < 3:
                continue

            data_output_path = data_output_base+is_roi+'/'
            if not os.path.isdir(data_output_path):
                os.makedirs(data_output_path)

            # copy file
            copy(fname, data_output_path)

            # general caracteristics
            flash_density = CG/area
            cell_ID_max_list.append(traj_ID[0])
            flash_density_max_list.append(np.max(flash_density))
            nflashes_max_list.append(CG[np.argmax(flash_density)])
            area_flash_max_list.append(area[np.argmax(flash_density)])
            time_flash_density_max_list.append(
                yyyymmddHHMM[np.argmax(flash_density)])
            flash_density_max_rank_list.append(
                RANKr[np.argmax(flash_density)])
            rank_max_list.append(np.max(RANKr))
            time_rank_max_list.append(yyyymmddHHMM[np.argmax(RANKr)])

            cell_ID_list = np.append(cell_ID_list, traj_ID)
            time_list = np.append(time_list, yyyymmddHHMM)
            lon_list = np.append(lon_list, lon)
            lat_list = np.append(lat_list, lat)
            flash_density_list = np.append(flash_density_list, flash_density)
            rank_flash_density_list = np.append(
                rank_flash_density_list, RANKr)
            area_list = np.append(area_list, area)
            nflash_list = np.append(nflash_list, CG)

            # Time series plots
            figfname = data_output_path+str(traj_ID[0])+'_flash_density.png'
            plot_timeseries(
                yyyymmddHHMM, [flash_density], [figfname], labelx='Time UTC',
                labely='Flash density [flashes/km2]',
                title=str(traj_ID[0])+' flash density')

            figfname = data_output_path+str(traj_ID[0])+'_area.png'
            plot_timeseries(
                yyyymmddHHMM, [area], [figfname], labelx='Time UTC',
                labely='area [km2]', title=str(traj_ID[0])+' cell area')

            figfname = data_output_path+str(traj_ID[0])+'_vel.png'
            plot_timeseries(
                yyyymmddHHMM, [vel_x, vel_y], [figfname], labelx='Time UTC',
                labely='Velocity [km/h]', labels=['x speed', 'y speed'],
                title=str(traj_ID[0])+' cell velocity')

            figfname = data_output_path+str(traj_ID[0])+'_det.png'
            plot_timeseries(
                yyyymmddHHMM, [det], [figfname], labelx='Time UTC',
                labely='Detection threshold [dBZ]',
                title=str(traj_ID[0])+' cell detection threshold')

            figfname = data_output_path+str(traj_ID[0])+'_rank.png'
            plot_timeseries(
                yyyymmddHHMM, [RANKr], [figfname], labelx='Time UTC',
                labely='Rank [-]', title=str(traj_ID[0])+' cell rank')

            figfname = data_output_path+str(traj_ID[0])+'_lightning.png'
            plot_timeseries(
                yyyymmddHHMM, [CG_n, CG_p, CG], [figfname], labelx='Time UTC',
                labely='N flash [-]', labels=['CG-', 'CG+', 'CG'],
                title=str(traj_ID[0])+' flashes in cell')

            figfname = data_output_path+str(traj_ID[0])+'_ET.png'
            plot_timeseries(
                yyyymmddHHMM, [ET45, ET45m, ET15, ET15m], [figfname],
                labelx='Time UTC', labely='Echo Top [km]',
                labels=['ET45', 'ET45m', 'ET15', 'ET15m'],
                title=str(traj_ID[0])+' Echo top')

            figfname = data_output_path+str(traj_ID[0])+'_VIL.png'
            plot_timeseries(
                yyyymmddHHMM, [VIL], [figfname], labelx='Time UTC',
                labely='VIL [Kg/m2]', labels=['VIL'],
                title=str(traj_ID[0])+' VIL')

            figfname = data_output_path+str(traj_ID[0])+'_maxH.png'
            plot_timeseries(
                yyyymmddHHMM, [maxH, maxHm], [figfname], labelx='Time UTC',
                labely='Max. Echo Height [Km]', labels=['maxH', 'maxHm'],
                title=str(traj_ID[0])+' Height of Max. Reflectivity')

            figfname = data_output_path+str(traj_ID[0])+'_POH.png'
            plot_timeseries(
                yyyymmddHHMM, [POH], [figfname], labelx='Time UTC',
                labely='POH [%]', labels=['POH'],
                title=str(traj_ID[0])+' Probability of Hail')

            # plot position
            # get time since start of cell in s
            td_vec = yyyymmddHHMM-yyyymmddHHMM[0]
            tt_s = np.empty(td_vec.size, dtype=float)
            for j, td in enumerate(td_vec):
                tt_s[j] = td.total_seconds()
            cb_label = (
                'Time since '+yyyymmddHHMM[0].strftime('%Y-%m-%d %H:%M') +
                ' [s]')
            figfname = data_output_path+str(traj_ID[0])+'_pos.png'
            figfname = plot_pos(
                lat, lon, tt_s, [figfname], cb_label=cb_label,
                titl=str(traj_ID[0])+' Cell Position')
            print('Plotted '+' '.join(figfname))

    fname = args.trtbase+'Santis_cell_scores.csv'
    write_trt_cell_scores(
        cell_ID_max_list, time_flash_density_max_list,
        flash_density_max_rank_list, nflashes_max_list, area_flash_max_list,
        flash_density_max_list, time_rank_max_list, rank_max_list, fname)

    fname = args.trtbase+'Santis_cell_euclid_lightning.csv'
    write_trt_cell_lightning(
        cell_ID_list, time_list, lon_list, lat_list, area_list,
        rank_flash_density_list, nflash_list, flash_density_list, fname)

    plot_scatter_comp(
        flash_density_list, rank_flash_density_list/10.,
        [args.trtbase+'hist_flash_density_rank'],
        labelx='flash density [flashes/km2]', labely='rank',
        titl='Flash density vs Rank', axis=None, metadata=None, dpi=72)


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
