#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_process_data_trt
================================================

This program activates the Pyrad processing to obtain
TRT cells trajectories and obtains post-processing products

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import os
import glob
from warnings import warn

import numpy as np

from pyrad.flow.flow_control import main as pyrad_main
from pyrad.io import get_fieldname_pyart, write_trt_cell_lightning
from pyrad.io import read_profile_ts, read_histogram_ts, read_quantiles_ts
from pyrad.io import read_trt_traj_data, read_thundertracking_info
from pyrad.graph import get_field_name, get_colobar_label
from pyrad.graph import _plot_time_range

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

    # keyword arguments
    parser.add_argument(
        '--days', type=str, default=None,
        help='Dates to process. Format YYYY-MM-DD. Coma separated')

    parser.add_argument(
        '--info_file', type=str,
        default='/store/msrad/radar/thundertracking/info/thundertracking_info.csv',
        help='configuration file path')

    parser.add_argument(
        '--years', type=str,
        default=None,
        help='Years to process. If None all years in file will be processed')

    parser.add_argument(
        '--max_rank', type=float,
        default=None,
        help='Max rank to process')

    parser.add_argument(
        '--trtbase', type=str,
        default='/store/msrad/radar/trt/',
        help='name of folder containing the TRT cell data')

    parser.add_argument(
        '--postproc', type=int,
        default=1,
        help='If true the data will be post-processed')

    parser.add_argument(
        '--radarbase', type=str,
        default='/store/msrad/radar/pyrad_products/thundertracking/',
        help='name of folder containing the radar data')

    parser.add_argument(
        '--cfgpath', type=str,
        default=os.path.expanduser('~')+'/pyrad/config/processing/',
        help='configuration file path')

    parser.add_argument(
        '--datatypes', type=str,
        default='RR,hydro,KDPc,dBZc,RhoHVc,TEMP,ZDRc',
        help='Name of the polarimetric moments to process. Coma separated')

    parser.add_argument(
        '--datasets', type=str,
        default='RR,hydro,KDPc,dBZc,RhoHVc,TEMP,ZDRc',
        help='Name of the directory containing the datasets')

    parser.add_argument(
        '--hres', type=float, default=250., help='Height resolution')

    parser.add_argument(
        '--path_structure', type=int, default=2,
        help='If true the data is at the cell center')

    args = parser.parse_args()

    print("====== PYRAD TRT data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD TRT data processing finished: ")

    print('config path: '+args.cfgpath)
    print('config file: '+args.proc_cfgfile)
    print('trt path: '+args.trtbase)
    print('radar data path: '+args.radarbase)

    cfgfile_proc = args.cfgpath+args.proc_cfgfile
    trajtype = 'trt'

    if args.days is not None:
        time_dir_list = args.days.split(',')
    else:
        # get the years to process
        years = None
        if args.years is not None:
            years = list(map(int, args.years.split(',')))

        _, max_rank, _, trt_time_start, trt_time_end = read_thundertracking_info(
            args.info_file)
        trt_times = np.append(trt_time_start, trt_time_end)

        trt_dates = np.array([], dtype=datetime.date)
        for trt_time in trt_times:
            trt_dates = np.append(trt_dates, trt_time.date())
        trt_dates = np.sort(np.unique(trt_dates))
        time_dir_list = []
        for rank, trt_date in zip(max_rank, trt_dates):
            if years is not None:
                if trt_date.year not in years:
                    continue

            if args.max_rank is not None:
                if rank > args.max_rank:
                    continue

            time_dir_list.append(trt_date.strftime("%Y-%m-%d"))


    if args.postproc:
        datatype_list = args.datatypes.split(',')
        dataset_list = args.datasets.split(',')

        if np.size(datatype_list) != np.size(dataset_list):
            warn(
                str(np.size(datatype_list))+' datatypes but ' +
                str(np.size(dataset_list)) +
                ' dataset directories. Their number must be equal')
            return

    # Find all TRT files in directory
    trt_list = []
    if args.days is not None:
        for time_dir in time_dir_list:
            trt_list.extend(glob.glob(
                args.trtbase+time_dir+'/TRTC_cell_plots/All/*.trt'))
            trt_list.extend(glob.glob(
                args.trtbase+time_dir+'/TRTC_cell_plots/Some/*.trt'))
    else:
        for time_dir in time_dir_list:
            trt_list.extend(glob.glob(
                args.trtbase+time_dir+'/TRTC_cell/*_tt.trt'))

    if len(trt_list) == 0:
        warn('No valid TRT files found in '+args.trtbase)
        return

    # Pyrad data processing
    trt_cell_id_list = []
    trt_file_list = []
    for fname in trt_list:
        print('processing TRT cell file '+fname)
        try:
            infostr = os.path.basename(fname).split('.')[0]
            infostr = infostr.replace('_tt', '')
            pyrad_main(
                cfgfile_proc, trajfile=fname, infostr=infostr,
                trajtype=trajtype)
            trt_cell_id_list.append(infostr)
            trt_file_list.append(fname)
        except:
            warn('Unable to process TRT cell file '+fname)

    if not args.postproc:
        return

    # plot time series and get altitude of graupel column
    if 'hydro' in datatype_list:
        cell_ID_list = np.asarray([], dtype=int)
        time_list = np.asarray([], dtype=datetime.datetime)
        lon_list = np.asarray([], dtype=float)
        lat_list = np.asarray([], dtype=float)
        area_list = np.asarray([], dtype=float)
        rank_list = np.asarray([], dtype=float)
        rm_hmin_list = np.ma.asarray([], dtype=float)
        rm_hmax_list = np.ma.asarray([], dtype=float)

    for i, trt_cell_id in enumerate(trt_cell_id_list):
        print('\n\nPost-processing cell: '+trt_cell_id)
        dt_str = trt_cell_id[0:12]
        dt_cell = datetime.datetime.strptime(dt_str, "%Y%m%d%H%M")
        time_dir = dt_cell.strftime("%Y-%m-%d")
        for j, datatype in enumerate(datatype_list):
            dataset = dataset_list[j]
            if args.path_structure == 1:
                file_base2 = args.radarbase+time_dir+'/'+dataset+'_trt_center_traj/'
            elif args.path_structure == 0:
                file_base2 = args.radarbase+time_dir+'/'+dataset+'_trt_traj/'
            elif args.path_structure == 2:
                file_base2 = args.radarbase+time_dir+'/trt_traj_tt/'

            field_name = get_fieldname_pyart(datatype)
            field_dict = get_metadata(field_name)
            titl = 'TRT cell '+trt_cell_id+'\n'+get_field_name(
                field_dict, field_name)

            # plot time-height
            if args.path_structure == 2:
                flist = glob.glob(
                    file_base2+'PROFILE_'+dataset+'/*_'+trt_cell_id+'_rhi_profile_*_' +
                    datatype+'_hres'+str(int(args.hres))+'.csv')
            else:
                flist = glob.glob(
                    file_base2+'PROFILE/*_'+trt_cell_id+'_rhi_profile_*_' +
                    datatype+'_hres'+str(int(args.hres))+'.csv')


            if not flist:
                warn('No profile files found in '+file_base2 +
                     'PROFILE/ for TRT cell ' +
                     trt_cell_id+' with resolution '+str(args.hres))
            else:
                if args.path_structure == 1:
                    labels = ['Mean', 'Min', 'Max']
                else:
                    labels = [
                        '50.0-percentile', '25.0-percentile', '75.0-percentile']
                    if datatype == 'RhoHVc':
                        labels = [
                            '80.0-percentile', '65.0-percentile',
                            '95.0-percentile']
                    elif datatype == 'hydro':
                        labels = [
                            'Mode', '2nd most common', '3rd most common',
                            '% points mode', '% points 2nd most common',
                            '% points 3rd most common']
                    elif datatype == 'entropy' or 'prop' in datatype:
                        labels = ['Mean', 'Min', 'Max']

                tbin_edges, hbin_edges, _, data_ma, start_time = (
                    read_profile_ts(flist, labels, hres=args.hres, t_res=None))

                basepath_out = os.path.dirname(flist[0])
                fname = (
                    basepath_out+'/'+trt_cell_id+'_trt_TIME_HEIGHT_' +
                    datatype+'_hres'+str(args.hres)+'.png')

                vmin = vmax = None
                if datatype == 'RhoHVc':
                    vmin = 0.95
                    vmax = 1.00

                xlabel = (
                    'time (s from '+start_time.strftime("%Y-%m-%d %H:%M:%S") +
                    ')')
                _plot_time_range(
                    tbin_edges, hbin_edges, data_ma, field_name, [fname],
                    titl=titl, xlabel=xlabel, ylabel='height (m MSL)',
                    figsize=[10, 8], vmin=vmin, vmax=vmax, dpi=72)

                print("----- plot to '%s'" % fname)

                # Get min and max altitude of graupel/hail area
                if datatype == 'hydro':
                    (traj_ID, yyyymmddHHMM, lon, lat, _, _, _, area, _, _, _,
                     RANKr, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
                     _) = read_trt_traj_data(trt_file_list[i])

                    hmin, hmax = get_graupel_column(
                        tbin_edges, hbin_edges, data_ma, start_time,
                        yyyymmddHHMM)

                    cell_ID_list = np.append(cell_ID_list, traj_ID)
                    time_list = np.append(time_list, yyyymmddHHMM)
                    lon_list = np.append(lon_list, lon)
                    lat_list = np.append(lat_list, lat)
                    area_list = np.append(area_list, area)
                    rank_list = np.append(rank_list, RANKr)
                    rm_hmin_list = np.ma.append(rm_hmin_list, hmin)
                    rm_hmax_list = np.ma.append(rm_hmax_list, hmax)

            # plot time-hist
            if args.path_structure == 2:
                flist = glob.glob(
                    file_base2+'HISTOGRAM_'+dataset+'/*_'+trt_cell_id+'_histogram_*_' +
                    datatype+'.csv')
            else:
                flist = glob.glob(
                    file_base2+'HISTOGRAM/*_'+trt_cell_id+'_histogram_*_' +
                    datatype+'.csv')

            if not flist:
                warn('No histogram files found in '+file_base2 +
                     'HISTOGRAM/ for TRT cell '+trt_cell_id)
            else:
                tbin_edges, bin_edges, data_ma, start_time = read_histogram_ts(
                    flist, datatype, t_res=None)

                basepath_out = os.path.dirname(flist[0])
                fname = (
                    basepath_out+'/'+trt_cell_id+'_trt_HISTOGRAM_'+datatype +
                    '.png')

                data_ma[data_ma == 0.] = np.ma.masked
                xlabel = (
                    'time (s from '+start_time.strftime("%Y-%m-%d %H:%M:%S") +
                    ')')
                _plot_time_range(
                    tbin_edges, bin_edges, data_ma, 'frequency_of_occurrence',
                    [fname], titl=titl, xlabel=xlabel,
                    ylabel=get_colobar_label(field_dict, field_name),
                    vmin=0., vmax=np.max(data_ma), figsize=[10, 8], dpi=72)

                print("----- plot to '%s'" % fname)

            # plot quantiles
            flist = glob.glob(
                file_base2+'QUANTILES/*_'+trt_cell_id+'_quantiles_*_' +
                datatype+'.csv')

            if not flist:
                warn('No quantiles files found in '+file_base2 +
                     'QUANTILES/ for TRT cell '+trt_cell_id)
                continue

            tbin_edges, qbin_edges, data_ma, start_time = read_quantiles_ts(
                flist, step=5., qmin=0., qmax=100., t_res=None)

            basepath_out = os.path.dirname(flist[0])
            fname = (
                basepath_out+'/'+trt_cell_id+'_trt_QUANTILES_'+datatype +
                '.png')

            vmin = vmax = None
            if datatype == 'RhoHVc':
                vmin = 0.95
                vmax = 1.00
            xlabel = (
                'time (s from '+start_time.strftime("%Y-%m-%d %H:%M:%S") +
                ')')
            _plot_time_range(
                tbin_edges, qbin_edges, data_ma, field_name, [fname],
                titl=titl, xlabel=xlabel, ylabel='Quantile', vmin=vmin,
                vmax=vmax, figsize=[10, 8], dpi=72)

            print("----- plot to '%s'" % fname)

    if 'hydro' in datatype_list:
        fname = args.trtbase+'cell_rimed_particles_column.csv'
        write_trt_cell_lightning(
            cell_ID_list, time_list, lon_list, lat_list, area_list,
            rank_list, rm_hmin_list, rm_hmax_list, fname)

        print("----- written to '%s'" % fname)


def get_graupel_column(tbin_edges, hbin_edges, data_ma, start_time,
                       yyyymmddHHMM):
    """
    Gets the minimum and maximum heigth of the graupel column

    Parameters
    ----------
    tbin_edges : 1D array of floats
        The time bin edges [s]
    hbin_edges : 1D array of floats
        The height bin edges [m MSL]
    data_ma : 2D array of ints
        Matrix containing the time-height hydrometeor classification
        information
    start_time : datetime object
        start time of the radar data
    yyyymmddHHMM : 1D array of datetime objects
        time steps of the TRT cell

    Returns
    -------
    hmin, hmax : 1D float arrays
        the minimum and maximum altitude of the rimed particles column

    """
    tbin_rights = tbin_edges[1:]
    hbin_lefts = hbin_edges[:-1]
    hbin_rights = hbin_edges[1:]
    H_lefts, _ = np.meshgrid(hbin_lefts, tbin_rights)
    H_lefts = np.ma.asarray(H_lefts)
    H_rights, _ = np.meshgrid(hbin_rights, tbin_rights)
    H_rights = np.ma.asarray(H_rights)
    mask = np.logical_not(np.logical_or(data_ma == 4, data_ma == 9))
    H_lefts[mask] = np.ma.masked
    H_rights[mask] = np.ma.masked
    hmin = np.ma.min(H_lefts, axis=1)
    hmax = np.ma.max(H_rights, axis=1)

    # All TRT cell time steps have a corresponding rimed column height value
    # Return the values
    if tbin_rights.size == yyyymmddHHMM.size:
        return hmin, hmax

    # Missing rimed height values. Determine those missing and put the data
    # in the right time step
    hmin_aux = np.ma.masked_all(yyyymmddHHMM.size, dtype=float)
    hmax_aux = np.ma.masked_all(yyyymmddHHMM.size, dtype=float)
    for k, dt_trt in enumerate(yyyymmddHHMM):
        for l, dt_rad in enumerate(tbin_rights):
            if dt_trt == start_time + datetime.timedelta(seconds=dt_rad):
                hmin_aux[k] = hmin[l]
                hmax_aux[k] = hmax[l]
                break
    return hmin_aux, hmax_aux


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
