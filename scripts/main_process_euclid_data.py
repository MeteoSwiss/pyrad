#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_process_euclid_data
================================================

This program reads EUCLID raw data and plots several features such as sources
position, histograms, etc.

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit

import numpy as np

from pyrad.io import read_meteorage
from pyrad.util import belongs_roi_indices
from pyrad.graph import plot_pos, plot_histogram

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
        help='Dates to process. Format YYYYMMDD')

    # keyword arguments
    parser.add_argument(
        '--basepath', type=str,
        default='/store/msrad/lightning/meteorage/',
        help='name of folder containing the EUCLID lightning data')

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

    print("====== EUCLID data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== EUCLID data processing finished: ")

    day_vec = []
    for day in args.days:
        day_vec.append(datetime.datetime.strptime(day, '%Y%m%d'))

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

    bin_edges_intens = np.arange(-100., 101., 1.)

    lat_all = np.asarray([], dtype=float)
    lon_all = np.asarray([], dtype=float)
    intens_all = np.asarray([], dtype=float)
    intra_all = np.asarray([], dtype=int)
    for day in day_vec:
        day_str = day.strftime('%y%m%d')
        fname = args.basepath+'THX/THX'+day.strftime('%y%j0000')+'.prd'

        print('Reading EUCLID data file '+fname)
        (stroke_time, lon, lat, intens, ns, mode, intra, ax, ki2, ecc, incl,
         sind) = read_meteorage(fname)

        print('N strokes: '+str(stroke_time.size))
        print('IC: '+str(intra[intra == 1].size))
        print('CG: '+str(intra[intra == 0].size))

        inds, is_roi = belongs_roi_indices(lat, lon, roi)

        if is_roi == 'None':
            print('No strokes in ROI')
            continue

        lon_roi = lon[inds]
        lat_roi = lat[inds]
        intens_roi = intens[inds]
        intra_roi = intra[inds]

        # add valid data to list
        lat_all = np.append(lat_all, lat_roi)
        lon_all = np.append(lon_all, lon_roi)
        intens_all = np.append(intens_all, intens_roi)
        intra_all = np.append(intra_all, intra_roi)

        print('N strokes: '+str(lon_roi.size))
        print('IC: '+str(intra_roi[intra_roi == 1].size))
        print('CG: '+str(intra_roi[intra_roi == 0].size))

        # plot position all strokes
        figfname = args.basepath+day_str+'_EUCLID_strokes_pos.png'
        figfname = plot_pos(
            lat_roi, lon_roi, intens_roi, [figfname],
            cb_label='stroke intensity [kA]',
            titl=day_str+' EUCLID stroke position')
        print('Plotted '+' '.join(figfname))

        # plot position IC
        figfname = args.basepath+day_str+'_EUCLID_IC_pos.png'
        figfname = plot_pos(
            lat_roi[intra_roi == 1], lon_roi[intra_roi == 1],
            intens_roi[intra_roi == 1], [figfname],
            cb_label='stroke intensity [kA]',
            titl=day_str+' EUCLID IC position')
        print('Plotted '+' '.join(figfname))

        # plot position CG
        lat_CG = lat_roi[intra_roi == 0]
        lon_CG = lon_roi[intra_roi == 0]
        intens_CG = intens_roi[intra_roi == 0]

        figfname = args.basepath+day_str+'_EUCLID_CG_pos.png'
        figfname = plot_pos(
            lat_CG, lon_CG, intens_CG, [figfname],
            cb_label='stroke intensity [kA]',
            titl=day_str+' EUCLID CG position')
        print('Plotted '+' '.join(figfname))

        # plot position CGp
        figfname = args.basepath+day_str+'_EUCLID_CGp_pos.png'
        figfname = plot_pos(
            lat_CG[intens_CG > 0.], lon_CG[intens_CG > 0.],
            intens_CG[intens_CG > 0.], [figfname],
            cb_label='stroke intensity [kA]',
            titl=day_str+' EUCLID CGp position')
        print('Plotted '+' '.join(figfname))

        # plot position CGn
        figfname = args.basepath+day_str+'_EUCLID_CGn_pos.png'
        figfname = plot_pos(
            lat_CG[intens_CG < 0.], lon_CG[intens_CG < 0.],
            intens_CG[intens_CG < 0.], [figfname],
            cb_label='stroke intensity [kA]',
            titl=day_str+' EUCLID CGn position')
        print('Plotted '+' '.join(figfname))

        # Plot histogram intensity all strokes
        fname_hist = args.basepath+day_str+'_EUCLID_strokes_hist_intens.png'
        fname_hist = plot_histogram(
            bin_edges_intens, intens_roi, [fname_hist], labelx='Intensity [kA]',
            titl=day_str+' EUCLID stroke intensity')
        print('Plotted '+' '.join(fname_hist))

        # Plot histogram intensity IC
        fname_hist = args.basepath+day_str+'_EUCLID_IC_hist_intens.png'
        fname_hist = plot_histogram(
            bin_edges_intens, intens_roi[intra_roi == 1], [fname_hist],
            labelx='Intensity [kA]',
            titl=day_str+' EUCLID IC intensity')
        print('Plotted '+' '.join(fname_hist))

        # Plot histogram intensity CG
        fname_hist = args.basepath+day_str+'_EUCLID_CG_hist_intens.png'
        fname_hist = plot_histogram(
            bin_edges_intens, intens_roi[intra_roi == 0], [fname_hist],
            labelx='Intensity [kA]',
            titl=day_str+' EUCLID CG intensity')

    lat_CG = lat_all[intra_all == 0]
    lon_CG = lon_all[intra_all == 0]
    intens_CG = intens_all[intra_all == 0]

    print('N strokes: '+str(lon_all.size))
    print('IC: '+str(intra_all[intra_all == 1].size))
    print('CG: '+str(lon_CG.size))
    print('CGp: '+str(lon_CG[intens_CG > 0].size))
    print('CGn: '+str(lon_CG[intens_CG < 0].size))

    # plot position all strokes
    figfname = args.basepath+'EUCLID_strokes_pos.png'
    figfname = plot_pos(
        lat_all, lon_all, intens_all, [figfname],
        cb_label='stroke intensity [kA]',
        titl=day_str+' EUCLID stroke position')
    print('Plotted '+' '.join(figfname))

    # plot position IC
    figfname = args.basepath+'EUCLID_IC_pos.png'
    figfname = plot_pos(
        lat_all[intra_all == 1], lon_all[intra_all == 1],
        intens_all[intra_all == 1], [figfname],
        cb_label='stroke intensity [kA]',
        titl='EUCLID IC position')
    print('Plotted '+' '.join(figfname))

    # plot position CG
    figfname = args.basepath+'EUCLID_CG_pos.png'
    figfname = plot_pos(
        lat_CG, lon_CG, intens_CG, [figfname],
        cb_label='stroke intensity [kA]',
        titl='EUCLID CG position')
    print('Plotted '+' '.join(figfname))

    # plot position CGp
    figfname = args.basepath+'EUCLID_CGp_pos.png'
    figfname = plot_pos(
        lat_CG[intens_CG > 0.], lon_CG[intens_CG > 0.],
        intens_CG[intens_CG > 0.], [figfname],
        cb_label='stroke intensity [kA]',
        titl='EUCLID CGp position')
    print('Plotted '+' '.join(figfname))

    # plot position CGn
    figfname = args.basepath+'EUCLID_CGn_pos.png'
    figfname = plot_pos(
        lat_CG[intens_CG < 0.], lon_CG[intens_CG < 0.],
        intens_CG[intens_CG < 0.], [figfname],
        cb_label='stroke intensity [kA]',
        titl='EUCLID CGn position')
    print('Plotted '+' '.join(figfname))

    # Plot histogram intensity all strokes
    fname_hist = args.basepath+'EUCLID_strokes_hist_intens.png'
    fname_hist = plot_histogram(
        bin_edges_intens, intens_all, [fname_hist], labelx='Intensity [kA]',
        titl='EUCLID stroke intensity')
    print('Plotted '+' '.join(fname_hist))

    # Plot histogram intensity IC
    fname_hist = args.basepath+'EUCLID_IC_hist_intens.png'
    fname_hist = plot_histogram(
        bin_edges_intens, intens_all[intra_all == 1], [fname_hist],
        labelx='Intensity [kA]',
        titl='EUCLID IC intensity')
    print('Plotted '+' '.join(fname_hist))

    # Plot histogram intensity CG
    fname_hist = args.basepath+'EUCLID_CG_hist_intens.png'
    fname_hist = plot_histogram(
        bin_edges_intens, intens_all[intra_all == 0], [fname_hist],
        labelx='Intensity [kA]',
        titl='EUCLID CG intensity')
    print('Plotted '+' '.join(fname_hist))


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
