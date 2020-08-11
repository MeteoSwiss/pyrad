#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_process_data_windmill_hist_all
================================================

This program puts together histograms of windmill radar returns

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

from pyrad.io import read_histogram, get_fieldname_pyart
from pyrad.io import write_histogram
from pyrad.graph import plot_histogram2, get_colobar_label

from pyart.config import get_metadata

print(__doc__)


def main():
    """
    """

    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    # keyword arguments
    parser.add_argument(
        '--database', type=str,
        default='/store/msrad/radar/pyrad_products/',
        help='base path to the radar data')

    parser.add_argument(
        '--datadirs', type=str,
        default=(
            'mals_sha_windmills_point_psr_filtered_WM1_20200304-20200311,'
            'mals_sha_windmills_point_psr_filtered_WM1_20200312-20200315,'
            'mals_sha_windmills_point_psr_filtered_WM1_20200316-20200320,'
            'mals_sha_windmills_point_psr_filtered_WM1_20200321-20200325'),
        help='directories containing data')

    parser.add_argument(
        '--datatypes', type=str,
        default='dBuZ,dBuZv,rcs_h,rcs_v,ZDRu,RhoHVu,uPhiDPu,Vu,Wu',
        help='Data types. Coma separated')

    args = parser.parse_args()

    print("====== PYRAD windmill data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD windmill data processing finished: ")

    datadirs = args.datadirs.split(',')
    datatypes = args.datatypes.split(',')

    # Read periods of processing
    for datatype in datatypes:
        first_read = False
        for datadir in datadirs:
            # Read data time series files
            flist = glob.glob(
                args.database+datadir+'/'+datatype +
                '_TS/TS/ts_POINT_MEASUREMENT_hist_'+datatype+'.csv')
            if not flist:
                continue

            hist_aux , bin_edges_aux = read_histogram(flist[0])
            if not first_read:
                hist = hist_aux
                bin_edges = bin_edges_aux
                first_read = True
                continue
            
            hist += hist_aux

        basepath = os.path.dirname(flist[0])+'/'

        # Histogram plots
        field_name = get_fieldname_pyart(datatype)
        field_dict = get_metadata(field_name)

        fname = args.database+'ts_POINT_MEASUREMENT_hist_'+datatype+'.png'

        bin_centers = bin_edges[:-1]+((bin_edges[1]-bin_edges[0])/2.)
        fname = plot_histogram2(
            bin_centers, hist, [fname],
            labelx=get_colobar_label(field_dict, field_name), titl=datatype)
        print('Plotted '+' '.join(fname))

        fname = args.database+'ts_POINT_MEASUREMENT_hist_'+datatype+'.csv'
        fname = write_histogram(bin_edges, hist, fname)
        print('Written '+fname)


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
