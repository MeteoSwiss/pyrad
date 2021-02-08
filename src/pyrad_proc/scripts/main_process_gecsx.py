#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
Pyrad: The MeteoSwiss Radar Processing framework
================================================

Welcome to Pyrad!

This program computes the radar visibility using the GECSX algorithm

References
    ----------
    Gabella, M., & Perona, G. (1998). Simulation of the Orographic Influence
    on Weather Radar Using a Geometric–Optics Approach, Journal of Atmospheric
    and Oceanic Technology, 15(6), 1485-1494.

To run the processing framework type:
    python main_process_data.py \
[config_file] --starttime [process_start_time] --endtime [process_end_time] \
--cfgpath [cfgpath] --gatherplots

cfgpath is an optional argument with default: \
'$HOME/pyrad/config/processing/'

if gatherplots is set to 1, all generated figures will be copied into a
new directory called "ALL_FIGURES" located in the output folder. This is
convenient since GECSX can generate many figures and they are placed by Pyrad
in separate folders.

There are two ways to use this program:
    1. By providing it with a set of valid radar scans, in this case the
    "datapath" entry in the main config file must be defined and a valid radar
    scan must be located in '<datapath>/<scanname>/<YYYY-MM-DD>/
    <YYYYMMDDHHMMSS00datatype>.<ext>', <scanname> can be any name, and datatype
    must correspond to the datatype entry for the GECSX dataset in the prod
    config file. You also need to provide starttime and endtime that include
    the timestamp of the radar scan
    2. Without radar data, by providing the following entries in the prod
    file for the GECSX dataset (you can choose any value, these are examples)
     rmax         FLOAT 50000.    # [m] maximum range
     azmin        FLOAT 0.        # [deg] minimum azimuth angle
     azmax        FLOAT 360.      # [deg] maximum azimuth angle
     anglestep    FLOAT 1.       # [deg] azimuth angle step
     range_resolution FLOAT 50. # [deg] range resolution
     antenna_elevations FLTARR 2 # deg
            0.7
            3.0
    as well as the following entries in the loc file (again choose any value)
    RadarPosition STRUCT 3
        latitude FLOAT 46.842473
        longitude FLOAT 6.918370
        altitude FLOAT 449.5

See the two examples pay_main_DX50.txt and pay_main_norad.txt in
$HOME/pyrad/config/gecsx/

Example:
    python main_process_gecsx.py pay_main_norad.txt
--cfgpath $HOME/pyrad/config/gecsx/ --gatherplots 1

    python main_process_gecsx.py pay_main_DX50.txt --starttime \
'20160101000000' --endtime '20170101001000' --cfgpath $HOME/pyrad/config/gecsx/
--gatherplots 1

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import os

from pyrad.flow.flow_control import main_gecsx

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
        '--starttime', type=str, default=None,
        help=('starting time of the data to be processed. ' +
              'Format ''YYYYMMDDhhmmss'''))
    parser.add_argument(
        '--endtime', type=str, default=None,
        help='end time of the data to be processed. Format ''YYYYMMDDhhmmss''')
    parser.add_argument(
        '--cfgpath', type=str,
        default=os.path.expanduser('~')+'/pyrad/config/processing/',
        help='configuration file path')
    parser.add_argument("-i", "--infostr", type=str,
                        help="Information string about the actual data "
                        "processing (e.g. 'RUN57'). This string is added "
                        "to the filenames of the product files.",
                        default="")
    parser.add_argument("-g", "--gatherplots", type=int,
                        help="If set to 1 will create a folder called ALL_FIGURES "
                        "in the output folder as defined by saveimgbasepath "
                        "in the main config file, and will copy all generated "
                        "figures in this folder (for convenience)",
                        default=1)
    args = parser.parse_args()

    print("====== PYRAD data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD data processing finished: ")

    print('config path: '+args.cfgpath)
    print('config file: '+args.proc_cfgfile)
    if args.starttime is not None:
        print('start time: '+args.starttime)
    else:
        print('start time not defined by user')
    if args.endtime is not None:
        print('end time: '+args.endtime)
    else:
        print('end time not defined by user')

    proc_starttime = None
    if args.starttime is not None:
        proc_starttime = datetime.datetime.strptime(
            args.starttime, '%Y%m%d%H%M%S')
    proc_endtime = None
    if args.endtime is not None:
        proc_endtime = datetime.datetime.strptime(args.endtime, '%Y%m%d%H%M%S')
    cfgfile_proc = args.cfgpath+args.proc_cfgfile
    gatherplots = args.gatherplots
    if args.infostr == 'None':
        infostr = ''
    else:
        infostr = args.infostr

    main_gecsx(cfgfile_proc, starttime=proc_starttime, endtime=proc_endtime,
               infostr=infostr, gather_plots=gatherplots)

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
