#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
Pyrad: The MeteoSwiss Radar Processing framework
================================================

Welcome to Pyrad!

This program processes and post-processes data over a time span

To run the processing framework type:
    python main_process_data.py \
[config_file] --starttime [process_start_time] --endtime [process_end_time] \
--postproc_cfgfile [postproc_config_file] --cfgpath [cfgpath]

If startime and endtime are not specified the program determines them from
the trajectory file or the last processed volume.
postproc_cfgfile is an optional argument with default: None
cfgpath is an optional argument with default: \
'$HOME/pyrad/config/processing/'
The trajectory file can be of type plane, lightning or proc_periods. If it is \
of type lightning the flash number can be specified

Example:
    python main_process_data.py 'paradiso_fvj_vol.txt' --starttime \
'20140523000000' --endtime '20140523001000' --postproc_cfgfile \
'paradiso_fvj_vol_postproc.txt' --cfgpath '$HOME/pyrad/config/processing/'

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import os

from pyrad.flow.flow_control import main as pyrad_main

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
        '--postproc_cfgfile', type=str, default=None,
        help='name of main post-processing configuration file')
    parser.add_argument(
        '--cfgpath', type=str,
        default=os.path.expanduser('~')+'/pyrad/config/processing/',
        help='configuration file path')
    parser.add_argument("-i", "--infostr", type=str,
                        help="Information string about the actual data "
                        "processing (e.g. 'RUN57'). This string is added "
                        "to the filenames of the product files.",
                        default="")
    parser.add_argument("-t", "--trajfile", type=str, default='',
                        help="Definition file of plane trajectory. "
                        "Configuration of scan sector, products, ...")
    parser.add_argument("--trajtype", type=str, default='plane',
                        help="Type of trajectory. "
                        "Can be either 'plane', 'lightning' or 'proc_periods'")
    parser.add_argument("--flashnr", type=int, default=0,
                        help="If type of trajectory is 'lightning', "
                        "flash number the data of which will be processed"
                        "0 means that all lightning data will be processed")
    parser.add_argument("--MULTIPROCESSING_DSET", type=int, default=0,
                        help="If 1 the generation of the datasets at the "
                        "same processing level will be parallelized")
    parser.add_argument("--MULTIPROCESSING_PROD", type=int, default=0,
                        help="If 1 the generation of the products of each "
                        "dataset will be parallelized")
    parser.add_argument("--PROFILE_MULTIPROCESSING", type=int, default=0,
                        help="If 1 the multiprocessing is profiled")

    args = parser.parse_args()

    print("====== PYRAD data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD data processing finished: ")

    print('config path: '+args.cfgpath)
    print('config file: '+args.proc_cfgfile)
    print('postproc config file: '+str(args.postproc_cfgfile))
    if args.starttime is not None:
        print('start time: '+args.starttime)
    else:
        print('start time not defined by user')
    if args.endtime is not None:
        print('end time: '+args.endtime)
    else:
        print('end time not defined by user')
    if args.MULTIPROCESSING_DSET:
        print('Dataset generation will be parallelized')
    if args.MULTIPROCESSING_PROD:
        print('Product generation will be parallelized')
    if args.PROFILE_MULTIPROCESSING:
        print('Parallel processing performance will be profiled')

    proc_starttime = None
    if args.starttime is not None:
        proc_starttime = datetime.datetime.strptime(
            args.starttime, '%Y%m%d%H%M%S')
    proc_endtime = None
    if args.endtime is not None:
        proc_endtime = datetime.datetime.strptime(args.endtime, '%Y%m%d%H%M%S')
    cfgfile_proc = args.cfgpath+args.proc_cfgfile

    if args.infostr == 'None':
        infostr = ''
    else:
        infostr = args.infostr

    pyrad_main(cfgfile_proc, starttime=proc_starttime, endtime=proc_endtime,
               trajfile=args.trajfile, infostr=infostr,
               trajtype=args.trajtype, flashnr=args.flashnr,
               MULTIPROCESSING_DSET=args.MULTIPROCESSING_DSET,
               MULTIPROCESSING_PROD=args.MULTIPROCESSING_PROD,
               PROFILE_MULTIPROCESSING=args.PROFILE_MULTIPROCESSING)

    if args.postproc_cfgfile is not None:
        cfgfile_postproc = args.cfgpath+args.postproc_cfgfile
        pyrad_main(cfgfile_postproc, starttime=proc_starttime,
                   endtime=proc_endtime, trajfile=args.trajfile,
                   infostr=infostr, trajtype=args.trajtype,
                   flashnr=args.flashnr,
                   MULTIPROCESSING_DSET=args.MULTIPROCESSING_DSET,
                   MULTIPROCESSING_PROD=args.MULTIPROCESSING_PROD,
                   PROFILE_MULTIPROCESSING=args.PROFILE_MULTIPROCESSING)


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
