#!/usr/bin/env python

"""
================================================
Pyrad: The MeteoSwiss Radar Processing framework
================================================

Welcome to Pyrad!

This program does the daily processing and post-processing over a period of \
time.

To run the processing framework type:
    python main_process_data_period.py \
[config_file] [process_start_date] [process_end_date] \
--starttime [process_start_time] --endtime [process_end_time] \
--postproc_cfgfile [postproc_config_file] --cfgpath [cfgpath]

starttime is an optional argument with default: '000000'
endtime is an optional argument with default: '235959'
postproc_cfgfile is an optional argument with default: None
cfgpath is an optional argument with default: \
'$HOME/pyrad/config/processing/'

Example:
    python main_process_data_period.py 'paradiso_fvj_vol.txt' '20140523' \
'20140525' --starttime '000000' --endtime '001000' \
--postproc_cfgfile 'mals_emm_vol_postproc.txt' \
--cfgpath '$HOME/pyrad/config/processing/'

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import os

from pyrad.flow import main as pyrad_main

print(__doc__)


def main():
    """
    """

    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    parser.add_argument(
        'proc_cfgfile', type=str, help='name of main configuration file')
    parser.add_argument(
        'startdate', type=str,
        help='starting date of the data to be processed. Format ''YYYYMMDD'' ')
    parser.add_argument(
        'enddate', type=str,
        help='end date of the data to be processed. Format ''YYYYMMDD'' ')

    # keyword arguments
    parser.add_argument(
        '--starttime', type=str, default='000000',
        help='starting date of the data to be processed. Format ''hhmmss'' ')
    parser.add_argument(
        '--endtime', type=str, default='235959',
        help='end date of the data to be processed. Format ''hhmmss'' ')

    parser.add_argument("-i", "--infostr",
                        help="Information string about the actual data "
                        "processing (e.g. 'RUN57'). This string is added "
                        "to the filenames of the product files.",
                        default="")

    parser.add_argument("--MULTIPROCESSING_DSET", type=int, default=0,
                        help="If 1 the generation of the datasets at the "
                        "same processing level will be parallelized")
    parser.add_argument("--MULTIPROCESSING_PROD", type=int, default=0,
                        help="If 1 the generation of the products of each "
                        "dataset will be parallelized")
    parser.add_argument("--PROFILE_MULTIPROCESSING", type=int, default=0,
                        help="If 1 the multiprocessing is profiled")

    parser.add_argument(
        '--postproc_cfgfile', type=str, default=None,
        help='name of main post-processing configuration file')
    parser.add_argument(
        '--cfgpath', type=str,
        default=os.path.expanduser('~')+'/pyrad/config/processing/',
        help='configuration file path')

    args = parser.parse_args()

    print("====== PYRAD data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD data processing finished: ")

    print('config path: '+args.cfgpath)
    print('config file: '+args.proc_cfgfile)
    print('postproc config file: '+str(args.postproc_cfgfile))
    print('start date: '+args.startdate)
    print('end date: '+args.enddate)
    print('start time each day: '+args.starttime)
    print('end time each day: '+args.endtime)
    if args.MULTIPROCESSING_DSET:
        print('Dataset generation will be parallelized')
    if args.MULTIPROCESSING_PROD:
        print('Product generation will be parallelized')
    if args.PROFILE_MULTIPROCESSING:
        print('Parallel processing performance will be profiled')

    proc_startdate = datetime.datetime.strptime(
        args.startdate, '%Y%m%d')
    proc_enddate = datetime.datetime.strptime(
        args.enddate, '%Y%m%d')
    proc_starttime = datetime.timedelta(
        hours=float(args.starttime[0:2]), minutes=float(args.starttime[2:4]),
        seconds=float(args.starttime[4:6]))
    proc_endtime = datetime.timedelta(
        hours=float(args.endtime[0:2]), minutes=float(args.endtime[2:4]),
        seconds=float(args.endtime[4:6]))

    cfgfile_proc = args.cfgpath+args.proc_cfgfile
    if args.postproc_cfgfile is not None:
        cfgfile_postproc = args.cfgpath+args.postproc_cfgfile

    ndays = (proc_enddate - proc_startdate).days + 1
    print('Number of days to process: '+str(ndays)+'\n\n')

    if args.infostr == 'None':
        infostr = ''
    else:
        infostr = args.infostr

    for day in range(ndays):
        current_date = proc_startdate + datetime.timedelta(days=day)
        proc_startdatetime = current_date + proc_starttime
        proc_enddatetime = current_date + proc_endtime
        try:
            pyrad_main(cfgfile_proc, starttime=proc_startdatetime,
                       endtime=proc_enddatetime, infostr=infostr,
                       MULTIPROCESSING_DSET=args.MULTIPROCESSING_DSET,
                       MULTIPROCESSING_PROD=args.MULTIPROCESSING_PROD,
                       PROFILE_MULTIPROCESSING=args.PROFILE_MULTIPROCESSING)
            if args.postproc_cfgfile is not None:
                pyrad_main(cfgfile_postproc, starttime=proc_startdatetime,
                           endtime=proc_enddatetime, infostr=infostr,
                           MULTIPROCESSING_DSET=args.MULTIPROCESSING_DSET,
                           MULTIPROCESSING_PROD=args.MULTIPROCESSING_PROD,
                           PROFILE_MULTIPROCESSING=args.PROFILE_MULTIPROCESSING)
        except ValueError:
            print(ValueError)


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
