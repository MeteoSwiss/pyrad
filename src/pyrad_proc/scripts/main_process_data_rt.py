#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
Pyrad: The MeteoSwiss Radar Processing framework
================================================

Welcome to Pyrad!

This program performs real time processing of the data

To run the processing framework type:
    python main_process_data_rt.py \
[config_files] --starttime [process_start_time] --endtime [process_end_time] \
--cfgpath [cfgpath] --proc_period [proc_period]

If startime or endtime are specified the program will start processing at the
specified time and end at the specified time. Otherwise the program ends when
the user interrupts it.
cfgpath is an optional argument with default: \
'$HOME/pyrad/config/processing/'
proc_period is the time that has to pass before attempting to restart the
processing in [s]
if proc_finish is not none it indicates the time the program is allowed to ran
berfore forcing it to end


Example:
    python main_process_data_rt.py 'paradiso_fvj_vol.txt' \
'paradiso_fvj_rhi.txt' --starttime '20140523000000' \
--endtime '20140523001000' --cfgpath '$HOME/pyrad/config/processing/' \
--proc_period 60 --proc_finish 120

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import os
import traceback
from warnings import warn

from pyrad.flow.flow_control import main_rt as pyrad_main

print(__doc__)


def main():
    """
    """

    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    # positional arguments
    parser.add_argument(
        'cfgfiles', nargs='+', type=str,
        help='name of main configuration file')

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

    parser.add_argument(
        '--proc_period', type=int, default=60,
        help='Period between processing rounds (s)')

    parser.add_argument(
        '--proc_finish', type=int, default=None,
        help='Processing time allowed before shutdown (s)')

    args = parser.parse_args()

    print("====== PYRAD data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD data processing finished: ")

    print('config path: '+args.cfgpath)
    cfgfile_list = []
    for ind, cfgfile in enumerate(args.cfgfiles):
        print('config file '+str(ind)+': '+cfgfile)
        cfgfile_list.append(args.cfgpath+cfgfile)
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

    end_proc = False
    while not end_proc:
        try:
            end_proc = pyrad_main(
                cfgfile_list, starttime=proc_starttime, endtime=proc_endtime,
                proc_period=args.proc_period, proc_finish=args.proc_finish)
        except:
            traceback.print_exc()
            if args.proc_finish is None:
                warn("An exception occurred. " +
                     "Restarting the real time processing")
            else:
                end_proc = True


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
