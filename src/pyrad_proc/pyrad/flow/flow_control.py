"""
pyrad.flow.flow_control
=======================

functions to control the Pyrad data processing flow

.. autosummary::
    :toctree: generated/

    main
    main_rt

"""
from __future__ import print_function
import warnings
from warnings import warn
import os
from datetime import datetime
from datetime import timedelta
import gc
import queue
import time

from pyart import version as pyart_version
from pyrad import version as pyrad_version

from .flow_aux import _warning_format, _initialize_listener
from .flow_aux import _create_cfg_dict, _create_datacfg_dict
from .flow_aux import _get_times_and_traj, _get_datatype_list
from .flow_aux import _get_datasets_list, _get_masterfile_list
from .flow_aux import _wait_for_files, _get_radars_data
from .flow_aux import _initialize_datasets
from .flow_aux import _process_datasets, _postprocess_datasets

from ..io.io_aux import get_datetime
from ..io.read_data_other import read_last_state
from ..io.write_data import write_last_state

ALLOW_USER_BREAK = False

try:
    from dask.diagnostics import Profiler, ResourceProfiler, CacheProfiler
    from dask.diagnostics import visualize
    from distributed import Client
    from bokeh.io import export_png
    _DASK_AVAILABLE = True
except ImportError:
    warn('dask not available: The processing will not be parallelized')
    _DASK_AVAILABLE = False


def main(cfgfile, starttime=None, endtime=None, trajfile="", trajtype='plane',
         flashnr=0, infostr="", MULTIPROCESSING_DSET=False,
         MULTIPROCESSING_PROD=False, PROFILE_MULTIPROCESSING=False):
    """
    Main flow control. Processes radar data off-line over a period of time
    given either by the user, a trajectory file, or determined by the last
    volume processed and the current time. Multiple radars can be processed
    simultaneously

    Parameters
    ----------
    cfgfile : str
        path of the main config file
    starttime, endtime : datetime object
        start and end time of the data to be processed
    trajfile : str
        path to file describing the trajectory
    trajtype : str
        type of trajectory file. Can be either 'plane', 'lightning' or
        'proc_periods'
    flashnr : int
        If larger than 0 will select a flash in a lightning trajectory file.
        If 0 the data corresponding to the trajectory of all flashes will be
        plotted
    infostr : str
        Information string about the actual data processing
        (e.g. 'RUN57'). This string is added to product files.
    MULTIPROCESSING_DSET : Bool
        If true the generation of datasets at the same processing level will
        be parallelized
    MULTIPROCESSING_PROD : Bool
        If true the generation of products from each dataset will be
        parallelized
    PROFILE_MULTIPROCESSING : Bool
        If true and code parallelized the multiprocessing is profiled

    """
    print("- PYRAD version: %s (compiled %s by %s)" %
          (pyrad_version.version, pyrad_version.compile_date_time,
           pyrad_version.username))
    print("- PYART version: " + pyart_version.version)

    # Define behaviour of warnings
    warnings.simplefilter('always')  # always print matching warnings
    # warnings.simplefilter('error')  # turn matching warnings into exceptions
    warnings.formatwarning = _warning_format  # define format

    if ALLOW_USER_BREAK:
        input_queue = _initialize_listener()

    if not _DASK_AVAILABLE:
        MULTIPROCESSING_DSET = False
        MULTIPROCESSING_PROD = False
        PROFILE_MULTIPROCESSING = False

    # check if multiprocessing profiling is necessary
    if not MULTIPROCESSING_DSET and not MULTIPROCESSING_PROD:
        PROFILE_MULTIPROCESSING = False
    elif MULTIPROCESSING_DSET and MULTIPROCESSING_PROD:
        PROFILE_MULTIPROCESSING = False

    if MULTIPROCESSING_DSET and MULTIPROCESSING_PROD:
        # necessary to launch tasks from tasks
        Client()

    if PROFILE_MULTIPROCESSING:
        prof = Profiler()
        rprof = ResourceProfiler()
        cprof = CacheProfiler()

        prof.register()
        rprof.register()
        cprof.register()

    cfg = _create_cfg_dict(cfgfile)
    datacfg = _create_datacfg_dict(cfg)

    starttimes, endtimes, traj = _get_times_and_traj(
        trajfile, starttime, endtime, cfg['ScanPeriod'],
        last_state_file=cfg['lastStateFile'], trajtype=trajtype,
        flashnr=flashnr)

    if infostr:
        print('- Info string : ' + infostr)

    # get data types and levels
    datatypesdescr_list = list()
    for i in range(1, cfg['NumRadars']+1):
        datatypesdescr_list.append(
            _get_datatype_list(cfg, radarnr='RADAR'+'{:03d}'.format(i)))

    dataset_levels = _get_datasets_list(cfg)

    masterfilelist, masterdatatypedescr, masterscan = _get_masterfile_list(
        datatypesdescr_list[0], starttimes, endtimes, datacfg,
        scan_list=datacfg['ScanList'])

    nvolumes = len(masterfilelist)
    if nvolumes == 0:
        raise ValueError(
            "ERROR: Could not find any valid volumes between " +
            starttimes[0].strftime('%Y-%m-%d %H:%M:%S') + " and " +
            endtimes[-1].strftime('%Y-%m-%d %H:%M:%S') + " for " +
            "master scan '" + str(masterscan) +
            "' and master data type '" + masterdatatypedescr +
            "'")
    print('- Number of volumes to process: ' + str(nvolumes))
    print('- Start time: ' + starttimes[0].strftime("%Y-%m-%d %H:%M:%S"))
    print('- end time: ' + endtimes[-1].strftime("%Y-%m-%d %H:%M:%S"))

    # initial processing of the datasets
    print('\n\n- Initializing datasets:')
    dscfg, traj = _initialize_datasets(
        dataset_levels, cfg, traj=traj, infostr=infostr)

    # process all data files in file list or until user interrupts processing
    for masterfile in masterfilelist:
        if ALLOW_USER_BREAK:
            # check if user has requested exit
            try:
                input_queue.get_nowait()
                warn('Program terminated by user')
                break
            except queue.Empty:
                pass

        print('\n- master file: ' + os.path.basename(masterfile))

        master_voltime = get_datetime(masterfile, masterdatatypedescr)

        radar_list = _get_radars_data(
            master_voltime, datatypesdescr_list, datacfg,
            num_radars=datacfg['NumRadars'])

        # process all data sets
        dscfg, traj = _process_datasets(
            dataset_levels, cfg, dscfg, radar_list, master_voltime, traj=traj,
            infostr=infostr, MULTIPROCESSING_DSET=MULTIPROCESSING_DSET,
            MULTIPROCESSING_PROD=MULTIPROCESSING_PROD)

        # delete variables
        del radar_list

        gc.collect()

    # post-processing of the datasets
    print('\n\n- Post-processing datasets:')
    dscfg, traj = _postprocess_datasets(
        dataset_levels, cfg, dscfg, traj=traj, infostr=infostr)

    if PROFILE_MULTIPROCESSING:
        prof.unregister()
        rprof.unregister()
        cprof.unregister()

        bokeh_plot = visualize([prof, rprof, cprof], show=False, save=False)

        profile_path = os.path.expanduser('~')+'/profiling/'
        if not os.path.isdir(profile_path):
            os.makedirs(profile_path)

        export_png(bokeh_plot, filename=(
            profile_path+datetime.utcnow().strftime('%Y%m%d%H%M%S') +
            '_profile.png'))

    print('- This is the end my friend! See you soon!')


def main_rt(cfgfile_list, starttime=None, endtime=None, infostr_list=None,
            proc_period=60, proc_finish=None):
    """
    main flow control. Processes radar data in real time. The start and end
    processing times can be determined by the user. This function is inteded
    for a single radar

    Parameters
    ----------
    cfgfile_list : list of str
        path of the main config files
    starttime, endtime : datetime object
        start and end time of the data to be processed
    infostr_list : list of str
        Information string about the actual data processing
        (e.g. 'RUN57'). This string is added to product files.
    proc_period : int
        period of time before starting a new processing round (seconds)
    cronjob_controlled : Boolean
        If True means that the program is started periodically from a cronjob
        and therefore finishes execution after processing
    proc_finish : int or None
        if set to a value the program will be forced to shut down after the
        value (in seconds) from start time has been exceeded

    Returns
    -------
    end_proc : Boolean
        If true the program has ended successfully

    """
    print("- PYRAD version: %s (compiled %s by %s)" %
          (pyrad_version.version, pyrad_version.compile_date_time,
           pyrad_version.username))
    print("- PYART version: " + pyart_version.version)

    # Define behaviour of warnings
    warnings.simplefilter('always')  # always print matching warnings
    # warnings.simplefilter('error')  # turn matching warnings into exceptions
    warnings.formatwarning = _warning_format  # define format

    # The processing will be allowed to run for a limited period
    if proc_finish is not None:
        startime_proc = datetime.utcnow()
        # for offline testing
        # startime_proc = startime_proc.replace(
        #     year=endtime.year, month=endtime.month, day=endtime.day)
        # startime_proc = startime_proc.replace(hour=10)

        endtime_proc = startime_proc+timedelta(seconds=proc_finish)

    if ALLOW_USER_BREAK:
        input_queue = _initialize_listener()

    cfg_list = []
    datacfg_list = []
    dscfg_list = []
    datatypesdescr_list_list = []
    dataset_levels_list = []
    last_processed_list = []

    for icfg, cfgfile in enumerate(cfgfile_list):
        cfg = _create_cfg_dict(cfgfile)
        if infostr_list is not None:
            infostr = infostr_list[icfg]
        else:
            infostr = ""
        datacfg = _create_datacfg_dict(cfg)

        if infostr:
            print('- Info string : ' + infostr)

        # find out last processed volume
        last_processed = read_last_state(cfg['lastStateFile'])
        if last_processed is None:
            print('- last processed volume unknown')
        else:
            print('- last processed volume: '+last_processed.strftime(
                '%Y%m%d%H%M%S'))
        last_processed_list.append(last_processed)

        # get data types and levels
        datatypesdescr_list = list()
        for i in range(1, cfg['NumRadars']+1):
            datatypesdescr_list.append(
                _get_datatype_list(cfg, radarnr='RADAR'+'{:03d}'.format(i)))

        dataset_levels = _get_datasets_list(cfg)

        # initial processing of the datasets
        print('\n\n- Initializing datasets:')
        dscfg, traj = _initialize_datasets(
            dataset_levels, cfg, infostr=infostr)

        cfg_list.append(cfg)
        datacfg_list.append(datacfg)
        dscfg_list.append(dscfg)
        datatypesdescr_list_list.append(datatypesdescr_list)
        dataset_levels_list.append(dataset_levels)

        # remove variables from memory
        del cfg
        del datacfg
        del dscfg
        del datatypesdescr_list
        del dataset_levels
        del last_processed
        del traj

        gc.collect()

    end_proc = False
    while not end_proc:
        if ALLOW_USER_BREAK:
            # check if user has requested exit
            try:
                user_input = input_queue.get_nowait()
                end_proc = user_input
                warn('Program terminated by user')
                break
            except queue.Empty:
                pass

        nowtime = datetime.utcnow()
        # for offline testing
        # nowtime = nowtime.replace(
        #     year=endtime.year, month=endtime.month, day=endtime.day)
        # nowtime = nowtime.replace(hour=10)

        # if processing end time exceeded finalize processing
        if proc_finish is not None:
            if nowtime >= endtime_proc:
                end_proc = True
                warn('Allowed processing time exceeded')
                break

        # end time has been set and current time older than end time
        # quit processing
        if endtime is not None:
            if nowtime > endtime:
                end_proc = True
                break

        # start time has been set. Check if current time has to be
        # processed. If not sleep until next proc_period
        if starttime is not None:
            if nowtime < starttime:
                time.sleep(proc_period)
                continue

        vol_processed = False
        for icfg, cfg in enumerate(cfg_list):
            if ALLOW_USER_BREAK:
                # check if user has requested exit
                try:
                    user_input = input_queue.get_nowait()
                    end_proc = user_input
                    warn('Program terminated by user')
                    break
                except queue.Empty:
                    pass

            datacfg = datacfg_list[icfg]
            dscfg = dscfg_list[icfg]
            datatypesdescr_list = datatypesdescr_list_list[icfg]
            dataset_levels = dataset_levels_list[icfg]
            last_processed = last_processed_list[icfg]
            if infostr_list is not None:
                infostr = infostr_list[icfg]
            else:
                infostr = ""

            # wait until new files are available
            masterfile, masterdatatypedescr, last_processed = _wait_for_files(
                nowtime, datacfg, datatypesdescr_list[0],
                last_processed=last_processed)
            if masterfile is None:
                last_processed_list[icfg] = last_processed
                if last_processed is not None:
                    write_last_state(last_processed, cfg['lastStateFile'])
                continue

            print('\n- master file: ' + os.path.basename(masterfile))
            master_voltime = get_datetime(masterfile, masterdatatypedescr)

            # get data of master radar
            radar_list = _get_radars_data(
                master_voltime, datatypesdescr_list, datacfg)

            # process all data sets
            dscfg, traj = _process_datasets(
                dataset_levels, cfg, dscfg, radar_list, master_voltime,
                infostr=infostr)

            last_processed_list[icfg] = master_voltime
            write_last_state(master_voltime, cfg['lastStateFile'])
            dscfg_list[icfg] = dscfg

            vol_processed = True

            # remove variables from memory
            del radar_list
            del cfg
            del datacfg
            del dscfg
            del datatypesdescr_list
            del dataset_levels
            del last_processed
            del traj

            gc.collect()

        nowtime_new = datetime.utcnow()
        # for offline testing
        # nowtime_new = nowtime_new.replace(
        #     year=endtime.year, month=endtime.month, day=endtime.day)
        # nowtime_new = nowtime_new.replace(hour=10)

        proc_time = (nowtime_new-nowtime).total_seconds()
        if vol_processed:
            print('Processing time %s s\n' % proc_time)

        # if processing end time exceeded finalize processing
        if proc_finish is not None:
            if nowtime_new >= endtime_proc:
                end_proc = True
                warn('Allowed processing time exceeded')
                break

        if proc_time < proc_period:
            time.sleep(proc_period-proc_time)

    # only do post processing if program properly terminated by user
    if end_proc:
        # post-processing of the datasets
        print('\n\n- Post-processing datasets:')
        for icfg, cfg in enumerate(cfg_list):
            dscfg = dscfg_list[icfg]
            dataset_levels = dataset_levels_list[icfg]
            if infostr_list is not None:
                infostr = infostr_list[icfg]
            else:
                infostr = ""
            dscfg, traj = _postprocess_datasets(
                dataset_levels, cfg, dscfg, infostr=None)

            # remove variables from memory
            del cfg
            del dscfg
            del dataset_levels
            del traj

            gc.collect()

    print('- This is the end my friend! See you soon!')

    return end_proc
