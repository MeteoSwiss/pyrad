"""
pyrad.flow.flow_aux
===================

Auxiliary functions to control the Pyrad data processing flow

.. autosummary::
    :toctree: generated/

    _initialize_listener
    _user_input_listener
    _get_times_and_traj
    _initialize_datasets
    _process_datasets
    _postprocess_datasets
    _wait_for_files
    _get_radars_data
    _generate_dataset
    _generate_prod
    _create_cfg_dict
    _create_datacfg_dict
    _create_dscfg_dict
    _create_prdcfg_dict
    _get_datatype_list
    _get_datasets_list
    _get_masterfile_list
    _add_dataset
    _warning_format

"""
from __future__ import print_function
import sys
from warnings import warn
import traceback
import os
from datetime import datetime
from datetime import timedelta
import inspect
import gc

import queue
import time
import threading
import glob
from copy import deepcopy
import numpy as np

try:
    from memory_profiler import profile as mprofile
    _MPROFILE_AVAILABLE = True
except ImportError:
    warn('Memory profiler not available')
    _MPROFILE_AVAILABLE = False

from pyrad import proc

from ..io.config import read_config
from ..io.read_data_radar import get_data
from ..io.io_aux import get_datetime, get_file_list, get_scan_list
from ..io.io_aux import get_dataset_fields, get_datatype_fields
from ..io.io_aux import get_new_rainbow_file_name, get_fieldname_pyart
from ..io.trajectory import Trajectory
from ..io.read_data_other import read_last_state, read_proc_periods

from ..proc.process_aux import get_process_func
from ..prod.product_aux import get_prodgen_func

try:
    import dask
except ImportError:
    warn('dask not available: The processing will not be parallelized')

PROFILE_LEVEL = 0


def profiler(level=1):
    """
    Function to be used as decorator for memory debugging. The function will
    be profiled or not according to its level respect to the global variable
    PROFILE_LEVEL

    Parameters
    ----------
    level : int
        profiling level

    Returns
    -------
    func or func wrapper : function
        The function or its wrapper for profiling

    """
    def profile_real_decorator(func):
        """
        real decorator

        Parameters
        ----------
        func : function
            function to profile

        Returns
        -------
        wrapper : function
            The function wrapper

        """
        def wrapper(*args, **kwargs):
            """
            wrapper

            Parameters
            ----------
            args, kwargs : arguments
                The arguments of the function

            Returns_initialize_datasets
            -------
            func : function
                The original function if no profiling has to be performed or
                the function decorated with the memory decorator

            """
            if not _MPROFILE_AVAILABLE:
                return func(*args, **kwargs)

            if ((PROFILE_LEVEL == 1 and level == 1) or
                    (PROFILE_LEVEL == 2 and (level in (1, 2))) or
                    PROFILE_LEVEL == 3):
                print('profiling '+str(func))
                func2 = mprofile(func)
                return func2(*args, **kwargs)
            return func(*args, **kwargs)
        return wrapper
    return profile_real_decorator


def _initialize_listener():
    """
    initialize the input listener

    Returns
    -------
    input_queue : queue object
        the queue object where to put the quit signal

    """
    input_queue = queue.Queue()
    pinput = threading.Thread(
        name='user_input_listener', target=_user_input_listener,
        daemon=True, args=(input_queue, ))
    # input_queue = mp.Queue()
    # pinput = mp.Process(
    #     name='user_input_listener', target=_user_input_listener,
    #      daemon=True, args=(input_queue, ))
    pinput.start()

    return input_queue


def _user_input_listener(input_queue):
    """
    Permanently listens to the keyword input until the user types "Return"

    Parameters
    ----------
    input_queue : queue object
        the queue object where to put the quit signal

    """
    print("Press Enter to quit: ")
    while True:
        user_input = sys.stdin.read(1)
        if '\n' in user_input or '\r' in user_input:
            warn('Exit requested by user')
            input_queue.put(True)
            break
        time.sleep(1)


@profiler(level=3)
def _get_times_and_traj(trajfile, starttime, endtime, scan_period,
                        last_state_file=None, trajtype='plane',
                        flashnr=0):
    """
    Gets the trajectory and the start time and end time if they have
    not been set

    Parameters
    ----------
    trajfile : str
        trajectory file
    starttime, endtime : datetime object or None
        the start and stop times of the processing
    scan_period : float
        the scan period in minutes
    last_state_file : str
        name of the file that stores the time of the last processed volume
    trajtype : str
        type of trajectory. Can be plane, lightning or proc_periods
    flashnr : int
        If type of trajectory is lightning, the flash number. 0 means all
        flash numbers included

    Returns
    -------
    starttimes, endtimes : array of datetime objects
        the start and end times of the data to process
    traj : trajectory object or None
        the trajectory object

    """
    if trajfile:
        if trajtype == 'proc_periods':
            print("- Processing periods file: " + trajfile)
            starttimes, endtimes = read_proc_periods(trajfile)
            traj = None
            if starttimes is None or endtimes is None:
                sys.exit(1)
        else:
            print("- Trajectory file: " + trajfile)
            try:
                traj = Trajectory(
                    trajfile, starttime=starttime, endtime=endtime,
                    trajtype=trajtype, flashnr=flashnr)
            except Exception as inst:
                warn(str(inst))
                sys.exit(1)

            # Derive start and end time (if not specified by arguments)
            if starttime is None:
                scan_min = scan_period * 2  # [min]
                starttimes = np.array(
                    [traj.get_start_time() - timedelta(minutes=scan_min)])
            else:
                starttimes = np.array([starttime])
            if endtime is None:
                scan_min = scan_period * 2  # [min]
                endtimes = np.array(
                    [traj.get_end_time() + timedelta(minutes=scan_min)])
            else:
                endtimes = np.array([endtime])
    else:
        traj = None
        if starttime is None:
            starttimes = None
        else:
            starttimes = np.array([starttime])
        if endtime is None:
            endtimes = None
        else:
            endtimes = np.array([endtime])

    # if start time is not defined and the file lastState exists and
    # contains a valid date start processing from the last valid date.
    # Otherwise start processing from yesterday at 00:00:00 UTC
    if starttimes is None and last_state_file is not None:
        filename = glob.glob(last_state_file)
        if not filename:
            nowtime = datetime.utcnow()
            starttimes = np.array([(nowtime - timedelta(days=1)).replace(
                hour=0, minute=0, second=0, microsecond=0)])
            warn('File '+last_state_file+' not found. ' +
                 'Start time set at ' +
                 starttimes[0].strftime('%Y-%m-%d %H:%M:%S'))
        else:
            starttime = read_last_state(last_state_file)
            if starttime is None:
                nowtime = datetime.utcnow()
                starttimes = np.array([(nowtime - timedelta(days=1)).replace(
                    hour=0, minute=0, second=0, microsecond=0)])
                warn('File '+last_state_file+' not valid. ' +
                     'Start time set at ' +
                     starttimes[0].strftime('%Y-%m-%d %H:%M:%S'))
            else:
                starttimes = np.array([starttime])

    if endtimes is None:
        endtimes = np.array([datetime.utcnow()])
        warn('End Time not defined. Set as ' +
             endtimes[0].strftime('%Y-%m-%d %H:%M:%S'))

    if endtimes[0] < starttimes[0]:
        raise ValueError(
            'Start time '+starttimes[0].strftime('%Y-%m-%d %H:%M:%S') +
            'older than end time ' +
            endtimes[0].strftime('%Y-%m-%d %H:%M:%S'))

    return starttimes, endtimes, traj


def _initialize_datasets(dataset_levels, cfg, traj=None, infostr=None):
    """
    Initializes datasets. Creates the data set configuration dictionary

    Parameters
    ----------
    dataset_levels : dict
        dictionary containing the list of data sets to be generated at each
        processing level
    cfg : dict
        processing configuration dictionary
    traj : trajectory object
        object containing the trajectory
    infostr : str
        Information string about the actual data processing
        (e.g. 'RUN57'). This string is added to product files.

    Returns
    -------
    dscfg : dict
        dictionary containing the configuration data for each dataset
    traj : trajectory object
        the modified trajectory object

    """
    dscfg = dict()
    for level in sorted(dataset_levels):
        print('-- Process level: '+level)
        for dataset in dataset_levels[level]:
            print('--- Processing dataset: '+dataset)
            dscfg.update({dataset: _create_dscfg_dict(cfg, dataset)})
            _, _, _, dscfg[dataset] = _generate_dataset(
                dataset, cfg, dscfg[dataset], proc_status=0,
                radar_list=None, voltime=None, trajectory=traj,
                runinfo=infostr)

            gc.collect()

    # manual garbage collection after initial processing
    gc.collect()

    return dscfg, traj


# @profile
@profiler(level=1)
def _process_datasets(dataset_levels, cfg, dscfg, radar_list, master_voltime,
                      traj=None, infostr=None, MULTIPROCESSING_DSET=False,
                      MULTIPROCESSING_PROD=False):
    """
    Processes the radar volumes for a particular time stamp.

    Parameters
    ----------
    dataset_levels : dict
        dictionary containing the list of data sets to be generated at each
        processing level
    cfg : dict
        processing configuration dictionary
    dscfg : dict
        dictionary containing the configuration data for each dataset
    radar_list : list of radar objects
        The radar objects to be processed
    master_voltime : datetime object
        the reference radar volume time
    traj : trajectory object
        and object containing the trajectory
    infostr : str
        Information string about the actual data processing
        (e.g. 'RUN57'). This string is added to product files.
    MULTIPROCESSING_DSET : Bool
        If true the generation of datasets at the same processing level will
        be parallelized
    MULTIPROCESSING_PROD : Bool
        If true the generation of products from each dataset will be
        parallelized

    Returns
    -------
    dscfg : dict
        the modified configuration dictionary
    traj : trajectory object
        the modified trajectory object

    """
    for level in sorted(dataset_levels):
        print('-- Process level: '+level)
        if MULTIPROCESSING_DSET:
            jobs = []
            make_global_list = []
            substitute_object_list = []
            fields_to_remove_list = []
            for dataset in dataset_levels[level]:
                print('--- Processing dataset: '+dataset)
                make_global_list.append(dscfg[dataset]['MAKE_GLOBAL'])
                substitute_object_list.append(
                    dscfg[dataset]['SUBSTITUTE_OBJECT'])
                fields_to_remove_list.append(
                    dscfg[dataset]['FIELDS_TO_REMOVE'])

                # delay the data hashing
                radar_list_aux = dask.delayed(radar_list)
                dscfg_aux = dask.delayed(dscfg[dataset])
                jobs.append(dask.delayed(_generate_dataset)(
                    dataset, cfg, dscfg_aux, proc_status=1,
                    radar_list=radar_list_aux, voltime=master_voltime,
                    trajectory=traj, runinfo=infostr,
                    MULTIPROCESSING_PROD=MULTIPROCESSING_PROD))

            try:
                jobs = dask.compute(*jobs)

                # add new dataset to radar object if necessary
                # update dataset config dictionary
                for i, job in enumerate(jobs):
                    dscfg[job[2]] = job[3]
                    _add_dataset(
                        job[0][0], radar_list, job[1],
                        make_global=make_global_list[i],
                        substitute_object=substitute_object_list[i],
                        fields_to_remove=fields_to_remove_list[i])

                del jobs
            except Exception as ee:
                warn(str(ee))
                traceback.print_exc()
        else:
            for dataset in dataset_levels[level]:
                print('--- Processing dataset: '+dataset)
                try:
                    new_dataset, ind_rad, _, dscfg[dataset] = (
                        _generate_dataset(
                            dataset, cfg, dscfg[dataset], proc_status=1,
                            radar_list=radar_list, voltime=master_voltime,
                            trajectory=traj, runinfo=infostr,
                            MULTIPROCESSING_PROD=MULTIPROCESSING_PROD))

                    # adds the first dataset generated to the object. Typically
                    # only one dataset is generated but gecsx generates two:
                    # The first one is a radar object and the second is a grid
                    # object.
                    _add_dataset(
                        new_dataset[0], radar_list, ind_rad,
                        make_global=dscfg[dataset]['MAKE_GLOBAL'],
                        substitute_object=dscfg[dataset]['SUBSTITUTE_OBJECT'],
                        fields_to_remove=dscfg[dataset]['FIELDS_TO_REMOVE'])

                    del new_dataset
                    gc.collect()
                except Exception as ee:
                    warn(str(ee))
                    traceback.print_exc()

    # manual garbage collection after processing each radar volume
    gc.collect()

    return dscfg, traj


def _postprocess_datasets(dataset_levels, cfg, dscfg, traj=None, infostr=None):
    """
    Processes the radar volumes for a particular time stamp.

    Parameters
    ----------
    dataset_levels : dict
        dictionary containing the list of data sets to be generated at each
        processing level
    cfg : dict
        processing configuration dictionary
    dscfg : dict
        dictionary containing the configuration data for each dataset
    traj : trajectory object
        and object containing the trajectory
    infostr : str
        Information string about the actual data processing
        (e.g. 'RUN57'). This string is added to product files.

    Returns
    -------
    dscfg : dict
        the modified configuration dictionary
    traj : trajectory object
        the modified trajectory object

    """
    for level in sorted(dataset_levels):
        print('-- Process level: '+level)
        for dataset in dataset_levels[level]:
            print('--- Processing dataset: '+dataset)
            _generate_dataset(
                dataset, cfg, dscfg[dataset], proc_status=2,
                radar_list=None, voltime=None, trajectory=traj,
                runinfo=infostr)

            gc.collect()

    # manual garbage collection after post-processing
    gc.collect()

    return dscfg, traj


def _wait_for_files(nowtime, datacfg, datatype_list, last_processed=None):
    """
    Waits for the master file and all files in a volume scan to be present
    returns the masterfile if the volume scan can be processed.

    Parameters
    ----------
    nowtime : datetime object
        the current time
    datacfg : dict
        dictionary containing the parameters to get the radar data
    last_processed : datetime or None
        The end time of the previously processed radar volume

    Returns
    -------
    masterfile : str or None
        name of the master file. None if the volume was not complete
    masterdatatypedescr : str
        the description of the master data type
    last_processed : datetime
        True of all scans found

    """
    endtime_loop = deepcopy(nowtime)

    nscans = 1
    if datacfg['ScanList'] is not None:
        nscans = len(datacfg['ScanList'][0])

    scan_min = datacfg['ScanPeriod'] * 2.  # [min]

    starttime_loop_default = endtime_loop - timedelta(minutes=scan_min)
    if last_processed is None:
        # last processed volume not known. Process last scan
        starttime_loop = starttime_loop_default
    elif last_processed > endtime_loop:
        warn('last processed volume too new. Reprocessing the data')
        starttime_loop = starttime_loop_default
        last_processed = starttime_loop_default
    elif (starttime_loop_default - last_processed
          > timedelta(minutes=scan_min*3)):
        warn('last processed volume too old. ' +
             'There may be missing processed volumes')
        starttime_loop = starttime_loop_default
        last_processed = starttime_loop_default
    else:
        # process from last processed scan
        starttime_loop = last_processed + timedelta(seconds=10)

    masterfilelist, masterdatatypedescr, _ = _get_masterfile_list(
        datatype_list, np.array([starttime_loop]), np.array([endtime_loop]),
        datacfg, scan_list=datacfg['ScanList'])

    nvolumes = len(masterfilelist)
    if nvolumes == 0:
        return None, None, last_processed

    # check if there are rainbow data types and how many
    nrainbow = 0
    datatype_rainbow = []
    for datatype_descr in datatype_list:
        _, datagroup, datatype, _, _ = get_datatype_fields(datatype_descr)
        if datagroup == 'RAINBOW':
            datatype_rainbow.append(datatype)
            nrainbow += 1

    if nscans == 1:
        if nrainbow < 2:
            return masterfilelist[0], masterdatatypedescr, last_processed

        # If more than one data type is of type rainbow we have to wait for
        # all data type files to be present
        masterfile = masterfilelist[0]
        rainbow_files = []
        for datatype in datatype_rainbow:
            rainbow_file = get_new_rainbow_file_name(
                masterfile, masterdatatypedescr, datatype)
            rainbow_files.append(rainbow_file)

        # allow 30 s for the transfer of all datatype files
        found_all = _wait_for_rainbow_datatypes(rainbow_files, period=30)
        if found_all:
            return masterfile, masterdatatypedescr, last_processed

        # if not all data types available skip the volume
        warn('Not all data types for master file: ' +
             os.path.basename(masterfile)+' arrived on time. ' +
             'The volume will be skipped')
        return (
            None, None, get_datetime(masterfilelist[0], masterdatatypedescr))

    # if there is more than one scan in the list wait until all files
    # for the first volume have arrived
    masterfile = masterfilelist[0]
    master_voltime = get_datetime(masterfile, masterdatatypedescr)
    wait_time = nowtime+timedelta(minutes=scan_min)
    found_all = False
    currenttime = deepcopy(nowtime)
    while currenttime <= wait_time or not found_all:
        currenttime = datetime.utcnow()
        # for offline testing
        # currenttime = currenttime.replace(
        #    year=nowtime.year, month=nowtime.month, day=nowtime.day)
        # currenttime = currenttime.replace(nowtime.hour)

        starttime_loop = master_voltime
        endtime_loop = master_voltime+timedelta(minutes=scan_min)
        filelist_vol = []
        for scan in datacfg['ScanList'][0]:
            filelist = get_file_list(
                masterdatatypedescr, np.array([starttime_loop]),
                np.array([endtime_loop]), datacfg, scan=scan)
            if not filelist:
                filelist_vol = []
                found_all = False
                break
            filelist_vol.append(filelist[0])
            found_all = True
        if found_all:
            if nrainbow < 2:
                return masterfile, masterdatatypedescr, last_processed
            break

    if not found_all:
        # if not all scans available skip the volume
        warn('Not all scans for master file: ' +
             os.path.basename(masterfile)+' arrived on time. ' +
             'The volume will be skipped')

        return None, None, get_datetime(masterfile, masterdatatypedescr)

    # If more than one data type is of type rainbow we have to wait
    # for all data type files for all scans to be present
    rainbow_files = []
    for file in filelist_vol:
        for datatype in datatype_rainbow:
            rainbow_file = get_new_rainbow_file_name(
                file, masterdatatypedescr, datatype)
            rainbow_files.append(rainbow_file)

    # allow 30 s for the transfer of all datatype files
    found_all = _wait_for_rainbow_datatypes(rainbow_files, period=30)
    if found_all:
        return masterfile, masterdatatypedescr, last_processed

    # if not all scans available skip the volume
    warn('Not all data types for all scans of master file: ' +
         os.path.basename(masterfile)+' arrived on time. ' +
         'The volume will be skipped')

    return None, None, get_datetime(masterfile, masterdatatypedescr)


def _wait_for_rainbow_datatypes(rainbow_files, period=30):
    """
    waits until the files for all rainbow data types are present.

    Parameters
    ----------
    rainbow_files : list of strings
        a list containing the names of all the rainbow files to wait for
    period : int
        the time it has to wait (s)

    Returns
    -------
    found_all : Boolean
        True if all files were present. False otherwise

    """
    currenttime = datetime.utcnow()
    # for offline testing
    # currenttime = currenttime.replace(
    #     year=2017, month=6, day=14)
    # startime_proc = currenttime.replace(10)

    wait_time = currenttime+timedelta(seconds=period)
    while currenttime <= wait_time:
        currenttime = datetime.utcnow()
        # for offline testing
        # currenttime = currenttime.replace(
        #    year=2017, month=6, day=14)
        # startime_proc = currenttime.replace(10)

        found_all = False
        for rainbow_file in rainbow_files:
            filename = glob.glob(rainbow_file)
            if not filename:
                found_all = False
                break
            found_all = True
        if found_all:
            return found_all

    return found_all


@profiler(level=2)
def _get_radars_data(master_voltime, datatypesdescr_list, datacfg,
                     num_radars=1):
    """
    Get the radars data.

    Parameters
    ----------
    master_voltime : datetime object
        reference time
    datatypesdescr_list : list of lists
        List of the raw data types to get from each radar
    datacfg : dict
        dictionary containing the parameters to get the radar data

    Returns
    -------
    radar_list : list
        a list containing the radar objects

    """
    # get data of master radar
    radar_list = list()
    radar_list.append(
        get_data(master_voltime, datatypesdescr_list[0], datacfg))

    if num_radars == 1:
        return radar_list

    # get data of rest of radars
    for i in range(1, num_radars):
        filelist_ref, datatypedescr_ref, _ = _get_masterfile_list(
            datatypesdescr_list[i],
            [master_voltime-timedelta(seconds=datacfg['TimeTol'])],
            [master_voltime+timedelta(seconds=datacfg['TimeTol'])],
            datacfg, scan_list=datacfg['ScanList'])

        nfiles_ref = len(filelist_ref)
        if nfiles_ref == 0:
            warn("ERROR: Could not find any valid volume for " +
                 " reference time " +
                 master_voltime.strftime('%Y-%m-%d %H:%M:%S') +
                 ' and radar RADAR'+'{:03d}'.format(i+1))
            radar_list.append(None)
        elif nfiles_ref == 1:
            voltime_ref = get_datetime(
                filelist_ref[0], datatypedescr_ref)
            radar_list.append(
                get_data(voltime_ref, datatypesdescr_list[i], datacfg))
        else:
            voltime_ref_list = []
            for j in range(nfiles_ref):
                voltime_ref_list.append(get_datetime(
                    filelist_ref[j], datatypedescr_ref))
            voltime_ref = min(
                voltime_ref_list, key=lambda x: abs(x-master_voltime))
            radar_list.append(
                get_data(voltime_ref, datatypesdescr_list[i], datacfg))

    return radar_list


@profiler(level=2)
def _generate_dataset(dsname, cfg, dscfg, proc_status=0, radar_list=None,
                      voltime=None, trajectory=None, runinfo=None,
                      MULTIPROCESSING_PROD=False):
    """
    generates new datasets

    Parameters
    ----------
    dsname : str
        name of the dataset being generated
    cfg : dict
        configuration data
    dscfg : dict
        dataset configuration data
    proc_status : int
        processing status 0: init 1: processing 2: final
    radar_list : list
        a list containing the radar objects
    voltime : datetime
        reference time of the radar(s)
    trajectory : trajectory object
        trajectory object
    runinfo : str
        string containing run info
    MULTIPROCESSING_PROD : Bool
        If true the generation of products from each dataset will be
        parallelized

    Returns
    -------
    new_dataset : list
        The list of new datasets generated. None otherwise
    ind_rad : int
        the index to the reference radar object
    dsname : str
        name of the dataset being generated
    dscfg : dict
        the modified dataset configuration dictionary


    """
    dscfg = deepcopy(dscfg)

    dscfg['timeinfo'] = voltime
    try:
        proc_ds_func, dsformat = get_process_func(
            dscfg['type'], dscfg['dsname'])
    except Exception as inst:
        warn(str(inst))
        raise

    if isinstance(proc_ds_func, str):
        proc_ds_func = getattr(proc, proc_ds_func)

    # Create dataset
    if 'trajectory' in inspect.getfullargspec(proc_ds_func).args:
        new_dataset, ind_rad = proc_ds_func(
            proc_status, dscfg, radar_list=radar_list, trajectory=trajectory)
    else:
        new_dataset, ind_rad = proc_ds_func(
            proc_status, dscfg, radar_list=radar_list)

    if new_dataset is None:
        return None, None, dsname, dscfg


    if not isinstance(dsformat, list):
        dsformat = [dsformat]
        new_dataset = [new_dataset]

    # Handles the case of hybrid e.g. GRID/VOL proc
    for dset, dsformat in zip(new_dataset, dsformat):
        try:
            prod_func = get_prodgen_func(
            dsformat, dscfg['dsname'], dscfg['type'])
        except Exception as inst:
            warn(str(inst))
            raise

        # create the data set products
        if 'products' in dscfg:
            if MULTIPROCESSING_PROD:
                jobs = []
                for product in dscfg['products']:
                    # delay the data hashing
                    new_dataset = dask.delayed(dset)
                    jobs.append(dask.delayed(_generate_prod)(
                        new_dataset, cfg, product, prod_func, dscfg['dsname'],
                        voltime, runinfo=runinfo))

                dask.compute(*jobs)

            else:
                for product in dscfg['products']:
                    _generate_prod(
                        dset, cfg, product, prod_func, dscfg['dsname'],
                        voltime, runinfo=runinfo)

                    gc.collect()
    return new_dataset, ind_rad, dsname, dscfg


@profiler(level=3)
def _generate_prod(dataset, cfg, prdname, prdfunc, dsname, voltime,
                   runinfo=None):
    """
    generates a product

    Parameters
    ----------
    dataset : object
        the dataset object
    cfg : dict
        configuration data
    prdname : str
        name of the product
    prdfunc : func
        name of the product processing function
    dsname : str
        name of the dataset
    voltime : datetime object
        reference time of the radar(s)
    runinfo : str
        string containing run info

    Returns
    -------
    error : bool
        False if the products could be generated

    """
    print('---- Processing product: ' + prdname)
    prdcfg = _create_prdcfg_dict(
        cfg, dsname, prdname, voltime, runinfo=runinfo)
    try:
        prdfunc(dataset, prdcfg)
        return False
    except Exception as inst:
        warn(str(inst))
        traceback.print_exc()
        return True


@profiler(level=3)
def _create_cfg_dict(cfgfile):
    """
    creates a configuration dictionary

    Parameters
    ----------
    cfgfile : str
        path of the main config file

    Returns
    -------
    cfg : dict
        dictionary containing the configuration data

    """
    cfg = dict({'configFile': cfgfile})
    try:
        print("- Main config file : %s" % cfgfile)
        cfg = read_config(cfg['configFile'], cfg=cfg)
        print("- Location config file : %s" % cfg['locationConfigFile'])
        cfg = read_config(cfg['locationConfigFile'], cfg=cfg)
        print("- Product config file : %s" % cfg['productConfigFile'])
        cfg = read_config(cfg['productConfigFile'], cfg=cfg)
    except Exception as inst:
        warn(str(inst))
        sys.exit(1)

    # check for mandatory config parameters
    param_must = ['name', 'configpath', 'saveimgbasepath', 'dataSetList']
    for param in param_must:
        if param not in cfg:
            raise Exception("ERROR config: Parameter '%s' undefined!" % param)

    # fill in defaults
    if 'NumRadars' not in cfg:
        cfg.update({'NumRadars': 1})
    if 'TimeTol' not in cfg:
        cfg.update({'TimeTol': 3600.})
    if 'ScanList' not in cfg:
        cfg.update({'ScanList': None})
    else:
        cfg.update({'ScanList': get_scan_list(cfg['ScanList'])})
    if 'lastStateFile' not in cfg:
        cfg.update({'lastStateFile': None})
    if 'datapath' not in cfg:
        cfg.update({'datapath': None})
    if 'satpath' not in cfg:
        cfg.update({'satpath': None})
    if 'path_convention' not in cfg:
        cfg.update({'path_convention': 'MCH'})
    if 'cosmopath' not in cfg:
        cfg.update({'cosmopath': None})
    if 'psrpath' not in cfg:
        cfg.update({'psrpath': None})
    if 'iqpath' not in cfg:
        cfg.update({'iqpath': None})
    if 'colocgatespath' not in cfg:
        cfg.update({'colocgatespath': None})
    if 'excessgatespath' not in cfg:
        cfg.update({'excessgatespath': None})
    if 'dempath' not in cfg:
        cfg.update({'dempath': None})
    if 'smnpath' not in cfg:
        cfg.update({'smnpath': None})
    if 'disdropath' not in cfg:
        cfg.update({'disdropath': None})
    if 'solarfluxpath' not in cfg:
        cfg.update({'solarfluxpath': None})
    if 'loadbasepath' not in cfg:
        cfg.update({'loadbasepath': None})
    if 'loadname' not in cfg:
        cfg.update({'loadname': None})
    if 'RadarName' not in cfg:
        cfg.update({'RadarName': None})
    if 'RadarRes' not in cfg:
        cfg.update({'RadarRes': None})
    if 'metranet_read_lib' not in cfg:
        cfg.update({'metranet_read_lib': 'C'})
    if 'ScanPeriod' not in cfg:
        warn('WARNING: Scan period not specified. ' +
             'Assumed default value 5 min')
        cfg.update({'ScanPeriod': 5})
    if 'CosmoRunFreq' not in cfg:
        warn('WARNING: COSMO run frequency not specified. ' +
             'Assumed default value 3h')
        cfg.update({'CosmoRunFreq': 3})
    if 'CosmoForecasted' not in cfg:
        warn('WARNING: Hours forecasted by COSMO not specified. ' +
             'Assumed default value 7h (including analysis)')
        cfg.update({'CosmoForecasted': 7})

    # Instrument parameters not in radar object attributes
    if 'lradomeh' not in cfg:
        cfg.update({
            'lradomeh': 0.*np.ones(cfg['NumRadars'], dtype=np.float32)})
    if 'lradomev' not in cfg:
        cfg.update({
            'lradomev': 0.*np.ones(cfg['NumRadars'], dtype=np.float32)})
    if 'lrxh' not in cfg:
        cfg.update({
            'lrxh': 0.*np.ones(cfg['NumRadars'], dtype=np.float32)})
    if 'lrxv' not in cfg:
        cfg.update({
            'lrxv': 0.*np.ones(cfg['NumRadars'], dtype=np.float32)})
    if 'ltxh' not in cfg:
        cfg.update({
            'ltxh': 0.*np.ones(cfg['NumRadars'], dtype=np.float32)})
    if 'ltxv' not in cfg:
        cfg.update({
            'ltxv': 0.*np.ones(cfg['NumRadars'], dtype=np.float32)})

    # Convert the following strings to string arrays
    strarr_list = [
        'datapath', 'cosmopath', 'dempath', 'loadbasepath', 'psrpath',
        'iqpath', 'satpath', 'loadname', 'RadarName', 'RadarRes', 'ScanList',
        'imgformat', 'frequency', 'radar_beam_width_h', 'radar_beam_width_v',
        'pulse_width', 'nyquist_velocity', 'AntennaGainH', 'AntennaGainV',
        'dBADUtodBmh', 'dBADUtodBmv', 'mflossh', 'mflossv', 'radconsth',
        'radconstv', 'txpwrh', 'txpwrv', 'attg']
    for param in strarr_list:
        if param in cfg and isinstance(cfg[param], str):
            cfg[param] = [cfg[param]]

    # Convert the following floats to float arrays
    fltarr_list = [
        'frequency', 'radar_beam_width_h', 'radar_beam_width_v',
        'pulse_width', 'nyquist_velocity', 'AntennaGainH', 'AntennaGainV',
        'dBADUtodBmh', 'dBADUtodBmv', 'mflossh', 'mflossv', 'radconsth',
        'radconstv', 'txpwrh', 'txpwrv', 'attg', 'lradomeh', 'lradomev',
        'lrxh', 'lrxv', 'ltxh', 'ltxv','mosotti_factor','refcorr','AzimTol']
    for param in fltarr_list:
        if param in cfg and isinstance(cfg[param], float):
            cfg[param] = [cfg[param]]

    # if specified in config, convert coordinates to arrays
    if 'RadarPosition' in cfg:
        fltarr_list = ['latitude', 'longitude', 'altitude']
        for param in fltarr_list:
            if isinstance(cfg['RadarPosition'][param], float):
                cfg['RadarPosition'][param] = [cfg['RadarPosition'][param]]

    return cfg


@profiler(level=3)
def _create_datacfg_dict(cfg):
    """
    creates a data configuration dictionary from a config dictionary

    Parameters
    ----------
    cfg : dict
        config dictionary

    Returns
    -------
    datacfg : dict
        data config dictionary

    """

    datacfg = dict({'datapath': cfg['datapath']})
    datacfg.update({'satpath': cfg['satpath']})
    datacfg.update({'psrpath': cfg['psrpath']})
    datacfg.update({'iqpath': cfg['iqpath']})
    datacfg.update({'ScanList': cfg['ScanList']})
    datacfg.update({'TimeTol': cfg['TimeTol']})
    datacfg.update({'NumRadars': cfg['NumRadars']})
    datacfg.update({'cosmopath': cfg['cosmopath']})
    datacfg.update({'dempath': cfg['dempath']})
    datacfg.update({'loadbasepath': cfg['loadbasepath']})
    datacfg.update({'loadname': cfg['loadname']})
    datacfg.update({'RadarName': cfg['RadarName']})
    datacfg.update({'RadarRes': cfg['RadarRes']})
    datacfg.update({'ScanPeriod': cfg['ScanPeriod']})
    datacfg.update({'CosmoRunFreq': int(cfg['CosmoRunFreq'])})
    datacfg.update({'CosmoForecasted': int(cfg['CosmoForecasted'])})
    datacfg.update({'path_convention': cfg['path_convention']})
    datacfg.update({'metranet_read_lib': cfg['metranet_read_lib']})

    # Modify size of radar or radar spectra object
    datacfg.update({'elmin': cfg.get('elmin', None)})
    datacfg.update({'elmax': cfg.get('elmax', None)})
    datacfg.update({'azmin': cfg.get('azmin', None)})
    datacfg.update({'azmax': cfg.get('azmax', None)})
    datacfg.update({'rmin': cfg.get('rmin', None)})
    datacfg.update({'rmax': cfg.get('rmax', None)})

    # Modify size of grid object
    datacfg.update({'latmin': cfg.get('latmin', None)})
    datacfg.update({'latmax': cfg.get('latmax', None)})
    datacfg.update({'lonmin': cfg.get('lonmin', None)})
    datacfg.update({'lonmax': cfg.get('lonmax', None)})
    datacfg.update({'altmin': cfg.get('altmin', None)})
    datacfg.update({'altmax': cfg.get('altmax', None)})
    datacfg.update({'nx': cfg.get('nx', None)})
    datacfg.update({'ny': cfg.get('ny', None)})
    datacfg.update({'nz': cfg.get('nz', None)})

    # variables to get the spectral data
    datacfg.update({'undo_txcorr': cfg.get('undo_txcorr', True)})
    datacfg.update({'fold': cfg.get('fold', True)})
    datacfg.update({'positive_away': cfg.get('positive_away', True)})
    datacfg.update({'cpi': cfg.get('cpi', 'low_prf')})
    datacfg.update({'ang_tol': cfg.get('ang_tol', 0.5)})

    # Radar position
    if 'RadarPosition' in cfg:
        datacfg.update({'RadarPosition': cfg['RadarPosition']})

    # Instrument parameters
    if 'frequency' in cfg:
        datacfg.update({'frequency': cfg['frequency']})
    if 'radar_beam_width_h' in cfg:
        datacfg.update({'radar_beam_width_h': cfg['radar_beam_width_h']})
    if 'radar_beam_width_v' in cfg:
        datacfg.update({'radar_beam_width_v': cfg['radar_beam_width_v']})
    if 'pulse_width' in cfg:
        datacfg.update({'pulse_width': cfg['pulse_width']})
    if 'nyquist_velocity' in cfg:
        datacfg.update({'nyquist_velocity': cfg['nyquist_velocity']})
    if 'AntennaGainH' in cfg:
        datacfg.update({'AntennaGainH': cfg['AntennaGainH']})
    if 'AntennaGainV' in cfg:
        datacfg.update({'AntennaGainV': cfg['AntennaGainV']})
    if 'AntennaGainV' in cfg:
        datacfg.update({'AntennaGainV': cfg['AntennaGainV']})
    if 'AntennaGainV' in cfg:
        datacfg.update({'AntennaGainV': cfg['AntennaGainV']})
    if 'altitude' in cfg:
        datacfg.update({'altitude': cfg['altitude']})
    if 'latitude' in cfg:
        datacfg.update({'latitude': cfg['latitude']})
    if 'longitude' in cfg:
        datacfg.update({'longitude': cfg['longitude']})

    # Radar calibration parameters
    if 'dBADUtodBmh' in cfg:
        datacfg.update({'dBADUtodBmh': cfg['dBADUtodBmh']})
    if 'dBADUtodBmv' in cfg:
        datacfg.update({'dBADUtodBmv': cfg['dBADUtodBmv']})
    if 'mflossh' in cfg:
        datacfg.update({'mflossh': cfg['mflossh']})
    if 'mflossv' in cfg:
        datacfg.update({'mflossv': cfg['mflossv']})
    if 'radconsth' in cfg:
        datacfg.update({'radconsth': cfg['radconsth']})
    if 'radconstv' in cfg:
        datacfg.update({'radconstv': cfg['radconstv']})
    if 'txpwrh' in cfg:
        datacfg.update({'txpwrh': cfg['txpwrh']})
    if 'txpwrv' in cfg:
        datacfg.update({'txpwrv': cfg['txpwrv']})
    if 'attg' in cfg:
        datacfg.update({'attg': cfg['attg']})
    if 'attg' in cfg:
        datacfg.update({'attg': cfg['attg']})
    if 'mosotti_factor' in cfg:
        datacfg.update({'mosotti_factor': cfg['mosotti_factor']})
    if 'refcorr' in cfg:
        datacfg.update({'refcorr': cfg['refcorr']})
    if 'AzimTol' in cfg:
        datacfg.update({'AzimTol': cfg['AzimTol']})
    return datacfg


@profiler(level=3)
def _create_dscfg_dict(cfg, dataset):
    """
    creates a dataset configuration dictionary

    Parameters
    ----------
    cfg : dict
        config dictionary
    dataset : str
        name of the dataset

    Returns
    -------
    dscfg : dict
        dataset config dictionary

    """
    dscfg = cfg[dataset]

    # Path related parameters
    dscfg.update({'configpath': cfg['configpath']})
    dscfg.update({'basepath': cfg['saveimgbasepath']})
    dscfg.update({'path_convention': cfg['path_convention']})
    dscfg.update({'procname': cfg['name']})
    dscfg.update({'dsname': dataset})
    dscfg.update({'solarfluxpath': cfg['solarfluxpath']})
    dscfg.update({'colocgatespath': cfg['colocgatespath']})
    dscfg.update({'excessgatespath': cfg['excessgatespath']})
    dscfg.update({'dempath': cfg['dempath']})
    dscfg.update({'cosmopath': cfg['cosmopath']})
    dscfg.update({'CosmoRunFreq': cfg['CosmoRunFreq']})
    dscfg.update({'CosmoForecasted': cfg['CosmoForecasted']})
    dscfg.update({'metranet_read_lib': cfg['metranet_read_lib']})
    dscfg.update({'lastStateFile': cfg['lastStateFile']})
    dscfg.update({'timeinfo': None})

    # Instrument parameters
    dscfg.update({'RadarName': cfg['RadarName']})
    dscfg.update({'ScanPeriod': cfg['ScanPeriod']})
    dscfg.update({'lrxh': cfg['lrxh']})
    dscfg.update({'lrxv': cfg['lrxv']})
    dscfg.update({'ltxh': cfg['ltxh']})
    dscfg.update({'ltxv': cfg['ltxv']})
    dscfg.update({'lradomeh': cfg['lradomeh']})
    dscfg.update({'lradomev': cfg['lradomev']})


    # PAR and ASR variable
    if 'par_azimuth_antenna' in cfg:
        dscfg.update({'par_azimuth_antenna': cfg['par_azimuth_antenna']})
    if 'par_elevation_antenna' in cfg:
        dscfg.update({'par_elevation_antenna': cfg['par_elevation_antenna']})
    if 'asr_highbeam_antenna' in cfg:
        dscfg.update({'asr_highbeam_antenna': cfg['asr_highbeam_antenna']})
    if 'asr_lowbeam_antenna' in cfg:
        dscfg.update({'asr_lowbeam_antenna': cfg['asr_lowbeam_antenna']})
    if 'target_radar_pos' in cfg:
        dscfg.update({'target_radar_pos': cfg['target_radar_pos']})

    # indicates the dataset has been initialized and aux data is available
    dscfg.update({'initialized': False})
    dscfg.update({'global_data': None})

    # Convert the following strings to string arrays
    strarr_list = ['datatype', 'FIELDS_TO_REMOVE']
    for param in strarr_list:
        if param in dscfg:
            if isinstance(dscfg[param], str):
                dscfg[param] = [dscfg[param]]

    # variables to make the data set available in the next level
    if 'MAKE_GLOBAL' not in dscfg:
        dscfg.update({'MAKE_GLOBAL': 0})
    if 'SUBSTITUTE_OBJECT' not in dscfg:
        dscfg.update({'SUBSTITUTE_OBJECT': 0})
    if 'FIELDS_TO_REMOVE' not in dscfg:
        dscfg.update({'FIELDS_TO_REMOVE': None})

    return dscfg


@profiler(level=3)
def _create_prdcfg_dict(cfg, dataset, product, voltime, runinfo=None):
    """
    creates a product configuration dictionary

    Parameters
    ----------
    cfg : dict
        config dictionary
    dataset : str
        name of the dataset used to create the product
    product : str
        name of the product
    voltime : datetime object
        time of the dataset

    Returns
    -------
    prdcfg : dict
        product config dictionary

    """

    # Ugly copying of dataset config parameters to product
    # config dict. Better: Make dataset config dict available to
    # the product generation.
    prdcfg = cfg[dataset]['products'][product]
    prdcfg.update({'procname': cfg['name']})
    prdcfg.update({'lastStateFile': cfg['lastStateFile']})
    prdcfg.update({'basepath': cfg['saveimgbasepath']})
    prdcfg.update({'smnpath': cfg['smnpath']})
    prdcfg.update({'disdropath': cfg['disdropath']})
    prdcfg.update({'cosmopath': cfg['cosmopath']})
    prdcfg.update({'ScanPeriod': cfg['ScanPeriod']})
    prdcfg.update({'imgformat': cfg['imgformat']})
    prdcfg.update({'RadarName': cfg['RadarName']})
    if 'ppiImageConfig' in cfg:
        prdcfg.update({'ppiImageConfig': cfg['ppiImageConfig']})
    if 'ppiMapImageConfig' in cfg:
        prdcfg.update({'ppiMapImageConfig': cfg['ppiMapImageConfig']})
    if 'rhiImageConfig' in cfg:
        prdcfg.update({'rhiImageConfig': cfg['rhiImageConfig']})
    if 'gridMapImageConfig' in cfg:
        prdcfg.update({'gridMapImageConfig': cfg['gridMapImageConfig']})
    if 'sunhitsImageConfig' in cfg:
        prdcfg.update({'sunhitsImageConfig': cfg['sunhitsImageConfig']})
    prdcfg.update({'dsname': dataset})
    prdcfg.update({'dstype': cfg[dataset]['type']})
    prdcfg.update({'prdname': product})
    prdcfg.update({'timeinfo': voltime})
    prdcfg.update({'runinfo': runinfo})
    if 'dssavename' in cfg[dataset]:
        prdcfg.update({'dssavename': cfg[dataset]['dssavename']})

    return prdcfg


@profiler(level=3)
def _get_datatype_list(cfg, radarnr='RADAR001'):
    """
    get list of unique input data types

    Parameters
    ----------
    cfg : dict
        config dictionary
    radarnr : str
        radar number identifier

    Returns
    -------
    datatypesdescr : list
        list of data type descriptors

    """
    datatypesdescr = set()

    for datasetdescr in cfg['dataSetList']:
        _, dataset = get_dataset_fields(datasetdescr)
        if 'datatype' not in cfg[dataset]:
            continue
        if isinstance(cfg[dataset]['datatype'], str):
            (radarnr_descr, datagroup, datatype_aux, dataset_save,
             product_save) = (
                 get_datatype_fields(cfg[dataset]['datatype']))
            if datagroup != 'PROC' and radarnr_descr == radarnr:
                if ((dataset_save is None) and (product_save is None)):
                    datatypesdescr.add(
                        radarnr_descr+":"+datagroup+":"+datatype_aux)
                elif ((dataset_save is not None) and (product_save is None)):
                    datatypesdescr.add(
                        radarnr_descr+":"+datagroup+":"+datatype_aux+"," +
                        dataset_save)
                else:
                    datatypesdescr.add(
                        radarnr_descr+":"+datagroup+":"+datatype_aux+"," +
                        dataset_save+","+product_save)
        else:
            for datatype in cfg[dataset]['datatype']:
                (radarnr_descr, datagroup, datatype_aux, dataset_save,
                 product_save) = (
                     get_datatype_fields(datatype))
                if datagroup != 'PROC' and radarnr_descr == radarnr:
                    if dataset_save is None and product_save is None:
                        datatypesdescr.add(
                            radarnr_descr+":"+datagroup+":"+datatype_aux)
                    elif dataset_save is not None and product_save is None:
                        datatypesdescr.add(
                            radarnr_descr+":"+datagroup+":"+datatype_aux+"," +
                            dataset_save)
                    else:
                        datatypesdescr.add(
                            radarnr_descr+":"+datagroup+":"+datatype_aux +
                            ","+dataset_save+","+product_save)

    datatypesdescr = list(datatypesdescr)

    return datatypesdescr


@profiler(level=3)
def _get_datasets_list(cfg):
    """
    get list of dataset at each processing level

    Parameters
    ----------
    cfg : dict
        config dictionary

    Returns
    -------
    dataset_levels : dict
        a dictionary containing the list of datasets at each processing level

    """
    dataset_levels = dict({'l00': list()})
    for datasetdescr in cfg['dataSetList']:
        proclevel, dataset = get_dataset_fields(datasetdescr)
        if proclevel in dataset_levels:
            dataset_levels[proclevel].append(dataset)
        else:
            dataset_levels.update({proclevel: [dataset]})

    return dataset_levels


@profiler(level=3)
def _get_masterfile_list(datatypesdescr, starttimes, endtimes, datacfg,
                         scan_list=None):
    """
    get master file list

    Parameters
    ----------
    datatypesdescr : list
        list of unique data type descriptors
    starttimes, endtimes : array of datetime objects
        start and end of processing periods
    datacfg : dict
        data configuration dictionary
    scan_list : list
        list of scans

    Returns
    -------
    masterfilelist : list
        the list of master files
    masterdatatypedescr : str
        the master data type descriptor

    """

    masterdatatypedescr = None
    masterscan = None
    for datatypedescr in datatypesdescr:
        radarnr, datagroup, _, _, _ = get_datatype_fields(datatypedescr)
        if (datagroup not in (
                'COSMO', 'RAD4ALPCOSMO', 'CFRADIALCOSMO', 'DEM', 'RAD4ALPDEM',
                'RAD4ALPHYDRO', 'RAD4ALPDOPPLER', 'RAD4ALPIQ', 'PSR',
                'PSRSPECTRA')):
            masterdatatypedescr = datatypedescr
            if scan_list is not None:
                masterscan = scan_list[int(radarnr[5:8])-1][0]
            break

    # if data type is not radar use dBZ as reference
    if masterdatatypedescr is None:
        for datatypedescr in datatypesdescr:
            radarnr, datagroup, _, _, _ = get_datatype_fields(datatypedescr)
            if datagroup in ('COSMO', 'DEM', 'PSR', 'PSRSPECTRA'):
                masterdatatypedescr = radarnr+':RAINBOW:dBZ'
                if scan_list is not None:
                    masterscan = scan_list[int(radarnr[5:8])-1][0]
                break
            if (datagroup in (
                    'RAD4ALPCOSMO', 'RAD4ALPDEM', 'RAD4ALPHYDRO',
                    'RAD4ALPDOPPLER', 'RAD4ALPIQ')):
                masterdatatypedescr = radarnr+':RAD4ALP:dBZ'
                if scan_list is not None:
                    masterscan = scan_list[int(radarnr[5:8])-1][0]
                break

    if masterdatatypedescr is None:
        warn('No data type can be used as master.'
             ' Add a dBZ data type to the config file')
        return [], None, None

    masterfilelist = get_file_list(
        masterdatatypedescr, starttimes, endtimes, datacfg,
        scan=masterscan)

    return masterfilelist, masterdatatypedescr, masterscan


@profiler(level=3)
def _add_dataset(new_dataset, radar_list, ind_rad, make_global=True,
                 substitute_object=False, fields_to_remove=None):
    """
    adds a new field to an existing radar object

    Parameters
    ----------
    new_dataset : dict
        dictionary with key radar_out containing the new fields
    radar : radar object
        the radar object containing the global data
    make_global : boolean
        if true a new field is added to the global data
    substitute_object : boolean
        if true the new object will substitute the previous one
    fields_to_remove : list or None
        List of fields to be removed from the object

    Returns
    -------
    0 if successful. None otherwise

    """
    if radar_list is None:
        return None

    if not make_global:
        return None

    if new_dataset is None:
        return None

    if 'radar_out' not in new_dataset:
        print('No radar_out field in new_dataset')
        return None

    if substitute_object:
        print('Substituting object')
        radar_list[ind_rad] = new_dataset['radar_out']
        return 0

    for field in new_dataset['radar_out'].fields:
        print('Adding field: '+field)
        radar_list[ind_rad].add_field(
            field, new_dataset['radar_out'].fields[field],
            replace_existing=True)

    if fields_to_remove is None:
        return 0

    for field in fields_to_remove:
        field_pyart = get_fieldname_pyart(field)
        print('Removing field: '+field_pyart)
        if field_pyart not in radar_list[ind_rad].fields:
            print('Field '+field_pyart+' not in radar object')
            continue
        del radar_list[ind_rad].fields[field_pyart]

    return 0


def _warning_format(message, category, filename, lineno, file=None,
                    line=None):
    return '%s (%s:%s)\n' % (message, filename, lineno)
