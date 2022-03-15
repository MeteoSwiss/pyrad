"""
pyrad.proc.process_traj
=============================

Trajectory functions. Functions to pass trajectory dataset data to
the product generation functions.

.. autosummary::
    :toctree: generated/

    process_trajectory
    process_traj_trt
    process_traj_trt_contour
    process_traj_lightning
    process_traj_atplane
    process_traj_antenna_pattern
    _get_ts_values_antenna_pattern
    _get_contour_trt
    _get_gates
    _get_gates_trt
    _get_gates_antenna_pattern
    _get_closest_bin
    _sample_out_of_sector
    TargetRadar
"""

from warnings import warn
import gc
from copy import deepcopy

import numpy as np
from netCDF4 import num2date, date2num

import pyart
from pyart.config import get_metadata
from pyart.core import Radar
from pyart.util import get_target_elevations

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..io.io_aux import get_field_unit, get_field_name
from ..io.timeseries import TimeSeries
from ..io.read_data_other import read_antenna_pattern

from ..util.stat_utils import quantiles_weighted
from ..util.radar_utils import belongs_roi_indices, find_nearest_gate


def process_trajectory(procstatus, dscfg, radar_list=None, trajectory=None):
    """
    Return trajectory

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
    radar_list : list of Radar objects
        Optional. list of radar objects
    trajectory : Trajectory object
        containing trajectory samples

    Returns
    -------
    new_dataset : Trajectory object
        radar object
    ind_rad : int
        None

    """

    if procstatus != 1:
        return None, None

    if not dscfg['initialized']:
        if trajectory is None:
            raise Exception("ERROR: Undefined trajectory for dataset '%s'"
                            % dscfg['dsname'])

        if radar_list is not None:
            for radar in radar_list:
                rad = trajectory.add_radar(radar)
                trajectory.calculate_velocities(rad)
        else:
            warn('ERROR: No valid radar found')
            return None, None

        dscfg['initialized'] = True
        return trajectory, None

    return None, None


def process_traj_trt(procstatus, dscfg, radar_list=None, trajectory=None):
    """
    Processes data according to TRT trajectory

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        time_tol : float. Dataset keyword
            tolerance between reference time of the radar volume and that of
            the TRT cell [s]. Default 100.
        alt_min, alt_max : float. Dataset keyword
            Minimum and maximum altitude of the data inside the TRT cell to
            retrieve [m MSL]. Default None
        cell_center : Bool. Dataset keyword
            If True only the range gate closest to the center of the cell is
            extracted. Default False
        latlon_tol : Float. Dataset keyword
            Tolerance in lat/lon when extracting data only from the center of
            the TRT cell. Default 0.01
    radar_list : list of Radar objects
        Optional. list of radar objects
    trajectory : Trajectory object
        containing trajectory samples

    Returns
    -------
    new_dataset : dictionary
        Dictionary containing radar_out, a radar object containing only data
        from inside the TRT cell
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # Process
    field_names = []
    datatypes = []
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_names.append(get_fieldname_pyart(datatype))
        datatypes.append(datatype)

    ind_rad = int(radarnr[5:8])-1
    if ((radar_list is None) or (radar_list[ind_rad] is None)):
        warn('ERROR: No valid radar found')
        return None, None

    # keep locally only field of interest in radar object
    radar = deepcopy(radar_list[ind_rad])
    radar.fields = dict()
    nfields_available = 0
    for field_name in field_names:
        if field_name not in radar_list[ind_rad].fields:
            warn("Datatype '%s' not available in radar data" % field_name)
            continue
        radar.add_field(field_name, radar_list[ind_rad].fields[field_name])
        nfields_available += 1

    if nfields_available == 0:
        warn("Fields not available in radar data")
        return None, None

    # get TRT cell corresponding to current radar volume
    time_tol = dscfg.get('TimeTol', 100.)
    alt_min = dscfg.get('alt_min', None)
    alt_max = dscfg.get('alt_max', None)
    cell_center = dscfg.get('cell_center', False)
    latlon_tol = dscfg.get('latlon_tol', 0.01)  # aprox. 1 km

    inds_ray, inds_rng, lat, lon, alt = _get_gates_trt(
        radar, trajectory, dscfg['timeinfo'], time_tol=time_tol,
        alt_min=alt_min, alt_max=alt_max, cell_center=cell_center,
        latlon_tol=latlon_tol)

    if inds_ray is None:
        return None, None

    # prepare new radar object output
    radar_roi = deepcopy(radar)

    radar_roi.range['data'] = radar.range['data'][inds_rng]
    radar_roi.ngates = inds_rng.size
    radar_roi.time['data'] = np.asarray([radar_roi.time['data'][0]])
    radar_roi.scan_type = 'roi'
    radar_roi.sweep_mode['data'] = np.array(['roi'])
    radar_roi.sweep_start_ray_index['data'] = np.array([0], dtype='int32')
    radar_roi.fixed_angle['data'] = np.array([], dtype='float64')
    radar_roi.sweep_number['data'] = np.array([0], dtype='int32')
    radar_roi.nsweeps = 1

    if radar.rays_are_indexed is not None:
        radar_roi.rays_are_indexed['data'] = np.array(
            [radar.rays_are_indexed['data'][0]])
    if radar.ray_angle_res is not None:
        radar_roi.ray_angle_res['data'] = np.array(
            [radar.ray_angle_res['data'][0]])

    radar_roi.sweep_end_ray_index['data'] = np.array([1], dtype='int32')
    radar_roi.rays_per_sweep = np.array([1], dtype='int32')
    radar_roi.azimuth['data'] = np.array([], dtype='float64')
    radar_roi.elevation['data'] = np.array([], dtype='float64')
    radar_roi.nrays = 1

    radar_roi.gate_longitude['data'] = np.empty(
        (radar_roi.nrays, radar_roi.ngates), dtype=float)
    radar_roi.gate_latitude['data'] = np.empty(
        (radar_roi.nrays, radar_roi.ngates), dtype=float)
    radar_roi.gate_altitude['data'] = np.empty(
        (radar_roi.nrays, radar_roi.ngates), dtype=float)

    radar_roi.gate_x['data'] = np.empty(
        (radar_roi.nrays, radar_roi.ngates), dtype=float)
    radar_roi.gate_y['data'] = np.empty(
        (radar_roi.nrays, radar_roi.ngates), dtype=float)
    radar_roi.gate_z['data'] = np.empty(
        (radar_roi.nrays, radar_roi.ngates), dtype=float)

    radar_roi.gate_longitude['data'][0, :] = lon
    radar_roi.gate_latitude['data'][0, :] = lat
    radar_roi.gate_altitude['data'][0, :] = alt

    radar_roi.gate_x['data'][0, :] = radar.gate_x['data'][inds_ray, inds_rng]
    radar_roi.gate_y['data'][0, :] = radar.gate_y['data'][inds_ray, inds_rng]
    radar_roi.gate_z['data'][0, :] = radar.gate_z['data'][inds_ray, inds_rng]

    radar_roi.fields = dict()
    for field_name in field_names:
        if field_name not in radar.fields:
            warn("Datatype '%s' not available in radar data" % field_name)
            continue

        field_dict = deepcopy(radar.fields[field_name])
        field_dict['data'] = np.ma.empty(
            (radar_roi.nrays, radar_roi.ngates), dtype=float)
        field_dict['data'][0, :] = radar.fields[field_name]['data'][
            inds_ray, inds_rng]
        radar_roi.add_field(field_name, field_dict)

    new_dataset = {'radar_out': radar_roi}

    return new_dataset, ind_rad


def process_traj_trt_contour(procstatus, dscfg, radar_list=None,
                             trajectory=None):
    """
    Gets the TRT cell contour corresponding to each radar volume

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        time_tol : float. Dataset keyword
            tolerance between reference time of the radar volume and that of
            the TRT cell [s]. Default 100.
    radar_list : list of Radar objects
        Optional. list of radar objects
    trajectory : Trajectory object
        containing trajectory samples

    Returns
    -------
    new_dataset : dict
        Dictionary containing radar_out and roi_dict. Radar out is the current
        radar object. roi_dict contains the positions defining the TRT cell
        contour
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # Process
    field_names = []
    datatypes = []
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_names.append(get_fieldname_pyart(datatype))
        datatypes.append(datatype)

    ind_rad = int(radarnr[5:8])-1
    if ((radar_list is None) or (radar_list[ind_rad] is None)):
        warn('ERROR: No valid radar found')
        return None, None

    # keep locally only field of interest in radar object
    radar = deepcopy(radar_list[ind_rad])
    radar.fields = dict()
    nfields_available = 0
    for field_name in field_names:
        if field_name not in radar_list[ind_rad].fields:
            warn("Datatype '%s' not available in radar data" % field_name)
            continue
        radar.add_field(field_name, radar_list[ind_rad].fields[field_name])
        nfields_available += 1

    if nfields_available == 0:
        warn("Fields not available in radar data")
        return None, None

    # get TRT cell corresponding to current radar volume
    time_tol = dscfg.get('TimeTol', 100.)

    roi_dict = _get_contour_trt(
        radar, trajectory, dscfg['timeinfo'], time_tol=time_tol)

    if roi_dict is None:
        return None, None

    new_dataset = {
        'radar_out': radar,
        'roi_dict': roi_dict}

    return new_dataset, ind_rad


def process_traj_lightning(procstatus, dscfg, radar_list=None,
                           trajectory=None):
    """
    Return time series according to lightning trajectory

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        data_is_log : dict. Dataset keyword
            Dictionary specifying for each field if it is in log (True) or
            linear units (False). Default False
        ang_tol : float. Dataset keyword
            Factor that multiplies the angle resolution. Used when determining
            the neighbouring rays. Default 1.2

    radar_list : list of Radar objects
        Optional. list of radar objects
    trajectory : Trajectory object
        containing trajectory samples

    Returns
    -------
    trajectory : Trajectory object
        Object holding time series
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        # first call: nothing to do
        return None, None

    if procstatus == 2:
        # last call: do the products
        if not dscfg['initialized']:
            warn('ERROR: No trajectory dataset available!')
            return None, None
        trajdict = dscfg['traj_atplane_dict']

        dataset = {
            'ts': trajdict['ts'],
            'final': 1
        }
        return dataset, trajdict['ind_rad']

    # Process
    field_names = []
    datatypes = []
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_names.append(get_fieldname_pyart(datatype))
        datatypes.append(datatype)

    ind_rad = int(radarnr[5:8])-1
    if ((radar_list is None) or (radar_list[ind_rad] is None)):
        warn('ERROR: No valid radar found')
        return None, None

    # keep locally only field of interest in radar object
    radar = deepcopy(radar_list[ind_rad])
    radar.fields = dict()
    nfields_available = 0
    for field_name in field_names:
        if field_name not in radar_list[ind_rad].fields:
            warn("Datatype '%s' not available in radar data" % field_name)
            continue
        radar.add_field(field_name, radar_list[ind_rad].fields[field_name])
        nfields_available += 1

    if nfields_available == 0:
        warn("Fields not available in radar data")
        return None, None

    ttask_start = radar.time['data'].min()
    dt_task_start = num2date(ttask_start, radar.time['units'],
                             radar.time['calendar'])

    if not dscfg['initialized']:
        # init
        if trajectory is None:
            raise Exception("ERROR: Undefined trajectory for dataset '%s'"
                            % dscfg['dsname'])

        rad_traj = trajectory.add_radar(radar)

        description = [
            "Description:",
            "Time series of a weather radar data type at the location",
            "of the lightning flash.",
            "The time samples where the flash was out of the weather radar",
            "sector are NOT included in this file.",
            "NaN (Not a number): No data detected at the flash location."
        ]

        ts_dict = dict()
        data_is_log = dict()
        for datatype, field_name in zip(datatypes, field_names):
            ts = TimeSeries(
                description, maxlength=trajectory.time_vector.size,
                datatype=datatype)

            unit = get_field_unit(datatype)
            name = get_field_name(datatype)
            # Append empty series: Note: sequence matters!
            ts.add_dataseries("#Flash", "", "", plot=False)
            ts.add_dataseries("Power", "", "dBm", plot=False)
            ts.add_dataseries("at_flash", name, unit, color='b')
            ts.add_dataseries("Mean", name, unit, color='r')
            ts.add_dataseries("Min", name, unit, color='k', linestyle=':')
            ts.add_dataseries("Max", name, unit, color='k', linestyle=':')
            ts.add_dataseries("#Valid", "", "", plot=False)

            ts_dict.update({field_name: ts})
            data_is_log.update({field_name: False})
            if 'data_is_log' in dscfg:
                if datatype in dscfg['data_is_log']:
                    data_is_log[field_name] = (dscfg['data_is_log'] != 0)
                else:
                    warn('Units type for data type '+datatype +
                         ' not specified. Assumed linear')

        trajdict = dict({
            'radar_traj': rad_traj,
            'radar_old': None,
            'radar_old2': None,
            'last_task_start_dt': None,
            'ind_rad': ind_rad,
            'ts_dict': ts_dict,
            'data_is_log': data_is_log})
        traj_ind = trajectory.get_samples_in_period(end=dt_task_start)

        dscfg['traj_atplane_dict'] = trajdict
        dscfg['initialized'] = True
    else:
        # init already done
        trajdict = dscfg['traj_atplane_dict']
        rad_traj = trajdict['radar_traj']
        traj_ind = trajectory.get_samples_in_period(
            start=trajdict['last_task_start_dt'], end=dt_task_start)
        ts_dict = trajdict['ts_dict']
        data_is_log = trajdict['data_is_log']

    traj_time_vec = date2num(trajectory.time_vector, radar.time['units'],
                             radar.time['calendar'])

    if np.size(traj_ind) == 0:
        warn('No trajectory samples within current period')

        trajdict['radar_old2'] = trajdict['radar_old']
        trajdict['radar_old'] = radar
        return None, None

    # User defined parameter
    ang_tol = dscfg.get('ang_tol', 1.2)

    az_list = []
    el_list = []
    rr_list = []
    tt_list = []
    for tind in np.nditer(traj_ind):
        az = rad_traj.azimuth_vec[tind]
        el = rad_traj.elevation_vec[tind]
        rr = rad_traj.range_vec[tind]
        tt = traj_time_vec[tind]
        dBm = trajectory.dBm[tind]
        flashnr = trajectory.flashnr_vec[tind]

        (radar_sel, traj_ray_ind, traj_rng_ind, cell_ray_inds,
         cell_rng_ind_min, cell_rng_ind_max) = _get_gates(
             radar, az, el, rr, tt, trajdict, ang_tol=ang_tol)

        if radar_sel is None:
            continue

        # Get data samples and compute statistics
        for field_name in field_names:
            if field_name not in radar_sel.fields:
                warn("Datatype '%s' not available in radar data" % field_name)
                continue
            rdata = radar_sel.fields[field_name]['data']
            val = rdata[traj_ray_ind, traj_rng_ind]
            cell_vals = rdata[
                cell_ray_inds, cell_rng_ind_min:cell_rng_ind_max+1]

            # Compute statistics and get number of valid data
            if data_is_log[field_name]:
                val_mean = np.ma.masked
                cell_vals_valid = cell_vals.compressed()
                if cell_vals_valid.size > 0:
                    vals_lin = np.ma.power(10., cell_vals_valid/10.)
                    val_mean = np.ma.mean(vals_lin)
                    val_mean = 10. * np.ma.log10(val_mean)
            else:
                val_mean = np.ma.mean(cell_vals)
            val_min = np.ma.min(cell_vals)
            val_max = np.ma.max(cell_vals)

            nvals_valid = np.count_nonzero(
                np.logical_not(np.ma.getmaskarray(cell_vals)))

            # Add to time series dict
            ts_dict[field_name].add_timesample(
                trajectory.time_vector[tind],
                (flashnr, dBm, val, val_mean, val_min, val_max, nvals_valid))

        az_list.append(az)
        el_list.append(el)
        rr_list.append(rr)
        tt_list.append(trajectory.time_vector[tind])
        # end loop over traj samples within period

    # output radar volume and flash coordinates respect to radar
    radar_sel = radar
    if trajdict['radar_old'] is not None:
        radar_sel = trajdict['radar_old']

    dataset = {
        'azi_traj': np.asarray(az_list),
        'ele_traj': np.asarray(el_list),
        'rng_traj': np.asarray(rr_list),
        'time_traj': np.asarray(tt_list),
        'radar': radar_sel,
        'final': 0
    }

    # update trajectory dictionary
    trajdict['last_task_start_dt'] = dt_task_start
    trajdict['radar_old2'] = trajdict['radar_old']
    trajdict['radar_old'] = radar

    return dataset, ind_rad


def process_traj_atplane(procstatus, dscfg, radar_list=None, trajectory=None):
    """
    Return time series according to trajectory

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        data_is_log : dict. Dataset keyword
            Dictionary specifying for each field if it is in log (True) or
            linear units (False). Default False
        ang_tol : float. Dataset keyword
            Factor that multiplies the angle resolution. Used when determining
            the neighbouring rays. Default 1.2
    radar_list : list of Radar objects
        Optional. list of radar objects
    trajectory : Trajectory object
        containing trajectory samples

    Returns
    -------
    trajectory : Trajectory object
        Object holding time series
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        # first call: nothing to do
        return None, None

    if procstatus == 2:
        # last call: do the products
        if not dscfg['initialized']:
            warn('ERROR: No trajectory dataset available!')
            return None, None
        trajdict = dscfg['traj_atplane_dict']

        dataset = {
            'ts_dict': trajdict['ts_dict'],
            'final': 1
        }
        return dataset, trajdict['ind_rad']

    # Process
    field_names = []
    datatypes = []
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_names.append(get_fieldname_pyart(datatype))
        datatypes.append(datatype)

    ind_rad = int(radarnr[5:8])-1
    if ((radar_list is None) or (radar_list[ind_rad] is None)):
        warn('ERROR: No valid radar found')
        return None, None

    # keep locally only field of interest in radar object
    radar = deepcopy(radar_list[ind_rad])
    radar.fields = dict()
    nfields_available = 0
    for field_name in field_names:
        if field_name not in radar_list[ind_rad].fields:
            warn("Datatype '%s' not available in radar data" % field_name)
            continue
        radar.add_field(field_name, radar_list[ind_rad].fields[field_name])
        nfields_available += 1

    if nfields_available == 0:
        warn("Fields not available in radar data")
        return None, None

    ttask_start = radar.time['data'].min()
    dt_task_start = num2date(ttask_start, radar.time['units'],
                             radar.time['calendar'])
    if not dscfg['initialized']:
        # init
        if trajectory is None:
            raise Exception("ERROR: Undefined trajectory for dataset '%s'"
                            % dscfg['dsname'])

        rad_traj = trajectory.add_radar(radar)

        description = [
            "Description:",
            "Time series of a weather radar data type at the location",
            "of the plane.",
            "The time samples where the plane was out of the weather radar",
            "sector are NOT included in this file.",
            "NaN (Not a number): No rain detected at the plane location."
        ]

        ts_dict = dict()
        data_is_log = dict()
        for datatype, field_name in zip(datatypes, field_names):
            ts = TimeSeries(
                description, maxlength=trajectory.time_vector.size,
                datatype=datatype)

            unit = get_field_unit(datatype)
            name = get_field_name(datatype)
            # Append empty series: Note: sequence matters!
            ts.add_dataseries("at plane", name, unit, color='b')
            ts.add_dataseries("Mean", name, unit, color='r')
            ts.add_dataseries("Min", name, unit, color='k', linestyle=':')
            ts.add_dataseries("Max", name, unit, color='k', linestyle=':')
            ts.add_dataseries("#Valid", "", "", plot=False)

            ts_dict.update({field_name: ts})
            data_is_log.update({field_name: False})
            if 'data_is_log' in dscfg:
                if datatype in dscfg['data_is_log']:
                    data_is_log[field_name] = (dscfg['data_is_log'] != 0)
                else:
                    warn('Units type for data type '+datatype +
                         ' not specified. Assumed linear')

        trajdict = dict({
            'radar_traj': rad_traj,
            'radar_old': None,
            'radar_old2': None,
            'last_task_start_dt': None,
            'ind_rad': ind_rad,
            'ts_dict': ts_dict,
            'data_is_log': data_is_log})
        traj_ind = trajectory.get_samples_in_period(end=dt_task_start)

        dscfg['traj_atplane_dict'] = trajdict
        dscfg['initialized'] = True
    else:
        # init already done
        trajdict = dscfg['traj_atplane_dict']
        rad_traj = trajdict['radar_traj']
        traj_ind = trajectory.get_samples_in_period(
            start=trajdict['last_task_start_dt'], end=dt_task_start)
        ts_dict = trajdict['ts_dict']
        data_is_log = trajdict['data_is_log']

    traj_time_vec = date2num(trajectory.time_vector, radar.time['units'],
                             radar.time['calendar'])

    if np.size(traj_ind) == 0:
        warn('No trajectory samples within current period')

        trajdict['radar_old2'] = trajdict['radar_old']
        trajdict['radar_old'] = radar
        return None, None

    # User defined parameter
    ang_tol = dscfg.get('ang_tol', 1.2)

    az_list = []
    el_list = []
    rr_list = []
    tt_list = []
    for tind in np.nditer(traj_ind):
        az = rad_traj.azimuth_vec[tind]
        el = rad_traj.elevation_vec[tind]
        rr = rad_traj.range_vec[tind]
        tt = traj_time_vec[tind]

        (radar_sel, traj_ray_ind, traj_rng_ind, cell_ray_inds,
         cell_rng_ind_min, cell_rng_ind_max) = _get_gates(
             radar, az, el, rr, tt, trajdict, ang_tol=ang_tol)

        if radar_sel is None:
            continue

        # Get data samples and compute statistics
        for field_name in field_names:
            if field_name not in radar_sel.fields:
                warn("Datatype '%s' not available in radar data" % field_name)
                continue
            rdata = radar_sel.fields[field_name]['data']
            val = rdata[traj_ray_ind, traj_rng_ind]
            cell_vals = rdata[
                cell_ray_inds, cell_rng_ind_min:cell_rng_ind_max+1]

            # Compute statistics and get number of valid data
            if data_is_log[field_name]:
                val_mean = np.ma.masked
                cell_vals_valid = cell_vals.compressed()
                if cell_vals_valid.size > 0:
                    vals_lin = np.ma.power(10., cell_vals_valid/10.)
                    val_mean = np.ma.mean(vals_lin)
                    val_mean = 10. * np.ma.log10(val_mean)
            else:
                val_mean = np.ma.mean(cell_vals)
            val_min = np.ma.min(cell_vals)
            val_max = np.ma.max(cell_vals)

            nvals_valid = np.count_nonzero(
                np.logical_not(np.ma.getmaskarray(cell_vals)))

            # Add to time series dict
            ts_dict[field_name].add_timesample(
                trajectory.time_vector[tind],
                (val, val_mean, val_min, val_max, nvals_valid))

        az_list.append(az)
        el_list.append(el)
        rr_list.append(rr)
        tt_list.append(trajectory.time_vector[tind])
        # end loop over traj samples within period

    # output radar volume and antenna coordinates respect to radar
    radar_sel = radar
    if trajdict['radar_old'] is not None:
        radar_sel = trajdict['radar_old']

    dataset = {
        'azi_traj': np.asarray(az_list),
        'ele_traj': np.asarray(el_list),
        'rng_traj': np.asarray(rr_list),
        'time_traj': np.asarray(tt_list),
        'radar': radar_sel,
        'final': 0
    }

    # update trajectory dictionary
    trajdict['last_task_start_dt'] = dt_task_start
    trajdict['radar_old2'] = trajdict['radar_old']
    trajdict['radar_old'] = radar

    return dataset, ind_rad


def process_traj_antenna_pattern(procstatus, dscfg, radar_list=None,
                                 trajectory=None):
    """
    Process a new array of data volumes considering a plane
    trajectory. As result a timeseries with the values
    transposed for a given antenna pattern is created.
    The result is created when the LAST flag is set.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        datatype : list of string. Dataset keyword
            The input data types
        antennaType : str. Dataset keyword
            Type of antenna of the radar we want to get the view from. Can
            be AZIMUTH, ELEVATION, LOWBEAM, HIGHBEAM
        par_azimuth_antenna : dict. Global ekyword
            Dictionary containing the parameters of the PAR azimuth antenna,
            i.e. name of the file with the antenna elevation pattern and fixed
            antenna angle
        par_elevation_antenna : dict. Global keyword
            Dictionary containing the parameters of the PAR elevation antenna,
            i.e. name of the file with the antenna azimuth pattern and fixed
            antenna angle
        asr_lowbeam_antenna : dict. Global keyword
            Dictionary containing the parameters of the ASR low beam antenna,
            i.e. name of the file with the antenna elevation pattern and fixed
            antenna angle
        asr_highbeam_antenna : dict. Global keyword
            Dictionary containing the parameters of the ASR high beam antenna,
            i.e. name of the file with the antenna elevation pattern and fixed
            antenna angle
        target_radar_pos : dict. Global keyword
            Dictionary containing the latitude, longitude and altitude of
            the radar we want to get the view from. If not specifying it will
            assume the radar is collocated
        range_all : Bool. Dataset keyword
            If the real radar and the synthetic radar are co-located and this
            parameter is true the statistics are going to be computed using
            all the data from range 0 to the position of the plane. Default
            False
        rhi_resolution : Bool. Dataset keyword
            Resolution of the synthetic RHI used to compute the data as viewed
            from the synthetic radar [deg]. Default 0.5
        max_altitude : float. Dataset keyword
            Max altitude of the data to use when computing the view from the
            synthetic radar [m MSL]. Default 12000.
        latlon_tol : float. Dataset keyword
            The tolerance in latitude and longitude to determine which
            synthetic radar gates are co-located with real radar gates [deg].
            Default 0.04
        alt_tol : float. Datset keyword
            The tolerance in altitude to determine which synthetic
            radar gates are co-located with real radar gates [m]. Default 1000.
        distance_upper_bound : float. Dataset keyword
            The maximum distance where to look for a neighbour when
            determining which synthetic radar gates are co-located with real
            radar gates [m]. Default 1000.
        use_cKDTree : Bool. Dataset keyword
            Which function to use to find co-located real radar gates with the
            synthetic radar. If True a function using cKDTree from
            scipy.spatial is used. This function uses parameter
            distance_upper_bound. If False a native implementation is used
            that takes as parameters latlon_tol and alt_tol. Default True.
        pattern_thres : float. Dataset keyword
            The minimum of the sum of the weights given to each value in order
            to consider the weighted quantile valid. It is related to the
            number of valid data points
        data_is_log : dict. Dataset keyword
            Dictionary specifying for each field if it is in log (True) or
            linear units (False). Default False
        use_nans : dict. Dataset keyword
            Dictionary specyfing whether the nans have to be used in the
            computation of the statistics for each field. Default False
        nan_value : dict. Dataset keyword
            Dictionary with the value to use to substitute the NaN values when
            computing the statistics of each field. Default 0
    radar_list : list of Radar objects
        Optional. list of radar objects
    trajectory : Trajectory object
        containing trajectory samples

    Returns
    -------
    trajectory : Trajectory object
        Object holding time series
    ind_rad : int
        radar index
    """

    if procstatus == 0:
        # first call: nothing to do
        return None, None

    if procstatus == 2:
        # last call: do the products
        if not dscfg['initialized']:
            warn('ERROR: No trajectory dataset available!')
            return None, None
        tadict = dscfg['traj_antenna_dict']
        dataset = {
            'ts_dict': tadict['ts_dict'],
            'final': 1
        }
        return dataset, tadict['ind_rad']

    # Process
    field_names = []
    datatypes = []
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_names.append(get_fieldname_pyart(datatype))
        datatypes.append(datatype)

    ind_rad = int(radarnr[5:8])-1
    if ((radar_list is None) or (radar_list[ind_rad] is None)):
        warn('ERROR: No valid radar found')
        return None, None

    # keep locally only field of interest in radar object
    radar = deepcopy(radar_list[ind_rad])
    radar.fields = dict()
    nfields_available = 0
    for field_name in field_names:
        if field_name not in radar_list[ind_rad].fields:
            warn("Datatype '%s' not available in radar data" % field_name)
            continue
        radar.add_field(field_name, radar_list[ind_rad].fields[field_name])
        nfields_available += 1

    if nfields_available == 0:
        warn("Fields not available in radar data")
        return None, None

    ttask_start = radar.time['data'].min()
    dt_task_start = num2date(ttask_start, radar.time['units'],
                             radar.time['calendar'])
    if not dscfg['initialized']:
        # === init ============================================================
        if trajectory is None:
            raise Exception("ERROR: Undefined trajectory for dataset '%s'"
                            % dscfg['dsname'])

        # Check config
        if 'antennaType' not in dscfg:
            raise Exception("ERROR: Undefined 'antennaType' for dataset '%s'"
                            % dscfg['dsname'])
        if 'configpath' not in dscfg:
            raise Exception("ERROR: Undefined 'configpath' for dataset '%s'"
                            % dscfg['dsname'])

        if 'target_radar_pos' not in dscfg:
            radar_antenna_atsameplace = True
            target_radar = None

            warn('No target radar position specified. ' +
                 'The radars are assumed colocated')
            rad_traj = trajectory.add_radar(radar)
        else:
            radar_antenna_atsameplace = False

            # create dummy radar object with target radar specs
            latitude = get_metadata('latitude')
            longitude = get_metadata('longitude')
            altitude = get_metadata('altitude')

            latitude['data'] = np.array(
                [dscfg['target_radar_pos']['latitude']], dtype='float64')
            longitude['data'] = np.array(
                [dscfg['target_radar_pos']['longitude']], dtype='float64')
            altitude['data'] = np.array(
                [dscfg['target_radar_pos']['altitude']], dtype='float64')

            target_radar = TargetRadar(latitude, longitude, altitude)
            rad_traj = trajectory.add_radar(target_radar)

        if dscfg['antennaType'] == 'AZIMUTH':
            is_azimuth_antenna = True
            info = 'parAzAnt'
            description = [
                "Antenna: PAR Azimuth antenna",
                "Description:",
                "Time series of a weather radar data type at the location",
                "of the plane weighted by the antenna pattern of the PAR",
                "antenna.",
                "The time samples where the plane was out of the weather",
                "radar sector are NOT included in this file.",
                "NaN (Not a number): No rain detected at the plane location."
            ]
            if 'par_azimuth_antenna' not in dscfg:
                raise Exception("ERROR: Undefined 'par_azimuth_antenna' for"
                                " dataset '%s'" % dscfg['dsname'])

            patternfile = dscfg['configpath'] + 'antenna/' \
                + dscfg['par_azimuth_antenna']['elPatternFile']
            fixed_angle = dscfg['par_azimuth_antenna']['fixed_angle']

        elif dscfg['antennaType'] == 'ELEVATION':
            is_azimuth_antenna = False
            info = 'parElAnt'
            description = [
                "Antenna: PAR Elevation antenna",
                "Description:",
                "Time series of a weather radar data type at the location",
                "of the plane weighted by the antenna pattern of the PAR",
                "antenna.",
                "The time samples where the plane was out of the weather ",
                "radar sector are NOT included in this file.",
                "NaN (Not a number): No rain detected at the plane location."
            ]

            if 'par_elevation_antenna' not in dscfg:
                raise Exception("ERROR: Undefined 'par_elevation_antenna' for"
                                " dataset '%s'" % dscfg['dsname'])

            patternfile = dscfg['configpath'] + 'antenna/' \
                + dscfg['par_elevation_antenna']['azPatternFile']
            fixed_angle = dscfg['par_elevation_antenna']['fixed_angle']

        elif dscfg['antennaType'] == 'LOWBEAM':
            is_azimuth_antenna = True
            info = 'asrLowBeamAnt'
            description = [
                "Antenna: ASR low beam antenna",
                "Description:",
                "Time series of a weather radar data type at the location",
                "of the plane weighted by the antenna pattern of the ASR",
                "antenna.",
                "The time samples where the plane was out of the weather",
                "radar sector are NOT included in this file.",
                "NaN (Not a number): No rain detected at the plane location."
            ]
            if 'asr_lowbeam_antenna' not in dscfg:
                raise Exception("ERROR: Undefined 'asr_lowbeam_antenna' for"
                                " dataset '%s'" % dscfg['dsname'])

            patternfile = dscfg['configpath'] + 'antenna/' \
                + dscfg['asr_lowbeam_antenna']['elPatternFile']
            fixed_angle = dscfg['asr_lowbeam_antenna']['fixed_angle']

        elif dscfg['antennaType'] == 'HIGHBEAM':
            is_azimuth_antenna = True
            info = 'asrHighBeamAnt'
            description = [
                "Antenna: ASR high beam antenna",
                "Description:",
                "Time series of a weather radar data type at the location",
                "of the plane weighted by the antenna pattern of the ASR",
                "antenna.",
                "The time samples where the plane was out of the weather",
                "radar sector are NOT included in this file.",
                "NaN (Not a number): No rain detected at the plane location."
            ]
            if 'asr_highbeam_antenna' not in dscfg:
                raise Exception("ERROR: Undefined 'asr_highbeam_antenna' for"
                                " dataset '%s'" % dscfg['dsname'])

            patternfile = dscfg['configpath'] + 'antenna/' \
                + dscfg['asr_highbeam_antenna']['elPatternFile']
            patternfile_low = dscfg['configpath'] + 'antenna/' \
                + dscfg['asr_lowbeam_antenna']['elPatternFile']
            fixed_angle = dscfg['asr_highbeam_antenna']['fixed_angle']
        else:
            raise Exception("ERROR: Unexpected antenna type '%s' for dataset"
                            " '%s'" % (dscfg['antennaType'], dscfg['dsname']))

        # Read dataset config parameters:
        weight_threshold = dscfg.get('pattern_thres', 0.)

        do_all_ranges = False
        if 'range_all' in dscfg:
            if dscfg['range_all'] != 0:
                do_all_ranges = True

        # Config parameters for processing when the weather radar and the
        # antenna are not at the same place:
        rhi_resolution = dscfg.get('rhi_resolution', 0.5)  # [deg]
        max_altitude = dscfg.get('max_altitude', 12000.)  # [m]
        latlon_tol = dscfg.get('latlon_tol', 0.04)  # [deg]
        alt_tol = dscfg.get('alt_tol', 1000.)  # [m]
        distance_upper_bound = dscfg.get('distance_upper_bound', 1000.)
        use_cKDTree = dscfg.get('use_cKDTree', True)

        # Get antenna pattern and make weight vector
        try:
            if info == 'asrHighBeamAnt':
                antpattern = read_antenna_pattern(
                    patternfile, linear=True, twoway=False)
                antpattern_low = read_antenna_pattern(
                    patternfile_low, linear=True, twoway=False)
                antpattern['attenuation'] *= antpattern_low['attenuation']
            else:
                antpattern = read_antenna_pattern(patternfile, linear=True,
                                                  twoway=True)
        except Exception as ee:
            warn(str(ee))
            raise

        pattern_angles = antpattern['angle'] + fixed_angle
        if not is_azimuth_antenna:
            pattern_angles[pattern_angles < 0] += 360.
        pattern_angles[pattern_angles >= 360.] -= 360.

        if radar_antenna_atsameplace:
            if is_azimuth_antenna:
                scan_angles = np.sort(np.unique(
                    radar.elevation['data'].round(decimals=1)))
            else:
                scan_angles = np.sort(np.unique(
                    radar.azimuth['data'].round(decimals=1)))
        else:
            scan_angles = np.arange(0, 90, rhi_resolution, dtype=float)

        weightvec = np.empty(scan_angles.size, dtype=float)
        for kk in range(scan_angles.size):
            ind = np.argmin(np.abs(pattern_angles - scan_angles[kk]))
            weightvec[kk] = antpattern['attenuation'][ind]

        ts_dict = dict()
        data_is_log = dict()
        use_nans = dict()
        nan_value = dict()
        for datatype, field_name in zip(datatypes, field_names):
            ts = TimeSeries(
                description, maxlength=trajectory.time_vector.size,
                datatype=datatype)

            unit = get_field_unit(datatype)
            name = get_field_name(datatype)
            # Quantiles of interest
            quantiles = [
                {"val": 0.1, "plot": False, "color": None, "ltype": None},
                {"val": 0.2, "plot": True, "color": 'k', "ltype": ':'},
                {"val": 0.3, "plot": False, "color": None, "ltype": None},
                {"val": 0.4, "plot": False, "color": None, "ltype": None},
                {"val": 0.5, "plot": True, "color": 'r', "ltype": None},
                {"val": 0.6, "plot": False, "color": None, "ltype": None},
                {"val": 0.7, "plot": False, "color": None, "ltype": None},
                {"val": 0.8, "plot": True, "color": 'k', "ltype": ':'},
                {"val": 0.9, "plot": False, "color": None, "ltype": None},
                {"val": 0.95, "plot": False, "color": None, "ltype": None}]

            ts.add_dataseries("Weighted average", name, unit, color='b')
            for qq in quantiles:
                label = "Quantile_%4.2f" % qq["val"]
                ts.add_dataseries(label, name, unit, plot=qq["plot"],
                                  color=qq["color"], linestyle=qq["ltype"])
            ts.add_dataseries("#Valid", "", "", plot=False)
            ts_dict.update({field_name: ts})

            data_is_log.update({field_name: False})
            if 'data_is_log' in dscfg:
                if datatype in dscfg['data_is_log']:
                    data_is_log[field_name] = (
                        dscfg['data_is_log'][datatype] != 0)
                else:
                    warn('Units type for data type '+datatype +
                         ' not specified. Assumed linear')

            use_nans.update({field_name: False})
            if 'use_nans' in dscfg:
                if datatype in dscfg['use_nans']:
                    use_nans[field_name] = (
                        dscfg['use_nans'][datatype] != 0)
                else:
                    warn('Use of nans not specified for data type '+datatype +
                         ' not specified. Assumed not used')

            nan_value.update({field_name: 0.})
            if 'nan_value' in dscfg:
                if datatype in dscfg['nan_value']:
                    nan_value[field_name] = dscfg['nan_value'][datatype]
                else:
                    warn('NaN value not specified for data type '+datatype +
                         ' not specified. Assumed 0')

        quants = np.array([ee['val'] for ee in quantiles])

        # Persistent data structure
        tadict = dict({
            'radar_traj': rad_traj,
            'radar_old': None,
            'radar_old2': None,
            'target_radar': target_radar,
            'last_task_start_dt': None,
            'ind_rad': ind_rad,
            'is_azimuth_antenna': is_azimuth_antenna,
            'info': info,
            'scan_angles': scan_angles,
            'radar_antenna_atsameplace': radar_antenna_atsameplace,
            'weightvec': weightvec,
            'quantiles': quants,
            'use_nans': use_nans,
            'nan_value': nan_value,
            'weight_threshold': weight_threshold,
            'do_all_ranges': do_all_ranges,
            'max_altitude': max_altitude,
            'latlon_tol': latlon_tol,
            'alt_tol': alt_tol,
            'distance_upper_bound': distance_upper_bound,
            'use_cKDTree': use_cKDTree,
            'data_is_log': data_is_log,
            'ts_dict': ts_dict})
        traj_ind = trajectory.get_samples_in_period(end=dt_task_start)

        dscfg['traj_antenna_dict'] = tadict
        dscfg['initialized'] = True
        # end init
    else:
        # init already done
        tadict = dscfg['traj_antenna_dict']

        traj_ind = trajectory.get_samples_in_period(
            start=tadict['last_task_start_dt'], end=dt_task_start)

    if not _get_ts_values_antenna_pattern(
            radar, trajectory, tadict, traj_ind, field_names):
        return None, None

    tadict['last_task_start_dt'] = dt_task_start
    tadict['radar_old2'] = tadict['radar_old']
    tadict['radar_old'] = radar

    # Collect garbage
    gc.collect()

    return None, None


def _get_ts_values_antenna_pattern(radar, trajectory, tadict, traj_ind,
                                   field_names):
    """
    Get the time series values of a trajectory using a synthetic antenna
    pattern

    Parameters
    ----------
    radar : radar object
        The radar volume with the data
    trajectory : trajectory object
        The plane trajectory
    tadict : dict
        A dictionary containing parameters useful for trajectory computation
    traj_ind : array
        The indices of trajectory data within the current radar volume time
    field_names : list of str
        list of names of the radar field

    Returns
    -------
    result : Bool
        A flag signaling whether radar data matching the trajectory was found

    """
    rad_traj = tadict['radar_traj']
    is_azimuth_antenna = tadict['is_azimuth_antenna']
    scan_angles = tadict['scan_angles']
    radar_antenna_atsameplace = tadict['radar_antenna_atsameplace']
    nan_value = tadict['nan_value']
    use_nans = tadict['use_nans']
    weight_threshold = tadict['weight_threshold']
    do_all_ranges = tadict['do_all_ranges']
    ts_dict = tadict['ts_dict']
    target_radar = tadict['target_radar']
    max_altitude = tadict['max_altitude']
    latlon_tol = tadict['latlon_tol']
    alt_tol = tadict['alt_tol']
    distance_upper_bound = tadict['distance_upper_bound']
    use_cKDTree = tadict['use_cKDTree']
    data_is_log = tadict['data_is_log']

    traj_time_vec = date2num(trajectory.time_vector, radar.time['units'],
                             radar.time['calendar'])

    if np.size(traj_ind) == 0:
        warn('No trajectory samples within current period')
        return False

    for tind in np.nditer(traj_ind):
        az = rad_traj.azimuth_vec[tind]
        el = rad_traj.elevation_vec[tind]
        rr = rad_traj.range_vec[tind]
        tt = traj_time_vec[tind]

        # Select radar object, find closest azimuth and elevation ray
        (radar_sel, ray_sel, rr_ind, el_vec_rnd, az_vec_rnd) = \
            _get_closest_bin(az, el, rr, tt, radar, tadict)

        # Check if traj sample is within scan sector
        if (_sample_out_of_sector(az, el, rr, radar_sel, ray_sel,
                                  rr_ind, el_vec_rnd, az_vec_rnd)):
            continue

        if radar_antenna_atsameplace:
            # ==============================================================
            # Radar and scanning antenna are at the SAME place
            # ==============================================================

            # ==============================================================
            # Get sample at bin

            if is_azimuth_antenna:
                angles = radar_sel.azimuth['data']
                angles_scan = radar_sel.elevation['data']
                ray_angle = radar_sel.azimuth['data'][ray_sel]
            else:
                angles = radar_sel.elevation['data']
                angles_scan = radar_sel.azimuth['data']
                ray_angle = radar_sel.elevation['data'][ray_sel]

            d_angle = np.abs(angles - ray_angle)
            ray_inds = np.where(d_angle < 0.09)[0]
            angles_sortind = np.argsort(angles_scan[ray_inds])

            ray_inds = ray_inds[angles_sortind]
            angles_sorted = angles_scan[ray_inds]

            # Set default values
            avg = None
            qvals = np.array([None] * tadict['quantiles'].size)
            nvals_valid = None

            if ((scan_angles.size != angles_sorted.size) or
                    (np.max(np.abs(scan_angles - angles_sorted)) > 0.1)):
                warn("Scan angle mismatch!")
                for field_name in field_names:
                    ts_dict[field_name].add_timesample(
                        trajectory.time_vector[tind],
                        np.concatenate([[avg], qvals, [nvals_valid]]))
                continue

            if do_all_ranges:
                rr_ind_min = 0
            else:
                rr_ind_min = rr_ind

            w_vec = tadict['weightvec']
            for field_name in field_names:
                if field_name not in radar_sel.fields:
                    warn("Datatype '%s' not available in radar data" %
                         field_name)
                    continue
                rdata = radar_sel.fields[field_name]['data']
                values = rdata[ray_inds, rr_ind_min:rr_ind+1]
                if use_nans[field_name]:
                    values_ma = np.ma.getmaskarray(values)
                    values[values_ma] = nan_value[field_name]
                try:
                    (avg, qvals, nvals_valid) = quantiles_weighted(
                        values,
                        weight_vector=w_vec,
                        quantiles=tadict['quantiles'],
                        weight_threshold=weight_threshold,
                        data_is_log=data_is_log[field_name])
                except Exception as ee:
                    warn(str(ee))
                    continue

                ts_dict[field_name].add_timesample(
                    trajectory.time_vector[tind],
                    np.concatenate([[avg], qvals, [nvals_valid]]))
        else:
            # ================================================================
            # Radar and scanning antenna are NOT at the same place
            # ================================================================
            ray_inds, rng_inds, w_inds = _get_gates_antenna_pattern(
                radar_sel, target_radar, az, rr, tt, scan_angles,
                alt_tol=alt_tol, latlon_tol=latlon_tol,
                max_altitude=max_altitude,
                distance_upper_bound=distance_upper_bound,
                use_cKDTree=use_cKDTree)

            w_vec = tadict['weightvec'][w_inds]

            for field_name in field_names:
                if field_name not in radar_sel.fields:
                    warn("Datatype '%s' not available in radar data" %
                         field_name)
                    continue
                rdata = radar_sel.fields[field_name]['data']
                values = rdata[ray_inds, rng_inds]
                if use_nans[field_name]:
                    values_ma = np.ma.getmaskarray(values)
                    values[values_ma] = nan_value[field_name]

                try:
                    (avg, qvals, nvals_valid) = quantiles_weighted(
                        values,
                        weight_vector=w_vec,
                        quantiles=tadict['quantiles'],
                        weight_threshold=weight_threshold,
                        data_is_log=data_is_log[field_name])
                except Exception as ee:
                    warn(str(ee))
                    continue

                ts_dict[field_name].add_timesample(
                    trajectory.time_vector[tind],
                    np.concatenate([[avg], qvals, [nvals_valid]]))
    return True


def _get_contour_trt(radar, trajectory, voltime, time_tol=100.):
    """
    Get the TRT cell contour corresponding to the current radar

    Parameters
    ----------
    radar : radar object
        The radar containing
    trajectory : trajectory object
        Object containing the TRT cell position and dimensions
    voltime : datetime object
        The radar volume reference time
    time_tol : float
        Time tolerance where to look for data [s]

    Returns
    -------
    lon_roi, lat_roi : array of floats
        The lat/lon defining the TRT cell contour

    """
    dt = np.empty(trajectory.time_vector.size, dtype=float)
    for i, time_traj in enumerate(trajectory.time_vector):
        dt[i] = np.abs((voltime - time_traj).total_seconds())
    if dt.min() > time_tol:
        warn('No TRT data for radar volume time {}'.format(voltime))
        return None
    ind = np.argmin(dt)

    lon_roi = trajectory.cell_contour[ind]['lon']
    lat_roi = trajectory.cell_contour[ind]['lat']
    lat_center = trajectory.wgs84_lat_deg[ind]
    lon_center = trajectory.wgs84_lon_deg[ind]

    roi_dict = {
        'lon': lon_roi,
        'lat': lat_roi,
        'lon_center': lon_center,
        'lat_center': lat_center}

    # extract the data within the ROI boundaries
    inds_ray, inds_rng = np.indices(np.shape(radar.gate_longitude['data']))

    mask = np.logical_and(
        np.logical_and(
            radar.gate_latitude['data'] >= roi_dict['lat'].min(),
            radar.gate_latitude['data'] <= roi_dict['lat'].max()),
        np.logical_and(
            radar.gate_longitude['data'] >= roi_dict['lon'].min(),
            radar.gate_longitude['data'] <= roi_dict['lon'].max()))

    if np.all(mask == 0):
        warn('No values within ROI')
        return None

    inds_ray = inds_ray[mask]
    inds_rng = inds_rng[mask]

    # extract the data inside the ROI
    lat = radar.gate_latitude['data'][mask]
    lon = radar.gate_longitude['data'][mask]
    inds, is_roi = belongs_roi_indices(lat, lon, roi_dict)

    if is_roi == 'None':
        warn('No values within ROI')
        return None

    return roi_dict


def _get_gates(radar, az, el, rr, tt, trajdict, ang_tol=1.2):
    """
    Find the gates of the radar object that have to be used to compute the
    data of a trajectory

    Parameters
    ----------
    radar : radar object
        The radar containing
    az, el, rr : floats
        The trajectory position respect to the radar
    tt : float
        the trajectory time respect to the beginning of the radar scan
    trajdict : dict
        Dictionary containing the trajectory parameters
    ang_tol : float
        Factor that multiplies the angle resolution. Used when determining
        the neighbouring rays

    Returns
    -------
    radar_sel : radar object
        The radar volume selected as closest to trajectory point
    ray_sel, rr_ind : ints
        ray and range indices of the radar gate closest to the trajectory
        position
    cell_ind : array of ints
        indices of the surrounding rays
    rr_min : int
        index of the minimum range of the surrounding gates
    rr_max : int
        index of the maximum range of the surrounding gates

    """
    # Find closest azimuth and elevation ray
    (radar_sel, ray_sel, rr_ind, el_vec_rnd, az_vec_rnd) = \
        _get_closest_bin(az, el, rr, tt, radar, trajdict)

    # Check if traj sample is within scan sector
    if (_sample_out_of_sector(az, el, rr, radar_sel, ray_sel,
                              rr_ind, el_vec_rnd, az_vec_rnd)):
        return None, None, None, None, None, None

    # Get indices of gates surrounding the cell (3x3 box)
    el_res = np.ma.median(np.ma.diff(el_vec_rnd))
    az_res = np.ma.median(np.ma.diff(az_vec_rnd))

    d_el = np.abs(radar_sel.elevation['data'] -
                  radar_sel.elevation['data'][ray_sel])
    d_az = np.abs(radar_sel.azimuth['data'] -
                  radar_sel.azimuth['data'][ray_sel])
    cell_ind = np.where((d_el < el_res*ang_tol) & (d_az < az_res*ang_tol))

    rr_min = rr_ind-1
    if rr_min < 0:
        rr_min = 0

    return radar_sel, ray_sel, rr_ind, cell_ind, rr_min, rr_ind+1


def _get_gates_trt(radar, trajectory, voltime, time_tol=100., alt_min=None,
                   alt_max=None, cell_center=False, latlon_tol=0.0005):
    """
    Find the gates of the radar object that belong to a TRT cell

    Parameters
    ----------
    radar : radar object
        The radar containing
    trajectory : trajectory object
        Object containing the TRT cell position and dimensions
    voltime : datetime object
        The radar volume reference time
    time_tol : float
        Time tolerance where to look for data [s]
    alt_min, alt_max : float
        Minimum and maximum altitude where to look for data [m]

    Returns
    -------
    inds_ray, inds_rng : array of ints
        The indices of the radar data inside the TRT cell

    """
    dt = np.empty(trajectory.time_vector.size, dtype=float)
    for i, time_traj in enumerate(trajectory.time_vector):
        dt[i] = np.abs((voltime - time_traj).total_seconds())
    if dt.min() > time_tol:
        warn('No TRT data for radar volume time')
        return None, None, None, None, None
    ind = np.argmin(dt)

    if not cell_center:
        lon_roi = trajectory.cell_contour[ind]['lon']
        lat_roi = trajectory.cell_contour[ind]['lat']

        roi_dict = {
            'lon': lon_roi,
            'lat': lat_roi,
            'alt_min': alt_min,
            'alt_max': alt_max}

        # extract the data within the ROI boundaries
        inds_ray, inds_rng = np.indices(np.shape(radar.gate_longitude['data']))

        mask = np.logical_and(
            np.logical_and(
                radar.gate_latitude['data'] >= roi_dict['lat'].min(),
                radar.gate_latitude['data'] <= roi_dict['lat'].max()),
            np.logical_and(
                radar.gate_longitude['data'] >= roi_dict['lon'].min(),
                radar.gate_longitude['data'] <= roi_dict['lon'].max()))

        if alt_min is not None:
            mask[radar.gate_altitude['data'] < alt_min] = 0
        if alt_max is not None:
            mask[radar.gate_altitude['data'] > alt_max] = 0

        if np.all(mask == 0):
            warn('No values within ROI')
            return None, None, None, None, None

        inds_ray = inds_ray[mask]
        inds_rng = inds_rng[mask]

        # extract the data inside the ROI
        lat = radar.gate_latitude['data'][mask]
        lon = radar.gate_longitude['data'][mask]
        inds, is_roi = belongs_roi_indices(lat, lon, roi_dict)

        if is_roi == 'None':
            warn('No values within ROI')
            return None, None, None, None, None

        inds_ray = inds_ray[inds]
        inds_rng = inds_rng[inds]

        lat = lat[inds]
        lon = lon[inds]
        alt = radar.gate_altitude['data'][inds_ray, inds_rng]
    else:
        radar_aux = deepcopy(radar)
        # transform radar into ppi over the required elevation
        if radar_aux.scan_type == 'rhi':
            target_elevations, el_tol = get_target_elevations(radar_aux)
            radar_ppi = cross_section_rhi(
                radar_aux, target_elevations, el_tol=el_tol)
        elif radar_aux.scan_type == 'ppi':
            radar_ppi = radar_aux
        else:
            warn('Error: unsupported scan type.')
            return None, None, None, None, None

        inds_ray = np.array([], dtype=int)
        inds_rng = np.array([], dtype=int)
        lat = np.array([])
        lon = np.array([])
        alt = np.array([])
        for sweep in range(radar_ppi.nsweeps):
            radar_aux = deepcopy(radar_ppi.extract_sweeps([sweep]))

            # find nearest gate to lat lon point
            ind_ray, ind_rng, _, _ = find_nearest_gate(
                radar_aux, trajectory.wgs84_lat_deg[ind],
                trajectory.wgs84_lon_deg[ind], latlon_tol=latlon_tol)

            if ind_ray is None:
                continue

            if alt_min is not None:
                if radar_aux.gate_altitude['data'][ind_ray][ind_rng] < alt_min:
                    continue
            if alt_max is not None:
                if radar_aux.gate_altitude['data'][ind_ray][ind_rng] > alt_max:
                    continue

            inds_ray = np.append(inds_ray, ind_ray)
            inds_rng = np.append(inds_rng, ind_rng)
            lat = np.append(
                lat, radar_aux.gate_latitude['data'][ind_ray][ind_rng])
            lon = np.append(
                lon, radar_aux.gate_longitude['data'][ind_ray][ind_rng])
            alt = np.append(
                alt, radar_aux.gate_altitude['data'][ind_ray][ind_rng])

        if inds_ray.size == 0:
            warn('No values in center of cell')
            return None, None, None, None, None

    return inds_ray, inds_rng, lat, lon, alt


def _get_gates_antenna_pattern(radar_sel, target_radar, az, rr, tt,
                               scan_angles, alt_tol=1000., latlon_tol=0.04,
                               max_altitude=12000.,
                               distance_upper_bound=1000., use_cKDTree=True):
    """
    Find the gates of the radar object that have to be used to compute the
    data of a trajectory as seen by another radar system

    Parameters
    ----------
    radar_sel : radar object
        The radar containing real data
    target_radar : radar object
        The virtual radar
    az, rr : floats
        The trajectory position respect to the radar
    tt : float
        the trajectory time respect to the beginning of the radar scan
    scan_angles : array
        The scan angles of the virtual radar object
    alt_tol : float
        The tolerance in altitude [m]
    latlon_tol : float
        The tolerance in latitude and longitude [deg]
    max_altitude : float
        The maximum altitude where to look for radar data [m MSL]
    distance_upper_bound : float
        The maximum distance where to look for neighbors [m]
    use_cKDTree : Bool
        If True the nearest neighbour to the synthetic radar is found using
        the cKDTree function of scipy.spatial. Otherwise it is found by Pyrad
        implementation

    Returns
    -------
    ray_ind, rng_ind : array of ints
        the indices of the radar data to use
    w_ind : array of ints
        The indices of the one-dimensional antenna pattern to use

    """
    # Create dummy Radar object with fix az and r, and all elevations
    n_rays = scan_angles.size

    r_time = {'data': np.array([tt * n_rays])}
    r_range = {'data': np.array([rr])}
    r_azimuth = get_metadata('azimuth')
    r_azimuth['data'] = np.ones(n_rays) * az
    r_elevation = {'data': scan_angles}
    r_sweep_number = {'data': np.array([0])}
    r_fields = {'colocated_gates': get_metadata('colocated_gates')}
    r_fields['colocated_gates']['data'] = 2*np.ma.ones((n_rays, 1), dtype=int)

    r_radar = Radar(r_time, r_range, r_fields, None, 'rhi',
                    target_radar.latitude, target_radar.longitude,
                    target_radar.altitude,
                    r_sweep_number, None, None, None, None,
                    r_azimuth, r_elevation)

    # flag regions with colocated usable data in r_radar
    r_ind_invalid = r_radar.gate_altitude['data'] > max_altitude
    r_radar.fields['colocated_gates']['data'][r_ind_invalid] = 1

    # flag regions with colocated usable data in radar_sel
    gate_coloc_radar_sel = pyart.util.intersection(
        radar_sel, r_radar, h_tol=alt_tol, latlon_tol=latlon_tol,
        vol_d_tol=None, vismin=None, hmin=None, hmax=max_altitude,
        rmin=None, rmax=None, elmin=None, elmax=None, azmin=None,
        azmax=None, visib_field=None,
        intersec_field='colocated_gates')
    radar_sel.add_field(
        'colocated_gates', gate_coloc_radar_sel, replace_existing=True)

    if use_cKDTree:
        colgates, r_radar_colg = pyart.util.colocated_gates2(
            r_radar, radar_sel, distance_upper_bound=distance_upper_bound)
    else:
        colgates, r_radar_colg = pyart.util.colocated_gates(
            r_radar, radar_sel, h_tol=alt_tol, latlon_tol=latlon_tol)

    w_ind = np.where(r_radar_colg['data'] > 1)[0]

    return colgates['rad2_ray_ind'], colgates['rad2_rng_ind'], w_ind


def _get_closest_bin(az, el, rr, tt, radar, tdict):
    """
    Get the radar bin closest to a certain trajectory position

    Parameters
    ----------
    az, el, rr : floats
        The trajectory position respect to the radar
    tt : float
        the trajectory time respect to the beginning of the radar scan
    radar : radar object
        the current radar object
    tdict : dict
        Dictionary containing trajectory parameters

    Returns
    -------
    radar_sel : radar object
        The selected radar (Current or one of the two previous ones)
    ray_sel, rr_ind : int
        The selected ray and range indices of the radar field
    el_vec_rnd, az_vec_rnd : array of floats
        The elevation and azimuth fields of the selected radar rounded to the
        first decimal

    """
    daz = np.abs(radar.azimuth['data'] - az)
    dele = np.abs(radar.elevation['data'] - el)
    dangle = np.sqrt(daz**2 + dele**2)
    ray_ind = np.argmin(dangle)
    dt = tt - radar.time['data'][ray_ind]

    if tdict['radar_old'] is not None:
        rad_old = tdict['radar_old']

        daz_old = np.abs(rad_old.azimuth['data'] - az)
        dele_old = np.abs(rad_old.elevation['data'] - el)
        dangle_old = np.sqrt(daz_old**2 + dele_old**2)
        ray_ind_old = np.argmin(dangle_old)
        dt_old = tt - rad_old.time['data'][ray_ind_old]

        # Find closest time
        if ((dt_old < 0.) and (tdict['radar_old2'] is not None)):
            rad_old2 = tdict['radar_old2']
            daz_old2 = np.abs(rad_old2.azimuth['data'] - az)
            dele_old2 = np.abs(rad_old2.elevation['data'] - el)
            dangle_old2 = np.sqrt(daz_old2**2 + dele_old2**2)
            ray_ind_old2 = np.argmin(dangle_old2)
            dt_old2 = tt - rad_old2.time['data'][ray_ind_old2]

            if np.abs(dt_old) < np.abs(dt_old2):
                radar_sel = rad_old
                ray_sel = ray_ind_old
            else:
                radar_sel = rad_old2
                ray_sel = ray_ind_old2
        else:
            if np.abs(dt) < np.abs(dt_old):
                radar_sel = radar
                ray_sel = ray_ind
            else:
                radar_sel = rad_old
                ray_sel = ray_ind_old
    else:
        radar_sel = radar
        ray_sel = ray_ind

    # Find closest range bin
    rr_ind = np.argmin(np.abs(radar_sel.range['data'] - rr))

    el_vec_rnd = np.sort(np.unique(
        radar_sel.elevation['data'].round(decimals=1)))
    az_vec_rnd = np.sort(np.unique(
        radar_sel.azimuth['data'].round(decimals=1)))

    return (radar_sel, ray_sel, rr_ind, el_vec_rnd, az_vec_rnd)


def _sample_out_of_sector(az, el, rr, radar_sel, ray_sel, rr_ind,
                          el_vec_rnd, az_vec_rnd):
    """
    Check if trajectory sample is within radar sector

    Parameters
    ----------
    az, el, rr : floats
        The trajectory position respect to the radar
    radar_sel : radar object
        The selected radar (Current or one of the two previous ones)
    ray_sel, rr_ind : int
        The selected ray and range indices of the radar field
    el_vec_rnd, az_vec_rnd : array of floats
        The elevation and azimuth fields of the selected radar rounded to the
        first decimal

    Returns
    -------
    result : bool
        False if the sample is out of sector. True otherwise

    """
    # Check if sample is within sector
    rr_min = radar_sel.range['data'][rr_ind]
    range_res = radar_sel.range['data'][1] - radar_sel.range['data'][0]
    if np.abs(rr_min - rr) > (2*range_res):
        # print("INFO: Trajectory sample out of range")
        return True

    if (((az_vec_rnd[0] - az) > 3.0) or ((az - az_vec_rnd[-1]) > 3.0)):
        # print("INFO: Trajectory sample out of azimuth angle")
        return True

    if (((el_vec_rnd[0] - el) > 3.0) or ((el - el_vec_rnd[-1]) > 3.0)):
        # print("INFO: Trajectory sample out of elevation angle")
        return True

    return False


class TargetRadar:
    """
    A class for dummy target radar object

    Attributes
    ----------
    latitude, longitude, altitude : float
        Position of the dummy radar

    """
    def __init__(self, latitude, longitude, altitude):
        self.latitude = latitude
        self.longitude = longitude
        self.altitude = altitude
