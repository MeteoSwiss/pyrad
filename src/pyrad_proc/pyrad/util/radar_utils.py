"""
pyrad.util.radar_utils
======================

Miscellaneous functions dealing with radar data

.. autosummary::
    :toctree: generated/

    get_data_along_rng
    get_data_along_azi
    get_data_along_ele
    get_ROI
    rainfall_accumulation
    time_series_statistics
    find_contiguous_times
    join_time_series
    get_range_bins_to_avg
    get_cercle_coords
    belongs_roi_indices
    find_ray_index
    find_rng_index
    find_ang_index
    find_nearest_gate
    find_colocated_indexes
    time_avg_range
    get_closest_solar_flux
    get_fixed_rng_data
    create_sun_hits_field
    create_sun_retrieval_field
    compute_quantiles
    compute_quantiles_from_hist
    compute_quantiles_sweep
    compute_histogram
    compute_histogram_sweep
    get_histogram_bins
    compute_2d_stats
    compute_1d_stats
    compute_2d_hist
    quantize_field
    compute_profile_stats
    project_to_vertical

"""
from warnings import warn
from copy import deepcopy
import datetime

import numpy as np
import scipy

try:
    import shapely
    _SHAPELY_AVAILABLE = True
except ImportError:
    warn('shapely not available')
    _SHAPELY_AVAILABLE = False

try:
    import pandas as pd
    _PANDAS_AVAILABLE = True
except ImportError:
    warn('Pandas not available')
    _PANDAS_AVAILABLE = False

import pyart

from .stat_utils import quantiles_weighted


def get_data_along_rng(radar, field_name, fix_elevations, fix_azimuths,
                       ang_tol=1., rmin=None, rmax=None, use_altitude=False):
    """
    Get data at particular (azimuths, elevations)

    Parameters
    ----------
    radar : radar object
        the radar object where the data is
    field_name : str
        name of the field to filter
    fix_elevations, fix_azimuths: list of floats or None
        List of elevations, azimuths couples [deg] if one of them is None
        all angles corresponding to one of the fixed angles specified will
        be gathered
    ang_tol : float
        Tolerance between the nominal angle and the radar angle [deg]
    rmin, rmax: float
        Min and Max range of the obtained data [m]
    use_altitude: bool
        If True instead of range the x_vals is gate altitude

    Returns
    -------
    xvals : list of float arrays
        The ranges of each azi, ele pair
    yvals : list of float arrays
        The values
    valid_azi, valid_ele : float arrays
        The azi, ele pairs

    """
    if rmin is None:
        rmin = 0.
    if rmax is None:
        rmax = np.max(radar.range['data'])

    xvals = []
    yvals = []
    valid_azi = []
    valid_ele = []
    if radar.scan_type in ('ppi', 'vertical_pointing'):
        if fix_elevations is None:
            warn('at least one fixed elevation has to be specified'
                 'when operating on PPI scans')
            return None, None, None, None
        if fix_azimuths is None:
            for ele in fix_elevations:
                ind_sweep = find_ang_index(
                    radar.fixed_angle['data'], ele, ang_tol=ang_tol)
                if ind_sweep is None:
                    warn('No elevation angle found '
                         'for fix_elevation {}'.format(ele))
                    continue
                if 'vertical_pointing':
                    new_dataset = deepcopy(radar)
                else:
                    new_dataset = radar.extract_sweeps([ind_sweep])

                for ind_azi, azi in enumerate(new_dataset.azimuth['data']):
                    if not use_altitude:
                        rng_mask = np.logical_and(
                            radar.range['data'] >= rmin,
                            radar.range['data'] <= rmax)

                        x = radar.range['data'][rng_mask]
                    else:
                        x = deepcopy(
                            new_dataset.gate_altitude['data'][ind_azi, :])
                        rng_mask = np.logical_and(x >= rmin, x <= rmax)
                        x = x[rng_mask]

                    yvals.append(
                        new_dataset.fields[field_name]['data'][
                            ind_azi, rng_mask])
                    xvals.append(x)
                    valid_azi.append(new_dataset.azimuth['data'][ind_azi])
                    valid_ele.append(new_dataset.elevation['data'][ind_azi])

            return xvals, yvals, valid_azi, valid_ele

        for ele, azi in zip(fix_elevations, fix_azimuths):
            ind_sweep = find_ang_index(
                radar.fixed_angle['data'], ele, ang_tol=ang_tol)
            if ind_sweep is None:
                warn('No elevation angle found for fix_elevation '+str(ele))
                continue
            new_dataset = radar.extract_sweeps([ind_sweep])

            try:
                dataset_line = pyart.util.cross_section_ppi(
                    new_dataset, [azi], az_tol=ang_tol)
            except EnvironmentError:
                warn(' No data found at azimuth '+str(azi) +
                     ' and elevation '+str(ele))
                continue
            if not use_altitude:
                rng_mask = np.logical_and(
                    radar.range['data'] >= rmin, radar.range['data'] <= rmax)

                x = radar.range['data'][rng_mask]
            else:
                x = deepcopy(dataset_line.gate_altitude['data'][0, :])
                rng_mask = np.logical_and(x >= rmin, x <= rmax)
                x = x[rng_mask]

            yvals.append(dataset_line.fields[field_name]['data'][0, rng_mask])
            xvals.append(x)
            valid_azi.append(dataset_line.azimuth['data'][0])
            valid_ele.append(dataset_line.elevation['data'][0])

        return xvals, yvals, valid_azi, valid_ele

    if fix_elevations is None or fix_azimuths is None:
        warn('at least one couple fixed elevation/fixed azimuth has to be'
             ' specified when operating on scans that are not PPI')
        return None, None, None, None

    for ele, azi in zip(fix_elevations, fix_azimuths):
        ind_sweep = find_ang_index(
            radar.fixed_angle['data'], azi, ang_tol=ang_tol)
        if ind_sweep is None:
            warn('No azimuth angle found for fix_azimuth '+str(azi))
            continue
        new_dataset = radar.extract_sweeps([ind_sweep])

        try:
            dataset_line = pyart.util.cross_section_rhi(
                new_dataset, [ele], el_tol=ang_tol)
        except EnvironmentError:
            warn(' No data found at azimuth '+str(azi) +
                 ' and elevation '+str(ele))
            continue
        if not use_altitude:
            rng_mask = np.logical_and(
                radar.range['data'] >= rmin, radar.range['data'] <= rmax)

            x = radar.range['data'][rng_mask]
        else:
            x = deepcopy(dataset_line.gate_altitude['data'][0, :])
            rng_mask = np.logical_and(x >= rmin, x <= rmax)
            x = x[rng_mask]

        yvals.append(
            dataset_line.fields[field_name]['data'][0, rng_mask])
        xvals.append(x)
        valid_azi.append(dataset_line.azimuth['data'][0])
        valid_ele.append(dataset_line.elevation['data'][0])

    return xvals, yvals, valid_azi, valid_ele


def get_data_along_azi(radar, field_name, fix_ranges, fix_elevations,
                       rng_tol=50., ang_tol=1., azi_start=None,
                       azi_stop=None):
    """
    Get data at particular (ranges, elevations)

    Parameters
    ----------
    radar : radar object
        the radar object where the data is
    field_name : str
        name of the field to filter
    fix_ranges, fix_elevations: list of floats
        List of ranges [m], elevations [deg] couples
    rng_tol : float
        Tolerance between the nominal range and the radar range [m]
    ang_tol : float
        Tolerance between the nominal angle and the radar angle [deg]
    azi_start, azi_stop: float
        Start and stop azimuth angle of the data [deg]

    Returns
    -------
    xvals : list of float arrays
        The ranges of each rng, ele pair
    yvals : list of float arrays
        The values
    valid_rng, valid_ele : float arrays
        The rng, ele pairs

    """
    if azi_start is None:
        azi_start = np.min(radar.azimuth['data'])
    if azi_stop is None:
        azi_stop = np.max(radar.azimuth['data'])

    yvals = []
    xvals = []
    valid_rng = []
    valid_ele = []
    for rng, ele in zip(fix_ranges, fix_elevations):
        ind_rng = find_rng_index(radar.range['data'], rng, rng_tol=rng_tol)
        if ind_rng is None:
            warn('No range gate found for fix_range '+str(rng))
            continue

        if radar.scan_type == 'ppi':
            ind_sweep = find_ang_index(
                radar.fixed_angle['data'], ele, ang_tol=ang_tol)
            if ind_sweep is None:
                warn('No elevation angle found for fix_elevation ' +
                     str(ele))
                continue
            new_dataset = radar.extract_sweeps([ind_sweep])
        else:
            try:
                new_dataset = pyart.util.cross_section_rhi(
                    radar, [ele], el_tol=ang_tol)
            except EnvironmentError:
                warn(
                    ' No data found at range '+str(rng) +
                    ' and elevation '+str(ele))
                continue
        if azi_start < azi_stop:
            azi_mask = np.logical_and(
                new_dataset.azimuth['data'] >= azi_start,
                new_dataset.azimuth['data'] <= azi_stop)
        else:
            azi_mask = np.logical_or(
                new_dataset.azimuth['data'] >= azi_start,
                new_dataset.azimuth['data'] <= azi_stop)
        yvals.append(
            new_dataset.fields[field_name]['data'][azi_mask, ind_rng])
        xvals.append(new_dataset.azimuth['data'][azi_mask])
        valid_rng.append(new_dataset.range['data'][ind_rng])
        valid_ele.append(new_dataset.elevation['data'][0])

    return xvals, yvals, valid_rng, valid_ele


def get_data_along_ele(radar, field_name, fix_ranges, fix_azimuths,
                       rng_tol=50., ang_tol=1., ele_min=None,
                       ele_max=None):
    """
    Get data at particular (ranges, azimuths)

    Parameters
    ----------
    radar : radar object
        the radar object where the data is
    field_name : str
        name of the field to filter
    fix_ranges, fix_azimuths: list of floats
        List of ranges [m], azimuths [deg] couples
    rng_tol : float
        Tolerance between the nominal range and the radar range [m]
    ang_tol : float
        Tolerance between the nominal angle and the radar angle [deg]
    ele_min, ele_max: float
        Min and max elevation angle [deg]

    Returns
    -------
    xvals : list of float arrays
        The ranges of each rng, ele pair
    yvals : list of float arrays
        The values
    valid_rng, valid_ele : float arrays
        The rng, ele pairs

    """
    if ele_min is None:
        ele_min = np.min(radar.elevation['data'])
    if ele_max is None:
        ele_max = np.max(radar.elevation['data'])

    yvals = []
    xvals = []
    valid_rng = []
    valid_azi = []
    for rng, azi in zip(fix_ranges, fix_azimuths):
        ind_rng = find_rng_index(radar.range['data'], rng, rng_tol=rng_tol)
        if ind_rng is None:
            warn('No range gate found for fix_range '+str(rng))
            continue

        if radar.scan_type == 'ppi':
            try:
                new_dataset = pyart.util.cross_section_ppi(
                    radar, [azi], az_tol=ang_tol)
            except EnvironmentError:
                warn(
                    ' No data found at range '+str(rng) +
                    ' and elevation '+str(azi))
                continue
        else:
            ind_sweep = find_ang_index(
                radar.fixed_angle['data'], azi, ang_tol=ang_tol)
            if ind_sweep is None:
                warn('No azimuth angle found for fix_azimuth '+str(azi))
                continue
            new_dataset = radar.extract_sweeps([ind_sweep])

        ele_mask = np.logical_and(
            new_dataset.elevation['data'] >= ele_min,
            new_dataset.elevation['data'] <= ele_max)
        yvals.append(
            new_dataset.fields[field_name]['data'][ele_mask, ind_rng])
        xvals.append(new_dataset.elevation['data'][ele_mask])
        valid_rng.append(new_dataset.range['data'][ind_rng])
        valid_azi.append(new_dataset.elevation['data'][0])

    return xvals, yvals, valid_rng, valid_azi


def get_ROI(radar, fieldname, sector):
    """
    filter out any data outside the region of interest defined by sector

    Parameters
    ----------
    radar : radar object
        the radar object where the data is
    fieldname : str
        name of the field to filter
    sector : dict
        a dictionary defining the region of interest

    Returns
    -------
    roi_flag : ndarray
        a field array with ones in gates that are in the Region of Interest

    """
    roi_flag = np.ma.ones((radar.nrays, radar.ngates), dtype=int)

    # check for altitude limits
    if sector['hmin'] is not None:
        roi_flag[radar.gate_altitude['data'] < sector['hmin']] = 0
    if sector['hmax'] is not None:
        roi_flag[radar.gate_altitude['data'] > sector['hmax']] = 0

    # check for range limits
    if sector['rmin'] is not None:
        roi_flag[:, radar.range['data'] < sector['rmin']] = 0
    if sector['rmax'] is not None:
        roi_flag[:, radar.range['data'] > sector['rmax']] = 0

    # check elevation angle limits
    if sector['elmin'] is not None:
        roi_flag[radar.elevation['data'] < sector['elmin'], :] = 0
    if sector['elmax'] is not None:
        roi_flag[radar.elevation['data'] > sector['elmax'], :] = 0

    # check min and max azimuth angle
    if sector['azmin'] is not None and sector['azmax'] is not None:
        if sector['azmin'] <= sector['azmax']:
            roi_flag[radar.azimuth['data'] < sector['azmin'], :] = 0
            roi_flag[radar.azimuth['data'] > sector['azmax'], :] = 0
        if sector['azmin'] > sector['azmax']:
            roi_flag[np.logical_and(
                radar.azimuth['data'] < sector['azmin'],
                radar.azimuth['data'] > sector['azmax']), :] = 0
    elif sector['azmin'] is not None:
        roi_flag[radar.azimuth['data'] < sector['azmin'], :] = 0
    elif sector['azmax'] is not None:
        roi_flag[radar.azimuth['data'] > sector['azmax'], :] = 0

    return roi_flag


def rainfall_accumulation(t_in_vec, val_in_vec, cum_time=3600.,
                          base_time=0., dropnan=False):
    """
    Computes the rainfall accumulation of a time series over a given period

    Parameters
    ----------
    t_in_vec : datetime array
        the input date and time array
    val_in_vec : float array
        the input values array [mm/h]
    cum_time : int
        accumulation time [s]
    base_time : int
        base time [s]
    dropnan : boolean
        if True remove NaN from the time series

    Returns
    -------
    t_out_vec : datetime array
        the output date and time array
    val_out_vec : float array
        the output values array
    np_vec : int array
        the number of samples at each period

    """
    # get the number of samples per interval
    t_out_vec, np_vec = time_series_statistics(
        t_in_vec, np.ones(len(val_in_vec), dtype=float), avg_time=cum_time,
        base_time=base_time, method='sum', dropnan=dropnan)

    np_vec[np.isnan(np_vec)] = 0
    np_vec = np_vec.astype(int)

    t_out_vec, val_out_vec = time_series_statistics(
        t_in_vec, val_in_vec, avg_time=cum_time, base_time=base_time,
        method='sum', dropnan=dropnan)

    t_sample = cum_time/np_vec  # find accumulation time of each sample
    val_out_vec *= (t_sample/3600.)  # conversion to mm in cum_time period

    val_out_vec = np.ma.asarray(val_out_vec)
    val_out_vec[np.isnan(val_out_vec)] = np.ma.masked

    return t_out_vec, val_out_vec, np_vec


def time_series_statistics(t_in_vec, val_in_vec, avg_time=3600,
                           base_time=1800, method='mean', dropnan=False):
    """
    Computes statistics over a time-averaged series. Only of package pandas is
    available otherwise returns None

    Parameters
    ----------
    t_in_vec : datetime array
        the input date and time array
    val_in_vec : float array
        the input values array
    avg_time : int
        averaging time [s]
    base_time : int
        base time [s]
    method : str
        statistical method
    dropnan : boolean
        if True remove NaN from the time series

    Returns
    -------
    t_out_vec : datetime array
        the output date and time array
    val_out_vec : float array
        the output values array

    """
    if not _PANDAS_AVAILABLE:
        warn('Pandas not available. Unable to compute time series statistics')
        return None, None

    df_in = pd.DataFrame(data=val_in_vec, index=pd.DatetimeIndex(t_in_vec))
    df_out = getattr(df_in.resample(
        str(avg_time)+'S', closed='right', label='right', base=base_time),
                     method)()
    if dropnan is True:
        df_out = df_out.dropna(how='any')
    t_out_vec = df_out.index.to_pydatetime()
    val_out_vec = df_out.values.flatten()

    return t_out_vec, val_out_vec


def find_contiguous_times(times, step=600):
    """
    Given and array of ordered times, find those contiguous according to
    a maximum time step

    Parameters
    ----------
    times : array of datetimes
        The array of times
    step : float
        The time step [s]

    Returns
    -------
    start_times, end_times : array of date times
        The start and end of each consecutive time period

    """
    run = []
    periods = []
    expect = None
    for time in times:
        if expect is None:
            run.append(time)
        elif time <= expect:
            run.append(time)
        else:
            run = [time]
            periods.append(run)
        expect = time+datetime.timedelta(seconds=step)

    if not periods:
        periods = [times]
    elif periods[0][0] != times[0]:
        periods.insert(0, [times[0]])

    print('number of consecutive periods: '+str(len(periods)))

    start_times = np.array([], dtype=datetime.datetime)
    end_times = np.array([], dtype=datetime.datetime)
    for period in periods:
        start_times = np.append(
            start_times, period[0]-datetime.timedelta(seconds=step))
        end_times = np.append(end_times, period[-1])

    return start_times, end_times


def join_time_series(t1, val1, t2, val2, dropnan=False):
    """
    joins time_series. Only of package pandas is available otherwise returns
    None.

    Parameters
    ----------
    t1 : datetime array
        time of first series
    val1 : float array
        value of first series
    t2 : datetime array
        time of second series
    val2 : float array
        value of second series
    dropnan : boolean
        if True remove NaN from the time series

    Returns
    -------
    t_out_vec : datetime array
        the resultant date time after joining the series
    val1_out_vec : float array
        value of first series
    val2_out_vec : float array
        value of second series

    """
    if not _PANDAS_AVAILABLE:
        warn('Pandas not available. Unable to join time series')
        return None, None, None

    df1 = pd.DataFrame(data=val1, index=pd.DatetimeIndex(t1))
    df2 = pd.DataFrame(data=val2, index=pd.DatetimeIndex(t2))
    df_out = pd.concat([df1, df2], join='outer', axis=1)
    if dropnan is True:
        df_out = df_out.dropna(how='any')
    t_out_vec = df_out.index.to_pydatetime()
    val1_out_vec = df_out.values[:, 0].flatten()
    val2_out_vec = df_out.values[:, 1].flatten()

    return t_out_vec, val1_out_vec, val2_out_vec


def get_range_bins_to_avg(rad1_rng, rad2_rng):
    """
    Compares the resolution of two radars and determines if and which radar
    has to be averaged and the length of the averaging window

    Parameters
    ----------
    rad1_rng : array
        the range of radar 1
    rad2_rng : datetime
        the range of radar 2

    Returns
    -------
    avg_rad1, avg_rad2 : Boolean
        Booleans specifying if the radar data has to be average in range
    avg_rad_lim : array with two elements
        the limits to the average (centered on each range gate)

    """
    rad1_res = rad1_rng[1]-rad1_rng[0]
    rad2_res = rad2_rng[1]-rad2_rng[0]
    res_ratio = rad1_res/rad2_res

    avg_rad1 = False
    avg_rad2 = False
    avg_rad_lim = None
    if res_ratio > 1.5:
        avg_rad2 = True
        nbins = int(res_ratio)
        if nbins % 2 == 0:
            avg_rad_lim = [-int(nbins/2)-1, int(nbins/2)]
        else:
            avg_rad_lim = [-int((nbins-1)/2), int((nbins-1)/2)]
    elif res_ratio < 1./1.5:
        avg_rad1 = True
        nbins = int(1./res_ratio)
        if nbins % 2 == 0:
            avg_rad_lim = [-int(nbins/2)-1, int(nbins/2)]
        else:
            avg_rad_lim = [-int((nbins-1)/2), int((nbins-1)/2)]

    return avg_rad1, avg_rad2, avg_rad_lim


def get_cercle_coords(x_centre, y_centre, radius=1000., resolution=16):
    """
    Get the points defining a cercle from the position of its centre and the
    radius

    Parameters
    ----------
    x_centre, y_centre : float
        The Cartesian coordinates of the center. Typically in m
    radius : float
        the cercle radius. Typically in m
    resolution : int
        The number of points used to define a quarter of the cercle

    Returns
    -------
    x_roi, y_roi : array of floats
        The position of the points defining the cercle

    """
    if not _SHAPELY_AVAILABLE:
        warn('shapely package not available. ' +
             'Unable to determine Region Of Interest')
        return None, None

    center_point = shapely.geometry.Point(x_centre, y_centre)
    cercle = center_point.buffer(radius, resolution=resolution)
    cercle_coords = list(cercle.exterior.coords)
    xy_roi = list(zip(*cercle_coords))
    x_roi = np.array(xy_roi[0])
    y_roi = np.array(xy_roi[1])

    return x_roi, y_roi


def belongs_roi_indices(lat, lon, roi, use_latlon=True):
    """
    Get the indices of points that belong to roi in a list of points

    Parameters
    ----------
    lat, lon : float arrays
        latitudes (or y) and longitudes(or x) to check
    roi : dict
        Dictionary describing the region of interest
    use_latlon : bool
        If True uses lat/lon as coordinates. Otherwise use Cartesian
        coordinates

    Returns
    -------
    inds : array of ints
        list of indices of points belonging to ROI
    is_roi : str
        Whether the list of points is within the region of interest.
        Can be 'All', 'None', 'Some'

    """
    if not _SHAPELY_AVAILABLE:
        warn('shapely package not available. ' +
             'Unable to determine if points belong to Region Of Interest')
        return np.asarray([]), 'None'

    lon_list = lon.flatten()
    lat_list = lat.flatten()

    if use_latlon:
        polygon = shapely.geometry.Polygon(list(zip(roi['lon'], roi['lat'])))
    else:
        polygon = shapely.geometry.Polygon(list(zip(roi['x'], roi['y'])))
    points = shapely.geometry.MultiPoint(list(zip(lon_list, lat_list)))

    inds = []
    if polygon.contains(points):
        warn('All points in the region of interest')
        is_roi = 'All'
        inds = np.indices(np.shape(lon))
    elif polygon.disjoint(points):
        warn('No points in the region of interest')
        is_roi = 'None'
    else:
        points_roi = points.intersection(polygon)
        if points_roi.geom_type == 'Point':
            ind = np.where(
                np.logical_and(lon == points_roi.x, lat == points_roi.y))
            if len(ind) == 1:
                ind = ind[0]
            inds.extend(ind)
        else:
            points_roi_list = list(points_roi)
            for point in points_roi_list:
                ind = np.where(np.logical_and(lon == point.x, lat == point.y))
                if len(ind) == 1:
                    ind = ind[0]
                inds.extend(ind)
        nroi = len(lat[inds])
        npoint = len(lat_list)
        warn(str(nroi)+' points out of '+str(npoint) +
             ' in the region of interest')
        is_roi = 'Some'

    return np.asarray(inds), is_roi


def find_ray_index(ele_vec, azi_vec, ele, azi, ele_tol=0., azi_tol=0.,
                   nearest='azi'):
    """
    Find the ray index corresponding to a particular elevation and azimuth

    Parameters
    ----------
    ele_vec, azi_vec : float arrays
        The elevation and azimuth data arrays where to look for
    ele, azi : floats
        The elevation and azimuth to search
    ele_tol, azi_tol : floats
        Tolerances [deg]
    nearest : str
        criteria to define wich ray to keep if multiple rays are within
        tolerance. azi: nearest azimuth, ele: nearest elevation

    Returns
    -------
    ind_ray : int
        The ray index

    """
    ind_ray = np.where(np.logical_and(
        np.logical_and(ele_vec <= ele+ele_tol, ele_vec >= ele-ele_tol),
        np.logical_and(azi_vec <= azi+azi_tol, azi_vec >= azi-azi_tol)))[0]

    if ind_ray.size == 0:
        return None
    if ind_ray.size == 1:
        return ind_ray[0]

    if nearest == 'azi':
        ind_min = np.argmin(np.abs(azi_vec[ind_ray]-azi))
    else:
        ind_min = np.argmin(np.abs(ele_vec[ind_ray]-ele))

    return ind_ray[ind_min]


def find_rng_index(rng_vec, rng, rng_tol=0.):
    """
    Find the range index corresponding to a particular range

    Parameters
    ----------
    rng_vec : float array
        The range data array where to look for
    rng : float
        The range to search
    rng_tol : float
        Tolerance [m]

    Returns
    -------
    ind_rng : int
        The range index

    """
    dist = np.abs(rng_vec-rng)
    ind_rng = np.argmin(dist)
    if dist[ind_rng] > rng_tol:
        return None

    return ind_rng


def find_ang_index(ang_vec, ang, ang_tol=0.):
    """
    Find the angle index corresponding to a particular fixed angle

    Parameters
    ----------
    ang_vec : float array
        The angle data array where to look for
    ang : float
        The angle to search
    ang_tol : float
        Tolerance [deg]

    Returns
    -------
    ind_ang : int
        The angle index

    """
    dist = np.abs(ang_vec-ang)
    ind_ang = np.argmin(dist)
    if dist[ind_ang] > ang_tol:
        return None

    return ind_ang


def find_nearest_gate(radar, lat, lon, latlon_tol=0.0005):
    """
    Find the radar gate closest to a lat,lon point

    Parameters
    ----------
    radar : radar object
        the radar object
    lat, lon : float
        The position of the point
    latlon_tol : float
        The tolerance around this point

    Returns
    -------
    ind_ray, ind_rng : int
        The ray and range index
    azi, rng : float
        the range and azimuth position of the gate

    """
    # find gates close to lat lon point
    inds_ray_aux, inds_rng_aux = np.where(np.logical_and(
        np.logical_and(
            radar.gate_latitude['data'] < lat+latlon_tol,
            radar.gate_latitude['data'] > lat-latlon_tol),
        np.logical_and(
            radar.gate_longitude['data'] < lon+latlon_tol,
            radar.gate_longitude['data'] > lon-latlon_tol)))

    if inds_ray_aux.size == 0:
        warn('No data found at point lat '+str(lat)+' +- ' +
             str(latlon_tol)+' lon '+str(lon)+' +- ' +
             str(latlon_tol)+' deg')

        return None, None, None, None

    # find closest latitude
    ind_min = np.argmin(np.abs(
        radar.gate_latitude['data'][inds_ray_aux, inds_rng_aux]-lat))
    ind_ray = inds_ray_aux[ind_min]
    ind_rng = inds_rng_aux[ind_min]

    azi = radar.azimuth['data'][ind_ray]
    rng = radar.range['data'][ind_rng]

    return ind_ray, ind_rng, azi, rng


def find_colocated_indexes(radar1, radar2, rad1_ele, rad1_azi, rad1_rng,
                           rad2_ele, rad2_azi, rad2_rng, ele_tol=0.5,
                           azi_tol=0.5, rng_tol=50.):
    """
    Given the theoretical elevation, azimuth and range of the co-located gates
    of two radars and a given tolerance returns the indices of the gates for
    the current radars

    Parameters
    ----------
    radar1, radar2 : radar objects
        the two radar objects
    rad1_ele, rad1_azi, rad1_rng : array of floats
        the radar coordinates of the radar1 gates
    rad2_ele, rad2_azi, rad2_rng : array of floats
        the radar coordinates of the radar2 gates
    ele_tol, azi_tol : floats
        azimuth and elevation angle tolerance [deg]
    rng_tol : float
        range Tolerance [m]

    Returns
    -------
    ind_ray_rad1, ind_rng_rad1, ind_ray_rad2, ind_rng_rad2 : array of ints
        the ray and range indexes of each radar gate

    """
    ngates = len(rad1_ele)
    ind_ray_rad1 = np.ma.masked_all(ngates, dtype=int)
    ind_rng_rad1 = np.ma.masked_all(ngates, dtype=int)
    ind_ray_rad2 = np.ma.masked_all(ngates, dtype=int)
    ind_rng_rad2 = np.ma.masked_all(ngates, dtype=int)
    for i in range(ngates):
        ind_ray_rad1_aux = find_ray_index(
            radar1.elevation['data'], radar1.azimuth['data'], rad1_ele[i],
            rad1_azi[i], ele_tol=ele_tol, azi_tol=azi_tol)
        if ind_ray_rad1_aux is None:
            continue
        ind_rng_rad1_aux = find_rng_index(
            radar1.range['data'], rad1_rng[i], rng_tol=rng_tol)
        if ind_rng_rad1_aux is None:
            continue

        ind_ray_rad2_aux = find_ray_index(
            radar2.elevation['data'], radar2.azimuth['data'], rad2_ele[i],
            rad2_azi[i], ele_tol=ele_tol, azi_tol=azi_tol)
        if ind_ray_rad2_aux is None:
            continue
        ind_rng_rad2_aux = find_rng_index(
            radar2.range['data'], rad2_rng[i], rng_tol=rng_tol)
        if ind_rng_rad2_aux is None:
            continue

        ind_ray_rad1[i] = ind_ray_rad1_aux
        ind_rng_rad1[i] = ind_rng_rad1_aux
        ind_ray_rad2[i] = ind_ray_rad2_aux
        ind_rng_rad2[i] = ind_rng_rad2_aux

    ind_ray_rad1 = ind_ray_rad1.compressed()
    ind_rng_rad1 = ind_rng_rad1.compressed()
    ind_ray_rad2 = ind_ray_rad2.compressed()
    ind_rng_rad2 = ind_rng_rad2.compressed()

    return ind_ray_rad1, ind_rng_rad1, ind_ray_rad2, ind_rng_rad2


def time_avg_range(timeinfo, avg_starttime, avg_endtime, period):
    """
    finds the new start and end time of an averaging

    Parameters
    ----------
    timeinfo : datetime
        the current volume time
    avg_starttime : datetime
        the current average start time
    avg_endtime: datetime
        the current average end time
    period: float
        the averaging period

    Returns
    -------
    new_starttime : datetime
        the new average start time
    new_endtime : datetime
        the new average end time

    """
    new_starttime = deepcopy(avg_starttime)
    new_endtime = deepcopy(avg_endtime)

    within_range = False
    while not within_range:
        if timeinfo > new_endtime:
            new_starttime += datetime.timedelta(seconds=period)
            new_endtime += datetime.timedelta(seconds=period)
        else:
            within_range = True

    return new_starttime, new_endtime


def get_closest_solar_flux(hit_datetime_list, flux_datetime_list,
                           flux_value_list):
    """
    finds the solar flux measurement closest to the sun hit

    Parameters
    ----------
    hit_datetime_list : datetime array
        the date and time of the sun hit
    flux_datetime_list : datetime array
        the date and time of the solar flux measurement
    flux_value_list: ndarray 1D
        the solar flux values

    Returns
    -------
    flux_datetime_closest_list : datetime array
        the date and time of the solar flux measurement closest to sun hit
    flux_value_closest_list : ndarray 1D
        the solar flux values closest to the sun hit time

    """
    flux_datetime_closest_list = list()
    flux_value_closest_list = np.ma.masked_all(len(hit_datetime_list))

    i = 0
    for hit_dt in hit_datetime_list:
        flux_datetime_closest = min(
            flux_datetime_list, key=lambda x: abs(x-hit_dt))
        flux_datetime_closest_list.append(flux_datetime_closest)

        # solar flux observation within 24h of sun hit
        time_diff = abs(flux_datetime_closest-hit_dt).total_seconds()
        if time_diff < 86400.:
            ind = flux_datetime_list.index(flux_datetime_closest)
            flux_value_closest_list[i] = flux_value_list[ind]
        else:
            warn('Nearest solar flux observation further than ' +
                 str(time_diff)+' s in time')
        i += 1

    return flux_datetime_closest_list, flux_value_closest_list


def get_fixed_rng_data(radar, field_names, fixed_rng, rng_tol=50.,
                       ele_min=None, ele_max=None, azi_min=None,
                       azi_max=None):
    """
    Creates a 2D-grid with (azi, ele) data at a fixed range

    Parameters
    ----------
    radar : radar object
        The radar object containing the data
    field_name : str
        The field name
    fixed_rng : float
        The fixed range [m]
    rng_tol : float
        The tolerance between the nominal range and the actual radar range [m]
    ele_min, ele_max, azi_min, azi_max : float or None
        The limits of the grid [deg]. If None the limits will be the limits
        of the radar volume

    Returns
    -------
    radar : radar object
        The radar object containing only the desired data

    """
    radar_aux = deepcopy(radar)

    ind_rng = find_rng_index(
        radar_aux.range['data'], fixed_rng, rng_tol=rng_tol)

    if ind_rng is None:
        warn('No range bin at range '+str(fixed_rng)+' with tolerance ' +
             str(rng_tol))
        return None, None, None

    # Determine angle limits
    if radar_aux.scan_type == 'ppi':
        if ele_min is None:
            ele_min = np.min(radar_aux.fixed_angle['data'])
        if ele_max is None:
            ele_max = np.max(radar_aux.fixed_angle['data'])
        if azi_min is None:
            azi_min = np.min(radar_aux.azimuth['data'])
        if azi_max is None:
            azi_max = np.max(radar_aux.azimuth['data'])
    else:
        if ele_min is None:
            ele_min = np.min(radar_aux.elevation['data'])
        if ele_max is None:
            ele_max = np.max(radar_aux.elevation['data'])
        if azi_min is None:
            azi_min = np.min(radar_aux.fixed_angle['data'])
        if azi_max is None:
            azi_max = np.max(radar_aux.fixed_angle['data'])

    if radar_aux.scan_type == 'ppi':
        # Get radar elevation angles within limits
        ele_vec = np.sort(radar_aux.fixed_angle['data'])
        ele_vec = ele_vec[
            np.logical_and(ele_vec >= ele_min, ele_vec <= ele_max)]
        if ele_vec is None:
            warn('No elevation angles between '+str(ele_min)+' and ' +
                 str(ele_max))
            return None, None, None

        # get sweeps corresponding to the desired elevation angles
        ind_sweeps = []
        for ele in ele_vec:
            ind_sweeps.append(
                np.where(radar_aux.fixed_angle['data'] == ele)[0][0])
        radar_aux = radar_aux.extract_sweeps(ind_sweeps)

        # Get indices of rays within limits
        if azi_min < azi_max:
            ind_rays = np.where(np.logical_and(
                radar_aux.azimuth['data'] >= azi_min,
                radar_aux.azimuth['data'] <= azi_max))[0]
        else:
            ind_rays = np.where(np.logical_or(
                radar_aux.azimuth['data'] >= azi_min,
                radar_aux.azimuth['data'] <= azi_max))[0]

    else:
        # Get radar azimuth angles within limits
        azi_vec = radar_aux.fixed_angle['data']
        if azi_min < azi_max:
            azi_vec = np.sort(azi_vec[
                np.logical_and(azi_vec >= azi_min, azi_vec <= azi_max)])
        else:
            azi_vec = azi_vec[
                np.logical_or(azi_vec >= azi_min, azi_vec <= azi_max)]
            azi_vec = np.append(
                np.sort(azi_vec[azi_vec >= azi_min]),
                np.sort(azi_vec[azi_vec < azi_min]))
        if azi_vec is None:
            warn('No azimuth angles between '+str(azi_min)+' and ' +
                 str(azi_max))
            return None, None, None

        # get sweeps corresponding to the desired azimuth angles
        ind_sweeps = []
        for azi in azi_vec:
            ind_sweeps.append(
                np.where(radar_aux.fixed_angle['data'] == azi)[0][0])
        radar_aux = radar_aux.extract_sweeps(ind_sweeps)

        # Get indices of rays within limits
        ind_rays = np.where(np.logical_and(
            radar_aux.elevation['data'] >= ele_min,
            radar_aux.elevation['data'] <= ele_max))[0]

    # get new sweep start index and stop index
    sweep_start_inds = deepcopy(radar_aux.sweep_start_ray_index['data'])
    sweep_end_inds = deepcopy(radar_aux.sweep_end_ray_index['data'])

    nrays = 0
    for j in range(radar_aux.nsweeps):
        # get azimuth indices for this elevation
        rays_in_sweep = np.size(
            ind_rays[np.logical_and(ind_rays >= sweep_start_inds[j],
                                    ind_rays <= sweep_end_inds[j])])
        radar_aux.rays_per_sweep['data'][j] = rays_in_sweep
        if j == 0:
            radar_aux.sweep_start_ray_index['data'][j] = 0
        else:
            radar_aux.sweep_start_ray_index['data'][j] = int(
                radar_aux.sweep_end_ray_index['data'][j-1]+1)
        radar_aux.sweep_end_ray_index['data'][j] = (
            radar_aux.sweep_start_ray_index['data'][j]+rays_in_sweep-1)
        nrays += rays_in_sweep

    # Get new fields
    for field_name in field_names:
        if field_name not in radar_aux.fields:
            warn('Field '+field_name+' not available')
            continue
        radar_aux.fields[field_name]['data'] = (
            radar_aux.fields[field_name]['data'][:, ind_rng])
        radar_aux.fields[field_name]['data'] = (
            radar_aux.fields[field_name]['data'][ind_rays, np.newaxis])

    # Update metadata
    radar_aux.time['data'] = radar_aux.time['data'][ind_rays]
    radar_aux.range['data'] = np.array([fixed_rng])
    radar_aux.azimuth['data'] = radar_aux.azimuth['data'][ind_rays]
    radar_aux.elevation['data'] = radar_aux.elevation['data'][ind_rays]
    radar_aux.init_gate_x_y_z()
    radar_aux.init_gate_longitude_latitude()
    radar_aux.init_gate_altitude()
    radar_aux.nrays = nrays
    radar_aux.ngates = 1

    return radar_aux


def create_sun_hits_field(rad_el, rad_az, sun_el, sun_az, data, imgcfg):
    """
    creates a sun hits field from the position and power of the sun hits

    Parameters
    ----------
    rad_el, rad_az, sun_el, sun_az : ndarray 1D
        azimuth and elevation of the radar and the sun respectively in degree
    data : masked ndarray 1D
        the sun hit data
    imgcfg: dict
        a dictionary specifying the ranges and resolution of the field to
        create

    Returns
    -------
    field : masked ndarray 2D
        the sun hit field

    """
    if data.compressed().size == 0:
        warn('No valid sun hits to plot.')
        return None

    azmin = imgcfg['azmin']
    azmax = imgcfg['azmax']
    elmin = imgcfg['elmin']
    elmax = imgcfg['elmax']
    azres = imgcfg['azres']
    elres = imgcfg['elres']

    mask = np.ma.getmaskarray(data)
    rad_el = rad_el[~mask]
    rad_az = rad_az[~mask]
    sun_el = sun_el[~mask]
    sun_az = sun_az[~mask]
    data = data[~mask]

    d_el = rad_el-sun_el
    d_az = (rad_az-sun_az)*np.cos(sun_el*np.pi/180.)

    npix_az = int((azmax-azmin)/azres)
    npix_el = int((elmax-elmin)/elres)
    field = np.ma.masked_all((npix_az, npix_el))

    ind_az = ((d_az+azmin)/azres).astype(int)
    ind_el = ((d_el+elmin)/elres).astype(int)

    field[ind_az, ind_el] = data

    return field


def create_sun_retrieval_field(par, field_name, imgcfg, lant=0.):
    """
    creates a sun retrieval field from the retrieval parameters

    Parameters
    ----------
    par : ndarray 1D
        the 5 retrieval parameters
    imgcfg: dict
        a dictionary specifying the ranges and resolution of the field to
        create

    Returns
    -------
    field : masked ndarray 2D
        the sun retrieval field

    """
    azmin = imgcfg['azmin']
    azmax = imgcfg['azmax']
    elmin = imgcfg['elmin']
    elmax = imgcfg['elmax']
    azres = imgcfg['azres']
    elres = imgcfg['elres']

    npix_az = int((azmax-azmin)/azres)
    npix_el = int((elmax-elmin)/elres)

    field = np.ma.masked_all((npix_az, npix_el))

    d_az = np.array(np.array(range(npix_az))*azres+azmin)
    d_el = np.array(np.array(range(npix_el))*elres+elmin)

    d_az_mat = np.broadcast_to(d_az.reshape(npix_az, 1), (npix_az, npix_el))
    d_el_mat = np.broadcast_to(d_el.reshape(1, npix_el), (npix_az, npix_el))

    field = (par[0]+par[1]*d_az_mat+par[2]*d_el_mat+par[3]*d_az_mat*d_az_mat +
             par[4]*d_el_mat*d_el_mat)
    if field_name in ('sun_est_power_h', 'sun_est_power_v'):
        # account for polarization of the antenna and scanning losses
        field += 3.+lant

    return field


def compute_quantiles(field, quantiles=None):
    """
    computes quantiles

    Parameters
    ----------
    field : ndarray 2D
        the radar field
    ray_start, ray_end : int
        starting and ending ray indexes
    quantiles: float array
        list of quantiles to compute

    Returns
    -------
    quantiles : float array
        list of quantiles
    values : float array
        values at each quantile

    """
    if quantiles is None:
        quantiles = [10., 20., 30., 40., 50., 60., 70., 80., 90., 95.]
        warn('No quantiles have been defined. Default ' + str(quantiles) +
             ' will be used')
    nquantiles = len(quantiles)
    values = np.ma.masked_all(nquantiles)

    data_valid = field.compressed()
    if np.size(data_valid) < 10:
        warn('Unable to compute quantiles. Not enough valid data')
        return quantiles, values

    for i in range(nquantiles):
        values[i] = np.percentile(data_valid, quantiles[i])

    return quantiles, values


def compute_quantiles_from_hist(bin_centers, hist, quantiles=None):
    """
    computes quantiles from histograms

    Parameters
    ----------
    bin_centers : ndarray 1D
        the bins
    hist : ndarray 1D
        the histogram
    quantiles: float array
        list of quantiles to compute

    Returns
    -------
    quantiles : float array
        list of quantiles
    values : float array
        values at each quantile

    """
    if quantiles is None:
        quantiles = [10., 20., 30., 40., 50., 60., 70., 80., 90.]
        warn('No quantiles have been defined. Default ' + str(quantiles) +
             ' will be used')
    nquantiles = len(quantiles)
    values = np.ma.masked_all(nquantiles)

    # check if all elements in histogram are masked values
    mask = np.ma.getmaskarray(hist)
    if mask.all():
        return quantiles, values

    np_t = np.ma.sum(hist)
    if np_t < 10:
        return quantiles, values

    freq = hist/np_t
    rel_freq = np.ma.cumsum(freq)

    percentiles = quantiles/100.
    for i in range(nquantiles):
        values[i] = bin_centers[rel_freq >= percentiles[i]][0]

    return quantiles, values


def compute_quantiles_sweep(field, ray_start, ray_end, quantiles=None):
    """
    computes quantiles of a particular sweep

    Parameters
    ----------
    field : ndarray 2D
        the radar field
    ray_start, ray_end : int
        starting and ending ray indexes
    quantiles: float array
        list of quantiles to compute

    Returns
    -------
    quantiles : float array
        list of quantiles
    values : float array
        values at each quantile

    """
    if quantiles is None:
        quantiles = [10., 20., 30., 40., 50., 60., 70., 80., 90.]
        warn('No quantiles have been defined. Default ' + str(quantiles) +
             ' will be used')
    nquantiles = len(quantiles)
    values = np.ma.masked_all(nquantiles)

    data_valid = field[ray_start:ray_end+1, :].compressed()
    if np.size(data_valid) < 10:
        warn('Unable to compute quantiles. Not enough valid data')
        return quantiles, values

    for i in range(nquantiles):
        values[i] = np.percentile(data_valid, quantiles[i])

    return quantiles, values


def compute_histogram(field, field_name, bin_edges=None, step=None,
                      vmin=None, vmax=None):
    """
    computes histogram of the data

    Parameters
    ----------
    field : ndarray 2D
        the radar field
    field_name: str or none
        name of the field
    bins_edges :ndarray 1D
        the bin edges
    step : float
        size of bin
    vmin, vmax : float
        The minimum and maximum value of the histogram

    Returns
    -------
    bin_edges : float array
        interval of each bin
    values : float array
        values at each bin

    """
    if bin_edges is None:
        if field_name is not None:
            bin_edges = get_histogram_bins(
                field_name, step=step, vmin=vmin, vmax=vmax)
        else:
            if vmin is None:
                vmin = np.ma.min(field)
            if vmax is None:
                vmax = np.ma.max(field)
            if step is None:
                step = (vmax-vmin)/100.
            bin_edges = np.arange(vmin, vmax+step, step)

    step_left = bin_edges[1]-bin_edges[0]
    step_right = bin_edges[-1]-bin_edges[-2]
    bin_center_left = bin_edges[0]+step_left/2.
    bin_center_right = bin_edges[-1]-step_right/2.
    values = np.ma.array(field).compressed()
    values[values < bin_center_left] = bin_center_left
    values[values > bin_center_right] = bin_center_right

    return bin_edges, values


def compute_histogram_sweep(field, ray_start, ray_end, field_name, step=None,
                            vmin=None, vmax=None):
    """
    computes histogram of the data in a particular sweep

    Parameters
    ----------
    field : ndarray 2D
        the radar field
    ray_start, ray_end : int
        starting and ending ray indexes
    field_name: str
        name of the field
    step : float
        size of bin
    vmin, vmax : float
        minimum and maximum values

    Returns
    -------
    bin_edges : float array
        interval of each bin
    values : float array
        values at each bin

    """
    bin_edges = get_histogram_bins(
        field_name, step=step, vmin=vmin, vmax=vmax)
    step_aux = bin_edges[1]-bin_edges[0]
    bin_centers = bin_edges[:-1]+step_aux/2.
    values = field[ray_start:ray_end+1, :].compressed()
    values[values < bin_centers[0]] = bin_centers[0]
    values[values > bin_centers[-1]] = bin_centers[-1]

    return bin_edges, values


def get_histogram_bins(field_name, step=None, vmin=None, vmax=None):
    """
    gets the histogram bins. If vmin or vmax are not define the range limits
    of the field as defined in the Py-ART config file are going to be used.

    Parameters
    ----------
    field_name: str
        name of the field
    step : float
        size of bin
    vmin, vmax : float
        The minimum and maximum value of the histogram

    Returns
    -------
    bin_edges : float array
        The bin edges

    """
    field_dict = pyart.config.get_metadata(field_name)
    if 'boundaries' in field_dict:
        return np.array(field_dict['boundaries'])

    vmin_aux, vmax_aux = pyart.config.get_field_limits(field_name)
    if vmin is None:
        vmin = vmin_aux
    if vmax is None:
        vmax = vmax_aux

    if step is None:
        step = (vmax-vmin)/50.
        warn('No step has been defined. Default '+str(step)+' will be used')

    return np.linspace(
        vmin-step/2., vmax+step/2., num=int((vmax-vmin)/step)+2)


def compute_2d_stats(field1, field2, field_name1, field_name2, step1=None,
                     step2=None):
    """
    computes a 2D histogram and statistics of the data

    Parameters
    ----------
    field1, field2 : ndarray 2D
        the two fields
    field_name1, field_nam2: str
        the name of the fields
    step1, step2 : float
        size of bin

    Returns
    -------
    hist_2d : array
        the histogram
    bin_edges1, bin_edges2 : float array
        The bin edges
    stats : dict
        a dictionary with statistics

    """
    if field1.size == 0 or field2.size == 0:
        warn('Unable to compute 2D histogram. Empty field')
        stats = {
            'npoints': 0,
            'meanbias': np.ma.asarray(np.ma.masked),
            'medianbias': np.ma.asarray(np.ma.masked),
            'quant25bias': np.ma.asarray(np.ma.masked),
            'quant75bias': np.ma.asarray(np.ma.masked),
            'modebias': np.ma.asarray(np.ma.masked),
            'corr': np.ma.asarray(np.ma.masked),
            'slope': np.ma.asarray(np.ma.masked),
            'intercep': np.ma.asarray(np.ma.masked),
            'intercep_slope_1': np.ma.asarray(np.ma.masked)
        }
        return None, None, None, stats

    hist_2d, bin_edges1, bin_edges2 = compute_2d_hist(
        field1, field2, field_name1, field_name2, step1=step1, step2=step2)
    step_aux1 = bin_edges1[1]-bin_edges1[0]
    bin_centers1 = bin_edges1[:-1]+step_aux1/2.
    step_aux2 = bin_edges2[1]-bin_edges2[0]
    bin_centers2 = bin_edges2[:-1]+step_aux2/2.
    npoints = len(field1)
    meanbias = 10.*np.ma.log10(
        np.ma.mean(np.ma.power(10., 0.1*field2)) /
        np.ma.mean(np.ma.power(10., 0.1*field1)))
    medianbias = np.ma.median(field2-field1)
    quant25bias = np.percentile((field2-field1).compressed(), 25.)
    quant75bias = np.percentile((field2-field1).compressed(), 75.)
    ind_max_val1, ind_max_val2 = np.where(hist_2d == np.ma.amax(hist_2d))
    modebias = bin_centers2[ind_max_val2[0]]-bin_centers1[ind_max_val1[0]]
    slope, intercep, corr, _, _ = scipy.stats.linregress(
        field1, y=field2)
    intercep_slope_1 = np.ma.mean(field2-field1)

    stats = {
        'npoints': npoints,
        'meanbias': np.ma.asarray(meanbias),
        'medianbias': np.ma.asarray(medianbias),
        'quant25bias': np.ma.asarray(quant25bias),
        'quant75bias': np.ma.asarray(quant75bias),
        'modebias': np.ma.asarray(modebias),
        'corr': np.ma.asarray(corr),
        'slope': np.ma.asarray(slope),
        'intercep': np.ma.asarray(intercep),
        'intercep_slope_1': np.ma.asarray(intercep_slope_1)
    }

    return hist_2d, bin_edges1, bin_edges2, stats


def compute_1d_stats(field1, field2):
    """
    returns statistics of data

    Parameters
    ----------
    field1, field2 : ndarray 1D
        the two fields to compare

    Returns
    -------
    stats : dict
        a dictionary with statistics

    """
    if field1.size == 0 or field2.size == 0:
        warn('Unable to compute statistics. Empty fields')
        stats = {
            'npoints': 0,
            'NB': np.ma.asarray(np.ma.masked),
            'corr': np.ma.asarray(np.ma.masked),
            'RMS': np.ma.asarray(np.ma.masked),
            'Nash': np.ma.asarray(np.ma.masked)
        }
        return stats

    npoints = len(field1)
    mean1 = np.ma.mean(field1)
    mean2 = np.ma.mean(field2)
    nb = mean2/mean1-1
    _, _, corr, _, _ = scipy.stats.linregress(field1, y=field2)
    rms = np.ma.sqrt(np.ma.sum(np.ma.power(field2-field1, 2.))/npoints)
    nash = (1.-np.ma.sum(np.ma.power(field2-field1, 2.)) /
            np.ma.sum(np.ma.power(field1-mean1, 2.)))

    stats = {
        'npoints': npoints,
        'NB': nb,
        'corr': corr,
        'RMS': rms,
        'Nash': nash
    }

    return stats


def compute_2d_hist(field1, field2, field_name1, field_name2, step1=None,
                    step2=None):
    """
    computes a 2D histogram of the data

    Parameters
    ----------
    field1, field2 : ndarray 2D
        the radar fields
    field_name1, field_name2 : str
        field names
    step1, step2 : float
        size of the bins

    Returns
    -------
    H : float array 2D
        The bi-dimensional histogram of samples x and y
    xedges, yedges : float array
        the bin edges along each dimension

    """
    bin_edges1 = get_histogram_bins(field_name1, step=step1)
    step_aux1 = bin_edges1[1]-bin_edges1[0]
    bin_centers1 = bin_edges1[:-1]+step_aux1/2.
    bin_edges2 = get_histogram_bins(field_name2, step=step2)
    step_aux2 = bin_edges2[1]-bin_edges2[0]
    bin_centers2 = bin_edges2[:-1]+step_aux2/2.

    field1[field1 < bin_centers1[0]] = bin_centers1[0]
    field1[field1 > bin_centers1[-1]] = bin_centers1[-1]

    field2[field2 < bin_centers2[0]] = bin_centers2[0]
    field2[field2 > bin_centers2[-1]] = bin_centers2[-1]

    fill_value = pyart.config.get_fillvalue()
    return np.histogram2d(
        field1.filled(fill_value=fill_value),
        field2.filled(fill_value=fill_value), bins=[bin_edges1, bin_edges2])


def quantize_field(field, field_name, step, vmin=None, vmax=None):
    """
    quantizes data

    Parameters
    ----------
    field : ndarray 2D
        the radar field
    field_name: str
        name of the field
    step : float
        size of bin
    vmin, vmax : float
        min and max values

    Returns
    -------
    fieldq : ndarray 2D
        The quantized field
    values : float array
        values at each bin

    """
    vmin_aux, vmax_aux = pyart.config.get_field_limits(field_name)
    if vmin is None:
        vmin = vmin_aux
    if vmax is None:
        vmax = vmax_aux
    field[field < vmin] = vmin
    field[field > vmax] = vmax
    fieldq = ((field+vmin)/step+1).astype(int)

    return fieldq.filled(fill_value=0)


def compute_profile_stats(field, gate_altitude, h_vec, h_res,
                          quantity='quantiles',
                          quantiles=np.array([0.25, 0.50, 0.75]),
                          nvalid_min=4, std_field=None, np_field=None,
                          make_linear=False, include_nans=False):
    """
    Compute statistics of vertical profile

    Parameters
    ----------
    field : ndarray
        the radar field
    gate_altitude: ndarray
        the altitude at each radar gate [m MSL]
    h_vec : 1D ndarray
        height vector [m MSL]
    h_res : float
        heigh resolution [m]
    quantity : str
        The quantity to compute. Can be
        ['quantiles', 'mode', 'regression_mean', 'mean'].
        If 'mean', the min, max, and average is computed.
    quantiles : 1D ndarray
        the quantiles to compute
    nvalid_min : int
        the minimum number of points to consider the stats valid
    std_field : ndarray
        the standard deviation of the regression at each range gate
    np_field : ndarray
        the number of points used to compute the regression at each range gate
    make_linear : Boolean
        If true the data is transformed into linear coordinates before taking
        the mean
    include_nans : Boolean
        If true NaN will be considered as zeros

    Returns
    -------
    vals : ndarray 2D
        The resultant statistics
    val_valid : ndarray 1D
        The number of points to compute the stats used at each height level

    """
    nh = h_vec.size

    if quantity == 'mean':
        vals = np.ma.empty((nh, 3), dtype=float)
        vals[:] = np.ma.masked
    elif quantity == 'mode':
        vals = np.ma.empty((nh, 6), dtype=float)
        vals[:, 0] = np.ma.masked
        vals[:, 2] = np.ma.masked
        vals[:, 4] = np.ma.masked
        vals[:, 1] = 0
        vals[:, 3] = 0
        vals[:, 5] = 0
    elif quantity == 'regression_mean':
        vals = np.ma.masked_all((nh, 2), dtype=float)
    else:
        vals = np.ma.masked_all((nh, quantiles.size), dtype=float)

    val_valid = np.zeros(nh, dtype=int)
    for i, h in enumerate(h_vec):
        data = field[np.logical_and(
            gate_altitude >= h-h_res/2., gate_altitude < h+h_res/2.)]
        if include_nans:
            data[np.ma.getmaskarray(data)] = 0.
        if data.size == 0:
            continue
        if quantity == 'mean':
            mask = np.ma.getmaskarray(data)
            nvalid = np.count_nonzero(np.logical_not(mask))
            if nvalid >= nvalid_min:
                if make_linear:
                    data = np.ma.power(10., 0.1*data)
                    vals[i, 0] = 10.*np.ma.log10(np.ma.mean(data))
                    vals[i, 1] = 10.*np.ma.log10(np.ma.min(data))
                    vals[i, 2] = 10.*np.ma.log10(np.ma.max(data))
                else:
                    vals[i, 0] = np.ma.mean(data)
                    vals[i, 1] = np.ma.min(data)
                    vals[i, 2] = np.ma.max(data)
                val_valid[i] = nvalid

        elif quantity == 'mode':
            mask = np.ma.getmaskarray(data)
            nvalid = np.count_nonzero(np.logical_not(mask))
            if nvalid >= nvalid_min:
                val_valid[i] = nvalid

                # get mode
                data = data.compressed()
                if data.size == 0:
                    continue
                mode, count = scipy.stats.mode(
                    data, axis=None, nan_policy='omit')
                vals[i, 0] = mode
                vals[i, 1] = count/nvalid*100.

                # get second most common
                data = np.ma.masked_where(data == mode, data).compressed()
                if data.size == 0:
                    continue
                mode, count = scipy.stats.mode(
                    data, axis=None, nan_policy='omit')
                vals[i, 2] = mode
                vals[i, 3] = count/nvalid*100.

                # get third most common
                data = np.ma.masked_where(data == mode, data).compressed()
                if data.size == 0:
                    continue
                mode, count = scipy.stats.mode(
                    data, axis=None, nan_policy='omit')
                vals[i, 4] = mode
                vals[i, 5] = count/nvalid*100.
        elif quantity == 'regression_mean':
            if std_field is None or np_field is None:
                warn('Unable to compute regression mean')
                return None, None
            data_std = std_field[np.logical_and(
                gate_altitude >= h-h_res/2., gate_altitude < h+h_res/2.)]
            data_np = np_field[np.logical_and(
                gate_altitude >= h-h_res/2., gate_altitude < h+h_res/2.)]

            val_valid[i] = np.sum(data_np)
            if val_valid[i] == 0.:
                continue

            data_var = np.ma.power(data_std, 2.)
            weights = (data_np-1)/(data_var+0.01)
            vals[i, 0] = np.ma.sum(weights*data)/np.ma.sum(weights)
            vals[i, 1] = np.ma.sqrt(
                np.ma.sum((data_np-1)*data_var)/np.ma.sum(data_np-1))
        else:
            _, quants, nvalid = quantiles_weighted(
                data, quantiles=quantiles)
            if nvalid is not None:
                if nvalid >= nvalid_min:
                    vals[i, :] = quants
                    val_valid[i] = nvalid

    return vals, val_valid


def project_to_vertical(data_in, data_height, grid_height, interp_kind='none',
                        fill_value=-9999.):
    """
    Projects radar data to a regular vertical grid

    Parameters
    ----------
    data_in : ndarray 1D
        the radar data to project
    data_height : ndarray 1D
        the height of each radar point
    grid_height : ndarray 1D
        the regular vertical grid to project to
    interp_kind : str
        The type of interpolation to use: 'none' or 'nearest'
    fill_value : float
        The fill value used for interpolation

    Returns
    -------
    data_out : ndarray 1D
        The projected data

    """
    if data_in.size == 0:
        data_out = np.ma.masked_all(grid_height.size)
        return data_out

    if interp_kind == 'none':
        hres = grid_height[1]-grid_height[0]
        data_out = np.ma.masked_all(grid_height.size)
        for ind_r, h in enumerate(grid_height):
            ind_h = find_rng_index(data_height, h, rng_tol=hres/2.)
            if ind_h is None:
                continue
            data_out[ind_r] = data_in[ind_h]
    elif interp_kind == 'nearest':
        data_filled = data_in.filled(fill_value=fill_value)
        f = scipy.interpolate.interp1d(
            data_height, data_filled, kind=interp_kind, bounds_error=False,
            fill_value=fill_value)
        data_out = np.ma.masked_values(f(grid_height), fill_value)
    else:
        valid = np.logical_not(np.ma.getmaskarray(data_in))
        height_valid = data_height[valid]
        data_valid = data_in[valid]
        f = scipy.interpolate.interp1d(
            height_valid, data_valid, kind=interp_kind, bounds_error=False,
            fill_value=fill_value)
        data_out = np.ma.masked_values(f(grid_height), fill_value)

    return data_out
