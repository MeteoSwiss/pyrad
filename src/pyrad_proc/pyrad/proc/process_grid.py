"""
pyrad.proc.process_grid
=======================

Functions to processes gridded data.

.. autosummary::
    :toctree: generated/

    process_raw_grid
    process_grid
    process_grid_point
    process_grid_time_stats
    process_grid_time_stats2
    process_grid_fields_diff
    process_grid_texture
    process_grid_mask
    process_normalize_luminosity
    process_pixel_filter

"""

from copy import deepcopy
from warnings import warn
import datetime
import numpy as np
import scipy
from netCDF4 import num2date

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..util.radar_utils import time_avg_range


def process_raw_grid(procstatus, dscfg, radar_list=None):
    """
    Dummy function that returns the initial input data set

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, _, _, _, _ = get_datatype_fields(datatypedescr)
        break
    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    new_dataset = {'radar_out': deepcopy(radar_list[ind_rad])}

    return new_dataset, ind_rad


def process_grid(procstatus, dscfg, radar_list=None):
    """
    Puts the radar data in a regular grid

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement
        gridconfig : dictionary. Dataset keyword
            Dictionary containing some or all of this keywords:
            xmin, xmax, ymin, ymax, zmin, zmax : floats
                minimum and maximum horizontal distance from grid origin [km]
                and minimum and maximum vertical distance from grid origin [m]
                Defaults -40, 40, -40, 40, 0., 10000.
            hres, vres : floats
                horizontal and vertical grid resolution [m]
                Defaults 1000., 500.
            latorig, lonorig, altorig : floats
                latitude and longitude of grid origin [deg] and altitude of
                grid origin [m MSL]
                Defaults the latitude, longitude and altitude of the radar
        wfunc : str. Dataset keyword
            the weighting function used to combine the radar gates close to a
            grid point. Possible values BARNES, BARNES2, CRESSMAN, NEAREST
            Default NEAREST
        roif_func : str. Dataset keyword
            the function used to compute the region of interest.
            Possible values: dist_beam, constant
        roi : float. Dataset keyword
            the (minimum) radius of the region of interest in m. Default half
            the largest resolution
        beamwidth : float. Dataset keyword
            the radar antenna beamwidth [deg]. If None that of the key
            radar_beam_width_h in attribute instrument_parameters of the radar
            object will be used. If the key or the attribute are not present
            a default 1 deg value will be used
        beam_spacing : float. Dataset keyword
            the beam spacing, i.e. the ray angle resolution [deg]. If None,
            that of the attribute ray_angle_res of the radar object will be
            used. If the attribute is None a default 1 deg value will be used

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the gridded data
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    field_names_aux = []
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_names_aux.append(get_fieldname_pyart(datatype))

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    # keep only fields present in radar object
    field_names = []
    nfields_available = 0
    for field_name in field_names_aux:
        if field_name not in radar.fields:
            warn('Field name '+field_name+' not available in radar object')
            continue
        field_names.append(field_name)
        nfields_available += 1

    if nfields_available == 0:
        warn("Fields not available in radar data")
        return None, None

    # default parameters
    xmin = -40.
    xmax = 40.
    ymin = -40.
    ymax = 40.
    zmin = 0.
    zmax = 10000.
    hres = 1000.
    vres = 500.
    lat = float(radar.latitude['data'])
    lon = float(radar.longitude['data'])
    alt = float(radar.altitude['data'])

    if 'gridConfig' in dscfg:
        if 'xmin' in dscfg['gridConfig']:
            xmin = dscfg['gridConfig']['xmin']
        if 'xmax' in dscfg['gridConfig']:
            xmax = dscfg['gridConfig']['xmax']
        if 'ymin' in dscfg['gridConfig']:
            ymin = dscfg['gridConfig']['ymin']
        if 'ymax' in dscfg['gridConfig']:
            ymax = dscfg['gridConfig']['ymax']
        if 'zmin' in dscfg['gridConfig']:
            zmin = dscfg['gridConfig']['zmin']
        if 'zmax' in dscfg['gridConfig']:
            zmax = dscfg['gridConfig']['zmax']
        if 'hres' in dscfg['gridConfig']:
            hres = dscfg['gridConfig']['hres']
        if 'vres' in dscfg['gridConfig']:
            vres = dscfg['gridConfig']['vres']
        if 'latorig' in dscfg['gridConfig']:
            lat = dscfg['gridConfig']['latorig']
        if 'lonorig' in dscfg['gridConfig']:
            lon = dscfg['gridConfig']['lonorig']
        if 'altorig' in dscfg['gridConfig']:
            alt = dscfg['gridConfig']['altorig']

    wfunc = dscfg.get('wfunc', 'NEAREST')
    roi_func = dscfg.get('roi_func', 'dist_beam')

    # number of grid points in cappi
    nz = int((zmax-zmin)/vres)+1
    ny = int((ymax-ymin)*1000./hres)+1
    nx = int((xmax-xmin)*1000./hres)+1

    min_radius = dscfg.get('roi', np.max([vres, hres])/2.)
    # parameters to determine the gates to use for each grid point
    beamwidth = dscfg.get('beamwidth', None)
    beam_spacing = dscfg.get('beam_spacing', None)

    if beamwidth is None:
        if (radar.instrument_parameters is not None and
                'radar_beam_width_h' in radar.instrument_parameters):
            beamwidth = radar.instrument_parameters[
                'radar_beam_width_h']['data'][0]
        else:
            warn('Unknown radar beamwidth. Default 1 deg will be used')
            beamwidth = 1

    if beam_spacing is None:
        if radar.ray_angle_res is not None:
            beam_spacing = radar.ray_angle_res['data'][0]
        else:
            warn('Unknown beam spacing. Default 1 deg will be used')
            beam_spacing = 1

    # cartesian mapping
    grid = pyart.map.grid_from_radars(
        (radar,), gridding_algo='map_to_grid',
        weighting_function=wfunc,
        roi_func=roi_func, h_factor=1.0, nb=beamwidth, bsp=beam_spacing,
        min_radius=min_radius, constant_roi=min_radius,
        grid_shape=(nz, ny, nx),
        grid_limits=((zmin, zmax), (ymin*1000., ymax*1000.),
                     (xmin*1000., xmax*1000.)),
        grid_origin=(lat, lon), grid_origin_alt=alt,
        fields=field_names)

    new_dataset = {'radar_out': grid}

    return new_dataset, ind_rad


def process_grid_point(procstatus, dscfg, radar_list=None):
    """
    Obtains the grid data at a point location.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement
        latlon : boolean. Dataset keyword
            if True position is obtained from latitude, longitude information,
            otherwise position is obtained from grid index (iz, iy, ix).
        lon : float. Dataset keyword
            the longitude [deg]. Use when latlon is True.
        lat : float. Dataset keyword
            the latitude [deg]. Use when latlon is True.
        alt : float. Dataset keyword
            altitude [m MSL]. Use when latlon is True.
        iz, iy, ix : int. Dataset keyword
            The grid indices. Use when latlon is False
        latlonTol : float. Dataset keyword
            latitude-longitude tolerance to determine which grid point to use
            [deg]
        altTol : float. Dataset keyword
            Altitude tolerance to determine which grid point to use [deg]

    radar_list : list of Radar objects
          Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the data and metadata at the point of interest
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        break
    field_name = get_fieldname_pyart(datatype)
    ind_rad = int(radarnr[5:8])-1

    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None, None

        # prepare for exit
        new_dataset = {
            'time': dscfg['global_data']['time'],
            'ref_time': dscfg['global_data']['ref_time'],
            'datatype': datatype,
            'point_coordinates_WGS84_lon_lat_alt': (
                dscfg['global_data']['point_coordinates_WGS84_lon_lat_alt']),
            'grid_points_iz_iy_ix': (
                dscfg['global_data']['grid_points_iz_iy_ix']),
            'final': True}

        return new_dataset, ind_rad

    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    grid = radar_list[ind_rad]

    if field_name not in grid.fields:
        warn('Unable to extract point measurement information. ' +
             'Field not available')
        return None, None

    if dscfg['latlon']:
        lon = dscfg['lon']
        lat = dscfg['lat']
        alt = dscfg.get('alt', 0.)
        latlon_tol = dscfg.get('latlonTol', 1.)
        alt_tol = dscfg.get('altTol', 100.)

        d_lon = np.min(np.abs(grid.point_longitude['data']-lon))
        if d_lon > latlon_tol:
            warn(' No grid point found for point (lat, lon, alt):(' +
                 str(lat)+', '+str(lon)+', '+str(alt) +
                 '). Minimum distance to longitude '+str(d_lon) +
                 ' larger than tolerance')
            return None, None
        d_lat = np.min(np.abs(grid.point_latitude['data']-lat))
        if d_lat > latlon_tol:
            warn(' No grid point found for point (lat, lon, alt):(' +
                 str(lat)+', '+str(lon)+', '+str(alt) +
                 '). Minimum distance to latitude '+str(d_lat) +
                 ' larger than tolerance')
            return None, None
        d_alt = np.min(np.abs(grid.point_altitude['data']-alt))
        if d_alt > alt_tol:
            warn(' No grid point found for point (lat, lon, alt):(' +
                 str(lat)+', '+str(lon)+', '+str(alt) +
                 '). Minimum distance to altitude '+str(d_alt) +
                 ' larger than tolerance')
            return None, None

        iz, iy, ix = np.unravel_index(
            np.argmin(
                np.abs(grid.point_longitude['data']-lon) +
                np.abs(grid.point_latitude['data']-lat)),
            grid.point_longitude['data'].shape)

        iz = np.argmin(np.abs(grid.point_altitude['data'][:, iy, ix]-alt))
    else:
        ix = dscfg['ix']
        iy = dscfg['iy']
        iz = dscfg['iz']

        lon = grid.point_longitude['data'][iz, iy, ix]
        lat = grid.point_latitude['data'][iz, iy, ix]
        alt = grid.point_altitude['data'][iz, iy, ix]

    val = grid.fields[field_name]['data'][iz, iy, ix]
    time = num2date(
        grid.time['data'][0], grid.time['units'], grid.time['calendar'])

    # initialize dataset
    if dscfg['initialized'] == 0:
        poi = {
            'point_coordinates_WGS84_lon_lat_alt': [lon, lat, alt],
            'grid_points_iz_iy_ix': [iz, iy, ix],
            'time': time,
            'ref_time': dscfg['timeinfo']}
        dscfg['global_data'] = poi
        dscfg['initialized'] = 1

    dscfg['global_data']['ref_time'] = dscfg['timeinfo']

    # prepare for exit
    new_dataset = dict()
    new_dataset.update({'value': val})
    new_dataset.update({'datatype': datatype})
    new_dataset.update({'time': time})
    new_dataset.update({'ref_time': dscfg['timeinfo']})
    new_dataset.update(
        {'point_coordinates_WGS84_lon_lat_alt': [lon, lat, alt]})
    new_dataset.update({'grid_points_iz_iy_ix': [iz, iy, ix]})
    new_dataset.update(
        {'used_coordinates_WGS84_lon_lat_alt': [
            grid.point_longitude['data'][iz, iy, ix],
            grid.point_latitude['data'][iz, iy, ix],
            grid.point_altitude['data'][iz, iy, ix]]})
    new_dataset.update({'final': False})

    return new_dataset, ind_rad


def process_grid_time_stats(procstatus, dscfg, radar_list=None):
    """
    computes the temporal statistics of a field

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        period : float. Dataset keyword
            the period to average [s]. If -1 the statistics are going to be
            performed over the entire data. Default 3600.
        start_average : float. Dataset keyword
            when to start the average [s from midnight UTC]. Default 0.
        lin_trans: int. Dataset keyword
            If 1 apply linear transformation before averaging
        use_nan : bool. Dataset keyword
            If true non valid data will be used
        nan_value : float. Dataset keyword
            The value of the non valid data. Default 0
        stat: string. Dataset keyword
            Statistic to compute: Can be mean, std, cov, min, max. Default
            mean
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_name = get_fieldname_pyart(datatype)
        break
    ind_rad = int(radarnr[5:8])-1

    start_average = dscfg.get('start_average', 0.)
    period = dscfg.get('period', 3600.)
    lin_trans = dscfg.get('lin_trans', 0)
    use_nan = dscfg.get('use_nan', 0)
    nan_value = dscfg.get('nan_value', 0.)
    stat = dscfg.get('stat', 'mean')

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        grid = radar_list[ind_rad]

        if field_name not in grid.fields:
            warn(field_name+' not available.')
            return None, None

        # Prepare auxiliary grid
        field_dict = deepcopy(grid.fields[field_name])
        if stat in ('mean', 'std', 'cov'):
            if lin_trans:
                field_dict['data'] = np.ma.power(10., 0.1*field_dict['data'])

            if use_nan:
                field_dict['data'] = np.ma.asarray(
                    field_dict['data'].filled(nan_value))

            if stat in ('std', 'cov'):
                sum2_dict = pyart.config.get_metadata('sum_squared')
                sum2_dict['data'] = field_dict['data']*field_dict['data']
        else:
            if use_nan:
                field_dict['data'] = np.ma.asarray(
                    field_dict['data'].filled(nan_value))

        npoints_dict = pyart.config.get_metadata('number_of_samples')
        npoints_dict['data'] = np.ma.asarray(
            np.logical_not(np.ma.getmaskarray(field_dict['data'])), dtype=int)

        grid_aux = deepcopy(grid)
        grid_aux.fields = dict()
        grid_aux.add_field(field_name, field_dict)
        grid_aux.add_field('number_of_samples', npoints_dict)

        if stat in ('std', 'cov'):
            grid_aux.add_field('sum_squared', sum2_dict)

        # first volume: initialize start and end time of averaging
        if dscfg['initialized'] == 0:
            avg_par = dict()
            if period != -1:
                date_00 = dscfg['timeinfo'].replace(
                    hour=0, minute=0, second=0, microsecond=0)

                avg_par.update(
                    {'starttime': date_00+datetime.timedelta(
                        seconds=start_average)})
                avg_par.update(
                    {'endtime': avg_par['starttime']+datetime.timedelta(
                        seconds=period)})
            else:
                avg_par.update({'starttime': dscfg['timeinfo']})
                avg_par.update({'endtime': dscfg['timeinfo']})

            avg_par.update({'timeinfo': dscfg['timeinfo']})
            dscfg['global_data'] = avg_par
            dscfg['initialized'] = 1

        if dscfg['initialized'] == 0:
            return None, None

        dscfg['global_data']['timeinfo'] = dscfg['timeinfo']
        # no grid object in global data: create it
        if 'grid_out' not in dscfg['global_data']:
            if period != -1:
                # get start and stop times of new grid object
                (dscfg['global_data']['starttime'],
                 dscfg['global_data']['endtime']) = (
                     time_avg_range(
                         dscfg['timeinfo'], dscfg['global_data']['starttime'],
                         dscfg['global_data']['endtime'], period))

                # check if volume time older than starttime
                if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
                    dscfg['global_data'].update({'grid_out': grid_aux})
            else:
                dscfg['global_data'].update({'grid_out': grid_aux})

            return None, None

        # still accumulating: add field_dict to global field_dict
        if (period == -1 or
                dscfg['timeinfo'] < dscfg['global_data']['endtime']):

            if period == -1:
                dscfg['global_data']['endtime'] = dscfg['timeinfo']

            dscfg['global_data']['grid_out'].fields['number_of_samples'][
                'data'] += npoints_dict['data']

            if stat in ('mean', 'std', 'cov'):
                masked_sum = np.ma.getmaskarray(
                    dscfg['global_data']['grid_out'].fields[
                        field_name]['data'])
                valid_sum = np.logical_and(
                    np.logical_not(masked_sum),
                    np.logical_not(np.ma.getmaskarray(field_dict['data'])))

                dscfg['global_data']['grid_out'].fields[
                    field_name]['data'][masked_sum] = (
                        field_dict['data'][masked_sum])

                dscfg['global_data']['grid_out'].fields[
                    field_name]['data'][valid_sum] += (
                        field_dict['data'][valid_sum])

                if stat in ('cov', 'std'):
                    dscfg['global_data']['grid_out'].fields[
                        'sum_squared']['data'][masked_sum] = (
                            field_dict['data'][masked_sum] *
                            field_dict['data'][masked_sum])

                    dscfg['global_data']['grid_out'].fields[
                        'sum_squared']['data'][valid_sum] += (
                            field_dict['data'][valid_sum] *
                            field_dict['data'][valid_sum])

            elif stat == 'max':
                dscfg['global_data']['grid_out'].fields[
                    field_name]['data'] = np.maximum(
                        dscfg['global_data']['grid_out'].fields[
                            field_name]['data'].filled(fill_value=-1.e300),
                        field_dict['data'].filled(fill_value=-1.e300))

                dscfg['global_data']['grid_out'].fields[
                    field_name]['data'] = np.ma.masked_values(
                        dscfg['global_data']['grid_out'].fields[
                            field_name]['data'], -1.e300)
            elif stat == 'min':
                dscfg['global_data']['grid_out'].fields[
                    field_name]['data'] = np.minimum(
                        dscfg['global_data']['grid_out'].fields[
                            field_name]['data'].filled(fill_value=1.e300),
                        field_dict['data'].filled(fill_value=1.e300))

                dscfg['global_data']['grid_out'].fields[
                    field_name]['data'] = np.ma.masked_values(
                        dscfg['global_data']['grid_out'].fields[
                            field_name]['data'], 1.e300)

            return None, None

        # we have reached the end of the accumulation period: do the averaging
        # and start a new object (only reachable if period != -1)
        if stat in ('mean', 'std', 'cov'):
            field_mean = (
                dscfg['global_data']['grid_out'].fields[field_name]['data'] /
                dscfg['global_data']['grid_out'].fields[
                    'number_of_samples']['data'])

            if stat == 'mean':
                if lin_trans:
                    dscfg['global_data']['grid_out'].fields[
                        field_name]['data'] = 10.*np.ma.log10(field_mean)
                else:
                    dscfg['global_data']['grid_out'].fields[
                        field_name]['data'] = field_mean
            elif stat in ('std', 'cov'):
                field_std = np.ma.sqrt(
                    dscfg['global_data']['grid_out'].fields[
                        'sum_squared']['data'] /
                    dscfg['global_data']['grid_out'].fields[
                        'number_of_samples']['data']-field_mean*field_mean)

                if stat == 'std':
                    if lin_trans:
                        dscfg['global_data']['grid_out'].fields[
                            field_name]['data'] = 10.*np.ma.log10(field_std)
                    else:
                        dscfg['global_data']['grid_out'].fields[
                            field_name]['data'] = field_std
                else:
                    if lin_trans:
                        dscfg['global_data']['grid_out'].fields[
                            field_name]['data'] = 10.*np.ma.log10(
                                field_std/field_mean)
                    else:
                        dscfg['global_data']['grid_out'].fields[
                            field_name]['data'] = field_std/field_mean

        new_dataset = {
            'radar_out': deepcopy(dscfg['global_data']['grid_out']),
            'timeinfo': dscfg['global_data']['endtime']}

        dscfg['global_data']['starttime'] += datetime.timedelta(
            seconds=period)
        dscfg['global_data']['endtime'] += datetime.timedelta(seconds=period)

        # remove old grid object from global_data dictionary
        dscfg['global_data'].pop('grid_out', None)

        # get start and stop times of new grid object
        dscfg['global_data']['starttime'], dscfg['global_data']['endtime'] = (
            time_avg_range(
                dscfg['timeinfo'], dscfg['global_data']['starttime'],
                dscfg['global_data']['endtime'], period))

        # check if volume time older than starttime
        if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
            dscfg['global_data'].update({'grid_out': grid_aux})

        return new_dataset, ind_rad

    # no more files to process if there is global data pack it up
    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None, None
        if 'grid_out' not in dscfg['global_data']:
            return None, None

        if stat in ('mean', 'std', 'cov'):
            field_mean = (
                dscfg['global_data']['grid_out'].fields[field_name]['data'] /
                dscfg['global_data']['grid_out'].fields[
                    'number_of_samples']['data'])

            if stat == 'mean':
                if lin_trans:
                    dscfg['global_data']['grid_out'].fields[
                        field_name]['data'] = 10.*np.ma.log10(field_mean)
                else:
                    dscfg['global_data']['grid_out'].fields[
                        field_name]['data'] = field_mean

            elif stat in ('std', 'cov'):
                field_std = np.ma.sqrt(
                    dscfg['global_data']['grid_out'].fields[
                        'sum_squared']['data'] /
                    dscfg['global_data']['grid_out'].fields[
                        'number_of_samples']['data']-field_mean*field_mean)
                if stat == 'std':
                    if lin_trans:
                        dscfg['global_data']['grid_out'].fields[
                            field_name]['data'] = 10.*np.ma.log10(field_std)
                    else:
                        dscfg['global_data']['grid_out'].fields[
                            field_name]['data'] = field_std
                else:
                    dscfg['global_data']['grid_out'].fields[
                        field_name]['data'] = field_std/field_mean

        new_dataset = {
            'radar_out': deepcopy(dscfg['global_data']['grid_out']),
            'timeinfo': dscfg['global_data']['endtime']}

        return new_dataset, ind_rad


def process_grid_time_stats2(procstatus, dscfg, radar_list=None):
    """
    computes temporal statistics of a field

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        period : float. Dataset keyword
            the period to average [s]. If -1 the statistics are going to be
            performed over the entire data. Default 3600.
        start_average : float. Dataset keyword
            when to start the average [s from midnight UTC]. Default 0.
        stat: string. Dataset keyword
            Statistic to compute: Can be median, mode, percentileXX
        use_nan : bool. Dataset keyword
            If true non valid data will be used
        nan_value : float. Dataset keyword
            The value of the non valid data. Default 0
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_name = get_fieldname_pyart(datatype)
        break
    ind_rad = int(radarnr[5:8])-1

    start_average = dscfg.get('start_average', 0.)
    period = dscfg.get('period', 3600.)
    use_nan = dscfg.get('use_nan', 0)
    nan_value = dscfg.get('nan_value', 0.)
    stat = dscfg.get('stat', 'median')
    if 'percentile' in stat:
        percentile = float(stat.replace('percentile', ''))

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        grid = radar_list[ind_rad]

        if field_name not in grid.fields:
            warn(field_name+' not available.')
            return None, None

        # prepare auxiliary radar
        field_dict = deepcopy(grid.fields[field_name])
        if use_nan:
            field_dict['data'] = np.ma.asarray(
                field_dict['data'].filled(nan_value))
        npoints_dict = pyart.config.get_metadata('number_of_samples')
        npoints_dict['data'] = np.ma.asarray(
            np.logical_not(np.ma.getmaskarray(field_dict['data'])), dtype=int)

        grid_aux = deepcopy(grid)
        grid_aux.fields = dict()
        grid_aux.add_field(field_name, field_dict)
        grid_aux.add_field('number_of_samples', npoints_dict)

        # first volume: initialize start and end time of averaging
        if dscfg['initialized'] == 0:
            avg_par = dict()
            if period != -1:
                date_00 = dscfg['timeinfo'].replace(
                    hour=0, minute=0, second=0, microsecond=0)

                avg_par.update(
                    {'starttime': date_00+datetime.timedelta(
                        seconds=start_average)})
                avg_par.update(
                    {'endtime': avg_par['starttime']+datetime.timedelta(
                        seconds=period)})
            else:
                avg_par.update({'starttime': dscfg['timeinfo']})
                avg_par.update({'endtime': dscfg['timeinfo']})

            avg_par.update({'timeinfo': dscfg['timeinfo']})
            dscfg['global_data'] = avg_par
            dscfg['initialized'] = 1

        if dscfg['initialized'] == 0:
            return None, None

        dscfg['global_data']['timeinfo'] = dscfg['timeinfo']
        # no grid object in global data: create it
        if 'grid_out' not in dscfg['global_data']:
            if period != -1:
                # get start and stop times of new grid object
                (dscfg['global_data']['starttime'],
                 dscfg['global_data']['endtime']) = (
                     time_avg_range(
                         dscfg['timeinfo'], dscfg['global_data']['starttime'],
                         dscfg['global_data']['endtime'], period))

                # check if volume time older than starttime
                if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
                    dscfg['global_data'].update({'grid_out': grid_aux})
                    dscfg['global_data'].update(
                        {'field_data': np.expand_dims(
                            grid_aux.fields[field_name]['data'], axis=0)})
            else:
                dscfg['global_data'].update({'grid_out': grid_aux})
                dscfg['global_data'].update(
                    {'field_data': np.expand_dims(
                        grid_aux.fields[field_name]['data'], axis=0)})

            return None, None

        # still accumulating: add field_dict to global field_dict
        if (period == -1 or
                dscfg['timeinfo'] < dscfg['global_data']['endtime']):

            if period == -1:
                dscfg['global_data']['endtime'] = dscfg['timeinfo']

            dscfg['global_data']['grid_out'].fields['number_of_samples'][
                'data'] += npoints_dict['data']

            dscfg['global_data']['field_data'] = np.ma.append(
                dscfg['global_data']['field_data'],
                np.expand_dims(field_dict['data'], axis=0), axis=0)

            return None, None

        # we have reached the end of the accumulation period: do the averaging
        # and start a new object (only reachable if period != -1)
        if stat == 'median':
            dscfg['global_data']['grid_out'].fields[
                field_name]['data'] = np.ma.median(
                    dscfg['global_data']['field_data'], axis=0)
        elif stat == 'mode':
            mode_data, _ = scipy.stats.mode(
                dscfg['global_data']['field_data'].filled(fill_value=np.nan),
                axis=0, nan_policy='omit')
            dscfg['global_data']['grid_out'].fields[field_name]['data'] = (
                np.ma.masked_invalid(np.squeeze(mode_data, axis=0)))
        elif 'percentile' in stat:
            percent_data = np.nanpercentile(
                dscfg['global_data']['field_data'].filled(fill_value=np.nan),
                percentile, axis=0)
            dscfg['global_data']['grid_out'].fields[field_name]['data'] = (
                np.ma.masked_invalid(percent_data))

        new_dataset = {
            'radar_out': deepcopy(dscfg['global_data']['grid_out']),
            'timeinfo': dscfg['global_data']['endtime']}

        dscfg['global_data']['starttime'] += datetime.timedelta(
            seconds=period)
        dscfg['global_data']['endtime'] += datetime.timedelta(seconds=period)

        # remove old grid object from global_data dictionary
        dscfg['global_data'].pop('grid_out', None)

        # get start and stop times of new grid object
        dscfg['global_data']['starttime'], dscfg['global_data']['endtime'] = (
            time_avg_range(
                dscfg['timeinfo'], dscfg['global_data']['starttime'],
                dscfg['global_data']['endtime'], period))

        # check if volume time older than starttime
        if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
            dscfg['global_data'].update({'grid_out': grid_aux})

        return new_dataset, ind_rad

    # no more files to process if there is global data pack it up
    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None, None
        if 'grid_out' not in dscfg['global_data']:
            return None, None

        if stat == 'median':
            dscfg['global_data']['grid_out'].fields[field_name]['data'] = (
                np.ma.median(dscfg['global_data']['field_data'], axis=0))
        elif stat == 'mode':
            mode_data, _ = scipy.stats.mode(
                dscfg['global_data']['field_data'].filled(fill_value=np.nan),
                axis=0, nan_policy='omit')
            dscfg['global_data']['grid_out'].fields[field_name]['data'] = (
                np.ma.masked_invalid(np.squeeze(mode_data, axis=0)))
        elif 'percentile' in stat:
            percent_data = np.nanpercentile(
                dscfg['global_data']['field_data'].filled(fill_value=np.nan),
                percentile, axis=0)
            dscfg['global_data']['grid_out'].fields[field_name]['data'] = (
                np.ma.masked_invalid(percent_data))

        new_dataset = {
            'radar_out': deepcopy(dscfg['global_data']['grid_out']),
            'timeinfo': dscfg['global_data']['endtime']}

        return new_dataset, ind_rad


def process_grid_fields_diff(procstatus, dscfg, radar_list=None):
    """
    Computes grid field differences

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

    Returns
    -------
    new_dataset : dict
        dictionary containing a radar object containing the field differences
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # create the list of data types for each radar
    if len(dscfg['datatype']) != 2:
        warn(
            'Two and only two fields are required to compute the field'
            ' differences')
        return None, None

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
    field_name_1 = get_fieldname_pyart(datatype)

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][1])
    field_name_2 = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None

    grid = radar_list[ind_rad]
    if ((field_name_1 not in grid.fields) or
            (field_name_2 not in grid.fields)):
        warn('Unable to compare fields '+field_name_1+'and '+field_name_2 +
             '. Fields missings')
        return None, None

    field_diff = pyart.config.get_metadata('fields_difference')
    field_diff['data'] = (
        grid.fields[field_name_1]['data'] -
        grid.fields[field_name_2]['data'])
    field_diff['long_name'] = field_name_1+' - '+field_name_2

    grid_diff = deepcopy(grid)
    grid_diff.fields = dict()
    grid_diff.add_field('fields_difference', field_diff)

    new_dataset = {'radar_out': grid_diff}

    return new_dataset, ind_rad


def process_grid_texture(procstatus, dscfg, radar_list=None):
    """
    Computes the 2D texture of a gridded field

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        xwind, ywind : int
            The size of the local window in the x and y axis. Default 7
        fill_value : float
            The value with which to fill masked data. Default np.NaN
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing a radar object containing the field differences
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # create the list of data types for each radar
    if len(dscfg['datatype']) != 1:
        warn('Texture can only be computed on one field')
        return None, None

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
    field_name = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None

    grid = radar_list[ind_rad]
    if field_name not in grid.fields:
        warn('Unable to compute field '+field_name+' texture. Field missing')
        return None, None

    xwind = dscfg.get('xwind', 7)
    ywind = dscfg.get('ywind', 7)
    fill_value = dscfg.get('fill_value', np.NaN)

    field_text = pyart.config.get_metadata('field_texture')
    field_text['data'] = np.ma.masked_all(
        (grid.nz, grid.ny, grid.nx),
        dtype=grid.fields[field_name]['data'].dtype)

    for level in range(grid.nz):
        field_array = deepcopy(grid.fields[field_name]['data'][level, :, :])
        field_array = field_array.filled(fill_value=fill_value)
        field_text['data'][level, :, :] = pyart.util.grid_texture_2d(
            field_array, xwind=xwind, ywind=ywind)

    grid_text = deepcopy(grid)
    grid_text.fields = dict()
    grid_text.add_field('field_texture', field_text)

    new_dataset = {'radar_out': grid_text}

    return new_dataset, ind_rad


def process_grid_mask(procstatus, dscfg, radar_list=None):
    """
    Mask data. Puts True if data is within thresholds, False if it is not.
    Thresholds can be min, max or both min and max

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration
    radar_list : list of Radar objects
        Optional. list of radar objects
    threshold_min : float or None
        Threshold used for the mask. Values below threshold are set to False.
        Above threshold are set to True. Default None.
    threshold_max : float or None
        Threshold used for the mask. Values above threshold are set to False.
        Below threshold are set to True. Default None.
    x_dir_ext, y_dir_ext : int
        Number of pixels by which to extend the mask on each side of the
        west-east direction and south-north direction

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
    field_name = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None

    grid = radar_list[ind_rad]

    if field_name not in grid.fields:
        warn('Unable to mask field '+field_name+' Field missing in grid')
        return None, None

    threshold_min = dscfg.get('threshold_min', None)
    threshold_max = dscfg.get('threshold_max', None)
    x_dir_ext = dscfg.get('x_dir_ext', 0)
    y_dir_ext = dscfg.get('y_dir_ext', 0)
    
    if threshold_min is None and threshold_max is None:
        warn('At least one threshold necessary')
        return None, None

    field_mask = pyart.config.get_metadata('field_mask')
    field_mask['data'] = np.ma.masked_all(
        (grid.nz, grid.ny, grid.nx), dtype=np.int8)
    field_mask['data'][:] = 0
    valid = np.logical_not(np.ma.getmaskarray(grid.fields[field_name]['data']))
    field_mask['data'][valid] = 1
    
    if threshold_min is not None and threshold_max is not None:
        field_mask['data'][
            (grid.fields[field_name]['data'] >= threshold_min) &
            (grid.fields[field_name]['data'] <= threshold_max)] = 2
        field_mask['long_name'] = '{} threshold {}-{}'.format(
            field_name, threshold_min, threshold_max)
    elif threshold_min is not None:
        field_mask['data'][grid.fields[field_name]['data'] >= threshold_min] = 2
        field_mask['long_name'] = '{} threshold_min {}'.format(
            field_name, threshold_min)
    elif threshold_max is not None:
        field_mask['data'][grid.fields[field_name]['data'] <= threshold_max] = 2
        field_mask['long_name'] = '{} threshold_max {}'.format(
            field_name, threshold_max)        

    if x_dir_ext > 0 or y_dir_ext > 0:
        ind_z, ind_y, ind_x = np.where(field_mask['data'] == 2)

        if ind_z.size > 0:
            for z, y, x in zip(ind_z, ind_y, ind_x):
                field_mask['data'][
                    0, y-y_dir_ext:y+y_dir_ext+1,
                    x-x_dir_ext:x+x_dir_ext+1] = 2

    grid_mask = deepcopy(grid)
    grid_mask.fields = dict()
    grid_mask.add_field('field_mask', field_mask)

    new_dataset = {'radar_out': grid_mask}

    return new_dataset, ind_rad


def process_normalize_luminosity(procstatus, dscfg, radar_list=None):
    """
    Normalize the data by the sinus of the sun elevation. The sun elevation is
    computed at the central pixel.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    radarnr, _, _, _, _ = get_datatype_fields(dscfg['datatype'][0])
    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None

    grid = radar_list[ind_rad]

    # get normalization factor
    lat_point = grid.point_latitude['data'][0, int(grid.ny/2), int(grid.nx/2)]
    lon_point = grid.point_longitude['data'][
        0, int(grid.ny/2), int(grid.nx/2)]

    el_sun, _ = pyart.correct.sun_position_pysolar(
        dscfg['timeinfo'], lat_point, lon_point)
    norm = np.sin(np.deg2rad(el_sun))

    grid_norm = deepcopy(grid)
    grid_norm.fields = dict()

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_name = get_fieldname_pyart(datatype)

        if field_name not in grid.fields:
            warn('Unable to normalize field '+field_name +
                 '. Field missing in grid')
            continue

        norm_field = pyart.config.get_metadata(field_name+'_norm')
        norm_field['data'] = grid.fields[field_name]['data']/norm

        grid_norm.add_field(field_name+'_norm', norm_field)

    new_dataset = {'radar_out': grid_norm}

    return new_dataset, ind_rad


def process_pixel_filter(procstatus, dscfg, radar_list=None):
    """
    Masks all pixels that are not of the class specified in keyword pixel_type

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        pixel_type : int or list of ints
            The type of pixels to keep: 0 No data, 1 Below threshold, 2 Above
            threshold. Default 2
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    mask_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'mask':
            mask_field = get_fieldname_pyart(datatype)
            break
    if mask_field is None:
        warn('mask field required to filter data')
        return None, None

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None

    grid = radar_list[ind_rad]

    if mask_field not in grid.fields:
        warn('Unable to filter data. Missing mask field')
        return None, None

    pixel_type = dscfg.get('pixel_type', 2)
    mask = np.ma.isin(
        grid.fields[mask_field]['data'], pixel_type, invert=True)

    grid_mask = deepcopy(grid)
    grid_mask.fields = dict()

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'mask':
            continue

        field_name = get_fieldname_pyart(datatype)
        if field_name not in grid.fields:
            warn('Unable to normalize field '+field_name +
                 '. Field missing in grid')
            continue

        field = deepcopy(grid.fields[field_name])
        field['data'] = np.ma.masked_where(mask, field['data'])

        grid_mask.add_field(field_name, field)

    new_dataset = {'radar_out': grid_mask}

    return new_dataset, ind_rad
