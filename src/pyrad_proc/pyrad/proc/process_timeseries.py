"""
pyrad.proc.process_timeseries
=============================

Functions to obtain time series of radar data

.. autosummary::
    :toctree: generated/

    process_point_measurement
    process_qvp
    process_rqvp
    process_evp
    process_svp
    process_time_height
    process_ts_along_coord

"""

from warnings import warn
import numpy as np
from netCDF4 import num2date

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart


def process_point_measurement(procstatus, dscfg, radar_list=None):
    """
    Obtains the radar data at a point location.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement
        single_point : boolean. Dataset keyword
            if True only one gate per radar volume is going to be kept.
            Otherwise all gates within the azimuth and elevation tolerance
            are going to be kept. This is useful to extract all data from
            fixed pointing scans. Default True
        latlon : boolean. Dataset keyword
            if True position is obtained from latitude, longitude information,
            otherwise position is obtained from antenna coordinates
            (range, azimuth, elevation).
        truealt : boolean. Dataset keyword
            if True the user input altitude is used to determine the point of
            interest.
            if False use the altitude at a given radar elevation ele over the
            point of interest.
        lon : float. Dataset keyword
            the longitude [deg]. Use when latlon is True.
        lat : float. Dataset keyword
            the latitude [deg]. Use when latlon is True.
        alt : float. Dataset keyword
            altitude [m MSL]. Use when latlon is True.
        ele : float. Dataset keyword
            radar elevation [deg]. Use when latlon is False or when latlon is
            True and truealt is False
        azi : float. Dataset keyword
            radar azimuth [deg]. Use when latlon is False
        rng : float. Dataset keyword
            range from radar [m]. Use when latlon is False
        AziTol : float. Dataset keyword
            azimuthal tolerance to determine which radar azimuth to use [deg]
        EleTol : float. Dataset keyword
            elevation tolerance to determine which radar elevation to use [deg]
        RngTol : float. Dataset keyword
            range tolerance to determine which radar bin to use [m]

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
            'antenna_coordinates_az_el_r': (
                dscfg['global_data']['antenna_coordinates_az_el_r']),
            'final': True}

        return new_dataset, ind_rad

    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if field_name not in radar.fields:
        warn('Unable to extract point measurement information. ' +
             'Field not available')
        return None, None

    projparams = dict()
    projparams.update({'proj': 'pyart_aeqd'})
    projparams.update({'lon_0': radar.longitude['data']})
    projparams.update({'lat_0': radar.latitude['data']})

    single_point = dscfg.get('single_point', True)
    if dscfg['latlon']:
        lon = dscfg['lon']
        lat = dscfg['lat']
        alt = dscfg['alt']
        x, y = pyart.core.geographic_to_cartesian(lon, lat, projparams)

        if not dscfg['truealt']:
            ke = 4./3.  # constant for effective radius
            a = 6378100.  # earth radius
            re = a * ke  # effective radius

            elrad = dscfg['ele'] * np.pi / 180.
            r_ground = np.sqrt(x ** 2. + y ** 2.)
            r = r_ground / np.cos(elrad)
            alt_radar = radar.altitude['data']+np.sqrt(
                r ** 2. + re ** 2. + 2. * r * re * np.sin(elrad)) - re
            alt_radar = alt_radar[0]
        else:
            alt_radar = dscfg['alt']

        r, az, el = pyart.core.cartesian_to_antenna(
            x, y, alt_radar-radar.altitude['data'])
        r = r[0]
        az = az[0]
        el = el[0]
    else:
        r = dscfg['rng']
        az = dscfg['azi']
        el = dscfg['ele']

        x, y, alt = pyart.core.antenna_to_cartesian(r/1000., az, el)
        lon, lat = pyart.core.cartesian_to_geographic(x, y, projparams)
        lon = lon[0]
        lat = lat[0]

    d_az = np.abs(radar.azimuth['data'] - az)
    if np.min(d_az) > dscfg['AziTol']:
        warn(' No radar bin found for point (az, el, r):(' +
             str(az)+', '+str(el)+', '+str(r) +
             '). Minimum distance to radar azimuth '+str(d_az) +
             ' larger than tolerance')
        return None, None

    d_el = np.abs(radar.elevation['data'] - el)
    if np.min(d_el) > dscfg['EleTol']:
        warn(' No radar bin found for point (az, el, r):(' +
             str(az)+', '+str(el)+', '+str(r) +
             '). Minimum distance to radar elevation '+str(d_el) +
             ' larger than tolerance')
        return None, None

    d_r = np.abs(radar.range['data'] - r)
    if np.min(d_r) > dscfg['RngTol']:
        warn(' No radar bin found for point (az, el, r):(' +
             str(az)+', '+str(el)+', '+str(r) +
             '). Minimum distance to radar range bin '+str(d_r) +
             ' larger than tolerance')
        return None, None

    if single_point:
        ind_ray = np.argmin(np.abs(radar.azimuth['data'] - az) +
                            np.abs(radar.elevation['data'] - el))
    else:
        ind_ray = np.where(np.logical_and(
            d_az <= dscfg['AziTol'], d_el <= dscfg['EleTol']))[0]
    ind_r = np.argmin(np.abs(radar.range['data'] - r))

    ant_coord = np.empty((3, ind_ray.size), np.float32)
    ant_coord[0, :] = radar.azimuth['data'][ind_ray]
    ant_coord[1, :] = radar.elevation['data'][ind_ray]
    ant_coord[2, :] = np.zeros(ind_ray.size)+radar.range['data'][ind_r]

    val = radar.fields[field_name]['data'].data[ind_ray, ind_r]
    time = num2date(radar.time['data'][ind_ray], radar.time['units'],
                    radar.time['calendar'])

    # initialize dataset
    if dscfg['initialized'] == 0:
        poi = {
            'point_coordinates_WGS84_lon_lat_alt': [lon, lat, alt],
            'antenna_coordinates_az_el_r': [az, el, r],
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
    new_dataset.update({'antenna_coordinates_az_el_r': [az, el, r]})
    new_dataset.update({'used_antenna_coordinates_az_el_r': ant_coord})
    new_dataset.update({'final': False})

    return new_dataset, ind_rad


def process_qvp(procstatus, dscfg, radar_list=None):
    """
    Computes quasi vertical profiles, by averaging over height levels
    PPI data.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement
        angle : int or float
            If the radar object contains a PPI volume, the sweep number to
            use, if it contains an RHI volume the elevation angle.
            Default 0.
        ang_tol : float
            If the radar object contains an RHI volume, the tolerance in the
            elevation angle for the conversion into PPI
        hmax : float
            The maximum height to plot [m]. Default 10000.
        hres : float
            The height resolution [m]. Default 50
        avg_type : str
            The type of averaging to perform. Can be either "mean" or "median"
            Default "mean"
        nvalid_min : int
            Minimum number of valid points to accept average. Default 30.
        interp_kind : str
            type of interpolation when projecting to vertical grid: 'none',
            or 'nearest', etc. Default 'none'
            'none' will select from all data points within the regular grid
            height bin the closest to the center of the bin.
            'nearest' will select the closest data point to the center of the
            height bin regardless if it is within the height bin or not.
            Data points can be masked values
            If another type of interpolation is selected masked values will be
            eliminated from the data points before the interpolation

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the QVP and a keyboard stating whether the
        processing has finished or not.
    ind_rad : int
        radar index

    Reference
    ---------
    Ryzhkov A., Zhang P., Reeves H., Kumjian M., Tschallener T., Trömel S.,
    Simmer C. 2016: Quasi-Vertical Profiles: A New Way to Look at Polarimetric
    Radar Data. JTECH vol. 33 pp 551-562

    """
    if procstatus == 0:
        return None, None

    if procstatus == 1:
        field_names = []
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            field_names.append(get_fieldname_pyart(datatype))

        ind_rad = int(radarnr[5:8])-1

        if (radar_list is None) or (radar_list[ind_rad] is None):
            warn('ERROR: No valid radar')
            return None, None

        radar = radar_list[ind_rad]

        # default parameters
        angle = dscfg.get('angle', 0)
        ang_tol = dscfg.get('ang_tol', 1.)
        hmax = dscfg.get('hmax', 10000.)
        hres = dscfg.get('hres', 50.)
        avg_type = dscfg.get('avg_type', 'mean')
        nvalid_min = dscfg.get('nvalid_min', 30)
        interp_kind = dscfg.get('interp_kind', 'none')

        # initialize dataset
        if not dscfg['initialized']:
            qvp = pyart.retrieve.compute_qvp(
                radar, field_names, ref_time=dscfg['timeinfo'],
                angle=angle, ang_tol=ang_tol, hmax=hmax, hres=hres,
                avg_type=avg_type, nvalid_min=nvalid_min,
                interp_kind=interp_kind, qvp=None)

            if qvp is None:
                warn('Unable to compute QVP')
                return None, None

            global_dict = dict()
            global_dict.update({'start_time': dscfg['timeinfo']})
            global_dict.update({'radar_out': qvp})
            dscfg['global_data'] = global_dict
            dscfg['initialized'] = 1
        else:
            qvp = pyart.retrieve.compute_qvp(
                radar, field_names, ref_time=dscfg['timeinfo'],
                angle=angle, ang_tol=ang_tol, hmax=hmax, hres=hres,
                avg_type=avg_type, nvalid_min=nvalid_min,
                interp_kind=interp_kind,
                qvp=dscfg['global_data']['radar_out'])

            if qvp is None:
                warn('Unable to compute QVP')
                return None, None

            dscfg['global_data']['radar_out'] = qvp

        new_dataset = dict()
        new_dataset.update({'radar_out': qvp})
        new_dataset.update({'radar_type': 'temporal'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad

    if procstatus == 2:
        if not dscfg['initialized']:
            return None, None

        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            break

        ind_rad = int(radarnr[5:8])-1

        qvp = dscfg['global_data']['radar_out']

        new_dataset = dict()
        new_dataset.update({'radar_out': qvp})
        new_dataset.update({'radar_type': 'final'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad


def process_rqvp(procstatus, dscfg, radar_list=None):
    """
    Computes range defined quasi vertical profiles, by averaging over height
    levels PPI data.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement
        hmax : float
            The maximum height to plot [m]. Default 10000.
        hres : float
            The height resolution [m]. Default 2.
        avg_type : str
            The type of averaging to perform. Can be either "mean" or "median"
            Default "mean"
        nvalid_min : int
            Minimum number of valid points to accept average. Default 30.
        interp_kind : str
            type of interpolation when projecting to vertical grid: 'none',
            or 'nearest', etc. Default 'nearest'
            'none' will select from all data points within the regular grid
            height bin the closest to the center of the bin.
            'nearest' will select the closest data point to the center of the
            height bin regardless if it is within the height bin or not.
            Data points can be masked values
            If another type of interpolation is selected masked values will be
            eliminated from the data points before the interpolation
        rmax : float
            ground range up to which the data is intended for use [m].
            Default 50000.
        weight_power : float
            Power p of the weighting function 1/abs(grng-(rmax-1))**p given to
            the data outside the desired range. -1 will set the weight to 0.
            Default 2.

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the QVP and a keyboard stating whether the
        processing has finished or not.
    ind_rad : int
        radar index

    Reference
    ---------
    Tobin D.M., Kumjian M.R. 2017: Polarimetric Radar and Surface-Based
    Precipitation-Type Observations of ice Pellet to Freezing Rain
    Transitions. Weather and Forecasting vol. 32 pp 2065-2082

    """
    if procstatus == 0:
        return None, None

    if procstatus == 1:
        field_names = []
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            field_names.append(get_fieldname_pyart(datatype))

        ind_rad = int(radarnr[5:8])-1

        if (radar_list is None) or (radar_list[ind_rad] is None):
            warn('ERROR: No valid radar')
            return None, None

        radar = radar_list[ind_rad]

        # default parameters
        hmax = dscfg.get('hmax', 10000.)
        hres = dscfg.get('hres', 2.)
        avg_type = dscfg.get('avg_type', 'mean')
        nvalid_min = dscfg.get('nvalid_min', 30)
        interp_kind = dscfg.get('interp_kind', 'nearest')
        rmax = dscfg.get('rmax', 50000.)
        weight_power = dscfg.get('weight_power', 2.)

        # initialize dataset
        if not dscfg['initialized']:
            qvp = pyart.retrieve.compute_rqvp(
                radar, field_names, ref_time=dscfg['timeinfo'],
                hmax=hmax, hres=hres, avg_type=avg_type,
                nvalid_min=nvalid_min, interp_kind=interp_kind, rmax=rmax,
                weight_power=weight_power, qvp=None)

            if qvp is None:
                warn('Unable to compute QVP')
                return None, None

            global_dict = dict()
            global_dict.update({'start_time': dscfg['timeinfo']})
            global_dict.update({'radar_out': qvp})
            dscfg['global_data'] = global_dict
            dscfg['initialized'] = 1
        else:
            qvp = pyart.retrieve.compute_rqvp(
                radar, field_names, ref_time=dscfg['timeinfo'],
                hmax=hmax, hres=hres, avg_type=avg_type,
                nvalid_min=nvalid_min, interp_kind=interp_kind, rmax=rmax,
                weight_power=weight_power,
                qvp=dscfg['global_data']['radar_out'])

            if qvp is None:
                warn('Unable to compute QVP')
                return None, None

            dscfg['global_data']['radar_out'] = qvp

        new_dataset = dict()
        new_dataset.update({'radar_out': qvp})
        new_dataset.update({'radar_type': 'temporal'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad

    if procstatus == 2:
        if not dscfg['initialized']:
            return None, None

        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            break

        ind_rad = int(radarnr[5:8])-1

        qvp = dscfg['global_data']['radar_out']

        new_dataset = dict()
        new_dataset.update({'radar_out': qvp})
        new_dataset.update({'radar_type': 'final'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad


def process_evp(procstatus, dscfg, radar_list=None):
    """
    Computes enhanced vertical profiles, by averaging over height levels
    PPI data.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement
        lat, lon : float
            latitude and longitude of the point of interest [deg]
        latlon_tol : float
            tolerance in latitude and longitude in deg. Default 0.0005
        delta_rng, delta_azi : float
            maximum range distance [m] and azimuth distance [degree] from the
            central point of the evp containing data to average. Default 5000.
            and 10.
        hmax : float
            The maximum height to plot [m]. Default 10000.
        hres : float
            The height resolution [m]. Default 250.
        avg_type : str
            The type of averaging to perform. Can be either "mean" or "median"
            Default "mean"
        nvalid_min : int
            Minimum number of valid points to consider the data valid when
            performing the averaging. Default 1
        interp_kind : str
            type of interpolation when projecting to vertical grid: 'none',
            or 'nearest', etc. Default 'none'.
            'none' will select from all data points within the regular grid
            height bin the closest to the center of the bin.
            'nearest' will select the closest data point to the center of the
            height bin regardless if it is within the height bin or not.
            Data points can be masked values
            If another type of interpolation is selected masked values will be
            eliminated from the data points before the interpolation

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the EVP and a keyboard stating whether the
        processing has finished or not.
    ind_rad : int
        radar index

    Reference
    ---------
    Kaltenboeck R., Ryzhkov A. 2016: A freezing rain storm explored with a
    C-band polarimetric weather radar using the QVP methodology.
    Meteorologische Zeitschrift vol. 26 pp 207-222

    """
    if procstatus == 0:
        return None, None

    if procstatus == 1:
        field_names = []
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            field_names.append(get_fieldname_pyart(datatype))

        ind_rad = int(radarnr[5:8])-1

        if (radar_list is None) or (radar_list[ind_rad] is None):
            warn('ERROR: No valid radar')
            return None, None

        radar = radar_list[ind_rad]

        # default parameters
        lon = dscfg['lon']
        lat = dscfg['lat']
        latlon_tol = dscfg.get('latlon_tol', 0.0005)
        delta_rng = dscfg.get('delta_rng', 15000.)
        delta_azi = dscfg.get('delta_azi', 10.)
        hmax = dscfg.get('hmax', 10000.)
        hres = dscfg.get('hres', 250.)
        avg_type = dscfg.get('avg_type', 'mean')
        nvalid_min = dscfg.get('nvalid_min', 1)
        interp_kind = dscfg.get('interp_kind', 'none')

        # initialize dataset
        if not dscfg['initialized']:
            qvp = pyart.retrieve.compute_evp(
                radar, field_names, lon, lat, ref_time=dscfg['timeinfo'],
                latlon_tol=latlon_tol, delta_rng=delta_rng,
                delta_azi=delta_azi, hmax=hmax, hres=hres, avg_type=avg_type,
                nvalid_min=nvalid_min, interp_kind=interp_kind, qvp=None)

            if qvp is None:
                warn('Unable to compute QVP')
                return None, None

            global_dict = dict()
            global_dict.update({'start_time': dscfg['timeinfo']})
            global_dict.update({'radar_out': qvp})
            dscfg['global_data'] = global_dict
            dscfg['initialized'] = 1
        else:
            qvp = pyart.retrieve.compute_evp(
                radar, field_names, lon, lat, ref_time=dscfg['timeinfo'],
                latlon_tol=latlon_tol, delta_rng=delta_rng,
                delta_azi=delta_azi, hmax=hmax, hres=hres, avg_type=avg_type,
                nvalid_min=nvalid_min, interp_kind=interp_kind,
                qvp=dscfg['global_data']['radar_out'])

            if qvp is None:
                warn('Unable to compute QVP')
                return None, None

        dscfg['global_data']['radar_out'] = qvp

        new_dataset = dict()
        new_dataset.update({'radar_out': qvp})
        new_dataset.update({'radar_type': 'temporal'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad

    if procstatus == 2:
        if not dscfg['initialized']:
            return None, None

        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            break

        ind_rad = int(radarnr[5:8])-1

        qvp = dscfg['global_data']['radar_out']

        new_dataset = dict()
        new_dataset.update({'radar_out': qvp})
        new_dataset.update({'radar_type': 'final'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad


def process_svp(procstatus, dscfg, radar_list=None):
    """
    Computes slanted vertical profiles, by averaging over height levels
    PPI data.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement
        angle : int or float
            If the radar object contains a PPI volume, the sweep number to
            use, if it contains an RHI volume the elevation angle.
            Default 0.
        ang_tol : float
            If the radar object contains an RHI volume, the tolerance in the
            elevation angle for the conversion into PPI. Default 1.
        lat, lon : float
            latitude and longitude of the point of interest [deg]
        latlon_tol : float
            tolerance in latitude and longitude in deg. Default 0.0005
        delta_rng, delta_azi : float
            maximum range distance [m] and azimuth distance [degree] from the
            central point of the svp containing data to average. Default 5000.
            and 10.
        hmax : float
            The maximum height to plot [m]. Default 10000.
        hres : float
            The height resolution [m]. Default 250.
        avg_type : str
            The type of averaging to perform. Can be either "mean" or "median"
            Default "mean"
        nvalid_min : int
            Minimum number of valid points to consider the data valid when
            performing the averaging. Default 1
        interp_kind : str
            type of interpolation when projecting to vertical grid: 'none',
            or 'nearest', etc. Default 'none'
            'none' will select from all data points within the regular grid
            height bin the closest to the center of the bin.
            'nearest' will select the closest data point to the center of the
            height bin regardless if it is within the height bin or not.
            Data points can be masked values
            If another type of interpolation is selected masked values will be
            eliminated from the data points before the interpolation

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the svp and a keyboard stating whether the
        processing has finished or not.
    ind_rad : int
        radar index

    Reference
    ---------
    Bukovcic P., Zrnic D., Zhang G. 2017: Winter Precipitation Liquid-Ice
    Phase Transitions Revealed with Polarimetric Radar and 2DVD Observations
    in Central Oklahoma. JTECH vol. 56 pp 1345-1363

    """
    if procstatus == 0:
        return None, None

    if procstatus == 1:
        field_names = []
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            field_names.append(get_fieldname_pyart(datatype))

        ind_rad = int(radarnr[5:8])-1

        if (radar_list is None) or (radar_list[ind_rad] is None):
            warn('ERROR: No valid radar')
            return None, None

        radar = radar_list[ind_rad]

        # default parameters
        angle = dscfg.get('angle', 0)
        ang_tol = dscfg.get('ang_tol', 1.)
        lon = dscfg['lon']
        lat = dscfg['lat']
        latlon_tol = dscfg.get('latlon_tol', 0.0005)
        delta_rng = dscfg.get('delta_rng', 30000.)
        delta_azi = dscfg.get('delta_azi', 10.)
        hmax = dscfg.get('hmax', 10000.)
        hres = dscfg.get('hres', 250.)
        avg_type = dscfg.get('avg_type', 'mean')
        nvalid_min = dscfg.get('nvalid_min', 1)
        interp_kind = dscfg.get('interp_kind', 'none')

        # initialize dataset
        if not dscfg['initialized']:
            qvp = pyart.retrieve.compute_svp(
                radar, field_names, lon, lat, angle,
                ref_time=dscfg['timeinfo'], ang_tol=ang_tol,
                latlon_tol=latlon_tol, delta_rng=delta_rng,
                delta_azi=delta_azi, hmax=hmax, hres=hres, avg_type=avg_type,
                nvalid_min=nvalid_min, interp_kind=interp_kind, qvp=None)

            if qvp is None:
                warn('Unable to compute QVP')
                return None, None

            global_dict = dict()
            global_dict.update({'start_time': dscfg['timeinfo']})
            global_dict.update({'radar_out': qvp})
            dscfg['global_data'] = global_dict
            dscfg['initialized'] = 1
        else:
            qvp = pyart.retrieve.compute_svp(
                radar, field_names, lon, lat, angle,
                ref_time=dscfg['timeinfo'], ang_tol=ang_tol,
                latlon_tol=latlon_tol, delta_rng=delta_rng,
                delta_azi=delta_azi, hmax=hmax, hres=hres, avg_type=avg_type,
                nvalid_min=nvalid_min, interp_kind=interp_kind,
                qvp=dscfg['global_data']['radar_out'])

            if qvp is None:
                warn('Unable to compute QVP')
                return None, None

        dscfg['global_data']['radar_out'] = qvp

        new_dataset = dict()
        new_dataset.update({'radar_out': qvp})
        new_dataset.update({'radar_type': 'temporal'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad

    if procstatus == 2:
        if not dscfg['initialized']:
            return None, None

        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            break

        ind_rad = int(radarnr[5:8])-1

        svp = dscfg['global_data']['radar_out']

        new_dataset = dict()
        new_dataset.update({'radar_out': svp})
        new_dataset.update({'radar_type': 'final'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad


def process_time_height(procstatus, dscfg, radar_list=None):
    """
    Produces time height radar objects at a point of interest defined by
    latitude and longitude. A time-height contains the evolution
    of the vertical structure of radar measurements above the location
    of interest.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement
        lat, lon : float
            latitude and longitude of the point of interest [deg]
        latlon_tol : float
            tolerance in latitude and longitude in deg. Default 0.0005
        hmax : float
            The maximum height to plot [m]. Default 10000.
        hres : float
            The height resolution [m]. Default 50
        interp_kind : str
            type of interpolation when projecting to vertical grid: 'none',
            or 'nearest', etc. Default 'none'
            'none' will select from all data points within the regular grid
            height bin the closest to the center of the bin.
            'nearest' will select the closest data point to the center of the
            height bin regardless if it is within the height bin or not.
            Data points can be masked values
            If another type of interpolation is selected masked values will be
            eliminated from the data points before the interpolation

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the QVP and a keyboard stating whether the
        processing has finished or not.
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        return None, None

    if procstatus == 1:
        field_names = []
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            field_names.append(get_fieldname_pyart(datatype))

        ind_rad = int(radarnr[5:8])-1

        if (radar_list is None) or (radar_list[ind_rad] is None):
            warn('ERROR: No valid radar')
            return None, None

        radar = radar_list[ind_rad]

        # default parameters
        lon = dscfg['lon']
        lat = dscfg['lat']
        latlon_tol = dscfg.get('latlon_tol', 0.0005)
        hmax = dscfg.get('hmax', 10000.)
        hres = dscfg.get('hres', 50.)
        interp_kind = dscfg.get('interp_kind', 'none')

        # initialize dataset
        if not dscfg['initialized']:
            qvp = pyart.retrieve.compute_vp(
                radar, field_names, lon, lat, ref_time=dscfg['timeinfo'],
                latlon_tol=latlon_tol, hmax=hmax, hres=hres,
                interp_kind=interp_kind, qvp=None)

            if qvp is None:
                warn('Unable to compute QVP')
                return None, None

            global_dict = dict()
            global_dict.update({'start_time': dscfg['timeinfo']})
            global_dict.update({'radar_out': qvp})
            dscfg['global_data'] = global_dict
            dscfg['initialized'] = 1
        else:
            qvp = pyart.retrieve.compute_vp(
                radar, field_names, lon, lat, ref_time=dscfg['timeinfo'],
                latlon_tol=latlon_tol, hmax=hmax, hres=hres,
                interp_kind=interp_kind,
                qvp=dscfg['global_data']['radar_out'])

            if qvp is None:
                warn('Unable to compute QVP')
                return None, None

        dscfg['global_data']['radar_out'] = qvp

        new_dataset = dict()
        new_dataset.update({'radar_out': qvp})
        new_dataset.update({'radar_type': 'temporal'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad

    if procstatus == 2:
        if not dscfg['initialized']:
            return None, None

        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            break

        ind_rad = int(radarnr[5:8])-1

        qvp = dscfg['global_data']['radar_out']

        new_dataset = dict()
        new_dataset.update({'radar_out': qvp})
        new_dataset.update({'radar_type': 'final'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad


def process_ts_along_coord(procstatus, dscfg, radar_list=None):
    """
    Produces time series along a particular antenna coordinate

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the time series
        mode : str
            coordinate to extract data along. Can be ALONG_AZI, ALONG_ELE or
            ALONG_RNG
        fixed_range, fixed_azimuth, fixed_elevation : float
            The fixed range [m], azimuth [deg] or elevation [deg] to extract.
            In each mode two of these parameters have to be defined. If they
            are not defined they default to 0.
        ang_tol, rng_tol : float
            The angle tolerance [deg] and range tolerance [m] around the fixed
            range or azimuth/elevation
        value_start, value_stop : float
            The minimum and maximum value at which the data along a coordinate
            start and stop

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the data and a keyboard stating whether the
        processing has finished or not.
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        return None, None

    if procstatus == 1:
        radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
        field_names = []
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            field_names.append(get_fieldname_pyart(datatype))

        ind_rad = int(radarnr[5:8])-1

        if (radar_list is None) or (radar_list[ind_rad] is None):
            warn('ERROR: No valid radar')
            return None, None

        radar = radar_list[ind_rad]

        mode = dscfg.get('mode', 'ALONG_AZI')
        if mode == 'ALONG_RNG':
            value_start = dscfg.get('value_start', 0.)
            value_stop = dscfg.get('value_stop', radar.range['data'][-1])
            ang_tol = dscfg.get('AngTol', 1.)
            rng_tol = dscfg.get('RngTol', 50.)
            fixed_range = dscfg.get('fixed_range', None)
            fixed_azimuth = dscfg.get('fixed_azimuth', 0.)
            fixed_elevation = dscfg.get('fixed_elevation', 0.)
        elif mode == 'ALONG_AZI':
            value_start = dscfg.get(
                'value_start', np.min(radar.azimuth['data']))
            value_stop = dscfg.get(
                'value_stop', np.max(radar.azimuth['data']))
            ang_tol = dscfg.get('AngTol', 1.)
            rng_tol = dscfg.get('RngTol', 50.)
            fixed_range = dscfg.get('fixed_range', 0.)
            fixed_azimuth = dscfg.get('fixed_azimuth', None)
            fixed_elevation = dscfg.get('fixed_elevation', 0.)
        elif mode == 'ALONG_ELE':
            value_start = dscfg.get(
                'value_start', np.min(radar.elevation['data']))
            value_stop = dscfg.get(
                'value_stop', np.max(radar.elevation['data']))
            ang_tol = dscfg.get('AngTol', 1.)
            rng_tol = dscfg.get('RngTol', 50.)
            fixed_range = dscfg.get('fixed_range', 0.)
            fixed_azimuth = dscfg.get('fixed_azimuth', 0.)
            fixed_elevation = dscfg.get('fixed_elevation', None)
        else:
            warn('Unknown plotting mode '+dscfg['mode'])
            return None, None

        # initialize dataset
        if not dscfg['initialized']:
            acoord = pyart.retrieve.compute_ts_along_coord(
                radar, field_names, mode=mode, fixed_range=fixed_range,
                fixed_azimuth=fixed_azimuth, fixed_elevation=fixed_elevation,
                ang_tol=ang_tol, rng_tol=rng_tol, value_start=value_start,
                value_stop=value_stop, ref_time=dscfg['timeinfo'],
                acoord=None)

            if acoord is None:
                warn('Unable to compute time series along coordinate')
                return None, None

            global_dict = dict()
            global_dict.update({'start_time': dscfg['timeinfo']})
            global_dict.update({'radar_out': acoord})
            dscfg['global_data'] = global_dict
            dscfg['initialized'] = 1
        else:
            acoord = pyart.retrieve.compute_ts_along_coord(
                radar, field_names, mode=mode, fixed_range=fixed_range,
                fixed_azimuth=fixed_azimuth, fixed_elevation=fixed_elevation,
                ang_tol=ang_tol, rng_tol=rng_tol, value_start=value_start,
                value_stop=value_stop, ref_time=dscfg['timeinfo'],
                acoord=dscfg['global_data']['radar_out'])

            if acoord is None:
                warn('Unable to compute time series along coordinate')
                return None, None

        dscfg['global_data']['radar_out'] = acoord

        new_dataset = dict()
        new_dataset.update({'radar_out': acoord})
        new_dataset.update({'radar_type': 'temporal'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad

    if procstatus == 2:
        if not dscfg['initialized']:
            return None, None

        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            break

        ind_rad = int(radarnr[5:8])-1

        acoord = dscfg['global_data']['radar_out']

        new_dataset = dict()
        new_dataset.update({'radar_out': acoord})
        new_dataset.update({'radar_type': 'final'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad
