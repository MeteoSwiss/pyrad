"""
pyrad.proc.process_Doppler
===========================

Functions for processing Doppler related parameters

.. autosummary::
    :toctree: generated/

    process_turbulence
    process_dealias_fourdd
    process_dealias_region_based
    process_dealias_unwrap_phase
    process_radial_velocity
    process_wind_vel
    process_windshear
    process_vad
    process_dda

"""

from copy import deepcopy
from warnings import warn

import numpy as np

import pyart

try:
    import pytda
    _PYTDA_AVAILABLE = True
except ImportError:
    _PYTDA_AVAILABLE = False

try:
    import pydda
    _PYDDA_AVAILABLE = True
except ImportError:
    _PYDDA_AVAILABLE = False

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart


def process_turbulence(procstatus, dscfg, radar_list=None):
    """
    Computes turbulence from the Doppler spectrum width and reflectivity using
    the PyTDA package

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The input data type
        radius : float. Dataset keyword
            Search radius for calculating Eddy Dissipation Rate (EDR).
            Default 2
        split_cut : Bool. Dataset keyword
            Set to True for split-cut volumes. Default False
        max_split_cut : Int. Dataset keyword
            Total number of tilts that are affected by split cuts. Only
            relevant if split_cut=True. Default 2
        xran, yran : float array. Dataset keyword
            Spatial range in X,Y to consider. Default [-100, 100] for both
            X and Y
        beamwidth : Float. Dataset keyword
            Radar beamwidth. Default None. If None it will be obtained from
            the radar object metadata. If cannot be obtained defaults to 1
            deg.
        compute_gate_pos : Bool. Dataset keyword
            If True the gate position is going to be computed in PyTDA.
            Otherwise the position from the radar object is used. Default
            False
        verbose : Bool. Dataset keyword
            True for verbose output. Default False

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """
    if not _PYTDA_AVAILABLE:
        warn('PyTDA package not available. Unable to compute turbulence')
        return None, None

    if procstatus != 1:
        return None, None

    width_field = None
    refl_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ('dBuZ', 'dBZ', 'dBZc', 'dBuZv', 'dBZv', 'dBZvc'):
            refl_field = get_fieldname_pyart(datatype)
        if datatype in ('W', 'Wv', 'Wu', 'Wvu'):
            width_field = get_fieldname_pyart(datatype)

    if width_field is None or refl_field is None:
        warn('Reflectivity and spectrum width fields required'
             ' to estimate turbulence')
        return None, None

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if width_field not in radar.fields or refl_field not in radar.fields:
        warn('Unable to compute turbulence. Missing data')
        return None, None

    # user defined parameters
    radius = dscfg.get('radius', 2.)
    split_cut = dscfg.get('split_cut', False)
    xran = dscfg.get('xran', [-100., 100.])
    yran = dscfg.get('yran', [-100., 100.])
    max_split_cut = dscfg.get('max_split_cut', 2)
    beamwidth = dscfg.get('beamwidth', None)
    verbose = dscfg.get('verbose', False)
    compute_gate_pos = dscfg.get('compute_gate_pos', False)

    if beamwidth is None:
        if (radar.instrument_parameters is not None and
                'radar_beam_width_h' in radar.instrument_parameters):
            beamwidth = radar.instrument_parameters[
                'radar_beam_width_h']['data'][0]
        else:
            warn('Unknown radar beamwidth. Default 1 deg will be used')
            beamwidth = 1

    rng_res = radar.range['data'][1]-radar.range['data'][0]/1000.

    radar_out = deepcopy(radar)
    radar_out.fields = dict()
    radar_out.add_field(refl_field, deepcopy(radar.fields[refl_field]))
    radar_out.add_field(width_field, deepcopy(radar.fields[width_field]))
    radar_out.fields[refl_field]['data'][
        np.ma.getmaskarray(radar_out.fields[refl_field]['data'])] = -32768
    radar_out.fields[width_field]['data'][
        np.ma.getmaskarray(radar_out.fields[width_field]['data'])] = -32768

    radar_out.fields[refl_field]['_FillValue'] = -32768
    radar_out.fields[width_field]['_FillValue'] = -32768

    if radar_out.scan_type == 'ppi':
        pytda.calc_turb_vol(
            radar_out, radius=radius, split_cut=split_cut, xran=xran,
            yran=yran, verbose=verbose, name_dz=refl_field,
            name_sw=width_field, turb_name='turbulence',
            max_split_cut=max_split_cut, use_ntda=True, beamwidth=beamwidth,
            gate_spacing=rng_res, compute_gate_pos=compute_gate_pos)
    elif radar_out.scan_type == 'rhi':
        pytda.calc_turb_rhi(
            radar_out, radius=radius, verbose=verbose, name_dz=refl_field,
            name_sw=width_field, turb_name='turbulence',
            use_ntda=True, beamwidth=beamwidth, gate_spacing=rng_res,
            compute_gate_pos=compute_gate_pos)
    else:
        warn('Radar volume of type '+radar_out.scan_type +
             '. Only volumes of type PPI or RHI are allowed')
        return None, None

    del radar_out.fields[refl_field]
    del radar_out.fields[width_field]

    radar_out.fields['turbulence']['data'] = np.ma.masked_values(
        radar_out.fields['turbulence']['data'], -32768.)

    # prepare for exit
    new_dataset = {'radar_out': radar_out}

    return new_dataset, ind_rad


def process_dealias_fourdd(procstatus, dscfg, radar_list=None):
    """
    Dealiases the Doppler velocity field using the 4DD technique
    from Curtis and Houze, 2001

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The input data type
        filt : int. Dataset keyword
            Flag controlling Bergen and Albers filter, 1 = yes, 0 = no.
        sign : int. Dataset keyword
            Sign convention which the radial velocities in the volume created
            from the sounding data will will. This should match the
            convention used in the radar data. A value of 1 represents when
            positive values velocities are towards the radar, -1 represents
            when negative velocities are towards the radar.


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

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
    vel_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if vel_field not in radar.fields:
        warn('Unable to correct Doppler aliasing. Missing data')
        return None, None

    corr_vel_field = 'dealiased_'+vel_field

    if not dscfg['initialized']:
        # Use phase unwraping method to obtain first guess

        # get user parameters
        interval_splits = dscfg.get('interval_splits', 3)
        skip_between_rays = dscfg.get('skip_between_rays', 100)
        skip_along_ray = dscfg.get('skip_along_ray', 100)
        centered = dscfg.get('centered', True)

        corr_vel_dict = pyart.correct.dealias_region_based(
            radar, ref_vel_field=None, interval_splits=interval_splits,
            interval_limits=None, skip_between_rays=skip_between_rays,
            skip_along_ray=skip_along_ray, centered=centered,
            nyquist_vel=None, check_nyquist_uniform=True, gatefilter=False,
            rays_wrap_around=None, keep_original=False, set_limits=False,
            vel_field=vel_field, corr_vel_field=corr_vel_field)

        dscfg['initialized'] = 1
    else:
        # get user parameters
        filt = dscfg.get('filt', 1)
        sign = dscfg.get('sign', 1)
        keep_mask = dscfg.get('keep_mask', True)

        last_radar = dscfg['global_data']

        corr_vel_dict = pyart.correct.dealias_fourdd(
            radar, last_radar=last_radar, sonde_profile=None,
            gatefilter=False, filt=filt, rsl_badval=131072.0,
            keep_original=False, set_limits=False, vel_field=vel_field,
            corr_vel_field=corr_vel_field, last_vel_field=corr_vel_field,
            debug=False, sign=sign)

        if keep_mask:
            mask = np.ma.getmaskarray(radar.fields[vel_field]['data'])
            corr_vel_dict['data'] = np.ma.masked_where(
                mask, corr_vel_dict['data'])

    # prepare for exit
    radar_out = deepcopy(radar)
    radar_out.fields = dict()
    radar_out.add_field(corr_vel_field, corr_vel_dict)
    new_dataset = {'radar_out': radar_out}

    # keep current corrected Doppler velocity field in memory
    dscfg['global_data'] = radar_out

    return new_dataset, ind_rad


def process_dealias_region_based(procstatus, dscfg, radar_list=None):
    """
    Dealiases the Doppler velocity field using a region based algorithm

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The input data type
        interval_splits : int, optional
            Number of segments to split the nyquist interval into when finding
            regions of similar velocity.  More splits creates a larger number
            of initial regions which takes longer to process but may result in
            better dealiasing.  The default value of 3 seems to be a good
            compromise between performance and artifact free dealiasing. This
            value is not used if the interval_limits parameter is not None.
        skip_between_rays, skip_along_ray : int, optional
            Maximum number of filtered gates to skip over when joining
            regions, gaps between region larger than this will not be
            connected. Parameters specify the maximum number of filtered gates
            between and along a ray. Set these parameters to 0 to disable
            unfolding across filtered gates.
        centered : bool, optional
            True to apply centering to each sweep after the dealiasing
            algorithm so that the average number of unfolding is near 0. False
            does not apply centering which may results in individual sweeps
            under or over folded by the nyquist interval.
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

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
    vel_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if vel_field not in radar.fields:
        warn('Unable to correct Doppler aliasing. Missing data')
        return None, None

    corr_vel_field = 'dealiased_'+vel_field

    # get user parameters
    interval_splits = dscfg.get('interval_splits', 3)
    skip_between_rays = dscfg.get('skip_between_rays', 100)
    skip_along_ray = dscfg.get('skip_along_ray', 100)
    centered = dscfg.get('centered', True)

    corr_vel_dict = pyart.correct.dealias_region_based(
        radar, ref_vel_field=None, interval_splits=interval_splits,
        interval_limits=None, skip_between_rays=skip_between_rays,
        skip_along_ray=skip_along_ray, centered=centered,
        nyquist_vel=None, check_nyquist_uniform=True, gatefilter=False,
        rays_wrap_around=None, keep_original=False, set_limits=False,
        vel_field=vel_field, corr_vel_field=corr_vel_field)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(corr_vel_field, corr_vel_dict)

    return new_dataset, ind_rad


def process_dealias_unwrap_phase(procstatus, dscfg, radar_list=None):
    """
    Dealiases the Doppler velocity field using multi-dimensional phase
    unwrapping

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The input data type
        unwrap_unit : {'ray', 'sweep', 'volume'}, optional
            Unit to unwrap independently.  'ray' will unwrap each ray
            individually, 'sweep' each sweep, and 'volume' will unwrap the
            entire volume in a single pass.  'sweep', the default, often gives
            superior results when the lower sweeps of the radar volume are
            contaminated by clutter. 'ray' does not use the gatefilter
            parameter and rays where gates ared masked will result in poor
            dealiasing for that ray.
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

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
    vel_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if vel_field not in radar.fields:
        warn('Unable to correct Doppler aliasing. Missing data')
        return None, None

    corr_vel_field = 'dealiased_'+vel_field

    # get user parameters
    unwrap_unit = dscfg.get('unwrap_unit', 'sweep')

    corr_vel_dict = pyart.correct.dealias_unwrap_phase(
        radar, unwrap_unit=unwrap_unit, nyquist_vel=None,
        check_nyquist_uniform=True, gatefilter=False,
        rays_wrap_around=None, keep_original=False, set_limits=False,
        vel_field=vel_field, corr_vel_field=corr_vel_field, skip_checks=False)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(corr_vel_field, corr_vel_dict)

    return new_dataset, ind_rad


def process_radial_velocity(procstatus, dscfg, radar_list=None):
    """
    Estimates the radial velocity respect to the radar from the wind velocity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The input data type
        latitude, longitude : float
            arbitrary coordinates [deg] from where to compute the radial
            velocity. If any of them is None it will be the radar position
        altitude : float
            arbitrary altitude [m MSL] from where to compute the radial
            velocity. If None it will be the radar altitude
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

    v_speed_field = None
    h_speed_field = None
    h_dir_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'wind_vel_v':
            v_speed_field = get_fieldname_pyart(datatype)
        if datatype == 'WIND_SPEED':
            h_speed_field = get_fieldname_pyart(datatype)
        if datatype == 'WIND_DIRECTION':
            h_dir_field = get_fieldname_pyart(datatype)

    if h_speed_field is None or h_dir_field is None:
        warn('Horizontal wind speed and direction fields required'
             ' to estimate radial velocity')
        return None, None

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if h_speed_field not in radar.fields or h_dir_field not in radar.fields:
        warn('Unable to estimate radial velocity. '
             'Missing horizontal wind')
        return None, None

    h_speed = radar.fields[h_speed_field]['data']
    h_dir = radar.fields[h_dir_field]['data']

    if v_speed_field is None or v_speed_field not in radar.fields:
        warn('Unknown vertical wind speed. Assumed 0')

        if v_speed_field is None:
            v_speed_field == 'vertical_wind_component'

        v_speed = np.ma.zeros((radar.nrays, radar.ngates))
    else:
        v_speed = radar.fields[v_speed_field]['data']

    # user defined parameters
    lat = dscfg.get('latitude', None)
    lon = dscfg.get('longitude', None)
    alt = dscfg.get('altitude', None)

    # get u and v wind components
    h_dir_rad = np.deg2rad(h_dir)
    speed_h_u = h_speed*np.sin(h_dir_rad)  # eastward component
    speed_h_v = h_speed*np.cos(h_dir_rad)  # northward component

    if lat is not None or lon is not None or alt is not None:
        # get antenna coordinates respect to new radar location
        if lat is None:
            lat = radar.latitude['data'][0]
        if lon is None:
            lon = radar.longitude['data'][0]
        if alt is None:
            alt = radar.altitude['data'][0]

        x, y = pyart.core.geographic_to_cartesian_aeqd(
            radar.gate_longitude['data'], radar.gate_latitude['data'], lon,
            lat)
        z = radar.gate_altitude['data'] - alt
        _, azimuths, elevations = pyart.core.cartesian_to_antenna(
            x, y, z)
        azi_2D_rad = np.deg2rad(azimuths)
        ele_2D_rad = np.deg2rad(elevations)
    else:
        azi_2D_rad = np.broadcast_to(
            np.deg2rad(radar.azimuth['data'])[:, np.newaxis],
            (radar.nrays, radar.ngates))
        ele_2D_rad = np.broadcast_to(
            np.deg2rad(radar.elevation['data'])[:, np.newaxis],
            (radar.nrays, radar.ngates))

    r_speed = pyart.config.get_metadata('velocity')

    # assuming no vertical velocity
    # r_speed['data'] = h_speed*np.cos(h_dir_rad-azi_2D_rad)*np.cos(ele_2D_rad)

    # with vertical velocity included
    r_speed['data'] = (
        (speed_h_u*np.sin(azi_2D_rad)+speed_h_v*np.cos(azi_2D_rad)) *
        np.cos(ele_2D_rad)+np.sin(ele_2D_rad)*v_speed)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field('velocity', r_speed)

    return new_dataset, ind_rad


def process_wind_vel(procstatus, dscfg, radar_list=None):
    """
    Estimates the horizontal or vertical component of the wind from the
    radial velocity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The input data type
        vert_proj : Boolean
            If true the vertical projection is computed. Otherwise the
            horizontal projection is computed
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

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
    vel_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if vel_field not in radar.fields:
        warn('Unable to retrieve wind speed. Missing data')
        return None, None

    vert_proj = dscfg.get('vert_proj', False)
    wind_field = 'azimuthal_horizontal_wind_component'
    if vert_proj:
        wind_field = 'vertical_wind_component'

    wind = pyart.retrieve.est_wind_vel(
        radar, vert_proj=vert_proj, vel_field=vel_field,
        wind_field=wind_field)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(wind_field, wind)

    return new_dataset, ind_rad


def process_windshear(procstatus, dscfg, radar_list=None):
    """
    Estimates the wind shear from the wind velocity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The input data type
        az_tol : float
            The tolerance in azimuth when looking for gates on top
            of the gate when computation is performed

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

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
    wind_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if wind_field not in radar.fields:
        warn('Unable to retrieve wind shear. Missing data')
        return None, None

    az_tol = dscfg.get('az_tol', 0.5)
    windshear_field = 'vertical_wind_shear'

    windshear = pyart.retrieve.est_vertical_windshear(
        radar, az_tol=az_tol, wind_field=wind_field,
        windshear_field=windshear_field)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(windshear_field, windshear)

    return new_dataset, ind_rad


def process_vad(procstatus, dscfg, radar_list=None):
    """
    Estimates vertical wind profile using the VAD (velocity Azimuth Display)
    technique

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The input data type
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

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
    vel_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if vel_field not in radar.fields:
        warn('Unable to retrieve wind speed. Missing data')
        return None, None

    # User defined parameters
    npoints_min = dscfg.get('npoints_min', 6)
    azi_spacing_max = dscfg.get('azi_spacing_max', 45.)
    vel_diff_max = dscfg.get('vel_diff_max', 10.)

    (u_vel_dict, v_vel_dict, w_vel_dict, vel_est_dict, vel_std_dict,
     vel_diff_dict) = pyart.retrieve.est_wind_profile(
         radar, npoints_min=npoints_min, azi_spacing_max=azi_spacing_max,
         vel_diff_max=vel_diff_max, rad_vel_field=vel_field,
         u_vel_field='eastward_wind_component',
         v_vel_field='northward_wind_component',
         w_vel_field='vertical_wind_component',
         vel_est_field='retrieved_velocity',
         vel_std_field='retrieved_velocity_std',
         vel_diff_field='velocity_difference')

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field('eastward_wind_component', u_vel_dict)
    new_dataset['radar_out'].add_field('northward_wind_component', v_vel_dict)
    new_dataset['radar_out'].add_field('vertical_wind_component', w_vel_dict)
    new_dataset['radar_out'].add_field('retrieved_velocity', vel_est_dict)
    new_dataset['radar_out'].add_field('retrieved_velocity_std', vel_std_dict)
    new_dataset['radar_out'].add_field('velocity_difference', vel_diff_dict)

    return new_dataset, ind_rad

def process_dda(procstatus, dscfg, radar_list=None):
    """
    Estimates horizontal wind speed and direction with a multi-doppler approach
    This method uses the python package pyDDA

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The input data type
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """
    if not _PYDDA_AVAILABLE:
        warn('PyDDA package not available. Unable to compute DDA')
        return None, None

    if procstatus != 1:
        return None, None

    import pdb; pdb.set_trace()
    # check how many radars are there
    radarnr_dict = dict()
    ind_radar_list = set()
    for datatypedescr in dscfg['datatype']:
        radarnr = datatypedescr.split(':')[0]
        radarnr_dict.update({radarnr: []})
        ind_radar_list.add(int(radarnr[5:8])-1)

    ind_radar_list = list(ind_radar_list)

    if len(radarnr_dict) < 2:
        warn('DDA requires data from at least two different radars')
        return None, None

    # create the list of data types for each radar
    radars = []
    for ind_radar in ind_radar_list:
        if radar_list[ind_radar_list] is None:
            warn('Unable to use DDA, one of the radar objects is empty')
            return None, None

        radars.append(radar_list[ind_radar_list])

    import pdb; pdb.set_trace()
    return new_dataset, ind_rad

