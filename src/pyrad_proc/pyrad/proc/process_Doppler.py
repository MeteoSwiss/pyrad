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
    process_windshear_lidar
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

from ..io.io_aux import (
    get_datatype_fields,
    get_fieldname_pyart,
    convert_pydda_to_pyart_grid,
)
from ..io import read_radiosounding_igra
from .process_grid import process_grid
from ..util import compute_average_vad


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
            The input data type, must contain,
            "dBuZ" or "dBZ" or "dBZc" or "dBuZv" or "dBZv" or "dBZvc" or "CNRc", and,
            "W" or "Wv" or "Wu" or "Wvu" or "WD" or "WDc"
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
        use_ntda : Bool. Dataset keyword
            Wether to use NCAR Turbulence Detection Algorithm (NTDA). Default
            True
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
        dictionary containing the output field "EDR"
    ind_rad : int
        radar index

    """
    if not _PYTDA_AVAILABLE:
        warn("PyTDA package not available. Unable to compute turbulence")
        return None, None

    if procstatus != 1:
        return None, None

    width_field = None
    refl_field = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("dBuZ", "dBZ", "dBZc", "dBuZv", "dBZv", "dBZvc", "CNRc"):
            refl_field = get_fieldname_pyart(datatype)
        if datatype in ("W", "Wv", "Wu", "Wvu", "WD", "WDc"):
            width_field = get_fieldname_pyart(datatype)

    if width_field is None or refl_field is None:
        warn(
            "Reflectivity and spectrum width fields required" " to estimate turbulence"
        )
        return None, None

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if width_field not in radar.fields or refl_field not in radar.fields:
        warn("Unable to compute turbulence. Missing data")
        return None, None

    # user defined parameters
    radius = dscfg.get("radius", 2.0)
    split_cut = dscfg.get("split_cut", False)
    xran = dscfg.get("xran", [-100.0, 100.0])
    yran = dscfg.get("yran", [-100.0, 100.0])
    max_split_cut = dscfg.get("max_split_cut", 2)
    use_ntda = dscfg.get("use_ntda", True)
    beamwidth = dscfg.get("beamwidth", None)
    verbose = dscfg.get("verbose", False)
    compute_gate_pos = dscfg.get("compute_gate_pos", False)

    if beamwidth is None:
        if (
            radar.instrument_parameters is not None
            and "radar_beam_width_h" in radar.instrument_parameters
        ):
            beamwidth = radar.instrument_parameters["radar_beam_width_h"]["data"][0]
        else:
            warn("Unknown radar beamwidth. Default 1 deg will be used")
            beamwidth = 1

    rng_res = (radar.range["data"][1] - radar.range["data"][0]) / 1000.0

    radar_out = deepcopy(radar)
    radar_out.fields = dict()
    radar_out.add_field(refl_field, deepcopy(radar.fields[refl_field]))
    radar_out.add_field(width_field, deepcopy(radar.fields[width_field]))
    radar_out.fields[refl_field]["data"][
        np.ma.getmaskarray(radar_out.fields[refl_field]["data"])
    ] = -32768
    radar_out.fields[width_field]["data"][
        np.ma.getmaskarray(radar_out.fields[width_field]["data"])
    ] = -32768

    radar_out.fields[refl_field]["_FillValue"] = -32768
    radar_out.fields[width_field]["_FillValue"] = -32768

    if radar_out.scan_type == "ppi":
        pytda.calc_turb_vol(
            radar_out,
            radius=radius,
            split_cut=split_cut,
            xran=xran,
            yran=yran,
            verbose=verbose,
            name_dz=refl_field,
            name_sw=width_field,
            turb_name="turbulence",
            max_split_cut=max_split_cut,
            use_ntda=use_ntda,
            beamwidth=beamwidth,
            gate_spacing=rng_res,
            compute_gate_pos=compute_gate_pos,
        )
    elif radar_out.scan_type == "rhi":
        pytda.calc_turb_rhi(
            radar_out,
            radius=radius,
            verbose=verbose,
            name_dz=refl_field,
            name_sw=width_field,
            turb_name="turbulence",
            use_ntda=use_ntda,
            beamwidth=beamwidth,
            gate_spacing=rng_res,
        )
    else:
        warn(
            "Radar volume of type "
            + radar_out.scan_type
            + ". Only volumes of type PPI or RHI are allowed"
        )
        return None, None

    del radar_out.fields[refl_field]
    del radar_out.fields[width_field]

    radar_out.fields["turbulence"]["data"] = np.ma.masked_values(
        radar_out.fields["turbulence"]["data"], -32768.0
    )

    # prepare for exit
    new_dataset = {"radar_out": radar_out}

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
            The input data type, must contain
            "V" or "Vc"
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
        dictionary containing the output field "dealV" or "dealVc" (if Vc was provided)
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg["datatype"][0])
    vel_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if vel_field not in radar.fields:
        warn("Unable to correct Doppler aliasing. Missing velocity field")
        warn(f"Radar available fields are {radar.fields}")
        return None, None

    corr_vel_field = "dealiased_" + vel_field

    if not dscfg["initialized"]:
        # Use phase unwraping method to obtain first guess

        # get user parameters
        interval_splits = dscfg.get("interval_splits", 3)
        skip_between_rays = dscfg.get("skip_between_rays", 100)
        skip_along_ray = dscfg.get("skip_along_ray", 100)
        centered = dscfg.get("centered", True)

        corr_vel_dict = pyart.correct.dealias_region_based(
            radar,
            ref_vel_field=None,
            interval_splits=interval_splits,
            interval_limits=None,
            skip_between_rays=skip_between_rays,
            skip_along_ray=skip_along_ray,
            centered=centered,
            nyquist_vel=None,
            check_nyquist_uniform=True,
            gatefilter=False,
            rays_wrap_around=None,
            keep_original=False,
            set_limits=False,
            vel_field=vel_field,
            corr_vel_field=corr_vel_field,
        )

        dscfg["initialized"] = 1
    else:
        # get user parameters
        filt = dscfg.get("filt", 1)
        sign = dscfg.get("sign", 1)
        keep_mask = dscfg.get("keep_mask", True)

        last_radar = dscfg["global_data"]

        corr_vel_dict = pyart.correct.dealias_fourdd(
            radar,
            last_radar=last_radar,
            sonde_profile=None,
            gatefilter=False,
            filt=filt,
            rsl_badval=131072.0,
            keep_original=False,
            set_limits=False,
            vel_field=vel_field,
            corr_vel_field=corr_vel_field,
            last_vel_field=corr_vel_field,
            debug=False,
            sign=sign,
        )

        if keep_mask:
            mask = np.ma.getmaskarray(radar.fields[vel_field]["data"])
            corr_vel_dict["data"] = np.ma.masked_where(mask, corr_vel_dict["data"])

    # prepare for exit
    radar_out = deepcopy(radar)
    radar_out.fields = dict()
    radar_out.add_field(corr_vel_field, corr_vel_dict)
    new_dataset = {"radar_out": radar_out}

    # keep current corrected Doppler velocity field in memory
    dscfg["global_data"] = radar_out

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
            The input data type, must contain,
            "V" or "Vc"
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
        nyquist_vel : float, optional
            Nyquist velocity of the aquired radar velocity.
            Usually this parameter is provided in the
            Radar object intrument_parameters. If it is not available it can
            be provided as a keyword here.
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "dealV" or "dealVc" (if Vc was provided)
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg["datatype"][0])
    vel_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if vel_field not in radar.fields:
        warn("Unable to correct Doppler aliasing. Missing velocity field")
        warn(f"Radar available fields are {radar.fields}")
        return None, None

    corr_vel_field = "dealiased_" + vel_field

    # get user parameters
    interval_splits = dscfg.get("interval_splits", 3)
    skip_between_rays = dscfg.get("skip_between_rays", 100)
    skip_along_ray = dscfg.get("skip_along_ray", 100)
    centered = dscfg.get("centered", True)
    nyquist_vel = dscfg.get("nyquist_vel", None)

    corr_vel_dict = pyart.correct.dealias_region_based(
        radar,
        ref_vel_field=None,
        interval_splits=interval_splits,
        interval_limits=None,
        skip_between_rays=skip_between_rays,
        skip_along_ray=skip_along_ray,
        centered=centered,
        nyquist_vel=nyquist_vel,
        check_nyquist_uniform=True,
        gatefilter=False,
        rays_wrap_around=None,
        keep_original=False,
        set_limits=False,
        vel_field=vel_field,
        corr_vel_field=corr_vel_field,
    )

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(radar)}
    new_dataset["radar_out"].fields = dict()
    new_dataset["radar_out"].add_field(corr_vel_field, corr_vel_dict)

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
            The input data type, must contain,
            "V" or "Vc"
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
        dictionary containing the output field "dealV" or "dealVc" (if Vc was provided)
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg["datatype"][0])
    vel_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if vel_field not in radar.fields:
        warn("Unable to correct Doppler aliasing. Missing data")
        return None, None

    corr_vel_field = "dealiased_" + vel_field

    # get user parameters
    unwrap_unit = dscfg.get("unwrap_unit", "sweep")

    corr_vel_dict = pyart.correct.dealias_unwrap_phase(
        radar,
        unwrap_unit=unwrap_unit,
        nyquist_vel=None,
        check_nyquist_uniform=True,
        gatefilter=False,
        rays_wrap_around=None,
        keep_original=False,
        set_limits=False,
        vel_field=vel_field,
        corr_vel_field=corr_vel_field,
        skip_checks=False,
    )

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(radar)}
    new_dataset["radar_out"].fields = dict()
    new_dataset["radar_out"].add_field(corr_vel_field, corr_vel_dict)

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
            The input data type, must contain
            WIND_SPEED, and,
            WIND_DIRECTION, and,
            wind_vel_v
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
        dictionary containing the output field "V"
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    v_speed_field = None
    h_speed_field = None
    h_dir_field = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == "wind_vel_v":
            v_speed_field = get_fieldname_pyart(datatype)
        if datatype == "WIND_SPEED":
            h_speed_field = get_fieldname_pyart(datatype)
        if datatype == "WIND_DIRECTION":
            h_dir_field = get_fieldname_pyart(datatype)

    if h_speed_field is None or h_dir_field is None:
        warn(
            "Horizontal wind speed and direction fields required"
            " to estimate radial velocity"
        )
        return None, None

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if h_speed_field not in radar.fields or h_dir_field not in radar.fields:
        warn("Unable to estimate radial velocity. " "Missing horizontal wind")
        return None, None

    h_speed = radar.fields[h_speed_field]["data"]
    h_dir = radar.fields[h_dir_field]["data"]

    if v_speed_field is None or v_speed_field not in radar.fields:
        warn("Unknown vertical wind speed. Assumed 0")

        if v_speed_field is None:
            v_speed_field == "vertical_wind_component"

        v_speed = np.ma.zeros((radar.nrays, radar.ngates))
    else:
        v_speed = radar.fields[v_speed_field]["data"]

    # user defined parameters
    lat = dscfg.get("latitude", None)
    lon = dscfg.get("longitude", None)
    alt = dscfg.get("altitude", None)

    # get u and v wind components
    h_dir_rad = np.deg2rad(h_dir)
    speed_h_u = h_speed * np.sin(h_dir_rad)  # eastward component
    speed_h_v = h_speed * np.cos(h_dir_rad)  # northward component

    if lat is not None or lon is not None or alt is not None:
        # get antenna coordinates respect to new radar location
        if lat is None:
            lat = radar.latitude["data"][0]
        if lon is None:
            lon = radar.longitude["data"][0]
        if alt is None:
            alt = radar.altitude["data"][0]

        x, y = pyart.core.geographic_to_cartesian_aeqd(
            radar.gate_longitude["data"], radar.gate_latitude["data"], lon, lat
        )
        z = radar.gate_altitude["data"] - alt
        _, azimuths, elevations = pyart.core.cartesian_to_antenna(x, y, z)
        azi_2D_rad = np.deg2rad(azimuths)
        ele_2D_rad = np.deg2rad(elevations)
    else:
        azi_2D_rad = np.broadcast_to(
            np.deg2rad(radar.azimuth["data"])[:, np.newaxis],
            (radar.nrays, radar.ngates),
        )
        ele_2D_rad = np.broadcast_to(
            np.deg2rad(radar.elevation["data"])[:, np.newaxis],
            (radar.nrays, radar.ngates),
        )

    r_speed = pyart.config.get_metadata("velocity")

    # assuming no vertical velocity
    # r_speed['data'] = h_speed*np.cos(h_dir_rad-azi_2D_rad)*np.cos(ele_2D_rad)

    # with vertical velocity included
    r_speed["data"] = (
        speed_h_u * np.sin(azi_2D_rad) + speed_h_v * np.cos(azi_2D_rad)
    ) * np.cos(ele_2D_rad) + np.sin(ele_2D_rad) * v_speed

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(radar)}
    new_dataset["radar_out"].fields = dict()
    new_dataset["radar_out"].add_field("velocity", r_speed)

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
            The input data type, must contain
            "V" or "Vc"
        vert_proj : Boolean
            If true the vertical projection is computed. Otherwise the
            horizontal projection is computed
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field
        "wind_vel_h_az", (if vert_proj is False), or,
        "wind_vel_v" (if vert_proj is True)
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg["datatype"][0])
    vel_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if vel_field not in radar.fields:
        warn("Unable to retrieve wind speed. Missing data")
        return None, None

    vert_proj = dscfg.get("vert_proj", False)
    wind_field = "azimuthal_horizontal_wind_component"
    if vert_proj:
        wind_field = "vertical_wind_component"

    wind = pyart.retrieve.est_wind_vel(
        radar, vert_proj=vert_proj, vel_field=vel_field, wind_field=wind_field
    )

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(radar)}
    new_dataset["radar_out"].fields = dict()
    new_dataset["radar_out"].add_field(wind_field, wind)

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
        dictionary containing the output field "windshear_v"
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg["datatype"][0])
    wind_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if wind_field not in radar.fields:
        warn("Unable to retrieve wind shear. Missing data")
        return None, None

    az_tol = dscfg.get("az_tol", 0.5)
    windshear_field = "vertical_wind_shear"

    windshear = pyart.retrieve.est_vertical_windshear(
        radar, az_tol=az_tol, wind_field=wind_field, windshear_field=windshear_field
    )

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(radar)}
    new_dataset["radar_out"].fields = dict()
    new_dataset["radar_out"].add_field(windshear_field, windshear)

    return new_dataset, ind_rad


def process_windshear_lidar(procstatus, dscfg, radar_list=None):
    """
    Estimates the wind shear from the wind velocity of lidar scans

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The input data type, must contain
            "V" or "Vc"
        az_tol : float
            The tolerance in azimuth when looking for gates on top
            of the gate when computation is performed

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "windshear_v"
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg["datatype"][0])
    wind_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if wind_field not in radar.fields:
        warn("Unable to retrieve wind shear. Missing data")
        return None, None

    az_tol = dscfg.get("az_tol", 0.5)
    windshear_field = "vertical_wind_shear"

    windshear = pyart.retrieve.est_vertical_windshear_lidar(
        radar, az_tol=az_tol, wind_field=wind_field, windshear_field=windshear_field
    )

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(radar)}
    new_dataset["radar_out"].fields = dict()
    new_dataset["radar_out"].add_field(windshear_field, windshear)

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
            The input data type, must contain
            "V" or "Vc"
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output fields
        "wind_vel_h_u", "wind_vel_h_v", "wind_vel_v",
        "estV", "stdV", and "diffV"
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg["datatype"][0])
    vel_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if vel_field not in radar.fields:
        warn("Unable to retrieve wind speed. Missing data")
        return None, None

    # User defined parameters
    npoints_min = dscfg.get("npoints_min", 6)
    azi_spacing_max = dscfg.get("azi_spacing_max", 45.0)
    vel_diff_max = dscfg.get("vel_diff_max", 10.0)

    (u_vel_dict, v_vel_dict, w_vel_dict, vel_est_dict, vel_std_dict, vel_diff_dict) = (
        pyart.retrieve.est_wind_profile(
            radar,
            npoints_min=npoints_min,
            azi_spacing_max=azi_spacing_max,
            vel_diff_max=vel_diff_max,
            rad_vel_field=vel_field,
            u_vel_field="eastward_wind_component",
            v_vel_field="northward_wind_component",
            w_vel_field="vertical_wind_component",
            vel_est_field="retrieved_velocity",
            vel_std_field="retrieved_velocity_std",
            vel_diff_field="velocity_difference",
        )
    )

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(radar)}
    new_dataset["radar_out"].fields = dict()
    new_dataset["radar_out"].add_field("eastward_wind_component", u_vel_dict)
    new_dataset["radar_out"].add_field("northward_wind_component", v_vel_dict)
    new_dataset["radar_out"].add_field("vertical_wind_component", w_vel_dict)
    new_dataset["radar_out"].add_field("retrieved_velocity", vel_est_dict)
    new_dataset["radar_out"].add_field("retrieved_velocity_std", vel_std_dict)
    new_dataset["radar_out"].add_field("velocity_difference", vel_diff_dict)

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
            The input data type, must contain
            "V" or "Vc", and,
            "dBuZ", "dBZ", or "dBZc"

        gridconfig : dictionary. Dataset keyword
            Dictionary containing some or all of this keywords:
            xmin, xmax, ymin, ymax, zmin, zmax : floats
                minimum and maximum horizontal distance from grid origin [km]
                and minimum and maximum vertical distance from grid origin [m]
                Defaults -40, 40, -40, 40, 0., 10000.
            latmin, latmax, lonmin, lonmax : floats
                minimum and maximum latitude and longitude [deg], if specified
                xmin, xmax, ymin, ymax will be ignored
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
        signs : list of integers
            The sign of the velocity field for every radar object.
            A value of 1 represents when
            positive values velocities are towards the radar, -1 represents
            when negative velocities are towards the radar.
        Co : float
            Weight for cost function related to observed radial velocities.
            Default: 1.
        Cm : float
            Weight for cost function related to the mass continuity equation.
            Default: 1500.
        Cx: float
            Smoothing coefficient for x-direction
        Cy: float
            Smoothing coefficient for y-direction
        Cz: float
            Smoothing coefficient for z-direction
        Cb: float
            Coefficient for sounding constraint
        Cv: float
            Weight for cost function related to vertical vorticity equation.
        Cmod: float
            Coefficient for model constraint
        Cpoint: float
            Coefficient for point constraint
        wind_tol: float
            Stop iterations after maximum change in winds is less than this
            value.
        frz : float
            The freezing level in meters. This is to tell PyDDA where to use
            ice particle fall speeds in the wind retrieval verus liquid.
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output fields
            "wind_vel_h_u", "wind_vel_h_v" and "wind_vel_v"
    ind_rad : int
        radar index

    """
    if not _PYDDA_AVAILABLE:
        warn("PyDDA and xarray package not available. Unable to compute wind fields")
        return None, None

    if procstatus != 1:
        return None, None

    if len(radar_list) < 2:
        warn("DDA requires data from at least two different radars")
        return None, None

    # check how many radars have to be compared and which datatype to use
    ind_radar_list = []
    field_name_list = []
    datatype_list = []
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_name = get_fieldname_pyart(datatype)
        ind_radar_list.append(int(radarnr[5:8]) - 1)
        datatype_list.append(datatype)
        field_name_list.append(field_name)
    ind_radar_list = np.array(ind_radar_list)
    datatype_list = np.array(datatype_list)
    field_name_list = np.array(field_name_list)

    # Get DDA numerical parameters
    Co = dscfg.get("Co", 1.0)
    Cm = dscfg.get("Cm", 1500.0)
    Cx = dscfg.get("Cx", 0.0)
    Cy = dscfg.get("Cy", 0.0)
    Cz = dscfg.get("Cz", 0.0)
    Cb = dscfg.get("Cb", 0.0)
    Cv = dscfg.get("Cv", 0.0)
    frz = dscfg.get("frz", 4500.0)
    Cmod = dscfg.get("Cmod", 0.0)
    Cpoint = dscfg.get("Cpoint", 0.0)
    wind_tol = dscfg.get("wind_tol", 0.1)

    signs = dscfg.get("signs", [1] * len(set(ind_radar_list)))
    sounding = dscfg.get("sounding", None)
    if sounding:
        dtime = pyart.graph.common.generate_radar_time_begin(radar_list[0])
        sounding_data = read_radiosounding_igra(sounding, dtime)
    else:
        sounding_data = None

    if sounding_data is not None:
        # Remove nan from souding
        sounding_data = sounding_data.dropna()
        # create wind profile from sounding
        wind_prof = pyart.core.HorizontalWindProfile(
            sounding_data["GPH"], sounding_data["WSPD"], sounding_data["WDIR"]
        )
    else:
        # Compute VAD for every radar to initialize DDA
        # z-vector for vad
        z_want = np.arange(
            dscfg["gridConfig"]["zmin"],
            dscfg["gridConfig"]["zmax"] + dscfg["gridConfig"]["vres"],
            dscfg["gridConfig"]["vres"],
        )

        wind_prof = compute_average_vad(
            radar_list,
            z_want,
            signs,
            dscfg["gridConfig"]["lonorig"],
            dscfg["gridConfig"]["latorig"],
        )

    # Get name of reflectivity and velocity fields for each radar
    refl_fields = []
    vel_fields = []
    for i, radar in enumerate(radar_list):
        # Find vel and refl fields
        field_names_rad = field_name_list[ind_radar_list == i]
        vel_field = field_names_rad[
            ["velocity" in vname for vname in field_name_list[ind_radar_list == i]]
        ][0]
        vel_fields.append(vel_field)
        refl_field = field_names_rad[
            ["reflectivity" in vname for vname in field_name_list[ind_radar_list == i]]
        ][0]
        refl_fields.append(refl_field)

    # Grid the variables
    grids = []
    for i, radar in enumerate(radar_list):
        # Now we grid
        dscfg_grid = deepcopy(dscfg)
        dscfg_grid["datatype"] = np.array(dscfg["datatype"])[ind_radar_list == i]
        grids.append(process_grid(1, dscfg_grid, radar_list)[0]["radar_out"])
    grids_pyDDA = []

    # Harmonize variables
    for i, grid in enumerate(grids):
        grid.fields["velocity"] = grid.fields.pop(vel_fields[i])
        if signs[i] == 1:
            grid.fields["velocity"]["data"] *= -1
        grid.fields["reflectivity"] = grid.fields.pop(refl_fields[i])
        grids_pyDDA.append(pydda.io.read_from_pyart_grid(grid))

    # DDA initialization
    grids_pyDDA[0] = pydda.initialization.make_wind_field_from_profile(
        grids_pyDDA[0], wind_prof, vel_field="velocity"
    )

    # Actual DDA computation
    new_grids, params = pydda.retrieval.get_dd_wind_field(
        grids_pyDDA,
        vel_name="velocity",
        refl_field="reflectivity",
        mask_outside_opt=True,
        engine="scipy",
        Co=Co,
        Cm=Cm,
        Cx=Cx,
        Cy=Cy,
        Cz=Cz,
        Cb=Cb,
        Cv=Cv,
        Cmod=Cmod,
        Cpoint=Cpoint,
        wind_tol=wind_tol,
        frz=frz,
    )

    for i, grid in enumerate(new_grids):
        # convert back to pyart grids
        new_grids[i] = convert_pydda_to_pyart_grid(grid)
        # add radar name
        new_grids[i].radar_name = grids[i].radar_name

    # pyDDA returns as many grid objects as input radars
    # the idea here is to merge these grid objects into one
    # and to replace the homogeneized vel and refl fields by the
    # original ones

    # create merged grid from first dda grid
    # First radar will have normal field names
    # Additional radars will have _radar1, radar_2, ..., radar_N
    # appended to their field names
    merged_grid = deepcopy(new_grids[0])
    for var in list(merged_grid.fields):
        if var in grids[0].fields:
            if var == "velocity":
                new_key = vel_fields[0]
            elif var == "reflectivity":
                new_key = refl_fields[0]
            else:
                new_key = var
            merged_grid.fields[new_key] = merged_grid.fields.pop(var)

    # add the other dda grids
    radar_metadata = []
    for i, additional_grid in enumerate(new_grids[1:]):
        for var in list(additional_grid.fields):
            if var in grids[i + 1].fields:
                if var == "velocity":
                    new_key = vel_fields[i + 1]
                elif var == "reflectivity":
                    new_key = refl_fields[i + 1]
                else:
                    new_key = var
                new_key += "_radar{:d}".format(i + 1)
                merged_grid.fields[new_key] = additional_grid.fields[var]

        # Get additional radar metadata
        mdata = {}
        mdata["radar_longitude"] = additional_grid.radar_longitude["data"]
        mdata["radar_latitude"] = additional_grid.radar_latitude["data"]
        # Try to find name of radar
        if additional_grid.radar_name["data"] != "":
            mdata["radar_name"] = additional_grid.radar_name["data"][0]
        elif dscfg["RadarName"][i + 1] != "":
            mdata["radar_name"] = dscfg["RadarName"][i + 1]
        else:
            mdata["radar_name"] = "Unknown"

        radar_metadata.append(mdata)

    # Rename DDA u, v, w fields to pyart names
    merged_grid.fields["eastward_wind_component"] = merged_grid.fields.pop("u")
    merged_grid.fields["northward_wind_component"] = merged_grid.fields.pop("v")
    merged_grid.fields["vertical_wind_component"] = merged_grid.fields.pop("w")

    # Add coordinates of additional radars in metadata
    merged_grid.metadata["additional_radars"] = radar_metadata

    new_dataset = {}
    new_dataset["radar_out"] = merged_grid
    return new_dataset, ind_radar_list
