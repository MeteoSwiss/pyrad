"""
pyrad.proc.process_spectra
==========================

Functions to processes spectral data.

.. autosummary::
    :toctree: generated/

    process_raw_spectra
    process_ifft
    process_spectra_point
    process_filter_0Doppler
    process_filter_srhohv
    process_filter_spectra_noise
    process_dealias_spectra
    process_spectra_ang_avg
    process_spectral_power
    process_spectral_noise
    process_spectral_phase
    process_spectral_reflectivity
    process_spectral_differential_reflectivity
    process_spectral_differential_phase
    process_spectral_rhohv
    process_pol_variables
    process_noise_power
    process_reflectivity
    process_differential_reflectivity
    process_differential_phase
    process_rhohv
    process_Doppler_velocity
    process_Doppler_width

"""

from copy import deepcopy
from warnings import warn
import numpy as np
from netCDF4 import num2date

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart


def process_raw_spectra(procstatus, dscfg, radar_list=None):
    """
    Dummy function that returns the initial input data set

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, _, _, _ = get_datatype_fields(datatypedescr)
        break
    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    new_dataset = {"radar_out": deepcopy(radar_list[ind_rad])}

    return new_dataset, ind_rad


def process_ifft(procstatus, dscfg, radar_list=None):
    """
    Compute the Doppler spectrum width from the spectral reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain,
            "ShhADU" or "ShhADUu", or,
            "SvvADU" or "SvvADUu", or,
            "sNADUh" or "sNADUv"
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output fields
        "IQhhADU" if "ShhADU" or "ShhADUu" were provided,
        "IQhvvDU" if "SvvADU" or "SvvADUu" were provided,
        "IQNADUh" if "sNADUh" was provided,
        "IQNADUv" if "sNADUv" was provided.
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg["datatype"][0])
    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    wind_params = dscfg.get("window", ["None"])
    if len(wind_params) == 1:
        window = wind_params[0]
        if window == "None":
            window = None
        else:
            try:
                window = float(window)
            except ValueError:
                pass
    else:
        window = wind_params
        for i in range(1, len(window)):
            window[i] = float(window[i])
        window = tuple(window)

    fields_in_list = []
    fields_out_list = []
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_name = get_fieldname_pyart(datatype)
        if field_name not in radar.fields:
            warn(field_name + " not in radar")
            continue
        if field_name in (
            "unfiltered_complex_spectra_hh_ADU",
            "complex_spectra_hh_ADU",
        ):
            fields_out_list.append("IQ_hh_ADU")
        elif field_name in (
            "unfiltered_complex_spectra_vv_ADU",
            "complex_spectra_vv_ADU",
        ):
            fields_out_list.append("IQ_vv_ADU")
        elif field_name == "spectral_noise_power_hh_ADU":
            fields_out_list.append("IQ_noise_power_hh_ADU")
        elif field_name == "spectral_noise_power_vv_ADU":
            fields_out_list.append("IQ_noise_power_vv_ADU")
        else:
            warn(field_name + " can not be inverse Fourier transformed")
        fields_in_list.append(field_name)

    radar_out = pyart.retrieve.compute_iq(
        radar, fields_in_list, fields_out_list, window=window
    )

    # prepare for exit
    new_dataset = {"radar_out": radar_out}

    return new_dataset, ind_rad


def process_spectra_point(procstatus, dscfg, radar_list=None):
    """
    Obtains the spectra or IQ data at a point location.

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
            (range, azimuth, elevation). Default False
        truealt : boolean. Dataset keyword
            if True the user input altitude is used to determine the point of
            interest.
            if False use the altitude at a given radar elevation ele over the
            point of interest. Default True
        lon : float. Dataset keyword
            the longitude [deg]. Use when latlon is True.
        lat : float. Dataset keyword
            the latitude [deg]. Use when latlon is True.
        alt : float. Dataset keyword
            altitude [m MSL]. Use when latlon is True. Default 0.
        ele : float. Dataset keyword
            radar elevation [deg]. Use when latlon is False or when latlon is
            True and truealt is False
        azi : float. Dataset keyword
            radar azimuth [deg]. Use when latlon is False
        rng : float. Dataset keyword
            range from radar [m]. Use when latlon is False
        AziTol : float. Dataset keyword
            azimuthal tolerance to determine which radar azimuth to use [deg].
            Default 0.5
        EleTol : float. Dataset keyword
            elevation tolerance to determine which radar elevation to use
            [deg]. Default 0.5
        RngTol : float. Dataset keyword
            range tolerance to determine which radar bin to use [m]. Default
            50.

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

    field_names = []
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_names.append(get_fieldname_pyart(datatype))

    ind_rad = int(radarnr[5:8]) - 1

    if procstatus == 2:
        if dscfg["initialized"] == 0:
            return None, None

        # prepare for exit
        new_dataset = {
            "radar_out": dscfg["global_data"]["psr_poi"],
            "point_coordinates_WGS84_lon_lat_alt": (
                dscfg["global_data"]["point_coordinates_WGS84_lon_lat_alt"]
            ),
            "antenna_coordinates_az_el_r": (
                dscfg["global_data"]["antenna_coordinates_az_el_r"]
            ),
            "final": True,
        }

        return new_dataset, ind_rad

    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid psr")
        return None, None
    psr = radar_list[ind_rad]

    projparams = dict()
    projparams.update({"proj": "pyart_aeqd"})
    projparams.update({"lon_0": psr.longitude["data"]})
    projparams.update({"lat_0": psr.latitude["data"]})

    truealt = dscfg.get("truealt", True)
    latlon = dscfg.get("latlon", False)

    single_point = dscfg.get("single_point", True)
    if latlon:
        lon = dscfg["lon"]
        lat = dscfg["lat"]
        alt = dscfg.get("alt", 0.0)
        dscfg.get("latlonTol", 1.0)
        dscfg.get("altTol", 100.0)

        x, y = pyart.core.geographic_to_cartesian(lon, lat, projparams)

        if not truealt:
            ke = 4.0 / 3.0  # constant for effective radius
            a = 6378100.0  # earth radius
            re = a * ke  # effective radius

            elrad = dscfg["ele"] * np.pi / 180.0
            r_ground = np.sqrt(x**2.0 + y**2.0)
            r = r_ground / np.cos(elrad)
            alt_psr = (
                psr.altitude["data"]
                + np.sqrt(r**2.0 + re**2.0 + 2.0 * r * re * np.sin(elrad))
                - re
            )
            alt_psr = alt_psr[0]
        else:
            alt_psr = alt

        r, az, el = pyart.core.cartesian_to_antenna(
            x, y, alt_psr - psr.altitude["data"]
        )
        r = r[0]
        az = az[0]
        el = el[0]
    else:
        r = dscfg["rng"]
        az = dscfg["azi"]
        el = dscfg["ele"]
        azi_tol = dscfg.get("AziTol", 0.5)
        ele_tol = dscfg.get("EleTol", 0.5)
        rng_tol = dscfg.get("RngTol", 50.0)

        x, y, alt = pyart.core.antenna_to_cartesian(r / 1000.0, az, el)
        lon, lat = pyart.core.cartesian_to_geographic(x, y, projparams)
        lon = lon[0]
        lat = lat[0]

    d_az = np.abs(psr.azimuth["data"] - az)
    if np.min(d_az) > azi_tol:
        warn(
            " No psr bin found for point (az, el, r):("
            + str(az)
            + ", "
            + str(el)
            + ", "
            + str(r)
            + "). Minimum distance to psr azimuth "
            + str(d_az)
            + " larger than tolerance"
        )
        return None, None

    d_el = np.abs(psr.elevation["data"] - el)
    if np.min(d_el) > ele_tol:
        warn(
            " No psr bin found for point (az, el, r):("
            + str(az)
            + ", "
            + str(el)
            + ", "
            + str(r)
            + "). Minimum distance to psr elevation "
            + str(d_el)
            + " larger than tolerance"
        )
        return None, None

    d_r = np.abs(psr.range["data"] - r)
    if np.min(d_r) > rng_tol:
        warn(
            " No psr bin found for point (az, el, r):("
            + str(az)
            + ", "
            + str(el)
            + ", "
            + str(r)
            + "). Minimum distance to psr range bin "
            + str(d_r)
            + " larger than tolerance"
        )
        return None, None

    if single_point:
        ind_ray = np.argmin(
            np.abs(psr.azimuth["data"] - az) + np.abs(psr.elevation["data"] - el)
        )
    else:
        ind_ray = np.where(np.logical_and(d_az <= azi_tol, d_el <= ele_tol))[0]
    ind_rng = np.argmin(np.abs(psr.range["data"] - r))
    nrays = ind_ray.size

    time_poi = num2date(
        psr.time["data"][ind_ray], psr.time["units"], psr.time["calendar"]
    )
    if nrays > 1:
        time_ref = time_poi[0]
    else:
        time_ref = time_poi
        time_poi = np.array([time_poi])

    # initialize dataset
    if not dscfg["initialized"]:
        psr_poi = deepcopy(psr)

        # prepare space for field
        psr_poi.fields = dict()
        for field_name in field_names:
            if field_name not in psr.fields:
                warn("Field " + field_name + " not in psr object")
                return None, None
            psr_poi.add_field(field_name, deepcopy(psr.fields[field_name]))
            psr_poi.fields[field_name]["data"] = np.array([])

        # fixed psr objects parameters
        psr_poi.range["data"] = np.array([r])
        psr_poi.ngates = 1

        psr_poi.time["units"] = pyart.io.make_time_unit_str(time_ref)
        psr_poi.time["data"] = np.array([])
        psr_poi.scan_type = "poi_time_series"
        psr_poi.sweep_number["data"] = np.array([], dtype=np.int32)
        psr_poi.nsweeps = 1
        psr_poi.sweep_mode["data"] = np.array(["poi_time_series"])
        psr_poi.rays_are_indexed = None
        psr_poi.ray_angle_res = None
        psr_poi.fixed_angle["data"] = np.array([az])

        # ray dependent psr objects parameters
        psr_poi.sweep_end_ray_index["data"] = np.array([-1], dtype="int32")
        psr_poi.rays_per_sweep["data"] = np.array([0], dtype="int32")
        psr_poi.azimuth["data"] = np.array([], dtype="float64")
        psr_poi.elevation["data"] = np.array([], dtype="float64")
        psr_poi.nrays = 0

        psr_poi.npulses["data"] = np.array([], dtype=np.int32)
        if psr_poi.Doppler_velocity is not None:
            psr_poi.Doppler_velocity["data"] = np.array([])
        if psr_poi.Doppler_frequency is not None:
            psr_poi.Doppler_frequency["data"] = np.array([])

        dscfg["global_data"] = {
            "psr_poi": psr_poi,
            "point_coordinates_WGS84_lon_lat_alt": [lon, lat, alt],
            "antenna_coordinates_az_el_r": [az, el, r],
        }

        dscfg["initialized"] = 1

    psr_poi = dscfg["global_data"]["psr_poi"]
    psr_poi.sweep_end_ray_index["data"][0] += nrays
    psr_poi.rays_per_sweep["data"][0] += nrays
    psr_poi.nrays += nrays
    psr_poi.azimuth["data"] = np.append(psr_poi.azimuth["data"], np.zeros(nrays) + az)
    psr_poi.elevation["data"] = np.append(
        psr_poi.elevation["data"], np.zeros(nrays) + el
    )
    start_time = num2date(0, psr_poi.time["units"], psr_poi.time["calendar"])
    for i in range(nrays):
        psr_poi.time["data"] = np.append(
            psr_poi.time["data"], (time_poi[i] - start_time).total_seconds()
        )

    psr_poi.gate_longitude["data"] = (
        np.ones((psr_poi.nrays, psr_poi.ngates), dtype="float64") * lon
    )
    psr_poi.gate_latitude["data"] = (
        np.ones((psr_poi.nrays, psr_poi.ngates), dtype="float64") * lat
    )
    psr_poi.gate_altitude["data"] = np.broadcast_to(
        alt, (psr_poi.nrays, psr_poi.ngates)
    )

    for field_name in field_names:
        dtype = psr_poi.fields[field_name]["data"].dtype
        if field_name not in psr.fields:
            warn("Field " + field_name + " not in psr object")
            poi_data = np.ma.masked_all((nrays, 1, psr.npulses_max), dtype=dtype)
        else:
            poi_data = psr.fields[field_name]["data"][ind_ray, ind_rng, :]
            poi_data = poi_data.reshape(nrays, 1, psr.npulses_max)

        # Put data in radar object
        if np.size(psr_poi.fields[field_name]["data"]) == 0:
            psr_poi.fields[field_name]["data"] = poi_data.reshape(
                nrays, 1, psr_poi.npulses_max
            )
        else:
            if psr_poi.npulses_max == psr.npulses_max:
                psr_poi.fields[field_name]["data"] = np.ma.append(
                    psr_poi.fields[field_name]["data"], poi_data, axis=0
                )
            elif psr.npulses_max < psr_poi.npulses_max:
                poi_data_aux = np.ma.masked_all(
                    (nrays, 1, psr_poi.npulses_max), dtype=dtype
                )
                for i in range(nrays):
                    poi_data_aux[i, 0, 0 : psr.npulses_max] = poi_data[i, 0, :]
                psr_poi.fields[field_name]["data"] = np.ma.append(
                    psr_poi.fields[field_name]["data"], poi_data_aux, axis=0
                )
            else:
                poi_data_aux = np.ma.masked_all(
                    (psr_poi.nrays, 1, psr.npulses_max), dtype=dtype
                )
                poi_data_aux[: psr_poi.nrays - nrays, :, 0 : psr_poi.npulses_max] = (
                    psr_poi.fields[field_name]["data"]
                )
                poi_data_aux[psr_poi.nrays - nrays :, :, :] = poi_data
                psr_poi.fields[field_name]["data"] = poi_data_aux

    psr_poi.npulses["data"] = np.append(
        psr_poi.npulses["data"], psr.npulses["data"][ind_ray]
    )
    if psr_poi.Doppler_velocity is not None:
        if np.size(psr_poi.Doppler_velocity["data"]) == 0:
            psr_poi.Doppler_velocity["data"] = psr.Doppler_velocity["data"][
                ind_ray, :
            ].reshape(nrays, psr_poi.npulses_max)
        else:
            Doppler_data = psr.Doppler_velocity["data"][ind_ray, :]
            Doppler_data = Doppler_data.reshape(nrays, psr.npulses_max)

            if psr_poi.npulses_max == psr.npulses_max:
                psr_poi.Doppler_velocity["data"] = np.ma.append(
                    psr_poi.Doppler_velocity["data"], Doppler_data, axis=0
                )
            elif psr.npulses_max < psr_poi.npulses_max:
                Doppler_aux = np.ma.masked_all((nrays, psr_poi.npulses_max))
                for i in range(nrays):
                    Doppler_aux[i, 0 : psr.npulses_max] = Doppler_data[i, :]
                psr_poi.Doppler_velocity["data"] = np.ma.append(
                    psr_poi.Doppler_velocity["data"], Doppler_aux, axis=0
                )
            else:
                Doppler_aux = np.ma.masked_all((psr_poi.nrays, psr.npulses_max))
                Doppler_aux[: psr_poi.nrays - nrays, 0 : psr_poi.npulses_max] = (
                    psr_poi.Doppler_velocity["data"]
                )
                Doppler_aux[psr_poi.nrays - nrays :, :] = Doppler_data
                psr_poi.Doppler_velocity["data"] = Doppler_aux

    if psr_poi.Doppler_frequency is not None:
        if np.size(psr_poi.Doppler_frequency["data"]) == 0:
            psr_poi.Doppler_frequency["data"] = psr.Doppler_frequency["data"][
                ind_ray, :
            ].reshape(nrays, psr_poi.npulses_max)
        else:
            Doppler_data = psr.Doppler_frequency["data"][ind_ray, :]
            Doppler_data = Doppler_data.reshape(nrays, psr.npulses_max)

            if psr_poi.npulses_max == psr.npulses_max:
                psr_poi.Doppler_frequency["data"] = np.ma.append(
                    psr_poi.Doppler_frequency["data"], Doppler_data, axis=0
                )
            elif psr.npulses_max < psr_poi.npulses_max:
                Doppler_aux = np.ma.masked_all((nrays, psr_poi.npulses_max))
                for i in range(nrays):
                    Doppler_aux[i, 0 : psr.npulses_max] = Doppler_data[i, :]
                psr_poi.Doppler_frequency["data"] = np.ma.append(
                    psr_poi.Doppler_frequency["data"], Doppler_aux, axis=0
                )
            else:
                Doppler_aux = np.ma.masked_all((psr_poi.nrays, psr.npulses_max))
                Doppler_aux[0 : psr_poi.nrays - nrays, 0 : psr_poi.npulses_max] = (
                    psr_poi.Doppler_frequency["data"]
                )
                Doppler_aux[psr_poi.nrays - nrays :, :] = Doppler_data
                psr_poi.Doppler_frequency["data"] = Doppler_aux

    psr_poi.npulses_max = max(psr_poi.npulses_max, psr.npulses_max)

    dscfg["global_data"]["psr_poi"] = psr_poi

    # prepare for exit
    new_dataset = {
        "radar_out": psr_poi,
        "point_coordinates_WGS84_lon_lat_alt": (
            dscfg["global_data"]["point_coordinates_WGS84_lon_lat_alt"]
        ),
        "antenna_coordinates_az_el_r": (
            dscfg["global_data"]["antenna_coordinates_az_el_r"]
        ),
        "final": False,
    }

    return new_dataset, ind_rad


def process_filter_0Doppler(procstatus, dscfg, radar_list=None):
    """
    Function to filter the 0-Doppler line bin and neighbours of the
    Doppler spectra

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, can be any of the spectral fields supported by pyrad

        filter_width : float
            The Doppler filter width. Default 0.
        filter_units : str
            Can be 'm/s' or 'Hz'. Default 'm/s'
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field, the names of the output fields
        is the same as the provided datatypes, except for unfiltered fields which are renamed in the following
        "dBuZ" => "dBZ"
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    field_name_list = []
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_name_list.append(get_fieldname_pyart(datatype))

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    filter_width = dscfg.get("filter_width", 0.0)
    filter_units = dscfg.get("filter_units", "m/s")

    if filter_units == "m/s":
        axis = psr.Doppler_velocity["data"]
    else:
        axis = psr.Doppler_frequency["data"]

    fields = dict()
    for field_name in field_name_list:
        if field_name not in psr.fields:
            warn("Unable to filter 0-Doppler. Missing field " + field_name)
            continue

        field_name_aux = field_name.replace("unfiltered_", "")
        field = pyart.config.get_metadata(field_name_aux)
        field["data"] = deepcopy(psr.fields[field_name]["data"])
        for ray in range(psr.nrays):
            ind = np.ma.where(
                np.logical_and(
                    axis[ray, :] >= -filter_width / 2.0,
                    axis[ray, :] <= filter_width / 2.0,
                )
            )
            field["data"][ray, :, ind] = np.ma.masked
        fields.update({field_name_aux: field})

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(psr)}
    new_dataset["radar_out"].fields = dict()
    for field_name in fields.keys():
        new_dataset["radar_out"].add_field(field_name, fields[field_name])

    return new_dataset, ind_rad


def process_filter_srhohv(procstatus, dscfg, radar_list=None):
    """
    Filter Doppler spectra as a function of spectral RhoHV

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain
            "sRhoHV" or "sRhoHVu",
            as well as any spectral field supported by pyrad
        sRhoHV_threshold : float
            Data with sRhoHV module above this threshold will be filtered.
            Default 1.
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field, the names of the output fields
        is the same as the provided datatypes, except for unfiltered fields which are renamed in the following
        "dBuZ" => "dBZ"
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    field_name_list = []
    sRhoHV_found = False
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("sRhoHV", "sRhoHVu") and not sRhoHV_found:
            sRhoHV_field = get_fieldname_pyart(datatype)
            sRhoHV_found = True
        else:
            field_name_list.append(get_fieldname_pyart(datatype))

    if not sRhoHV_found:
        warn("sRhoHV field is required for sRhoHV filtering")
        return None, None

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    if sRhoHV_field not in psr.fields:
        warn("Unable to obtain apply sRhoHV filter. Missing field " + sRhoHV_field)
        return None, None

    sRhoHV_threshold = dscfg.get("sRhoHV_threshold", 0.9)
    sRhoHV = psr.fields[sRhoHV_field]["data"]

    fields = dict()
    for field_name in field_name_list:
        if field_name not in psr.fields:
            warn("Unable to filter according to sRhoHV. Missing field " + field_name)
            continue

        field_name_aux = field_name.replace("unfiltered_", "")
        field = pyart.config.get_metadata(field_name_aux)
        field["data"] = deepcopy(psr.fields[field_name]["data"])
        field["data"][np.ma.abs(sRhoHV) <= sRhoHV_threshold] = np.ma.masked
        fields.update({field_name_aux: field})

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(psr)}
    new_dataset["radar_out"].fields = dict()
    for field_name in fields.keys():
        new_dataset["radar_out"].add_field(field_name, fields[field_name])

    return new_dataset, ind_rad


def process_filter_spectra_noise(procstatus, dscfg, radar_list=None):
    """
    Filter the noise of the Doppler spectra by clipping any data below
    the noise level plus a margin

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain,
            "ShhADU" or "SvvADU" or "ShhADUu" or "SvvADUu", and,
            "sNADUh" or "sNADUv"
        clipping_level : float
            The clipping level [dB above noise level]. Default 10.
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field, the names of the output fields
        is the same as the provided datatypes
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    field_name_list = []
    signal_found = False
    noise_found = False
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("ShhADU", "SvvADU", "ShhADUu", "SvvADUu") and not signal_found:
            signal_field = get_fieldname_pyart(datatype)
            signal_found = True
        elif datatype in ("sNADUh", "sNADUv") and not noise_found:
            noise_field = get_fieldname_pyart(datatype)
            noise_found = True
        else:
            field_name_list.append(get_fieldname_pyart(datatype))

    if not signal_found or not noise_found:
        warn("Signal and noise fields are required for noise filtering")
        return None, None

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    if signal_field not in psr.fields or noise_field not in psr.fields:
        warn("Unable to obtain apply spectral noise filter. Missing fields")
        return None, None

    clipping_level = dscfg.get("clipping_level", 10.0)

    # get Doppler bins below clipping level
    clip_pwr = psr.fields[noise_field]["data"] * np.power(10.0, 0.1 * clipping_level)

    s_pwr = pyart.retrieve.compute_spectral_power(
        psr, units="ADU", signal_field=signal_field, noise_field=noise_field
    )

    mask = np.ma.less_equal(s_pwr["data"], clip_pwr)

    # filter data
    new_dataset = {"radar_out": deepcopy(psr)}
    new_dataset["radar_out"].fields = dict()
    for field_name in field_name_list:
        if field_name not in psr.fields:
            warn("Unable to filter field " + field_name)
            continue

        new_dataset["radar_out"].add_field(field_name, psr.fields[field_name])
        new_dataset["radar_out"].fields[field_name]["data"][mask] = np.ma.masked

    return new_dataset, ind_rad


def process_dealias_spectra(procstatus, dscfg, radar_list=None):
    """
    Performs a dealiasing of spectra data, assuming at most one fold

    The method is quite simple and works in the following way at every
    radar gate

    - aliasing check
        check if there is no noise either on the left of the right of the spectrum
    - left/right tail computation
        identify left and right tail of the aliased spectrum
    - detect direction of aliasing
        compute Doppler mean velocity with right or left shift of the spectrum,
        perform the shift which minimizes the difference in Doppler velocity to
        previous range bin

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain,
            "ShhADU" or "SvvADU" or "ShhADUu" or "SvvADUu", and,
            "sNADUh" or "sNADUv"

    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output fields. The output fields are the same as the original
        ones but are dealiased.
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    field_name_list = []
    signal_found = False
    noise_found = False
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("ShhADU", "SvvADU", "ShhADUu", "SvvADUu") and not signal_found:
            signal_field = get_fieldname_pyart(datatype)
            signal_found = True
        elif datatype in ("sNADUh", "sNADUv") and not noise_found:
            noise_field = get_fieldname_pyart(datatype)
            noise_found = True
        else:
            field_name_list.append(get_fieldname_pyart(datatype))

    if not signal_found or not noise_found:
        warn("Signal and noise fields are required for noise filtering")
        return None, None

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    if signal_field not in psr.fields or noise_field not in psr.fields:
        warn("Unable to obtain apply spectral noise filter. Missing fields")
        return None, None

    s_pwr = pyart.retrieve.compute_spectral_power(
        psr, units="ADU", signal_field=signal_field, noise_field=noise_field
    )

    # Add spectral power field
    psr.fields["spectral_power"] = s_pwr

    # Dealias dataset
    new_dataset = pyart.retrieve.spectra.dealias_spectra(
        psr, "spectral_power", field_name_list
    )

    return new_dataset, ind_rad


def process_spectra_ang_avg(procstatus, dscfg, radar_list=None):
    """
    Function to average the spectra over the rays. This function is
    intended mainly for vertically pointing scans. The function assumes
    the volume is composed of a single sweep, it averages over the number
    of rays specified by the user and produces a single ray output.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types,
            any spectral datatype supported by pyrad
        navg : int
            Number of spectra to average. If -1 all spectra will be averaged.
            Default -1.
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the same output fields as the provided datatypes
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    field_name_list = []
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_name_list.append(get_fieldname_pyart(datatype))

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    navg = dscfg.get("navg", -1)

    # keep only fields of interest
    psr_aux = deepcopy(psr)
    psr_aux.fields = dict()
    for field_name in field_name_list:
        if field_name not in psr.fields:
            warn("Field " + field_name + " missing")
            continue
        psr_aux.add_field(field_name, psr.fields[field_name])

    psr_aux = pyart.util.interpol_spectra(psr_aux)

    if navg == -1:
        navg = psr.nrays
    elif navg > psr.nrays:
        warn(
            "Number of rays " + str(psr.nrays) + " smaller than number of "
            "desired spectra to average " + str(navg)
        )
        navg = psr.nrays

    for field_name in psr_aux.fields.keys():
        data_mean = np.ma.mean(
            psr_aux.fields[field_name]["data"][0:navg, :, 0 : psr_aux.npulses_max],
            axis=0,
        )
        psr_aux.fields[field_name]["data"] = np.ma.masked_all(
            (1, psr_aux.ngates, psr_aux.npulses_max),
            dtype=psr_aux.fields[field_name]["data"].dtype,
        )
        psr_aux.fields[field_name]["data"][0, :, :] = data_mean

    psr_aux.time["data"] = np.array([psr_aux.time["data"][int(navg / 2)]])
    psr_aux.azimuth["data"] = np.array([0], dtype=np.float32)
    psr_aux.elevation["data"] = np.array([psr_aux.elevation["data"][int(navg / 2)]])
    psr_aux.nrays = 1
    psr_aux.sweep_end_ray_index["data"] = np.array([0.0], dtype=np.int32)

    psr_aux.init_rays_per_sweep()
    psr_aux.init_gate_x_y_z()
    psr_aux.init_gate_longitude_latitude()
    psr_aux.init_gate_altitude()

    if psr_aux.Doppler_velocity is not None:
        psr_aux.Doppler_velocity["data"] = np.ma.expand_dims(
            psr_aux.Doppler_velocity["data"][0, :], axis=0
        )
    if psr_aux.Doppler_frequency is not None:
        psr_aux.Doppler_frequency["data"] = np.ma.expand_dims(
            psr_aux.Doppler_frequency["data"][0, :], axis=0
        )

    # prepare for exit
    new_dataset = {"radar_out": psr_aux}

    return new_dataset, ind_rad


def process_spectral_power(procstatus, dscfg, radar_list=None):
    """
    Computes the spectral power

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain,
            "ShhADU" or "SvvADU" or "ShhADUu" or "SvvADUu", and,
            "sNADUh", or "sNADUv"
        units : str
            The units of the returned signal. Can be 'ADU', 'dBADU' or 'dBm'
        subtract_noise : Bool
            If True noise will be subtracted from the signal
        smooth_window : int or None
            Size of the moving Gaussian smoothing window. If none no smoothing
            will be applied
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field
        "sPhhADU" or "sPhhADUu", or
        "sPvvADU" or "sPvvADUu", or
        "sPhhdBADU" or "sPhhdBADUu", or
        "sPvvdBADU" or "sPvvdBADUu", or
        "sPhhAdBm" or "sPhhdBmu", or
        "sPvvdBm" or "sPvvdBmu",
        depending on which input datatype and units were provided
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    noise_field = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("ShhADU", "SvvADU", "ShhADUu", "SvvADUu"):
            signal_field = get_fieldname_pyart(datatype)
        elif datatype in ("sNADUh", "sNADUv"):
            noise_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    if signal_field not in psr.fields:
        warn("Unable to obtain spectral signal power. Missing field " + signal_field)
        return None, None

    units = dscfg.get("units", "dBADU")
    subtract_noise = dscfg.get("subtract_noise", False)
    smooth_window = dscfg.get("smooth_window", None)

    s_pwr = pyart.retrieve.compute_spectral_power(
        psr,
        units=units,
        subtract_noise=subtract_noise,
        smooth_window=smooth_window,
        signal_field=signal_field,
        noise_field=noise_field,
    )

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(psr)}
    new_dataset["radar_out"].fields = dict()
    new_dataset["radar_out"].add_field(s_pwr["standard_name"], s_pwr)

    return new_dataset, ind_rad


def process_spectral_noise(procstatus, dscfg, radar_list=None):
    """
    Computes the spectral noise

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain,
            "ShhADU" or "SvvADU" or "ShhADUu" or "SvvADUu"
        units : str
            The units of the returned signal. Can be 'ADU', 'dBADU' or 'dBm'
        navg : int
            Number of spectra averaged
        rmin : int
            Range from which the data is used to estimate the noise
        nnoise_min : int
            Minimum number of samples to consider the estimated noise power
            valid
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field
            "sNADUh" or
            "sNADUv" or
            "sNdBADUh" or
            "sNdBADUv" or
            "sNdBmh" or
            "sNdBmv"
            depending on which input datatype and units were provided
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("ShhADU", "SvvADU", "ShhADUu", "SvvADUu"):
            signal_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    if signal_field not in psr.fields:
        warn("Unable to obtain spectral noise power. Missing field " + signal_field)
        return None, None

    units = dscfg.get("units", "ADU")
    navg = dscfg.get("navg", 1)
    rmin = dscfg.get("rmin", 0.0)
    nnoise_min = dscfg.get("nnoise_min", 100)

    s_pwr = pyart.retrieve.compute_spectral_noise(
        psr,
        units=units,
        navg=navg,
        rmin=rmin,
        nnoise_min=nnoise_min,
        signal_field=signal_field,
    )

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(psr)}
    new_dataset["radar_out"].fields = dict()
    new_dataset["radar_out"].add_field(s_pwr["standard_name"], s_pwr)

    return new_dataset, ind_rad


def process_spectral_phase(procstatus, dscfg, radar_list=None):
    """
    Computes the spectral phase

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain,
            "ShhADU" or "SvvADU" or "ShhADUu" or "SvvADUu"
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field
        "SPhasehh" if "ShhADU" was provided as input
        "SPhasehhu" if "ShhADUu" was provided as input
        "SPhasevv" if "SvvADU" was provided as input
        "SPhasevvu" if "SvvADUu" was provided as input

    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("ShhADU", "SvvADU", "ShhADUu", "SvvADUu"):
            signal_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    if signal_field not in psr.fields:
        warn("Unable to obtain spectral phase. Missing field " + signal_field)
        return None, None

    s_phase = pyart.retrieve.compute_spectral_phase(psr, signal_field=signal_field)

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(psr)}
    new_dataset["radar_out"].fields = dict()
    new_dataset["radar_out"].add_field(s_phase["standard_name"], s_phase)

    return new_dataset, ind_rad


def process_spectral_reflectivity(procstatus, dscfg, radar_list=None):
    """
    Computes spectral reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain,
            either a combination of signal and noise
            "ShhADU" or "SvvADU" or "ShhADUu" or "SvvADUu", and,
            "sNADUh" or "sNADUv", or
            the power signal
            "sPhhADU" or "sPvvADU" or "sPhhADUu" or "sPvvADUu"
        subtract_noise : Bool
            If True noise will be subtracted from the signal
        smooth_window : int or None
            Size of the moving Gaussian smoothing window. If none no smoothing
            will be applied
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field
        "sdBZ" if "ShhADU" (or "sPhhADU") was provided as input
        "sdBuZ" if "ShhADUu" (or "sPhhADUu") was provided as input
        "sdBZv" if "SvvADU" (or "sPvvADU") was provided as input
        "sdBuZv" if "SvvADUu" (or "sPvvADUu") was provided as input
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None
    noise_field = None
    signal_field = None
    pwr_field = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("ShhADU", "SvvADU", "ShhADUu", "SvvADUu"):
            signal_field = get_fieldname_pyart(datatype)
        elif datatype in ("sNADUh", "sNADUv"):
            noise_field = get_fieldname_pyart(datatype)
        elif datatype in ("sPhhADU", "sPvvADU", "sPhhADUu", "sPvvADUu"):
            pwr_field = get_fieldname_pyart(datatype)

    if pwr_field is None and signal_field is None:
        warn("Either signal or power fields must be specified")
        return None, None

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    compute_power = True
    if pwr_field is not None:
        compute_power = False

    if compute_power and signal_field not in psr.fields:
        warn("Unable to obtain spectral reflectivity. Missing field " + signal_field)
        return None, None

    if not compute_power and pwr_field not in psr.fields:
        warn("Unable to obtain spectral reflectivity. Missing field " + pwr_field)
        return None, None

    subtract_noise = dscfg.get("subtract_noise", False)
    smooth_window = dscfg.get("smooth_window", None)

    sdBZ = pyart.retrieve.compute_spectral_reflectivity(
        psr,
        compute_power=compute_power,
        subtract_noise=subtract_noise,
        smooth_window=smooth_window,
        pwr_field=pwr_field,
        signal_field=signal_field,
        noise_field=noise_field,
    )

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(psr)}
    new_dataset["radar_out"].fields = dict()
    new_dataset["radar_out"].add_field(sdBZ["standard_name"], sdBZ)

    return new_dataset, ind_rad


def process_spectral_differential_reflectivity(procstatus, dscfg, radar_list=None):
    """
    Computes spectral differential reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain,
            either a combination of signal and noise
            ("ShhADU" and "SvvADU") or ("ShhADUu" and "SvvADUu"), and,
            ("sNADUh" and "sNADUv"), or
            the power signal
            ("sPhhADU" and "sPvvADU") or ("sPhhADUu" and "sPvvADUu")
        subtract_noise : Bool
            If True noise will be subtracted from the signal
        smooth_window : int or None
            Size of the moving Gaussian smoothing window. If none no smoothing
            will be applied
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field
        "sZDR" if "ShhADU" and "SvvADU" were provided as input
        "sZDRu" if "ShhADUu" and "SvvADUu" were provided as input
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    noise_h_field = None
    noise_v_field = None
    signal_h_field = None
    signal_v_field = None
    pwr_h_field = None
    pwr_v_field = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("ShhADU", "ShhADUu"):
            signal_h_field = get_fieldname_pyart(datatype)
        elif datatype in ("SvvADU", "SvvADUu"):
            signal_v_field = get_fieldname_pyart(datatype)
        elif datatype == "sNADUh":
            noise_h_field = get_fieldname_pyart(datatype)
        elif datatype == "sNADUv":
            noise_v_field = get_fieldname_pyart(datatype)
        elif datatype in ("sPhhADU", "sPhhADUu"):
            pwr_h_field = get_fieldname_pyart(datatype)
        elif datatype in ("sPvvADU", "sPvvADUu"):
            pwr_v_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    compute_power = True
    if pwr_h_field is not None and pwr_v_field is not None:
        compute_power = False

    if compute_power and (
        signal_h_field not in psr.fields or signal_v_field not in psr.fields
    ):
        warn("Unable to obtain spectral differential reflectivity. " + "Missing fields")
        return None, None
    if not compute_power and (
        pwr_h_field not in psr.fields or pwr_v_field not in psr.fields
    ):
        warn("Unable to obtain spectral differential reflectivity. " + "Missing fields")
        return None, None

    subtract_noise = dscfg.get("subtract_noise", False)
    smooth_window = dscfg.get("smooth_window", None)

    sZDR = pyart.retrieve.compute_spectral_differential_reflectivity(
        psr,
        compute_power=compute_power,
        subtract_noise=subtract_noise,
        smooth_window=smooth_window,
        pwr_h_field=pwr_h_field,
        pwr_v_field=pwr_v_field,
        signal_h_field=signal_h_field,
        signal_v_field=signal_v_field,
        noise_h_field=noise_h_field,
        noise_v_field=noise_v_field,
    )

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(psr)}
    new_dataset["radar_out"].fields = dict()
    new_dataset["radar_out"].add_field(sZDR["standard_name"], sZDR)

    return new_dataset, ind_rad


def process_spectral_differential_phase(procstatus, dscfg, radar_list=None):
    """
    Computes the spectral differential phase

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain,
            either a combination of signal and noise
            ("ShhADU" and "SvvADU") or ("ShhADUu" and "SvvADUu"), and,
            "sRhoHV" or "sRhoHVu" (Optional)
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output fields
        "sPhiDP" if "ShhADU" and "SvvADU" were provided
        "sPhiDPu" if "ShhADUu" and "SvvADUu" were provided
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    signal_h_field = None
    signal_v_field = None
    srhohv_field = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("ShhADU", "ShhADUu"):
            signal_h_field = get_fieldname_pyart(datatype)
        elif datatype in ("SvvADU", "SvvADUu"):
            signal_v_field = get_fieldname_pyart(datatype)
        elif datatype in ("sRhoHV", "sRhoHVu"):
            srhohv_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    use_rhohv = False
    if srhohv_field is not None:
        use_rhohv = True

    if not use_rhohv and (
        signal_h_field not in psr.fields or signal_v_field not in psr.fields
    ):
        warn("Unable to obtain spectral signal differential phase. " + "Missing fields")
        return None, None
    if use_rhohv and srhohv_field not in psr.fields:
        warn("Unable to obtain spectral signal differential phase. " + "Missing fields")
        return None, None

    sPhiDP = pyart.retrieve.compute_spectral_differential_phase(
        psr,
        use_rhohv=use_rhohv,
        srhohv_field=srhohv_field,
        signal_h_field=signal_h_field,
        signal_v_field=signal_v_field,
    )

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(psr)}
    new_dataset["radar_out"].fields = dict()
    new_dataset["radar_out"].add_field(sPhiDP["standard_name"], sPhiDP)

    return new_dataset, ind_rad


def process_spectral_rhohv(procstatus, dscfg, radar_list=None):
    """
    Computes the spectral RhoHV

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain,
            ("ShhADU" and "SvvADU") or ("ShhADUu" and "SvvADUu"), and,
            ("sNADUh" and "sNADUv")
        subtract_noise : Bool
            If True noise will be subtracted from the signal
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output fields
        "sRhoHV" if "ShhADU" and "SvvADU" were provided
        "sRhoHVu" if "ShhADUu" and "SvvADUu" were provided
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    noise_h_field = None
    noise_v_field = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("ShhADU", "ShhADUu"):
            signal_h_field = get_fieldname_pyart(datatype)
        elif datatype in ("SvvADU", "SvvADUu"):
            signal_v_field = get_fieldname_pyart(datatype)
        elif datatype == "sNADUh":
            noise_h_field = get_fieldname_pyart(datatype)
        elif datatype == "sNADUv":
            noise_v_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    if signal_h_field not in psr.fields or signal_v_field not in psr.fields:
        warn("Unable to obtain spectral RhoHV. " + "Missing fields")
        return None, None

    subtract_noise = dscfg.get("subtract_noise", False)

    sRhoHV = pyart.retrieve.compute_spectral_rhohv(
        psr,
        subtract_noise=subtract_noise,
        signal_h_field=signal_h_field,
        signal_v_field=signal_v_field,
        noise_h_field=noise_h_field,
        noise_v_field=noise_v_field,
    )

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(psr)}
    new_dataset["radar_out"].fields = dict()
    new_dataset["radar_out"].add_field(sRhoHV["standard_name"], sRhoHV)

    return new_dataset, ind_rad


def process_pol_variables(procstatus, dscfg, radar_list=None):
    """
    Computes the polarimetric variables from the complex spectra

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain,
            either a combination of signal and noise
            ("ShhADU" and "SvvADU") or ("ShhADUu" and "SvvADUu"), and,
            ("sNADUh" and "sNADUv"), or
            the power signal
            ("sPhhADU" and "sPvvADU") or ("sPhhADUu" and "sPvvADUu"), and
            ("sRhoHV" or "sRhoHVu")

        subtract_noise : Bool
            If True noise will be subtracted from the signal. Default False
        smooth_window : int or None
            Size of the moving Gaussian smoothing window. If none no smoothing
            will be applied. Default None
        variables : list of str
            list of variables to compute. Default dBZ
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the all outputs fields, that correspond to the
        specified "variables" keyword
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    noise_h_field = None
    noise_v_field = None
    signal_h_field = None
    signal_v_field = None
    pwr_h_field = None
    pwr_v_field = None
    srhohv_field = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("ShhADU", "ShhADUu"):
            signal_h_field = get_fieldname_pyart(datatype)
        elif datatype in ("SvvADU", "SvvADUu"):
            signal_v_field = get_fieldname_pyart(datatype)
        elif datatype == "sNADUh":
            noise_h_field = get_fieldname_pyart(datatype)
        elif datatype == "sNADUv":
            noise_v_field = get_fieldname_pyart(datatype)
        elif datatype in ("sPhhADU", "sPhhADUu"):
            pwr_h_field = get_fieldname_pyart(datatype)
        elif datatype in ("sPvvADU", "sPvvADUu"):
            pwr_v_field = get_fieldname_pyart(datatype)
        elif datatype in ("sRhoHV", "sRhoHVu"):
            srhohv_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    use_pwr = False
    if pwr_h_field is not None or pwr_v_field is not None or srhohv_field is not None:
        use_pwr = True

    if not use_pwr and (
        signal_h_field not in psr.fields and signal_v_field not in psr.fields
    ):
        warn("Unable to obtain polarimetric variables. Missing fields")
        return None, None

    if use_pwr and (
        pwr_h_field not in psr.fields
        and pwr_h_field not in psr.fields
        and srhohv_field not in psr.fields
    ):
        warn("Unable to obtain polarimetric variables. Missing fields")
        return None, None

    subtract_noise = dscfg.get("subtract_noise", False)
    smooth_window = dscfg.get("smooth_window", None)
    variables = dscfg.get("variables", ["dBZ"])

    fields_list = []
    for variable in variables:
        fields_list.append(get_fieldname_pyart(variable))

    radar = pyart.retrieve.compute_pol_variables(
        psr,
        fields_list,
        use_pwr=use_pwr,
        subtract_noise=subtract_noise,
        smooth_window=smooth_window,
        srhohv_field=srhohv_field,
        pwr_h_field=pwr_h_field,
        pwr_v_field=pwr_v_field,
        signal_h_field=signal_h_field,
        signal_v_field=signal_v_field,
        noise_h_field=noise_h_field,
        noise_v_field=noise_v_field,
    )

    # prepare for exit
    new_dataset = {"radar_out": radar}

    return new_dataset, ind_rad


def process_noise_power(procstatus, dscfg, radar_list=None):
    """
    Computes the noise power from the spectra

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain
            "ShhADU" or "SvvADU" or "ShhADUu" or "SvvADUu"
        units : str
            The units of the returned signal. Can be 'ADU', 'dBADU' or 'dBm'
        navg : int
            Number of spectra averaged
        rmin : int
            Range from which the data is used to estimate the noise
        nnoise_min : int
            Minimum number of samples to consider the estimated noise power
            valid
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field
        "NdBADUh" or "NdBADUv", or
        "NdBmh" or "NdBmv", or
        "Nh" or "Nv"
        depending on which input datatype and units were provided
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("ShhADU", "SvvADU", "ShhADUu", "SvvADUu"):
            signal_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    if signal_field not in psr.fields:
        warn("Unable to obtain spectral noise power. Missing field " + signal_field)
        return None, None

    units = dscfg.get("units", "ADU")
    navg = dscfg.get("navg", 1)
    rmin = dscfg.get("rmin", 0.0)
    nnoise_min = dscfg.get("nnoise_min", 100)

    noise = pyart.retrieve.compute_noise_power(
        psr,
        units=units,
        navg=navg,
        rmin=rmin,
        nnoise_min=nnoise_min,
        signal_field=signal_field,
    )

    # prepare for exit
    new_dataset = {"radar_out": pyart.util.radar_from_spectra(psr)}
    new_dataset["radar_out"].add_field(noise["standard_name"], noise)

    return new_dataset, ind_rad


def process_reflectivity(procstatus, dscfg, radar_list=None):
    """
    Computes reflectivity from the spectral reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain
            "sdBZ" or sdBZv" or "sdBuZ" or "sdBuZv"
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "dBZ", "dBZv",
        "dBuZ" or "dBuZv" depending on the provided input datatype
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("sdBZ", "sdBZv", "sdBuZ", "sdBuZv"):
            sdBZ_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    if sdBZ_field not in psr.fields:
        warn("Unable to obtain reflectivity. " + "Missing field " + sdBZ_field)
        return None, None

    dBZ = pyart.retrieve.compute_reflectivity(psr, sdBZ_field=sdBZ_field)

    reflectivity_field = "reflectivity"
    if datatype in ("sdBZv", "sdBuZv"):
        reflectivity_field += "vv"

    if datatype in ("sdBuZ", "sdBuZv"):
        reflectivity_field = "unfiltered_" + reflectivity_field

    # prepare for exit
    new_dataset = {"radar_out": pyart.util.radar_from_spectra(psr)}
    new_dataset["radar_out"].add_field(reflectivity_field, dBZ)

    return new_dataset, ind_rad


def process_differential_reflectivity(procstatus, dscfg, radar_list=None):
    """
    Computes differential reflectivity from the horizontal and vertical
    spectral reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain,
            "sdBZ" and "sdBZv", or
            "sdBuZ" and "sdBZuv"
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output fields
        "ZDR" or "ZDRu" depending on the specified input datatype
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("sdBZ", "sdBuZ"):
            sdBZ_field = get_fieldname_pyart(datatype)
        elif datatype in ("sdBZv", "sdBuZv"):
            sdBZv_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    if sdBZ_field not in psr.fields or sdBZv_field not in psr.fields:
        warn("Unable to obtain differential reflectivity. " + "Missing fields.")
        return None, None

    zdr = pyart.retrieve.compute_differential_reflectivity(
        psr, sdBZ_field=sdBZ_field, sdBZv_field=sdBZv_field
    )

    zdr_field = "differential_reflectivity"
    if "unfiltered" in sdBZ_field:
        zdr_field = "unfiltered_" + zdr_field

    # prepare for exit
    new_dataset = {"radar_out": pyart.util.radar_from_spectra(psr)}
    new_dataset["radar_out"].add_field(zdr_field, zdr)

    return new_dataset, ind_rad


def process_differential_phase(procstatus, dscfg, radar_list=None):
    """
    Computes the differential phase from the spectral differential phase and
    the spectral reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "uPhiDPu" or
        "uPhiDP" depending on the provided datatypes
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("sdBZ", "sdBZv", "sdBuZ", "sdBuZv"):
            sdBZ_field = get_fieldname_pyart(datatype)
        elif datatype in ("sPhiDP", "sPhiDPu"):
            sPhiDP_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    if sdBZ_field not in psr.fields or sPhiDP_field not in psr.fields:
        warn("Unable to obtain PhiDP. Missing fields")
        return None, None

    uphidp = pyart.retrieve.compute_differential_phase(
        psr, sdBZ_field=sdBZ_field, sPhiDP_field=sPhiDP_field
    )

    uphidp_field = "uncorrected_differential_phase"
    if "unfiltered" in sPhiDP_field:
        uphidp_field = "uncorrected_unfiltered_differential_phase"

    # prepare for exit
    new_dataset = {"radar_out": pyart.util.radar_from_spectra(psr)}
    new_dataset["radar_out"].add_field(uphidp_field, uphidp)

    return new_dataset, ind_rad


def process_rhohv(procstatus, dscfg, radar_list=None):
    """
    Computes RhoHV from the complex spectras

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain,
            either a combination of signal and noise
            ("ShhADU" and "SvvADU") or ("ShhADUu" and "SvvADUu"), and,
            ("sNADUh" and "sNADUv"), or
            the power signal
            ("sPhhADU" and "sPvvADU") or ("sPhhADUu" and "sPvvADUu"), and
            optionally ("sRhoHV" or "sRhoHVu")

        subtract_noise : Bool
            If True noise will be subtracted from the signal
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "RhoHV"
        or "uRhoHV" depending on the provided datatypes

    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    noise_h_field = None
    noise_v_field = None
    signal_h_field = None
    signal_v_field = None
    pwr_h_field = None
    pwr_v_field = None
    srhohv_field = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("ShhADU", "ShhADUu"):
            signal_h_field = get_fieldname_pyart(datatype)
        elif datatype in ("SvvADU", "SvvADUu"):
            signal_v_field = get_fieldname_pyart(datatype)
        elif datatype == "sNADUh":
            noise_h_field = get_fieldname_pyart(datatype)
        elif datatype == "sNADUv":
            noise_v_field = get_fieldname_pyart(datatype)
        elif datatype in ("sPhhADU", "sPhhADUu"):
            pwr_h_field = get_fieldname_pyart(datatype)
        elif datatype in ("sPvvADU", "sPvvADUu"):
            pwr_v_field = get_fieldname_pyart(datatype)
        elif datatype in ("sRhoHV", "sRhoHVu"):
            srhohv_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    if srhohv_field is not None:
        use_rhohv = True
        rhohv_field = "cross_correlation_ratio"
        if "unfiltered" in srhohv_field:
            rhohv_field = "unfiltered_cross_correlation_ratio"
    else:
        use_rhohv = False
        rhohv_field = "cross_correlation_ratio"
        if "unfiltered" in signal_h_field:
            rhohv_field = "unfiltered_cross_correlation_ratio"

    if not use_rhohv and (
        signal_h_field not in psr.fields or signal_v_field not in psr.fields
    ):
        warn("Unable to obtain RhoHV. Missing fields")
        return None, None

    if use_rhohv and (
        srhohv_field not in psr.fields
        or pwr_h_field not in psr.fields
        or pwr_v_field not in psr.fields
    ):
        warn("Unable to obtain RhoHV. Missing fields")
        return None, None

    subtract_noise = dscfg.get("subtract_noise", False)

    rhohv = pyart.retrieve.compute_rhohv(
        psr,
        use_rhohv=use_rhohv,
        subtract_noise=subtract_noise,
        srhohv_field=srhohv_field,
        pwr_h_field=pwr_h_field,
        pwr_v_field=pwr_v_field,
        signal_h_field=signal_h_field,
        signal_v_field=signal_v_field,
        noise_h_field=noise_h_field,
        noise_v_field=noise_v_field,
    )

    # prepare for exit
    new_dataset = {"radar_out": pyart.util.radar_from_spectra(psr)}
    new_dataset["radar_out"].add_field(rhohv_field, rhohv)

    return new_dataset, ind_rad


def process_Doppler_velocity(procstatus, dscfg, radar_list=None):
    """
    Compute the Doppler velocity from the spectral reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain,
            "sdBZ" or "sdBZv" or "sdBuZ" or "sdBuZv"
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field
        "V" if "sdBZ" was provided
        "Vv" if "sdBZv" was provided
        "Vu" if "sdBuZ" was provided
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("sdBZ", "sdBZv", "sdBuZ", "sdBuZv"):
            sdBZ_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    if sdBZ_field not in psr.fields:
        warn("Unable to obtain Doppler velocity. " + "Missing field " + sdBZ_field)
        return None, None

    vel = pyart.retrieve.compute_Doppler_velocity(psr, sdBZ_field=sdBZ_field)

    vel_field = "velocity"
    if datatype in ("sdBZv", "sdBuZv"):
        vel_field += "_vv"
    if datatype in ("sdBuZ", "sdBuZv"):
        vel_field = "unfiltered_" + vel_field

    # prepare for exit
    new_dataset = {"radar_out": pyart.util.radar_from_spectra(psr)}
    new_dataset["radar_out"].add_field(vel_field, vel)

    return new_dataset, ind_rad


def process_Doppler_width(procstatus, dscfg, radar_list=None):
    """
    Compute the Doppler spectrum width from the spectral reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain,
            "sdBZ" or "sdBZv" or "sdBuZ" or "sdBuZv"
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field
        "W" if "sdBZ" was provided
        "Wv" if "sdBZv" was provided
        "Wu" if "sdBuZ" was provided
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("sdBZ", "sdBZv", "sdBuZ", "sdBuZv"):
            sdBZ_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    psr = radar_list[ind_rad]

    if sdBZ_field not in psr.fields:
        warn(
            "Unable to obtain Doppler spectrum width. " + "Missing field " + sdBZ_field
        )
        return None, None

    width = pyart.retrieve.compute_Doppler_width(psr, sdBZ_field=sdBZ_field)

    width_field = "spectrum_width"
    if datatype in ("sdBZv", "sdBuZv"):
        width_field += "_vv"
    if datatype in ("sdBuZ", "sdBuZv"):
        width_field = "unfiltered_" + width_field

    # prepare for exit
    new_dataset = {"radar_out": pyart.util.radar_from_spectra(psr)}
    new_dataset["radar_out"].add_field(width_field, width)

    return new_dataset, ind_rad
