"""
pyrad.proc.process_iq
=====================

Functions to processes IQ data.

.. autosummary::
    :toctree: generated/

    process_raw_iq
    process_pol_variables_iq
    process_reflectivity_iq
    process_st1_iq
    process_st2_iq
    process_wbn_iq
    process_differential_reflectivity_iq
    process_mean_phase_iq
    process_differential_phase_iq
    process_rhohv_iq
    process_Doppler_velocity_iq
    process_Doppler_width_iq
    process_fft

"""

from copy import deepcopy
from warnings import warn

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart


def process_raw_iq(procstatus, dscfg, radar_list=None):
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


def process_pol_variables_iq(procstatus, dscfg, radar_list=None):
    """
    Computes the polarimetric variables from the IQ data

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain
            "IQhhADU", "IQvvADU", "IQNADUh" and "IQNADUv"
        subtract_noise : Bool
            If True noise will be subtracted from the signal
        lag : int
            The time lag to use in the estimators
        direction : str
            The convention used in the Doppler mean field. Can be
            negative_away or negative_towards
        variables : list of str
            list of variables to compute. Default dBZ
        phase_offset : float. Dataset keyword
            The system differential phase offset to remove
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output fields corresponding to the specified
        "variables"
        ""
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    noise_h_field = None
    noise_v_field = None
    signal_h_field = None
    signal_v_field = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == "IQhhADU":
            signal_h_field = get_fieldname_pyart(datatype)
        elif datatype == "IQvvADU":
            signal_v_field = get_fieldname_pyart(datatype)
        elif datatype == "IQNADUh":
            noise_h_field = get_fieldname_pyart(datatype)
        elif datatype == "IQNADUv":
            noise_v_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    subtract_noise = dscfg.get("subtract_noise", False)
    lag = dscfg.get("lag", None)
    direction = dscfg.get("direction", "negative_away")
    variables = dscfg.get("variables", ["dBZ"])
    phase_offset = dscfg.get("phase_offset", 0.0)

    fields_list = []
    for variable in variables:
        fields_list.append(get_fieldname_pyart(variable))

    radar = pyart.retrieve.compute_pol_variables_iq(
        radar,
        fields_list,
        subtract_noise=subtract_noise,
        lag=lag,
        direction=direction,
        phase_offset=phase_offset,
        signal_h_field=signal_h_field,
        signal_v_field=signal_v_field,
        noise_h_field=noise_h_field,
        noise_v_field=noise_v_field,
    )

    # prepare for exit
    new_dataset = {"radar_out": radar}

    return new_dataset, ind_rad


def process_reflectivity_iq(procstatus, dscfg, radar_list=None):
    """
    Computes reflectivity from the IQ data

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain
            "IQhhADU" or "IQvvADU", and,
            "IQNADUh" or "IQNADUv"
        subtract_noise : Bool
            If True noise will be subtracted from the signal
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field
        "dBZ" if "IQhhADU" and "IQNADUh" are specified
        "dBZv" if "IQvvADU" and "IQNADUv" are specified
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    noise_field = None
    signal_field = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("IQhhADU", "IQvvADU"):
            signal_field = get_fieldname_pyart(datatype)
        elif datatype in ("IQNADUh", "IQNADUv"):
            noise_field = get_fieldname_pyart(datatype)

    if signal_field is None:
        warn("Signal field must be specified")
        return None, None

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if signal_field not in radar.fields:
        warn("Unable to obtain reflectivity. Missing field " + signal_field)
        return None, None

    subtract_noise = dscfg.get("subtract_noise", False)

    dBZ = pyart.retrieve.compute_reflectivity_iq(
        radar,
        subtract_noise=subtract_noise,
        signal_field=signal_field,
        noise_field=noise_field,
    )

    reflectivity_field = "reflectivity"
    if signal_field in ("IQ_vv_ADU",):
        reflectivity_field += "_vv"

    # prepare for exit
    new_dataset = {"radar_out": pyart.util.radar_from_spectra(radar)}
    new_dataset["radar_out"].add_field(reflectivity_field, dBZ)

    return new_dataset, ind_rad


def process_st1_iq(procstatus, dscfg, radar_list=None):
    """
    Computes the statistical test one lag fluctuation from the horizontal or
    vertical IQ data

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain
            "IQhhADU" or "IQvvADU"
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "ST1" (stat_test_lag1)
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg["datatype"][0])
    signal_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if signal_field not in radar.fields:
        warn("Unable to obtain ST1. Missing fields")
        return None, None

    st1 = pyart.retrieve.compute_st1_iq(radar, signal_field=signal_field)

    st1_field = "stat_test_lag1"
    if signal_field == "IQ_vv_ADU":
        st1_field += "_vv"

    # prepare for exit
    new_dataset = {"radar_out": pyart.util.radar_from_spectra(radar)}
    new_dataset["radar_out"].add_field(st1_field, st1)

    return new_dataset, ind_rad


def process_st2_iq(procstatus, dscfg, radar_list=None):
    """
    Computes the statistical test two lag fluctuation from the horizontal or
    vertical IQ data

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain
            "IQhhADU" or "IQvvADU"
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "ST2"
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg["datatype"][0])
    signal_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if signal_field not in radar.fields:
        warn("Unable to obtain ST2. Missing fields")
        return None, None

    st2 = pyart.retrieve.compute_st2_iq(radar, signal_field=signal_field)

    st2_field = "stat_test_lag2"
    if signal_field == "IQ_vv_ADU":
        st2_field += "_vv"

    # prepare for exit
    new_dataset = {"radar_out": pyart.util.radar_from_spectra(radar)}
    new_dataset["radar_out"].add_field(st2_field, st2)

    return new_dataset, ind_rad


def process_wbn_iq(procstatus, dscfg, radar_list=None):
    """
    Computes the wide band noise from the horizontal or vertical IQ data

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain
            "IQhhADU" or "IQvvADU"
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "WBN" (wide-band noise)
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg["datatype"][0])
    signal_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if signal_field not in radar.fields:
        warn("Unable to obtain WBN. Missing fields")
        return None, None

    wbn = pyart.retrieve.compute_wbn_iq(radar, signal_field=signal_field)

    wbn_field = "wide_band_noise"
    if signal_field == "IQ_vv_ADU":
        wbn_field += "_vv"

    # prepare for exit
    new_dataset = {"radar_out": pyart.util.radar_from_spectra(radar)}
    new_dataset["radar_out"].add_field(wbn_field, wbn)

    return new_dataset, ind_rad


def process_differential_reflectivity_iq(procstatus, dscfg, radar_list=None):
    """
    Computes differential reflectivity from the horizontal and vertical
    IQ data

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain
            "IQhhADU", "IQvvADU", "IQNADUh" and "IQNADUv"
        subtract_noise : Bool
            If True noise will be subtracted from the signal
        lag : int
            The time lag to use in the estimators
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "ZDR"
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    noise_h_field = None
    noise_v_field = None
    signal_h_field = None
    signal_v_field = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == "IQhhADU":
            signal_h_field = get_fieldname_pyart(datatype)
        elif datatype == "IQvvADU":
            signal_v_field = get_fieldname_pyart(datatype)
        elif datatype == "IQNADUh":
            noise_h_field = get_fieldname_pyart(datatype)
        elif datatype == "IQNADUv":
            noise_v_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if signal_h_field not in radar.fields or signal_v_field not in radar.fields:
        warn("Unable to obtain spectral differential reflectivity. " + "Missing fields")
        return None, None

    subtract_noise = dscfg.get("subtract_noise", False)
    lag = dscfg.get("lag", 0)

    zdr = pyart.retrieve.compute_differential_reflectivity_iq(
        radar,
        subtract_noise=subtract_noise,
        lag=lag,
        signal_h_field=signal_h_field,
        signal_v_field=signal_v_field,
        noise_h_field=noise_h_field,
        noise_v_field=noise_v_field,
    )

    # prepare for exit
    new_dataset = {"radar_out": pyart.util.radar_from_spectra(radar)}
    new_dataset["radar_out"].add_field("differential_reflectivity", zdr)

    return new_dataset, ind_rad


def process_mean_phase_iq(procstatus, dscfg, radar_list=None):
    """
    Computes the mean phase from the horizontal or vertical IQ data

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain
            "IQhhADU" or "IQvvADU"
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "MPH" (mean_phase)
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg["datatype"][0])
    signal_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if signal_field not in radar.fields:
        warn("Unable to obtain MPH. Missing fields")
        return None, None

    mph = pyart.retrieve.compute_mean_phase_iq(radar, signal_field=signal_field)

    mean_phase_field = "mean_phase"
    if signal_field == "IQ_vv_ADU":
        mean_phase_field += "_vv"

    # prepare for exit
    new_dataset = {"radar_out": pyart.util.radar_from_spectra(radar)}
    new_dataset["radar_out"].add_field(mean_phase_field, mph)

    return new_dataset, ind_rad


def process_differential_phase_iq(procstatus, dscfg, radar_list=None):
    """
    Computes the differential phase from the horizontal and vertical IQ data

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain
            "IQhhADU" and "IQvvADU"
        phase_offset : float. Dataset keyword
            The system differential phase offset to remove
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field uPhiDP
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    signal_h_field = None
    signal_v_field = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == "IQhhADU":
            signal_h_field = get_fieldname_pyart(datatype)
        elif datatype == "IQvvADU":
            signal_v_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if signal_h_field not in radar.fields or signal_v_field not in radar.fields:
        warn("Unable to obtain PhiDP. Missing fields")
        return None, None

    phase_offset = dscfg.get("phase_offset", 0.0)

    uphidp = pyart.retrieve.compute_differential_phase_iq(
        radar,
        phase_offset=phase_offset,
        signal_h_field=signal_h_field,
        signal_v_field=signal_v_field,
    )

    # prepare for exit
    new_dataset = {"radar_out": pyart.util.radar_from_spectra(radar)}
    new_dataset["radar_out"].add_field("uncorrected_differential_phase", uphidp)

    return new_dataset, ind_rad


def process_rhohv_iq(procstatus, dscfg, radar_list=None):
    """
    Computes RhoHV from the horizontal and vertical IQ data

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain
            "IQhhADU", "IQvvADU", "IQNADUh" and "IQNADUv"
        subtract_noise : Bool
            If True noise will be subtracted from the signal
        lag : int
            Time lag used in the computation
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "RhoHV"
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    noise_h_field = None
    noise_v_field = None
    signal_h_field = None
    signal_v_field = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == "IQhhADU":
            signal_h_field = get_fieldname_pyart(datatype)
        elif datatype == "IQvvADU":
            signal_v_field = get_fieldname_pyart(datatype)
        elif datatype == "IQNADUh":
            noise_h_field = get_fieldname_pyart(datatype)
        elif datatype == "IQNADUv":
            noise_v_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if signal_h_field not in radar.fields or signal_v_field not in radar.fields:
        warn("Unable to obtain RhoHV. Missing fields")
        return None, None

    subtract_noise = dscfg.get("subtract_noise", False)
    lag = dscfg.get("lag", 0)

    rhohv = pyart.retrieve.compute_rhohv_iq(
        radar,
        subtract_noise=subtract_noise,
        lag=lag,
        signal_h_field=signal_h_field,
        signal_v_field=signal_v_field,
        noise_h_field=noise_h_field,
        noise_v_field=noise_v_field,
    )

    # prepare for exit
    new_dataset = {"radar_out": pyart.util.radar_from_spectra(radar)}
    new_dataset["radar_out"].add_field("cross_correlation_ratio", rhohv)

    return new_dataset, ind_rad


def process_Doppler_velocity_iq(procstatus, dscfg, radar_list=None):
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
            The input data types, must contain
            "IQhhADU" or "IQvvADU"
        direction : str
            The convention used in the Doppler mean field. Can be
            negative_away or negative_towards
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "V"
        (if IQhhADU was provided) or "Vv" (if IQvvADU was provided)
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("IQhhADU", "IQvvADU"):
            signal_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if signal_field not in radar.fields:
        warn("Unable to obtain Doppler velocity. " + "Missing field " + signal_field)
        return None, None

    direction = dscfg.get("direction", "negative_away")

    vel = pyart.retrieve.compute_Doppler_velocity_iq(
        radar, signal_field=signal_field, direction=direction
    )

    vel_field = "velocity"
    if signal_field == "IQ_vv_ADU":
        vel_field += "_vv"

    # prepare for exit
    new_dataset = {"radar_out": pyart.util.radar_from_spectra(radar)}
    new_dataset["radar_out"].add_field(vel_field, vel)

    return new_dataset, ind_rad


def process_Doppler_width_iq(procstatus, dscfg, radar_list=None):
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
            The input data types, must contain
            "IQhhADU" or "IQvvADU"
        subtract_noise : Bool
            If True noise will be subtracted from the signals
        lag : int
            Time lag used in the denominator of the computation
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "W" (if IQhhADU was provided),
        or "Wv" (if IQvvADU was provided)
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    noise_field = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("IQhhADU", "IQvvADU"):
            signal_field = get_fieldname_pyart(datatype)
        elif datatype in ("IQNADUh", "IQNADUv"):
            noise_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if signal_field not in radar.fields:
        warn(
            "Unable to obtain Doppler spectrum width. "
            + "Missing field "
            + signal_field
        )
        return None, None

    dscfg.get("subtract_noise", False)
    lag = dscfg.get("lag", 1)

    width = pyart.retrieve.compute_Doppler_width_iq(
        radar,
        subtract_noise=True,
        signal_field=signal_field,
        noise_field=noise_field,
        lag=lag,
    )

    width_field = "spectrum_width"
    if signal_field == "IQ_vv_ADU":
        width_field += "_vv"

    # prepare for exit
    new_dataset = {"radar_out": pyart.util.radar_from_spectra(radar)}
    new_dataset["radar_out"].add_field(width_field, width)

    return new_dataset, ind_rad


def process_fft(procstatus, dscfg, radar_list=None):
    """
    Compute the Doppler spectra form the IQ data with a Fourier transform

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain,
            "IQNdBADUh" and/or "IQNdBADUv" and/or
            "IQNADUh" and/or "IQNADUv" (see new_dataset below)
        window : list of str
            Parameters of the window used to obtain the spectra. The
            parameters are the ones corresponding to function
            scipy.signal.windows.get_window. It can also be ['None'].
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output fields
        "ShhADUu" (unfiltered_complex_spectra_hh_ADU) if IQNdBADUh was provided,
        "SvvADUu" (unfiltered_complex_spectra_vv_ADU) if IQNdBADUv was provided,
        "sNADUh" (spectral_noise_power_hh_ADU) if IQNADUh was provided,
        "sNADUv" (spectral_noise_power_vv_ADU) if IQNADUv was provided,
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
        if field_name == "IQ_hh_ADU":
            fields_out_list.append("unfiltered_complex_spectra_hh_ADU")
        elif field_name == "IQ_vv_ADU":
            fields_out_list.append("unfiltered_complex_spectra_vv_ADU")
        elif field_name == "IQ_noise_power_hh_ADU":
            fields_out_list.append("spectral_noise_power_hh_ADU")
        elif field_name == "IQ_noiseADU_vv":
            fields_out_list.append("spectral_noise_power_vv_ADU")
        else:
            warn(field_name + " can not be Fourier transformed")
        fields_in_list.append(field_name)

    radar_out = pyart.retrieve.compute_spectra(
        radar, fields_in_list, fields_out_list, window=window
    )

    # prepare for exit
    new_dataset = {"radar_out": radar_out}

    return new_dataset, ind_rad
