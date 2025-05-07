"""
pyrad.proc.process_monitoring
=============================

Functions for monitoring of the polarimetric variables

.. autosummary::
    :toctree: generated/

    process_selfconsistency_kdp_phidp
    process_selfconsistency_bias
    process_selfconsistency_bias2
    process_estimate_phidp0
    process_rhohv_rain
    process_zdr_precip
    process_zdr_snow
    process_monitoring

"""

from copy import deepcopy
from warnings import warn
import numpy as np

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..io.read_data_other import read_selfconsistency
from ..io.read_data_sensor import read_fzl_igra
from ..io.read_data_radar import interpol_field

from ..util.radar_utils import get_histogram_bins
from ..util.stat_utils import ratio_bootstrapping


def process_selfconsistency_kdp_phidp(procstatus, dscfg, radar_list=None):
    """
    Computes specific differential phase and differential phase in rain using
    the selfconsistency between Zdr, Zh and KDP

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of strings. Dataset keyword
            The input data types, must contain
            "dBZ" or "dBZc", and,
            "ZDR" or "ZDRc", and,
            "PhiDP" or "PhiDPc", and,
            "uRhoHV" or "RhoHV" or "RhoHVc", and,
            "TEMP" (Optional), and,
            "H_ISO0" (Optional), and,
            "hydro" (Optional, only used if filter_rain)
        parametrization : str
            The type of parametrization for the self-consistency curves. Can
            be 'None', 'Gourley', 'Wolfensberger', 'Louf', 'Gorgucci' or
            'Vaccarono'
            'None' will use tables from config files. Default 'None'.
        rsmooth : float. Dataset keyword
            length of the smoothing window [m]. Default 2000.
        min_rhohv : float. Dataset keyword
            minimum valid RhoHV. Default 0.92
        filter_rain : Bool. Dataset keyword
            If True the hydrometeor classification is used to filter out gates
            that are not rain. Default True
        max_phidp : float. Dataset keyword
            maximum valid PhiDP [deg]. Default 20.
        ml_thickness : float. Dataset keyword
            assumed melting layer thickness [m]. Default 700.
        fzl : float. Dataset keyword
            The default freezing level height. It will be used if no
            temperature field name is specified or the temperature field is
            not in the radar object. Default 2000.
        sounding : str. Dataset keyword
            The nearest radiosounding WMO code (5 int digits). It will be used to
            compute the freezing level, if no temperature field name is specified,
            if the temperature field isin the radar object or if no freezing_level
            is explicitely defined.
        frequency : float. Dataset keyword
            the radar frequency [Hz]. If None that of the key
            frequency in attribute instrument_parameters of the radar
            object will be used. If the key or the attribute are not present
            the selfconsistency will not be computed
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output fields KDP and PhiDP
        ""
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    temp = None
    iso0 = None
    hydro = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == "dBZc":
            refl = "corrected_reflectivity"
        if datatype == "dBZ":
            refl = "reflectivity"
        if datatype == "ZDRc":
            zdr = "corrected_differential_reflectivity"
        if datatype == "ZDR":
            zdr = "differential_reflectivity"
        if datatype == "PhiDPc":
            phidp = "corrected_differential_phase"
        if datatype == "PhiDP":
            phidp = "differential_phase"
        if datatype == "TEMP":
            temp = "temperature"
        if datatype == "H_ISO0":
            iso0 = "height_over_iso0"
        if datatype == "hydro":
            hydro = "radar_echo_classification"
        if datatype == "RhoHV":
            rhohv = "cross_correlation_ratio"
        if datatype == "RhoHVc":
            rhohv = "corrected_cross_correlation_ratio"
        if datatype == "uRhoHV":
            rhohv = "uncorrected_cross_correlation_ratio"

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if (
        (refl not in radar.fields)
        or (zdr not in radar.fields)
        or (phidp not in radar.fields)
        or (rhohv not in radar.fields)
    ):
        warn("Unable to estimate PhiDP/KDP using selfconsistency. " + "Missing data")
        return None, None

    fzl = None
    temp_ref = "fixed_fzl"
    # determine which freezing level reference
    if temp is not None:
        if temp in radar.fields:
            temp_ref = "temperature"
        else:
            warn(
                "COSMO temperature field not available. "
                + "Using fixed freezing level height to determine liquid phase"
            )
    elif iso0 is not None:
        if iso0 in radar.fields:
            temp_ref = "height_over_iso0"
        else:
            warn(
                "Height over iso0 field not available. "
                + "Using fixed freezing level height to determine liquid phase"
            )
    else:
        # determine freezing level height if necessary
        temp_ref = "fixed_fzl"
        if "fzl" in dscfg:
            fzl = dscfg["fzl"]
        elif "sounding" in dscfg:
            sounding_code = dscfg["sounding"]
            t0 = pyart.util.datetime_utils.datetime_from_radar(radar)
            fzl = read_fzl_igra(sounding_code, t0, dscfg=dscfg)
        else:
            warn("Freezing level height not defined. Using default " + str(fzl) + " m")
            fzl = 2000

    # get self-consistency parametrization or curves
    parametrization = dscfg.get("parametrization", "None")
    if dscfg["initialized"] == 0:
        # get frequency band
        freq = dscfg.get("frequency", None)
        if freq is None:
            if (
                radar.instrument_parameters is not None
                and "frequency" in radar.instrument_parameters
            ):
                freq = radar.instrument_parameters["frequency"]["data"][0]
            else:
                warn(
                    "Unable to estimate Zh bias using "
                    + "self-consistency. Unknown radar frequency"
                )
                return None, None
        freq_band = pyart.retrieve.get_freq_band(freq)

        if parametrization == "None":
            # find unique elevations
            el_vec = np.unique(
                (10.0 * np.round(radar.elevation["data"], decimals=1)).astype(int)
            )
            zdr_kdpzh_list = list()
            el_list = list()
            for el in el_vec:
                if "selfconsistencypath" in dscfg:
                    selfcons_dir = dscfg["selfconsistencypath"]
                else:
                    selfcons_dir = dscfg["configpath"] + "selfconsistency/"
                fname = (
                    selfcons_dir
                    + "selfconsistency_zdr_zhkdp_"
                    + freq_band
                    + "band_temp10_elev"
                    + "{:03d}".format(el)
                    + "_mu05.txt"
                )
                zdr_kdpzh_table = read_selfconsistency(fname)
                if zdr_kdpzh_table is not None:
                    zdr_kdpzh_list.append(zdr_kdpzh_table)
                    el_list.append((el / 10.0).astype(int))
            if not el_list:
                warn(
                    "Unable to retrieve PhiDP and KDP using self-consistency. "
                    + "No selfconsistency files for the radar elevations."
                )

                return None, None
            zdr_kdpzh_dict = {
                "zdr_kdpzh": zdr_kdpzh_list,
                "elev": el_list,
                "freq_band": freq_band,
            }
        else:
            zdr_kdpzh_dict = {"zdr_kdpzh": None, "elev": None, "freq_band": freq_band}

        dscfg["global_data"] = zdr_kdpzh_dict
        dscfg["initialized"] = 1

    if dscfg["initialized"] == 1:
        # get user defined values
        rsmooth = dscfg.get("rsmooth", 2000.0)
        min_rhohv = dscfg.get("min_rhohv", 0.92)
        filter_rain = dscfg.get("filter_rain", True)
        max_phidp = dscfg.get("max_phidp", 20.0)
        ml_thickness = dscfg.get("ml_thickness", 700.0)

        kdpsim_field = "specific_differential_phase"
        phidpsim_field = "differential_phase"
        r_res = radar.range["data"][1] - radar.range["data"][0]
        smooth_wind_len = int(rsmooth / r_res)

        kdpsim, phidpsim = pyart.correct.selfconsistency_kdp_phidp(
            radar,
            dscfg["global_data"],
            min_rhohv=min_rhohv,
            filter_rain=filter_rain,
            max_phidp=max_phidp,
            smooth_wind_len=smooth_wind_len,
            doc=15,
            fzl=fzl,
            thickness=ml_thickness,
            parametrization=parametrization,
            refl_field=refl,
            phidp_field=phidp,
            zdr_field=zdr,
            temp_field=temp,
            iso0_field=iso0,
            hydro_field=hydro,
            rhohv_field=rhohv,
            kdpsim_field=kdpsim_field,
            phidpsim_field=phidpsim_field,
            temp_ref=temp_ref,
        )

        # prepare for exit
        new_dataset = {"radar_out": deepcopy(radar)}
        new_dataset["radar_out"].fields = dict()

        new_dataset["radar_out"].add_field(kdpsim_field, kdpsim)
        new_dataset["radar_out"].add_field(phidpsim_field, phidpsim)

        return new_dataset, ind_rad


def process_selfconsistency_bias(procstatus, dscfg, radar_list=None):
    """
    Estimates the reflectivity bias by means of the selfconsistency
    algorithm by Gourley

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain
            "dBZ" or "dBZc", and,
            "ZDR" or "ZDRc", and,
            "PhiDP" or "PhiDPc", and,
            "uRhoHV" or "RhoHV" or "RhoHVc", and,
            "TEMP" (Optional), and,
            "H_ISO0" (Optional), and,
            "hydro" (Optional, only used if filter_rain)
        parametrization : str
            The type of parametrization for the self-consistency curves. Can
            be 'None', 'Gourley', 'Wolfensberger', 'Louf', 'Gorgucci' or
            'Vaccarono'
            'None' will use tables from config files. Default 'None'.
        fzl : float. Dataset keyword
            Default freezing level height. Default 2000.
        sounding : str. Dataset keyword
            The nearest radiosounding WMO code (5 int digits). It will be used to
            compute the freezing level, if no temperature field name is specified,
            if the temperature field isin the radar object or if no freezing_level
            is explicitely defined.
        rsmooth : float. Dataset keyword
            length of the smoothing window [m]. Default 2000.
        min_rhohv : float. Dataset keyword
            minimum valid RhoHV. Default 0.92
        filter_rain : Bool. Dataset keyword
            If True the hydrometeor classification is used to filter out gates
            that are not rain. Default True
        max_phidp : float. Dataset keyword
            maximum valid PhiDP [deg]. Default 20.
        ml_thickness : float. Dataset keyword
            Melting layer thickness [m]. Default 700.
        rcell : float. Dataset keyword
            length of continuous precipitation to consider the precipitation
            cell a valid phidp segment [m]. Default 15000.
        dphidp_min : float. Dataset keyword
            minimum phase shift [deg]. Default 2.
        dphidp_max : float. Dataset keyword
            maximum phase shift [deg]. Default 16.
        frequency : float. Dataset keyword
            the radar frequency [Hz]. If None that of the key
            frequency in attribute instrument_parameters of the radar
            object will be used. If the key or the attribute are not present
            the selfconsistency will not be computed
        check_wet_radome : Bool. Dataset keyword
            if True the average reflectivity of the closest gates to the radar
            is going to be check to find out whether there is rain over the
            radome. If there is rain no bias will be computed. Default True.
        wet_radome_refl : Float. Dataset keyword
            Average reflectivity [dBZ] of the gates close to the radar to
            consider the radome as wet. Default 25.
        wet_radome_rng_min, wet_radome_rng_max : Float. Dataset keyword
            Min and max range [m] of the disk around the radar used to compute
            the average reflectivity to determine whether the radome is wet.
            Default 2000 and 4000.
        wet_radome_ngates_min : int
            Minimum number of valid gates to consider that the radome is wet.
            Default 180
        valid_gates_only : Bool
            If True the reflectivity bias obtained for each valid ray is going
            to be assigned only to gates of the segment used. That will give
            more weight to longer segments when computing the total bias.
            Default False
        keep_points : Bool
            If True the ZDR, ZH and KDP of the gates used in the self-
            consistency algorithm are going to be stored for further analysis.
            Default False
        rkdp : float
            The length of the window used to compute KDP with the single
            window least square method [m]. Default 6000.
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "dBZ_bias" (refl. bias)
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        return None, None

    keep_points = dscfg.get("keep_points", False)
    if procstatus == 2:
        if not keep_points or dscfg["initialized"] == 0:
            return None, None

        return (
            {"selfconsistency_points": dscfg["global_data"]["selfconsistency_points"]},
            None,
        )

    temp = None
    iso0 = None
    hydro = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == "dBZc":
            refl = "corrected_reflectivity"
        if datatype == "dBZ":
            refl = "reflectivity"
        if datatype == "ZDRc":
            zdr = "corrected_differential_reflectivity"
        if datatype == "ZDR":
            zdr = "differential_reflectivity"
        if datatype == "PhiDPc":
            phidp = "corrected_differential_phase"
        if datatype == "PhiDP":
            phidp = "differential_phase"
        if datatype == "TEMP":
            temp = "temperature"
        if datatype == "H_ISO0":
            iso0 = "height_over_iso0"
        if datatype == "hydro":
            hydro = "radar_echo_classification"
        if datatype == "RhoHV":
            rhohv = "cross_correlation_ratio"
        if datatype == "RhoHVc":
            rhohv = "corrected_cross_correlation_ratio"
        if datatype == "uRhoHV":
            rhohv = "uncorrected_cross_correlation_ratio"

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if (
        (refl not in radar.fields)
        or (zdr not in radar.fields)
        or (phidp not in radar.fields)
        or (rhohv not in radar.fields)
    ):
        warn(
            "Unable to estimate reflectivity bias using selfconsistency. "
            + "Missing data"
        )
        return None, None

    fzl = None
    temp_ref = "fixed_fzl"
    # determine which freezing level reference
    if temp is not None:
        if temp in radar.fields:
            temp_ref = "temperature"
        else:
            warn(
                "COSMO temperature field not available. "
                + "Using fixed freezing level height to determine liquid phase"
            )
    elif iso0 is not None:
        if iso0 in radar.fields:
            temp_ref = "height_over_iso0"
        else:
            warn(
                "Height over iso0 field not available. "
                + "Using fixed freezing level height to determine liquid phase"
            )
    else:
        # determine freezing level height if necessary
        temp_ref = "fixed_fzl"
        if "fzl" in dscfg:
            fzl = dscfg["fzl"]
        elif "sounding" in dscfg:
            sounding_code = dscfg["sounding"]
            t0 = pyart.util.datetime_utils.datetime_from_radar(radar)
            fzl = read_fzl_igra(sounding_code, t0, dscfg=dscfg)
        else:
            warn("Freezing level height not defined. Using default " + str(fzl) + " m")
            fzl = 2000

    # get self-consistency parametrization or curves
    parametrization = dscfg.get("parametrization", "None")
    if dscfg["initialized"] == 0:
        # get frequency band
        freq = dscfg.get("frequency", None)
        if freq is None:
            if (
                radar.instrument_parameters is not None
                and "frequency" in radar.instrument_parameters
            ):
                freq = radar.instrument_parameters["frequency"]["data"][0]
            else:
                warn(
                    "Unable to estimate Zh bias using "
                    + "self-consistency. Unknown radar frequency"
                )
                return None, None
        freq_band = pyart.retrieve.get_freq_band(freq)

        if parametrization == "None":
            # find unique elevations
            el_vec = np.unique(
                (10.0 * np.round(radar.elevation["data"], decimals=1)).astype(int)
            )
            zdr_kdpzh_list = list()
            el_list = list()
            for el in el_vec:
                if "selfconsistencypath" in dscfg:
                    selfcons_dir = dscfg["selfconsistencypath"]
                else:
                    selfcons_dir = dscfg["configpath"] + "selfconsistency/"
                fname = (
                    selfcons_dir
                    + "selfconsistency_zdr_zhkdp_"
                    + freq_band
                    + "band_temp10_elev"
                    + "{:03d}".format(el)
                    + "_mu05.txt"
                )
                zdr_kdpzh_table = read_selfconsistency(fname)
                if zdr_kdpzh_table is not None:
                    zdr_kdpzh_list.append(zdr_kdpzh_table)
                    el_list.append((el / 10.0).astype(int))
            if not el_list:
                warn(
                    "Unable to retrieve Zh bias using self-consistency. "
                    + "No selfconsistency files for the radar elevations."
                )

                return None, None
            zdr_kdpzh_dict = {
                "zdr_kdpzh": zdr_kdpzh_list,
                "elev": el_list,
                "freq_band": freq_band,
            }
        else:
            zdr_kdpzh_dict = {"zdr_kdpzh": None, "elev": None, "freq_band": freq_band}

        dscfg["global_data"] = {"zdr_kdpzh_dict": zdr_kdpzh_dict}
        if keep_points:
            dscfg["global_data"].update(
                {
                    "selfconsistency_points": {
                        "zdr": [],
                        "kdp": [],
                        "zh": [],
                        "timeinfo": dscfg["timeinfo"],
                        "parametrization": parametrization,
                        "zdr_kdpzh_dict": zdr_kdpzh_dict,
                    }
                }
            )
        dscfg["initialized"] = 1

    if dscfg["initialized"] == 1:
        # get user defined values
        rsmooth = dscfg.get("rsmooth", 2000.0)
        min_rhohv = dscfg.get("min_rhohv", 0.92)
        filter_rain = dscfg.get("filter_rain", True)
        max_phidp = dscfg.get("max_phidp", 20.0)
        ml_thickness = dscfg.get("ml_thickness", 700.0)
        rcell = dscfg.get("rcell", 15000.0)
        dphidp_min = dscfg.get("dphidp_min", 2.0)
        dphidp_max = dscfg.get("dphidp_max", 16.0)
        check_wet_radome = dscfg.get("check_wet_radome", True)
        wet_radome_refl = dscfg.get("wet_radome_refl", 25.0)
        wet_radome_rng_min = dscfg.get("wet_radome_rng_min", 2000.0)
        wet_radome_rng_max = dscfg.get("wet_radome_rng_max", 4000.0)
        dscfg.get("wet_radome_ngates_min", 180)
        valid_gates_only = dscfg.get("valid_gates_only", False)
        rkdp = dscfg.get("rkdp", 6000.0)

        r_res = radar.range["data"][1] - radar.range["data"][0]
        smooth_wind_len = int(rsmooth / r_res)
        kdp_wind_len = int(rkdp / r_res)
        min_rcons = int(rcell / r_res)
        wet_radome_len_min = int(wet_radome_rng_min / r_res)
        wet_radome_len_max = int(wet_radome_rng_max / r_res)

        refl_bias, selfconsistency_dict = pyart.correct.selfconsistency_bias(
            radar,
            dscfg["global_data"]["zdr_kdpzh_dict"],
            min_rhohv=min_rhohv,
            filter_rain=filter_rain,
            max_phidp=max_phidp,
            smooth_wind_len=smooth_wind_len,
            doc=15,
            fzl=fzl,
            thickness=ml_thickness,
            min_rcons=min_rcons,
            dphidp_min=dphidp_min,
            dphidp_max=dphidp_max,
            parametrization=parametrization,
            refl_field=refl,
            phidp_field=phidp,
            zdr_field=zdr,
            temp_field=temp,
            iso0_field=iso0,
            hydro_field=hydro,
            rhohv_field=rhohv,
            temp_ref=temp_ref,
            check_wet_radome=check_wet_radome,
            wet_radome_refl=wet_radome_refl,
            wet_radome_len_min=wet_radome_len_min,
            wet_radome_len_max=wet_radome_len_max,
            valid_gates_only=valid_gates_only,
            keep_points=keep_points,
            kdp_wind_len=kdp_wind_len,
        )

        if keep_points:
            if selfconsistency_dict is not None:
                dscfg["global_data"]["selfconsistency_points"]["zdr"].extend(
                    selfconsistency_dict["zdr"]
                )
                dscfg["global_data"]["selfconsistency_points"]["zh"].extend(
                    selfconsistency_dict["zh"]
                )
                dscfg["global_data"]["selfconsistency_points"]["kdp"].extend(
                    selfconsistency_dict["kdp"]
                )

        # prepare for exit
        new_dataset = {"radar_out": deepcopy(radar)}
        new_dataset["radar_out"].fields = dict()

        new_dataset["radar_out"].add_field("reflectivity_bias", refl_bias)

        return new_dataset, ind_rad


def process_selfconsistency_bias2(procstatus, dscfg, radar_list=None):
    """
    Estimates the reflectivity bias by means of the selfconsistency
    algorithm by Gourley

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain
            "dBZ" or "dBZc", and,
            "ZDR" or "ZDRc", and,
            "PhiDP" or "PhiDPc", and,
            "uRhoHV" or "RhoHV" or "RhoHVc", and,
            "TEMP" (Optional), and,
            "H_ISO0" (Optional), and,
            "hydro" (Optional, only used if filter_rain)
        parametrization : str
            The type of parametrization for the self-consistency curves. Can
            be 'None', 'Gourley', 'Wolfensberger', 'Louf', 'Gorgucci' or
            'Vaccarono'
            'None' will use tables from config files. Default 'None'.
        fzl : float. Dataset keyword
            Default freezing level height. Default 2000.
        sounding : str. Dataset keyword
            The nearest radiosounding WMO code (5 int digits). It will be used to
            compute the freezing level, if no temperature field name is specified,
            if the temperature field isin the radar object or if no freezing_level
            is explicitely defined.
        rsmooth : float. Dataset keyword
            length of the smoothing window [m]. Default 2000.
        min_rhohv : float. Dataset keyword
            minimum valid RhoHV. Default 0.92
        filter_rain : Bool. Dataset keyword
            If True the hydrometeor classification is used to filter out gates
            that are not rain. Default True
        max_phidp : float. Dataset keyword
            maximum valid PhiDP [deg]. Default 20.
        ml_thickness : float. Dataset keyword
            Melting layer thickness [m]. Default 700.
        rcell : float. Dataset keyword
            length of continuous precipitation to consider the precipitation
            cell a valid phidp segment [m]. Default 15000.
        frequency : float. Dataset keyword
            the radar frequency [Hz]. If None that of the key
            frequency in attribute instrument_parameters of the radar
            object will be used. If the key or the attribute are not present
            the selfconsistency will not be computed
        check_wet_radome : Bool. Dataset keyword
            if True the average reflectivity of the closest gates to the radar
            is going to be check to find out whether there is rain over the
            radome. If there is rain no bias will be computed. Default True.
        wet_radome_refl : Float. Dataset keyword
            Average reflectivity [dBZ] of the gates close to the radar to
            consider the radome as wet. Default 25.
        wet_radome_rng_min, wet_radome_rng_max : Float. Dataset keyword
            Min and max range [m] of the disk around the radar used to compute
            the average reflectivity to determine whether the radome is wet.
            Default 2000 and 4000.
        wet_radome_ngates_min : int
            Minimum number of valid gates to consider that the radome is wet.
            Default 180
        keep_points : Bool
            If True the ZDR, ZH and KDP of the gates used in the self-
            consistency algorithm are going to be stored for further analysis.
            Default False
        bias_per_gate : Bool
            If True the bias per gate will be computed
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "dBZ_bias" (refl. bias)
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        return None, None

    keep_points = dscfg.get("keep_points", False)
    bias_type = dscfg.get("bias_type", "cumulative")
    provide_confidence = dscfg.get("provide_confidence", False)
    nsamples_confidence = dscfg.get("nsamples_confidence", 1000)
    if procstatus == 2:
        if dscfg["initialized"] == 0:
            return None, None

        dataset = None
        if bias_type == "cumulative":
            kdp_obs = np.ma.array(dscfg["global_data"]["kdp_data_dict"]["kdp_obs"])
            kdp_sim = np.ma.array(dscfg["global_data"]["kdp_data_dict"]["kdp_sim"])
            reflectivity_bias = {
                "value": 10.0 * np.ma.log10(np.ma.sum(kdp_sim) / np.ma.sum(kdp_obs)),
                "npoints": kdp_obs.size,
                "timeinfo": dscfg["global_data"]["kdp_data_dict"]["timeinfo"],
                "bias_type": "cumulative",
            }

            if provide_confidence:
                samples = ratio_bootstrapping(
                    kdp_sim, kdp_obs, nsamples=nsamples_confidence
                )
                reflectivity_bias.update({"samples": 10.0 * np.ma.log10(samples)})
            dataset = {"reflectivity_bias": reflectivity_bias}

        if keep_points:
            if dataset is None:
                dataset = {
                    "selfconsistency_points": dscfg["global_data"][
                        "selfconsistency_points"
                    ]
                }
            else:
                dataset.update(
                    {
                        "selfconsistency_points": dscfg["global_data"][
                            "selfconsistency_points"
                        ]
                    }
                )

        return dataset, None

    temp = None
    iso0 = None
    hydro = None
    phidp = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == "dBZc":
            refl = "corrected_reflectivity"
        if datatype == "dBZ":
            refl = "reflectivity"
        if datatype == "ZDRc":
            zdr = "corrected_differential_reflectivity"
        if datatype == "ZDR":
            zdr = "differential_reflectivity"
        if datatype == "PhiDPc":
            phidp = "corrected_differential_phase"
        if datatype == "PhiDP":
            phidp = "differential_phase"
        if datatype == "KDPc":
            kdp = "corrected_specific_differential_phase"
        if datatype == "KDP":
            kdp = "specific_differential_phase"
        if datatype == "TEMP":
            temp = "temperature"
        if datatype == "H_ISO0":
            iso0 = "height_over_iso0"
        if datatype == "hydro":
            hydro = "radar_echo_classification"
        if datatype == "RhoHV":
            rhohv = "cross_correlation_ratio"
        if datatype == "RhoHVc":
            rhohv = "corrected_cross_correlation_ratio"
        if datatype == "uRhoHV":
            rhohv = "uncorrected_cross_correlation_ratio"

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if (
        (refl not in radar.fields)
        or (zdr not in radar.fields)
        or (kdp not in radar.fields)
        or (rhohv not in radar.fields)
    ):
        warn(
            "Unable to estimate reflectivity bias using selfconsistency. "
            + "Missing data"
        )
        return None, None

    fzl = None
    temp_ref = "fixed_fzl"
    # determine which freezing level reference
    if temp is not None:
        if temp in radar.fields:
            temp_ref = "temperature"
        else:
            warn(
                "COSMO temperature field not available. "
                + "Using fixed freezing level height to determine liquid phase"
            )
    elif iso0 is not None:
        if iso0 in radar.fields:
            temp_ref = "height_over_iso0"
        else:
            warn(
                "Height over iso0 field not available. "
                + "Using fixed freezing level height to determine liquid phase"
            )
    else:
        # determine freezing level height if necessary
        temp_ref = "fixed_fzl"
        if "fzl" in dscfg:
            fzl = dscfg["fzl"]
        elif "sounding" in dscfg:
            sounding_code = dscfg["sounding"]
            t0 = pyart.util.datetime_utils.datetime_from_radar(radar)
            fzl = read_fzl_igra(sounding_code, t0, dscfg=dscfg)
        else:
            warn("Freezing level height not defined. Using default " + str(fzl) + " m")
            fzl = 2000

    # get self-consistency parametrization or curves
    parametrization = dscfg.get("parametrization", "None")
    if dscfg["initialized"] == 0:
        # get frequency band
        freq = dscfg.get("frequency", None)
        if freq is None:
            if (
                radar.instrument_parameters is not None
                and "frequency" in radar.instrument_parameters
            ):
                freq = radar.instrument_parameters["frequency"]["data"][0]
            else:
                warn(
                    "Unable to estimate Zh bias using "
                    + "self-consistency. Unknown radar frequency"
                )
                return None, None
        freq_band = pyart.retrieve.get_freq_band(freq)

        if parametrization == "None":
            # find unique elevations
            el_vec = np.unique(
                (10.0 * np.round(radar.elevation["data"], decimals=1)).astype(int)
            )
            zdr_kdpzh_list = list()
            el_list = list()
            for el in el_vec:
                fname = (
                    dscfg["configpath"]
                    + "selfconsistency/"
                    + "selfconsistency_zdr_zhkdp_"
                    + freq_band
                    + "band_temp10_elev"
                    + "{:03d}".format(el)
                    + "_mu05.txt"
                )
                zdr_kdpzh_table = read_selfconsistency(fname)
                if zdr_kdpzh_table is not None:
                    zdr_kdpzh_list.append(zdr_kdpzh_table)
                    el_list.append((el / 10.0).astype(int))
            if not el_list:
                warn(
                    "Unable to retrieve Zh bias using self-consistency. "
                    + "No selfconsistency files for the radar elevations."
                )

                return None, None

            zdr_kdpzh_dict = {
                "zdr_kdpzh": zdr_kdpzh_list,
                "elev": el_list,
                "freq_band": freq_band,
            }
        else:
            zdr_kdpzh_dict = {"zdr_kdpzh": None, "elev": None, "freq_band": freq_band}

        dscfg["global_data"] = {"zdr_kdpzh_dict": zdr_kdpzh_dict}
        dscfg["global_data"].update(
            {
                "kdp_data_dict": {
                    "kdp_obs": [],
                    "kdp_sim": [],
                    "timeinfo": dscfg["timeinfo"],
                }
            }
        )
        if keep_points:
            dscfg["global_data"].update(
                {
                    "selfconsistency_points": {
                        "zdr": [],
                        "kdp": [],
                        "zh": [],
                        "timeinfo": dscfg["timeinfo"],
                        "parametrization": parametrization,
                        "zdr_kdpzh_dict": zdr_kdpzh_dict,
                    }
                }
            )
        dscfg["initialized"] = 1

    if dscfg["initialized"] == 1:
        # get user defined values
        rsmooth = dscfg.get("rsmooth", 2000.0)
        min_rhohv = dscfg.get("min_rhohv", 0.92)
        min_zdr = dscfg.get("min_zdr", 0.2)
        filter_rain = dscfg.get("filter_rain", True)
        max_phidp = dscfg.get("max_phidp", 20.0)
        ml_thickness = dscfg.get("ml_thickness", 700.0)
        rcell = dscfg.get("rcell", 15000.0)
        check_wet_radome = dscfg.get("check_wet_radome", True)
        wet_radome_refl = dscfg.get("wet_radome_refl", 25.0)
        wet_radome_rng_min = dscfg.get("wet_radome_rng_min", 2000.0)
        wet_radome_rng_max = dscfg.get("wet_radome_rng_max", 4000.0)
        dscfg.get("wet_radome_ngates_min", 180)
        bias_per_gate = dscfg.get("bias_per_gate", False)

        r_res = radar.range["data"][1] - radar.range["data"][0]
        smooth_wind_len = int(rsmooth / r_res)
        min_rcons = int(rcell / r_res)
        wet_radome_len_min = int(wet_radome_rng_min / r_res)
        wet_radome_len_max = int(wet_radome_rng_max / r_res)

        kdp_data_dict, refl_bias, selfcons_dict = pyart.correct.selfconsistency_bias2(
            radar,
            dscfg["global_data"]["zdr_kdpzh_dict"],
            min_rhohv=min_rhohv,
            min_zdr=min_zdr,
            filter_rain=filter_rain,
            max_phidp=max_phidp,
            smooth_wind_len=smooth_wind_len,
            doc=15,
            fzl=fzl,
            thickness=ml_thickness,
            min_rcons=min_rcons,
            parametrization=parametrization,
            refl_field=refl,
            phidp_field=phidp,
            zdr_field=zdr,
            temp_field=temp,
            iso0_field=iso0,
            hydro_field=hydro,
            rhohv_field=rhohv,
            kdp_field=kdp,
            temp_ref=temp_ref,
            check_wet_radome=check_wet_radome,
            wet_radome_refl=wet_radome_refl,
            wet_radome_len_min=wet_radome_len_min,
            wet_radome_len_max=wet_radome_len_max,
            keep_points=keep_points,
            bias_per_gate=bias_per_gate,
        )

        if keep_points:
            if selfcons_dict is not None:
                dscfg["global_data"]["selfconsistency_points"]["zdr"].extend(
                    selfcons_dict["zdr"]
                )
                dscfg["global_data"]["selfconsistency_points"]["zh"].extend(
                    selfcons_dict["zh"]
                )
                dscfg["global_data"]["selfconsistency_points"]["kdp"].extend(
                    selfcons_dict["kdp"]
                )

        if kdp_data_dict is not None:
            dscfg["global_data"]["kdp_data_dict"]["kdp_sim"].extend(
                kdp_data_dict["kdp_sim"]
            )
            dscfg["global_data"]["kdp_data_dict"]["kdp_obs"].extend(
                kdp_data_dict["kdp_obs"]
            )

        # prepare for exit
        dataset = None
        if bias_type == "instant":
            reflectivity_bias = {
                "value": np.ma.masked,
                "npoints": 0,
                "timeinfo": dscfg["timeinfo"],
                "bias_type": "instant",
            }

            if kdp_data_dict is not None:
                kdp_obs = np.ma.array(kdp_data_dict["kdp_obs"])
                kdp_sim = np.ma.array(kdp_data_dict["kdp_sim"])

                reflectivity_bias["value"] = 10.0 * np.ma.log10(
                    np.ma.sum(kdp_sim) / np.ma.sum(kdp_obs)
                )
                reflectivity_bias["npoints"] = kdp_obs.size

                if provide_confidence:
                    samples = ratio_bootstrapping(
                        kdp_sim, kdp_obs, iter=nsamples_confidence
                    )
                    reflectivity_bias.update({"samples": 10.0 * np.ma.log10(samples)})

            dataset = {"reflectivity_bias": reflectivity_bias}

        if bias_per_gate:
            if refl_bias is not None:
                if dataset is None:
                    dataset = {"radar_out": deepcopy(radar)}
                else:
                    dataset.update({"radar_out": deepcopy(radar)})
                dataset["radar_out"].fields = dict()
                dataset["radar_out"].add_field("reflectivity_bias", refl_bias)

        return dataset, ind_rad


def process_estimate_phidp0(procstatus, dscfg, radar_list=None):
    """
    estimates the system differential phase offset at each ray

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain
            "dBZ" or "dBZc", and,
            "PhiDP" or "PhiDPc" or "uPhiDP"
        rmin : float. Dataset keyword
            The minimum range where to look for valid data [m]
        rmax : float. Dataset keyword
            The maximum range where to look for valid data [m]
        rcell : float. Dataset keyword
            The length of a continuous cell to consider it valid precip [m]
        Zmin : float. Dataset keyword
            The minimum reflectivity [dBZ]
        Zmax : float. Dataset keyword
            The maximum reflectivity [dBZ]
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output fields "PhiDP0" (system diff. phase) and
        "PhiDP0_bin" (first gate diff. phase)
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == "dBZ":
            refl_field = "reflectivity"
        if datatype == "dBZc":
            refl_field = "corrected_reflectivity"
        if datatype == "PhiDP":
            psidp_field = "differential_phase"
        if datatype == "PhiDPc":
            psidp_field = "corrected_differential_phase"
        if datatype == "uPhiDP":
            psidp_field = "uncorrected_differential_phase"

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if (refl_field not in radar.fields) or (psidp_field not in radar.fields):
        warn("Unable to estimate PhiDP system offset. Missing data")
        return None, None

    ind_rmin = np.where(radar.range["data"] > dscfg["rmin"])[0][0]
    ind_rmax = np.where(radar.range["data"] < dscfg["rmax"])[0][-1]
    r_res = radar.range["data"][1] - radar.range["data"][0]
    min_rcons = int(dscfg["rcell"] / r_res)

    phidp0, first_gates = pyart.correct.det_sys_phase_ray(
        radar,
        ind_rmin=ind_rmin,
        ind_rmax=ind_rmax,
        min_rcons=min_rcons,
        zmin=dscfg["Zmin"],
        zmax=dscfg["Zmax"],
        phidp_field=psidp_field,
        refl_field=refl_field,
    )

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(radar)}
    new_dataset["radar_out"].fields = dict()

    new_dataset["radar_out"].add_field("system_differential_phase", phidp0)
    new_dataset["radar_out"].add_field("first_gate_differential_phase", first_gates)

    return new_dataset, ind_rad


def process_rhohv_rain(procstatus, dscfg, radar_list=None):
    """
    Keeps only suitable data to evaluate the 80 percentile of RhoHV in rain

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain
            "RhoHV" or "RhoHVc" or "uRhoHV", and,
            "dBZ" or "dBZc", and
            "TEMP" (Optional), or
            "H_ISO0" (Optional)
        rmin : float. Dataset keyword
            minimum range where to look for rain [m]. Default 1000.
        rmax : float. Dataset keyword
            maximum range where to look for rain [m]. Default 50000.
        Zmin : float. Dataset keyword
            minimum reflectivity to consider the bin as precipitation [dBZ].
            Default 20.
        Zmax : float. Dataset keyword
            maximum reflectivity to consider the bin as precipitation [dBZ]
            Default 40.
        ml_thickness : float. Dataset keyword
            assumed thickness of the melting layer. Default 700.
        fzl : float. Dataset keyword
            The default freezing level height. It will be used if no
            temperature field name is specified or the temperature field is
            not in the radar object. Default 2000.
        sounding : str. Dataset keyword
            The nearest radiosounding WMO code (5 int digits). It will be used to
            compute the freezing level, if no temperature field name is specified,
            if the temperature field isin the radar object or if no freezing_level
            is explicitely defined.
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "RhoHV_rain" (RhoHV in rain)
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    temp_field = None
    iso0_field = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == "RhoHV":
            rhohv_field = "cross_correlation_ratio"
        if datatype == "RhoHVc":
            rhohv_field = "corrected_cross_correlation_ratio"
        if datatype == "uRhoHV":
            rhohv_field = "uncorrected_cross_correlation_ratio"
        if datatype == "dBZc":
            refl_field = "corrected_reflectivity"
        if datatype == "dBZ":
            refl_field = "reflectivity"
        if datatype == "TEMP":
            temp_field = "temperature"
        if datatype == "H_ISO0":
            iso0_field = "height_over_iso0"

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if (refl_field not in radar.fields) or (rhohv_field not in radar.fields):
        warn("Unable to estimate RhoHV in rain. Missing data")
        return None, None

    fzl = None
    temp_ref = "fixed_fzl"
    # determine which freezing level reference
    if temp_field is not None:
        if temp_field in radar.fields:
            temp_ref = "temperature"
        else:
            warn(
                "COSMO temperature field not available. "
                + "Using fixed freezing level height to determine liquid phase"
            )
    elif iso0_field is not None:
        if iso0_field in radar.fields:
            temp_ref = "height_over_iso0"
        else:
            warn(
                "Height over iso0 field not available. "
                + "Using fixed freezing level height to determine liquid phase"
            )
    else:
        # determine freezing level height if necessary
        temp_ref = "fixed_fzl"
        if "fzl" in dscfg:
            fzl = dscfg["fzl"]
        elif "sounding" in dscfg:
            sounding_code = dscfg["sounding"]
            t0 = pyart.util.datetime_utils.datetime_from_radar(radar)
            fzl = read_fzl_igra(sounding_code, t0, dscfg=dscfg)
        else:
            warn("Freezing level height not defined. Using default " + str(fzl) + " m")
            fzl = 2000

    # default values
    rmin = 1000.0
    rmax = 50000.0
    zmin = 20.0
    zmax = 40.0
    thickness = 700.0

    # user defined values
    if "rmin" in dscfg:
        rmin = dscfg["rmin"]
    if "rmax" in dscfg:
        rmax = dscfg["rmax"]
    if "Zmin" in dscfg:
        zmin = dscfg["Zmin"]
    if "Zmax" in dscfg:
        zmax = dscfg["Zmax"]
    if "ml_thickness" in dscfg:
        thickness = dscfg["ml_thickness"]

    ind_rmin = np.where(radar.range["data"] > rmin)[0][0]
    ind_rmax = np.where(radar.range["data"] < rmax)[0][-1]

    rhohv_rain = pyart.correct.est_rhohv_rain(
        radar,
        ind_rmin=ind_rmin,
        ind_rmax=ind_rmax,
        zmin=zmin,
        zmax=zmax,
        thickness=thickness,
        doc=15,
        fzl=fzl,
        rhohv_field=rhohv_field,
        temp_field=temp_field,
        iso0_field=iso0_field,
        refl_field=refl_field,
        temp_ref=temp_ref,
    )

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(radar)}
    new_dataset["radar_out"].fields = dict()
    new_dataset["radar_out"].add_field("cross_correlation_ratio_in_rain", rhohv_rain)

    return new_dataset, ind_rad


def process_zdr_precip(procstatus, dscfg, radar_list=None):
    """
    Keeps only suitable data to evaluate the differential reflectivity in
    moderate rain or precipitation (for vertical scans)

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of strings. Dataset keyword
            The input data types, must contain
            "dBZ" or "dBZc", and,
            "ZDR" or "ZDRc", and,
            "PhiDP" or "PhiDPc", and,
            "uRhoHV" or "RhoHV" or "RhoHVc", and,
            "TEMP" (Optional), and,
            "H_ISO0" (Optional)
        ml_filter : boolean. Dataset keyword
            indicates if a filter on data in and above the melting layer is
            applied. Default True.
        rmin : float. Dataset keyword
            minimum range where to look for rain [m]. Default 1000.
        rmax : float. Dataset keyword
            maximum range where to look for rain [m]. Default 50000.
        Zmin : float. Dataset keyword
            minimum reflectivity to consider the bin as precipitation [dBZ].
            Default 20.
        Zmax : float. Dataset keyword
            maximum reflectivity to consider the bin as precipitation [dBZ]
            Default 22.
        RhoHVmin : float. Dataset keyword
            minimum RhoHV to consider the bin as precipitation
            Default 0.97
        PhiDPmax : float. Dataset keyword
            maximum PhiDP to consider the bin as precipitation [deg]
            Default 10.
        elmax : float. Dataset keyword
            maximum elevation angle where to look for precipitation [deg]
            Default None.
        ml_thickness : float. Dataset keyword
            assumed thickness of the melting layer. Default 700.
        fzl : float. Dataset keyword
            The default freezing level height. It will be used if no
            temperature field name is specified or the temperature field is
            not in the radar object. Default 2000.
        sounding : str. Dataset keyword
            The nearest radiosounding WMO code (5 int digits). It will be used to
            compute the freezing level, if no temperature field name is specified,
            if the temperature field isin the radar object or if no freezing_level
            is explicitely defined.
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "ZDR_prec"
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    temp_field = None
    iso0_field = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == "ZDR":
            zdr_field = "differential_reflectivity"
        if datatype == "ZDRc":
            zdr_field = "corrected_differential_reflectivity"
        if datatype == "PhiDP":
            phidp_field = "differential_phase"
        if datatype == "PhiDPc":
            phidp_field = "corrected_differential_phase"
        if datatype == "RhoHV":
            rhohv_field = "cross_correlation_ratio"
        if datatype == "RhoHVc":
            rhohv_field = "corrected_cross_correlation_ratio"
        if datatype == "uRhoHV":
            rhohv_field = "uncorrected_cross_correlation_ratio"
        if datatype == "dBZc":
            refl_field = "corrected_reflectivity"
        if datatype == "dBZ":
            refl_field = "reflectivity"
        if datatype == "TEMP":
            temp_field = "temperature"
        if datatype == "H_ISO0":
            iso0_field = "height_over_iso0"

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if (
        (refl_field not in radar.fields)
        or (rhohv_field not in radar.fields)
        or (zdr_field not in radar.fields)
        or (phidp_field not in radar.fields)
    ):
        warn("Unable to estimate ZDR in rain. Missing data")
        return None, None

    # if data in and above the melting layer has to be filtered determine the
    # field to use
    fzl = None
    temp_ref = "fixed_fzl"
    ml_filter = True
    if "ml_filter" in dscfg:
        ml_filter = dscfg["ml_filter"]

    if ml_filter:
        # determine which freezing level reference
        if temp_field is not None:
            if temp_field in radar.fields:
                temp_ref = "temperature"
            else:
                warn(
                    "COSMO temperature field not available. "
                    + "Using fixed freezing level height to determine liquid phase"
                )
        elif iso0_field is not None:
            if iso0_field in radar.fields:
                temp_ref = "height_over_iso0"
            else:
                warn(
                    "Height over iso0 field not available. "
                    + "Using fixed freezing level height to determine liquid phase"
                )
        else:
            # determine freezing level height if necessary
            temp_ref = "fixed_fzl"
            if "fzl" in dscfg:
                fzl = dscfg["fzl"]
            elif "sounding" in dscfg:
                sounding_code = dscfg["sounding"]
                t0 = pyart.util.datetime_utils.datetime_from_radar(radar)
                fzl = read_fzl_igra(sounding_code, t0, dscfg=dscfg)
            else:
                warn(
                    "Freezing level height not defined. Using default "
                    + str(fzl)
                    + " m"
                )
            fzl = 2000
    else:
        temp_ref = None

    # default values
    rmin = 1000.0
    rmax = 50000.0
    zmin = 20.0
    zmax = 22.0
    rhohvmin = 0.97
    phidpmax = 10.0
    elmax = None
    thickness = 700.0

    # user defined values
    if "rmin" in dscfg:
        rmin = dscfg["rmin"]
    if "rmax" in dscfg:
        rmax = dscfg["rmax"]
    if "Zmin" in dscfg:
        zmin = dscfg["Zmin"]
    if "Zmax" in dscfg:
        zmax = dscfg["Zmax"]
    if "RhoHVmin" in dscfg:
        rhohvmin = dscfg["RhoHVmin"]
    if "PhiDPmax" in dscfg:
        phidpmax = dscfg["PhiDPmax"]
    if "elmax" in dscfg:
        elmax = dscfg["elmax"]
    if "ml_thickness" in dscfg:
        thickness = dscfg["ml_thickness"]

    ind_rmin = np.where(radar.range["data"] > rmin)[0][0]
    ind_rmax = np.where(radar.range["data"] < rmax)[0][-1]

    zdr_precip = pyart.correct.est_zdr_precip(
        radar,
        ind_rmin=ind_rmin,
        ind_rmax=ind_rmax,
        zmin=zmin,
        zmax=zmax,
        rhohvmin=rhohvmin,
        phidpmax=phidpmax,
        elmax=elmax,
        thickness=thickness,
        doc=15,
        fzl=fzl,
        zdr_field=zdr_field,
        rhohv_field=rhohv_field,
        phidp_field=phidp_field,
        temp_field=temp_field,
        iso0_field=iso0_field,
        refl_field=refl_field,
        temp_ref=temp_ref,
    )

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(radar)}
    new_dataset["radar_out"].fields = dict()

    new_dataset["radar_out"].add_field(
        "differential_reflectivity_in_precipitation", zdr_precip
    )

    return new_dataset, ind_rad


def process_zdr_snow(procstatus, dscfg, radar_list=None):
    """
    Keeps only suitable data to evaluate the differential reflectivity in
    snow

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of strings. Dataset keyword
            The input data types, must contain
            "dBZ" or "dBZc", and,
            "ZDR" or "ZDRc", and,
            "PhiDP" or "PhiDPc", and,
            "uRhoHV" or "RhoHV" or "RhoHVc", and,
            "hydro", and,
            "TEMP" (Optional), and,
            "SNRh" or "SNRv" (Optional, used to filter with SNRmin and SNRmax)
        rmin : float. Dataset keyword
            minimum range where to look for rain [m]. Default 1000.
        rmax : float. Dataset keyword
            maximum range where to look for rain [m]. Default 50000.
        Zmin : float. Dataset keyword
            minimum reflectivity to consider the bin as snow [dBZ].
            Default 0.
        Zmax : float. Dataset keyword
            maximum reflectivity to consider the bin as snow [dBZ]
            Default 30.
        SNRmin : float. Dataset keyword
            minimum SNR to consider the bin as snow [dB].
            Default 10.
        SNRmax : float. Dataset keyword
            maximum SNR to consider the bin as snow [dB]
            Default 50.
        RhoHVmin : float. Dataset keyword
            minimum RhoHV to consider the bin as snow
            Default 0.97
        PhiDPmax : float. Dataset keyword
            maximum PhiDP to consider the bin as snow [deg]
            Default 10.
        elmax : float. Dataset keyword
            maximum elevation angle where to look for snow [deg]
            Default None.
        KDPmax : float. Dataset keyword
            maximum KDP to consider the bin as snow [deg]
            Default None
        TEMPmin : float. Dataset keyword
            minimum temperature to consider the bin as snow [deg C].
            Default None
        TEMPmax : float. Dataset keyword
            maximum temperature to consider the bin as snow [deg C]
            Default None
        hydroclass : list of ints. Dataset keyword
            list of hydrometeor classes to keep for the analysis
            Default [2] (dry snow)
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "ZDR_snow"
    ind_rad : int
        radar index
    """

    if procstatus != 1:
        return None, None

    temp_field = None
    kdp_field = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == "ZDR":
            zdr_field = "differential_reflectivity"
        if datatype == "ZDRc":
            zdr_field = "corrected_differential_reflectivity"
        if datatype == "PhiDP":
            phidp_field = "differential_phase"
        if datatype == "PhiDPc":
            phidp_field = "corrected_differential_phase"
        if datatype == "RhoHV":
            rhohv_field = "cross_correlation_ratio"
        if datatype == "RhoHVc":
            rhohv_field = "corrected_cross_correlation_ratio"
        if datatype == "uRhoHV":
            rhohv_field = "uncorrected_cross_correlation_ratio"
        if datatype == "dBZc":
            refl_field = "corrected_reflectivity"
        if datatype == "dBZ":
            refl_field = "reflectivity"
        if datatype == "TEMP":
            temp_field = "temperature"
        if datatype == "PhiDP":
            kdp_field = "specific_differential_phase"
        if datatype == "PhiDPc":
            kdp_field = "corrected_specific_differential_phase"
        if datatype == "SNRh":
            snr_field = "signal_to_noise_ratio_hh"
        if datatype == "SNRv":
            snr_field = "signal_to_noise_ratio_vv"
        if datatype == "hydro":
            hydro_field = "radar_echo_classification"

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if (
        (refl_field not in radar.fields)
        or (rhohv_field not in radar.fields)
        or (zdr_field not in radar.fields)
        or (phidp_field not in radar.fields)
        or (hydro_field not in radar.fields)
    ):
        warn("Unable to estimate ZDR in snow. Missing data")
        return None, None

    # User defined values
    rmin = dscfg.get("rmin", 1000.0)
    rmax = dscfg.get("rmax", 50000.0)
    zmin = dscfg.get("Zmin", 0.0)
    zmax = dscfg.get("Zmax", 30.0)
    snrmin = dscfg.get("SNRmin", 10.0)
    snrmax = dscfg.get("SNRmax", 50.0)
    rhohvmin = dscfg.get("RhoHVmin", 0.97)
    phidpmax = dscfg.get("PhiDPmax", 10.0)
    elmax = dscfg.get("elmax", None)
    kdpmax = dscfg.get("KDPmax", None)
    tempmin = dscfg.get("TEMPmin", None)
    tempmax = dscfg.get("TEMPmax", None)
    hydroclass = dscfg.get("hydroclass", [2])

    ind_rmin = np.where(radar.range["data"] > rmin)[0][0]
    ind_rmax = np.where(radar.range["data"] < rmax)[0][-1]

    zdr_snow = pyart.correct.est_zdr_snow(
        radar,
        ind_rmin=ind_rmin,
        ind_rmax=ind_rmax,
        zmin=zmin,
        zmax=zmax,
        snrmin=snrmin,
        snrmax=snrmax,
        rhohvmin=rhohvmin,
        kept_values=hydroclass,
        phidpmax=phidpmax,
        kdpmax=kdpmax,
        tempmin=tempmin,
        tempmax=tempmax,
        elmax=elmax,
        zdr_field=zdr_field,
        rhohv_field=rhohv_field,
        phidp_field=phidp_field,
        temp_field=temp_field,
        snr_field=snr_field,
        hydro_field=hydro_field,
        kdp_field=kdp_field,
        refl_field=refl_field,
    )

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(radar)}
    new_dataset["radar_out"].fields = dict()

    new_dataset["radar_out"].add_field("differential_reflectivity_in_snow", zdr_snow)

    return new_dataset, ind_rad


def process_monitoring(procstatus, dscfg, radar_list=None):
    """
    computes monitoring statistics

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            Arbitrary datatype supported by pyrad
        step : float. Dataset keyword
            The width of the histogram bin. Default is None. In that case the
            default step in function get_histogram_bins is used
        max_rays : int. Dataset keyword
            The maximum number of rays per sweep used when computing the
            histogram. If set above 0 the number of rays per sweep will be
            checked and if above max_rays the last rays of the sweep will be
            removed
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : Radar
        radar object containing histogram data
    ind_rad : int
        radar index

    """

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        for datatypedescr in dscfg["datatype"]:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            field_name = get_fieldname_pyart(datatype)
            break
        ind_rad = int(radarnr[5:8]) - 1
        if radar_list[ind_rad] is None:
            warn("No valid radar")
            return None, None
        radar = radar_list[ind_rad]

        if field_name not in radar.fields:
            warn(field_name + " not available.")
            return None, None

        step = dscfg.get("step", None)
        max_rays = dscfg.get("max_rays", 0)

        bin_edges = get_histogram_bins(field_name, step=step)
        nbins = len(bin_edges) - 1
        step = bin_edges[1] - bin_edges[0]
        bin_centers = bin_edges[:-1] + step / 2.0

        radar_aux = deepcopy(radar)
        if max_rays > 0:
            # remove excess of rays
            ind_above_max = np.where(radar.rays_per_sweep["data"] > max_rays)[0]
            if ind_above_max.size > 0:
                radar_aux.rays_per_sweep["data"][ind_above_max] = max_rays
                for ind in ind_above_max:
                    excess_rays = radar.rays_per_sweep["data"][ind] - max_rays
                    radar_aux.sweep_end_ray_index["data"][ind:] -= excess_rays
                    if ind < radar.nsweeps - 1:
                        radar_aux.sweep_start_ray_index["data"][ind + 1 :] = (
                            radar_aux.sweep_end_ray_index["data"][ind:-1] + 1
                        )

                radar_aux.nrays = np.sum(radar_aux.rays_per_sweep["data"])
                radar_aux.fields[field_name]["data"] = np.ma.masked_all(
                    (radar_aux.nrays, radar_aux.ngates),
                    dtype=radar.fields[field_name]["data"].dtype,
                )
                radar_aux.azimuth["data"] = np.empty(
                    radar_aux.nrays, dtype=radar.azimuth["data"].dtype
                )
                radar_aux.elevation["data"] = np.empty(
                    radar_aux.nrays, dtype=radar.elevation["data"].dtype
                )
                radar_aux.time["data"] = np.empty(
                    radar_aux.nrays, dtype=radar.time["data"].dtype
                )

                for sweep in range(radar.nsweeps):
                    ind_start_old = radar.sweep_start_ray_index["data"][sweep]
                    nrays_sweep = radar_aux.rays_per_sweep["data"][sweep]
                    ind_start_new = radar_aux.sweep_start_ray_index["data"][sweep]
                    ind_end_new = radar_aux.sweep_end_ray_index["data"][sweep]

                    radar_aux.fields[field_name]["data"][
                        ind_start_new : ind_end_new + 1, :
                    ] = radar.fields[field_name]["data"][
                        ind_start_old : ind_start_old + nrays_sweep, :
                    ]
                    radar_aux.azimuth["data"][ind_start_new : ind_end_new + 1] = (
                        radar.azimuth["data"][
                            ind_start_old : ind_start_old + nrays_sweep
                        ]
                    )
                    radar_aux.elevation["data"][ind_start_new : ind_end_new + 1] = (
                        radar.elevation["data"][
                            ind_start_old : ind_start_old + nrays_sweep
                        ]
                    )
                    radar_aux.time["data"][ind_start_new : ind_end_new + 1] = (
                        radar.time["data"][ind_start_old : ind_start_old + nrays_sweep]
                    )

        radar_hist = deepcopy(radar_aux)
        radar_hist.fields = dict()
        radar_hist.range["data"] = bin_centers
        radar_hist.ngates = nbins

        field_dict = pyart.config.get_metadata(field_name)
        field_dict["data"] = np.ma.zeros((radar_aux.nrays, nbins), dtype=int)

        field = deepcopy(radar_aux.fields[field_name]["data"])

        # put gates with values off limits to limit
        mask = np.ma.getmaskarray(field)
        ind = np.where(np.logical_and(mask is False, field < bin_centers[0]))
        field[ind] = bin_centers[0]

        ind = np.where(np.logical_and(mask is False, field > bin_centers[-1]))
        field[ind] = bin_centers[-1]

        for ray in range(radar_aux.nrays):
            field_dict["data"][ray, :], bin_edges = np.histogram(
                field[ray, :].compressed(), bins=bin_edges
            )

        radar_hist.add_field(field_name, field_dict)
        start_time = pyart.graph.common.generate_radar_time_begin(radar_hist)
        # keep histogram in Memory or add to existing histogram
        if dscfg["initialized"] == 0:
            dscfg["global_data"] = {"hist_obj": radar_hist, "timeinfo": start_time}
            dscfg["initialized"] = 1
        else:
            field_interp = interpol_field(
                dscfg["global_data"]["hist_obj"], radar_hist, field_name, fill_value=0
            )
            dscfg["global_data"]["hist_obj"].fields[field_name]["data"] += (
                field_interp["data"].filled(fill_value=0)
            ).astype("int64")

        #    dscfg['global_data']['timeinfo'] = dscfg['timeinfo']

        dataset = dict()
        dataset.update({"hist_obj": radar_hist})
        dataset.update({"hist_type": "instant"})
        dataset.update({"timeinfo": start_time})

        return dataset, ind_rad

    if procstatus == 2:
        if dscfg["initialized"] == 0:
            return None, None

        for datatypedescr in dscfg["datatype"]:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            field_name = get_fieldname_pyart(datatype)
            break
        ind_rad = int(radarnr[5:8]) - 1

        dataset = dict()
        dataset.update({"hist_obj": dscfg["global_data"]["hist_obj"]})
        dataset.update({"hist_type": "cumulative"})
        dataset.update({"timeinfo": dscfg["global_data"]["timeinfo"]})

        return dataset, ind_rad
