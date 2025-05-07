"""
pyrad.io.io_aux
===============

Auxiliary functions for reading/writing files

.. autosummary::
    :toctree: generated/

    get_rad4alp_prod_fname
    map_hydro
    mf_sname_to_wmo_number
    map_Doppler
    get_save_dir
    make_filename
    generate_field_name_str
    get_datatype_metranet
    get_datatype_knmi
    get_datatype_skyecho
    get_datatype_odim
    get_fieldname_pyart
    get_fieldname_icon
    get_field_unit
    get_field_name
    get_scan_files_to_merge
    get_scan_files_to_merge_s3
    get_file_list
    get_file_list_s3
    get_rad4alp_dir
    get_rad4alp_grid_dir
    get_trtfile_list
    get_scan_list
    get_new_rainbow_file_name
    get_datatype_fields
    get_dataset_fields
    get_datetime
    find_raw_icon_file
    find_icon_file
    find_hzt_file
    find_iso0_file
    find_iso0_grib_file
    find_rad4alpicon_file
    find_pyradicon_file
    _get_datetime
    find_date_in_file_name
    convert_pydda_to_pyart_grid

"""

import os
import glob
import re
import datetime
import fnmatch

from warnings import warn
from copy import deepcopy
import numpy as np
import pyart

from pyart.config import get_metadata

try:
    from pyart.aux_io import get_sweep_time_coverage
except ImportError:
    warn("Could not get the get_sweep_time_coverage from pyart")
    warn("You are likely using ARM Py-ART and not the MCH fork")
    warn("You will not be able to read skyecho data")

try:
    import boto3

    _BOTO3_AVAILABLE = True
except ImportError:
    warn("boto3 is not installed, data cannot be retrieved from S3 buckets!")
    _BOTO3_AVAILABLE = False


pyrad_to_pyart_keys_dict = {
    "dBZ": "reflectivity",
    "dBZ_flag": "reflectivity_flag",
    "dBZ_MF": "reflectivity",
    "Zn": "normalized_reflectivity",
    "Zlin": "linear_reflectivity",
    "dBuZ": "unfiltered_reflectivity",
    "dBZc": "corrected_reflectivity",
    "dBuZc": "corrected_unfiltered_reflectivity",
    "dBZv": "reflectivity_vv",
    "dBZvc": "corrected_reflectivity_vv",
    "dBuZv": "unfiltered_reflectivity_vv",
    "dBuZvc": "corrected_unfiltered_reflectivity_vv",
    "dBZhv": "reflectivity_hv",
    "dBZvh": "reflectivity_vh",
    "dBZ_bias": "reflectivity_bias",
    "eta_h": "volumetric_reflectivity",
    "eta_v": "volumetric_reflectivity_vv",
    "rcs_h": "radar_cross_section_hh",
    "rcs_v": "radar_cross_section_vv",
    "VPRcorr": "vpr_correction",
    "VPRFEATURES": "vpr_features",
    "sigma_zh": "sigma_zh",
    "ZDR": "differential_reflectivity",
    "ZDRu": "unfiltered_differential_reflectivity",
    "ZDRc": "corrected_differential_reflectivity",
    "ZDRuc": "corrected_unfiltered_differential_reflectivity",
    "ZDR_prec": "differential_reflectivity_in_precipitation",
    "ZDR_snow": "differential_reflectivity_in_snow",
    "ZDR_col": "differential_reflectivity_column_height",
    "LDRhv": "linear_depolarization_ratio_hv",
    "LDRvh": "linear_depolarization_ratio_vh",
    "dBm": "signal_power_hh",
    "dBmv": "signal_power_vv",
    "Nh": "noisedBZ_hh",
    "Nv": "noisedBZ_vv",
    "NdBADUh": "noisedBADU_hh",
    "NdBADUv": "noisedBADU_vv",
    "NdBmh": "noisedBm_hh",
    "NdBmv": "noisedBm_vv",
    "NADUh": "noiseADU_hh",
    "NADUv": "noiseADU_vv",
    "noise_pos_h": "noise_pos_h",
    "noise_pos_v": "noise_pos_v",
    "Nclip_h": "noise_clipping_level_hh_dBZ",
    "Nclip_v": "noise_clipping_level_vv_dBZ",
    "WBN": "wide_band_noise",
    "WBNc": "corrected_wide_band_noise",
    "ST1": "stat_test_lag1",
    "ST1c": "corrected_stat_test_lag1",
    "ST2": "stat_test_lag2",
    "ST2c": "corrected_stat_test_lag2",
    "TXh": "transmitted_signal_power_h",
    "TXv": "transmitted_signal_power_v",
    "SNRh": "signal_to_noise_ratio_hh",
    "SNRv": "signal_to_noise_ratio_vv",
    "CCORh": "clutter_correction_ratio_hh",
    "CCORv": "clutter_correction_ratio_vv",
    "dBm_sun_hit": "sun_hit_power_h",
    "dBmv_sun_hit": "sun_hit_power_v",
    "ZDR_sun_hit": "sun_hit_differential_reflectivity",
    "dBm_sun_est": "sun_est_power_h",
    "dBmv_sun_est": "sun_est_power_v",
    "ZDR_sun_est": "sun_est_differential_reflectivity",
    "sun_pos_h": "sun_hit_h",
    "sun_pos_v": "sun_hit_v",
    "sun_pos_zdr": "sun_hit_zdr",
    "RhoHV": "cross_correlation_ratio",
    "RhoHVu": "unfiltered_cross_correlation_ratio",
    "uRhoHV": "uncorrected_cross_correlation_ratio",
    "RhoHVc": "corrected_cross_correlation_ratio",
    "RhoHV_rain": "cross_correlation_ratio_in_rain",
    "RhoHV_theo": "theoretical_cross_correlation_ratio",
    "L": "logarithmic_cross_correlation_ratio",
    "CDR": "circular_depolarization_ratio",
    "PhiDP": "differential_phase",
    "uPhiDPu": "uncorrected_unfiltered_differential_phase",
    "uPhiDP": "uncorrected_differential_phase",
    "PhiDPc": "corrected_differential_phase",
    "PhiDP0": "system_differential_phase",
    "PhiDP0_bin": "first_gate_differential_phase",
    "KDP": "specific_differential_phase",
    "uKDP": "uncorrected_specific_differential_phase",
    "KDPc": "corrected_specific_differential_phase",
    "MPH": "mean_phase",
    "MPHc": "corrected_mean_phase",
    "V": "velocity",
    "Vv": "velocity_vv",
    "Vhv": "velocity_hv",
    "Vvh": "velocity_vh",
    "Vu": "unfiltered_velocity",
    "dealV": "dealiased_velocity",
    "Vc": "corrected_velocity",
    "dealVc": "dealiased_corrected_velocity",
    "estV": "retrieved_velocity",
    "stdV": "retrieved_velocity_std",
    "diffV": "velocity_difference",
    "W": "spectrum_width",
    "Wv": "spectrum_width_vv",
    "Whv": "spectrum_width_hv",
    "Wvh": "spectrum_width_vh",
    "Wu": "unfiltered_spectrum_width",
    "Wc": "corrected_spectrum_width",
    "wind_vel_h_az": "azimuthal_horizontal_wind_component",
    "wind_vel_v": "vertical_wind_component",
    "wind_vel_h_u": "eastward_wind_component",
    "wind_vel_h_v": "northward_wind_component",
    "windshear_v": "vertical_wind_shear",
    "WIND_SPEED": "wind_speed",
    "WIND_DIRECTION": "wind_direction",
    "EDR": "turbulence",
    "Ah": "specific_attenuation",
    "Ahc": "corrected_specific_attenuation",
    "PIA": "path_integrated_attenuation",
    "PIAc": "corrected_path_integrated_attenuation",
    "Adp": "specific_differential_attenuation",
    "Adpc": "corrected_specific_differential_attenuation",
    "PIDA": "path_integrated_differential_attenuation",
    "PIDAc": "corrected_path_integrated_differential_attenuation",
    "TEMP": "temperature",
    "ISO0": "iso0",
    "H_ISO0": "height_over_iso0",
    "H_ISO0c": "corrected_height_over_iso0",
    "HZT": "iso0_height",
    "HZTc": "corrected_iso0_height",
    "icon_index": "icon_index",
    "hzt_index": "hzt_index",
    "ml": "melting_layer",
    "VIS": "visibility",
    "HGHT": "height",
    "vis": "visibility",
    "visibility": "visibility",
    "visibility_polar": "visibility_polar",
    "terrain_altitude": "terrain_altitude",
    "bent_terrain_altitude": "bent_terrain_altitude",
    "terrain_slope": "terrain_slope",
    "terrain_aspect": "terrain_aspect",
    "elevation_angle": "elevation_angle",
    "min_vis_elevation": "min_vis_elevation",
    "min_vis_altitude": "min_vis_altitude",
    "min_vis_height_above_ground": "min_vis_height_above_ground",
    "min_rad_vis_height_above_ground": "min_rad_vis_height_above_ground",
    "incident_angle": "incident_angle",
    "sigma_0": "sigma_0",
    "effective_area": "effective_area",
    "dBm_clutter": "dBm_clutter",
    "dBZ_clutter": "dBZ_clutter",
    "echoID": "radar_echo_id",
    "CLT": "clutter_exit_code",
    "occurrence": "occurrence",
    "freq_occu": "frequency_of_occurrence",
    "RR": "radar_estimated_rain_rate",
    "RR_MP": "Marshall_Palmer_radar_estimated_rain_rate",
    "RR_Z": "radar_reflectivity_estimated_rain_rate",
    "RR_KDP": "radar_kdp_estimated_rain_rate",
    "RR_flag": "radar_estimated_rain_rate_flag",
    "RRc": "corrected_radar_estimated_rain_rate",
    "Raccu": "rainfall_accumulation",
    "RaccuMF": "rainfall_accumulation",
    "QIMF": "signal_quality_index",
    "QI": "signal_quality_index",
    "AF": "adjustment_factor",
    "radar_R_rel": "radar_rainrate_relation",
    "prec_type": "precipitation_type",
    "hydro": "radar_echo_classification",
    "hydroMF": "radar_echo_classification_MF",
    "hydroc": "corrected_radar_echo_classification",
    "confidence": "hydroclass_confidence",
    "entropy": "hydroclass_entropy",
    "propAG": "proportion_AG",
    "propCR": "proportion_CR",
    "propLR": "proportion_LR",
    "propRP": "proportion_RP",
    "propRN": "proportion_RN",
    "propVI": "proportion_VI",
    "propWS": "proportion_WS",
    "propMH": "proportion_MH",
    "propIH": "proportion_IH",
    "probAG": "probability_AG",
    "probCR": "probability_CR",
    "probLR": "probability_LR",
    "probRP": "probability_RP",
    "probRN": "probability_RN",
    "probVI": "probability_VI",
    "probWS": "probability_WS",
    "probMH": "probability_MH",
    "probIH": "probability_IH",
    "time_avg_flag": "time_avg_flag",
    "colocated_gates": "colocated_gates",
    "nsamples": "number_of_samples",
    "bird_density": "bird_density",
    "std": "standard_deviation",
    "sum": "sum",
    "sum2": "sum_squared",
    "diff": "fields_difference",
    "mask": "field_mask",
    "texture": "field_texture",
    "ShhADU": "complex_spectra_hh_ADU",
    "ShhADUu": "unfiltered_complex_spectra_hh_ADU",
    "SvvADU": "complex_spectra_vv_ADU",
    "SvvADUu": "unfiltered_complex_spectra_vv_ADU",
    "sPhhADU": "spectral_power_hh_ADU",
    "sPhhADUu": "unfiltered_spectral_power_hh_ADU",
    "sPvvADU": "spectral_power_vv_ADU",
    "sPvvADUu": "unfiltered_spectral_power_vv_ADU",
    "sPhhdBADU": "spectral_power_hh_dBADU",
    "sPhhdBADUu": "unfiltered_spectral_power_hh_dBADU",
    "sPvvdBADU": "spectral_power_vv_dBADU",
    "sPvvdBADUu": "unfiltered_spectral_power_vv_dBADU",
    "sPhhdBm": "spectral_power_hh_dBm",
    "sPhhdBmu": "unfiltered_spectral_power_hh_dBm",
    "sPvvdBm": "spectral_power_vv_dBm",
    "sPvvdBmu": "unfiltered_spectral_power_vv_dBm",
    "sNh": "spectral_noise_power_hh_dBZ",
    "sNv": "spectral_noise_power_vv_dBZ",
    "sNdBADUh": "spectral_noise_power_hh_dBADU",
    "sNdBADUv": "spectral_noise_power_vv_dBADU",
    "sNdBmh": "spectral_noise_power_hh_dBm",
    "sNdBmv": "spectral_noise_power_vv_dBm",
    "sNADUh": "spectral_noise_power_hh_ADU",
    "sNADUv": "spectral_noise_power_vv_ADU",
    "sPhasehh": "spectral_phase_hh",
    "sPhasehhu": "unfiltered_spectral_phase_hh",
    "sPhasevv": "spectral_phase_vv",
    "sPhasevvu": "unfiltered_spectral_phase_vv",
    "sdBZ": "spectral_reflectivity_hh",
    "sdBuZ": "unfiltered_spectral_reflectivity_hh",
    "sdBZv": "spectral_reflectivity_vv",
    "sdBuZv": "unfiltered_spectral_reflectivity_vv",
    "sZDR": "spectral_differential_reflectivity",
    "sZDRu": "unfiltered_spectral_differential_reflectivity",
    "sPhiDP": "spectral_differential_phase",
    "sPhiDPu": "unfiltered_spectral_differential_phase",
    "sRhoHV": "spectral_copolar_correlation_coefficient",
    "sRhoHVu": "unfiltered_spectral_copolar_correlation_coefficient",
    "IQhhADU": "IQ_hh_ADU",
    "IQvvADU": "IQ_vv_ADU",
    "IQNh": "IQ_noise_power_hh_dBZ",
    "IQNv": "IQ_noise_power_vv_dBZ",
    "IQNdBADUh": "IQ_noise_power_hh_dBADU",
    "IQNdBADUv": "IQ_noise_power_vv_dBADU",
    "IQNdBmh": "IQ_noise_power_hh_dBm",
    "IQNdBmv": "IQ_noise_power_vv_dBm",
    "IQNADUh": "IQ_noise_power_hh_ADU",
    "IQNADUv": "IQ_noise_power_vv_ADU",
    "POH": "probability_of_hail",
    "VIL": "vertically_integrated_liquid",
    "ETOP15": "echo_top_15dBZ",
    "ETOP20": "echo_top_20dBZ",
    "ETOP45": "echo_top_45dBZ",
    "ETOP50": "echo_top_50dBZ",
    "MAXECHO": "maximum_echo",
    "HMAXECHO": "maximum_echo_height",
    "AZC01": "rainfall_accumulation",
    "AZC03": "rainfall_accumulation",
    "AZC06": "rainfall_accumulation",
    "aZC01": "rainfall_accumulation",
    "aZC03": "rainfall_accumulation",
    "aZC06": "rainfall_accumulation",
    "CPC0005": "radar_estimated_rain_rate",
    "CPC0060": "rainfall_accumulation",
    "CPC0180": "rainfall_accumulation",
    "CPC0360": "rainfall_accumulation",
    "CPC0720": "rainfall_accumulation",
    "CPC1440": "rainfall_accumulation",
    "CPC2880": "rainfall_accumulation",
    "CPC4320": "rainfall_accumulation",
    "CPCH0005": "radar_estimated_rain_rate",
    "CPCH0060": "rainfall_accumulation",
    "CPCH0180": "rainfall_accumulation",
    "CPCH0360": "rainfall_accumulation",
    "CPCH0720": "rainfall_accumulation",
    "CPCH1440": "rainfall_accumulation",
    "CPCH2880": "rainfall_accumulation",
    "CPCH4320": "rainfall_accumulation",
    "nowpal60_P60": "rainfall_accumulation",
    "nowpal90_P90": "rainfall_accumulation",
    "nowpal180_P180": "rainfall_accumulation",
    "nowpal360_P360": "rainfall_accumulation",
    "nowpal720_P720": "rainfall_accumulation",
    "nowpal90_P30": "rainfall_accumulation",
    "nowpal90_P30_F60": "rainfall_accumulation",
    "nowpal90_F60": "rainfall_accumulation",
    "nowpal180_P60": "rainfall_accumulation",
    "nowpal180_P60_F120": "rainfall_accumulation",
    "nowpal180_F120": "rainfall_accumulation",
    "nowpal360_P120": "rainfall_accumulation",
    "nowpal360_P120_F240": "rainfall_accumulation",
    "nowpal360_F240": "rainfall_accumulation",
    "nowpal720_P360": "rainfall_accumulation",
    "nowpal720_P360_F360": "rainfall_accumulation",
    "nowpal720_F360": "rainfall_accumulation",
    "dACC": "rainfall_accumulation",
    "dACCH": "rainfall_accumulation",
    "dARC": "rainfall_accumulation",
    "RZC": "radar_estimated_rain_rate",
    "R1F": "radar_estimated_rain_rate",
    "rZC": "radar_estimated_rain_rate",
    "RZF": "radar_estimated_rain_rate",
    "dRZC": "radar_estimated_rain_rate",
    "BZC": "probability_of_hail",
    "dBZC": "probability_of_hail",
    "MZC": "maximum_expected_severe_hail_size",
    "dMZC": "maximum_expected_severe_hail_size",
    "GZC": "probability_of_hail",
    "dGZC": "probability_of_hail",
    "CZC": "maximum_echo",
    "dCZC": "maximum_echo",
    "HZC": "maximum_echo_height",
    "EZC15": "echo_top_15dBz",
    "EZC20": "echo_top_20dBz",
    "EZC45": "echo_top_45dBz",
    "EZC50": "echo_top_50dBz",
    "dEZC15": "echo_top_15dBz",
    "dEZC20": "echo_top_20dBz",
    "dEZC45": "echo_top_45dBz",
    "dEZC50": "echo_top_50dBz",
    "LZC": "vertically_integrated_liquid",
    "dLZC": "vertically_integrated_liquid",
    "OZC01": "reflectivity",
    "OZC02": "reflectivity",
    "OZC03": "reflectivity",
    "OZC04": "reflectivity",
    "OZC05": "reflectivity",
    "OZC06": "reflectivity",
    "OZC07": "reflectivity",
    "OZC08": "reflectivity",
    "OZC09": "reflectivity",
    "OZC10": "reflectivity",
    "OZC11": "reflectivity",
    "OZC12": "reflectivity",
    "OZC13": "reflectivity",
    "OZC14": "reflectivity",
    "OZC15": "reflectivity",
    "OZC16": "reflectivity",
    "OZC17": "reflectivity",
    "OZC18": "reflectivity",
    "ff": "wind_speed",
    "dd": "wind_direction",
    "u": "eastward_wind_component",
    "v": "northward_wind_component",
    "w": "vertical_wind_component",
    "width": "height_resolution",
    "gap": "gap",
    "dbz": "bird_reflectivity",
    "eta": "volumetric_reflectivity",
    "dens": "bird_density",
    "n": "number_of_samples_velocity",
    "n_dbz": "number_of_samples_reflectivity",
    "sd_vvp": "retrieved_velocity_std",
    "DBZH": "reflectivity",
    "n_all": "number_of_samples_velocity_all",
    "n_dbz_all": "number_of_samples_reflectivity_all",
    "VOL2BIRD_CLASS": "vol2bird_echo_classification",
    "VOL2BIRD_WEATHER": "vol2bird_weather",
    "VOL2BIRD_BACKGROUND": "vol2bird_background",
    "VOL2BIRD_BIOLOGY": "vol2bird_biology",
    "wind_vel_rad": "radial_wind_speed",
    "wind_vel_rad_filtered": "corrected_radial_wind_speed",
    "wind_vel_rad_ci": "radial_wind_speed_ci",
    "wind_vel_rad_status": "radial_wind_speed_status",
    "WD": "doppler_spectrum_width",
    "WDc": "corrected_doppler_spectrum_width",
    "WD_err": "doppler_spectrum_mean_error",
    "atmos_type": "atmospherical_structures_type",
    "beta_rel": "relative_beta",
    "beta_abs": "absolute_beta",
    "CNR": "cnr",
    "CNRc": "corrected_cnr",
    "HRV": "HRV",
    "VIS006": "VIS006",
    "VIS008": "VIS008",
    "IR_016": "IR_016",
    "IR_039": "IR_039",
    "WV_062": "WV_062",
    "WV_073": "WV_073",
    "IR_087": "IR_087",
    "IR_097": "IR_097",
    "IR_108": "IR_108",
    "IR_120": "IR_120",
    "IR_134": "IR_134",
    "CTH": "CTH",
    "HRV_norm": "HRV_norm",
    "VIS006_norm": "VIS006_norm",
    "VIS008_norm": "VIS008_norm",
    "IR_016_norm": "IR_016_norm",
    "SNR": "signal_to_noise_ratio",
    "VEL": "VEL",
    "RMS": "RMS",
    "LDR": "LDR",
    "NPK": "NPK",
    "SNRgc": "SNRgc",
    "VELgc": "VELgc",
    "RMSgc": "RMSgc",
    "LDRgc": "LDRgc",
    "NPKgc": "NPKgc",
    "SNRg": "SNRg",
    "VELg": "VELg",
    "RMSg": "RMSg",
    "LDRg": "LDRg",
    "NPKg": "NPKg",
    "SNRplank": "SNRplank",
    "VELplank": "VELplank",
    "RMSplank": "RMSplank",
    "LDRplank": "LDRplank",
    "NPKplank": "NPKplank",
    "SNRrain": "SNRrain",
    "VELrain": "VELrain",
    "RMSrain": "RMSrain",
    "LDRrain": "LDRrain",
    "NPKrain": "NPKrain",
    "SNRcl": "SNRcl",
    "VELcl": "VELcl",
    "RMScl": "RMScl",
    "LDRcl": "LDRcl",
    "NPKcl": "NPKcl",
    "SNRice": "SNRice",
    "VELice": "VELice",
    "RMSice": "RMSice",
    "LDRice": "LDRice",
    "NPKice": "NPKice",
    "RHO": "RHO",
    "DPS": "DPS",
    "LDRnormal": "LDRnormal",
    "RHOwav": "RHOwav",
    "DPSwav": "DPSwav",
    "SKWg": "SKWg",
    "HSDco": "HSDco",
    "HSDcx": "HSDcx",
    "Ze": "Ze",
    "Zg": "Zg",
    "Z": "Z",
    "RRcr": "RR",
    "LWCcr": "LWC",
    "TEMPcr": "TEMP",
    "ISDRco": "ISDRco",
    "ISDRcx": "ISDRcx",
    "SNRcx": "SNRcx",
    "SNRCorFaCo": "SNRCorFaCo",
    "avgdBZ": "avg_reflectivity",
    "NdBZ": "npoints_reflectivity",
    "quant05dBZ": "quant05_reflectivity",
    "quant10dBZ": "quant10_reflectivity",
    "quant20dBZ": "quant20_reflectivity",
    "quant50dBZ": "quant50_reflectivity",
    "quant80dBZ": "quant80_reflectivity",
    "quant90dBZ": "quant90_reflectivity",
    "quant95dBZ": "quant95_reflectivity",
    "avgRR": "avg_radar_estimated_rain_rate",
    "NRR": "npoints_radar_estimated_rain_rate",
    "quant05RR": "quant05_radar_estimated_rain_rate",
    "quant10RR": "quant10_radar_estimated_rain_rate",
    "quant20RR": "quant20_radar_estimated_rain_rate",
    "quant50RR": "quant50_radar_estimated_rain_rate",
    "quant80RR": "quant80_radar_estimated_rain_rate",
    "quant90RR": "quant90_radar_estimated_rain_rate",
    "quant95RR": "quant95_radar_estimated_rain_rate",
    "avgV": "avg_velocity",
    "NV": "npoints_velocity",
    "quant05V": "quant05_velocity",
    "quant10V": "quant10_velocity",
    "quant20V": "quant20_velocity",
    "quant50V": "quant50_velocity",
    "quant80V": "quant80_velocity",
    "quant90V": "quant90_velocity",
    "quant95V": "quant95_velocity",
    "avgVc": "avg_corrected_velocity",
    "NVc": "npoints_corrected_velocity",
    "quant05Vc": "quant05_corrected_velocity",
    "quant10Vc": "quant10_corrected_velocity",
    "quant20Vc": "quant20_corrected_velocity",
    "quant50Vc": "quant50_corrected_velocity",
    "quant80Vc": "quant80_corrected_velocity",
    "quant90Vc": "quant90_corrected_velocity",
    "quant95Vc": "quant95_corrected_velocity",
    "avgdealV": "avg_dealiased_velocity",
    "NdealV": "npoints_dealiased_velocity",
    "quant05dealV": "quant05_dealiased_velocity",
    "quant10dealV": "quant10_dealiased_velocity",
    "quant20dealV": "quant20_dealiased_velocity",
    "quant50dealV": "quant50_dealiased_velocity",
    "quant80dealV": "quant80_dealiased_velocity",
    "quant90dealV": "quant90_dealiased_velocity",
    "quant95dealV": "quant95_dealiased_velocity",
}


def get_rad4alp_prod_fname(datatype):
    """
    Given a datatype find the corresponding start and termination of the
    METRANET product file name

    Parameters
    ----------
    datatype : str
        the data type

    Returns
    -------
    acronym : str
        The start of the METRANET file name
    termination : str
        The end of the METRANET file name

    """
    # datatype missing:
    # NHC (Hail forecast)
    # NZC (Thunderstorm forecast based on TRT)

    termination = ".???"
    # Polar products
    if datatype == "hydro":
        acronym = "YM"
    elif datatype == "dealV":
        acronym = "DV"

    # rainfall accumulation products
    elif datatype == "AZC01":
        acronym = "AZC"  # Rain rate accu with local bias corrected
        termination = ".801"
    elif datatype == "AZC03":
        acronym = "AZC"  # Rain rate accu with local bias corrected
        termination = ".803"
    elif datatype == "AZC06":
        acronym = "AZC"  # Rain rate accu with local bias corrected
        termination = ".806"
    elif datatype == "aZC01":
        acronym = "aZC"  # Rain rate accu with local bias not corrected
        termination = ".801"
    elif datatype == "aZC03":
        acronym = "aZC"  # Rain rate accu with local bias not corrected
        termination = ".803"
    elif datatype == "aZC06":
        acronym = "aZC"  # Rain rate accu with local bias not corrected
        termination = ".806"

    # CPC
    elif datatype == "CPC0005":
        acronym = "CPC"
        termination = "_00005.801.gif"
    elif datatype == "CPC0060":
        acronym = "CPC"
        termination = "_00060.801.gif"
    elif datatype == "CPC0180":
        acronym = "CPC"
        termination = "_00180.801.gif"
    elif datatype == "CPC0360":
        acronym = "CPC"
        termination = "_00360.801.gif"
    elif datatype == "CPC0720":
        acronym = "CPC"
        termination = "_00720.801.gif"
    elif datatype == "CPC1440":
        acronym = "CPC"
        termination = "_01440.801.gif"
    elif datatype == "CPC2880":
        acronym = "CPC"
        termination = "_02880.801.gif"
    elif datatype == "CPC4320":
        acronym = "CPC"
        termination = "_04320.801.gif"

    elif datatype == "CPCH0005":
        acronym = "CPC"
        termination = "_00005.801.gif"
    elif datatype == "CPCH0060":
        acronym = "CPC"
        termination = "_00060.801.gif"
    elif datatype == "CPCH0180":
        acronym = "CPC"
        termination = "_00180.801.gif"
    elif datatype == "CPCH0360":
        acronym = "CPC"
        termination = "_00360.801.gif"
    elif datatype == "CPCH0720":
        acronym = "CPC"
        termination = "_00720.801.gif"
    elif datatype == "CPCH1440":
        acronym = "CPC"
        termination = "_01440.801.gif"
    elif datatype == "CPCH2880":
        acronym = "CPC"
        termination = "_02880.801.gif"
    elif datatype == "CPCH4320":
        acronym = "CPC"
        termination = "_04320.801.gif"

    elif datatype == "nowpal60_P60":
        acronym = "AZC"
        termination = ".accu_0060_RZC_60"
    elif datatype == "nowpal90_P90":
        acronym = "AZC"
        termination = ".accu_0090_RZC_90"
    elif datatype == "nowpal180_P180":
        acronym = "AZC"
        termination = ".accu_0180_RZC_180"
    elif datatype == "nowpal360_P360":
        acronym = "AZC"
        termination = ".accu_0360_RZC_360"
    elif datatype == "nowpal720_P720":
        acronym = "AZC"
        termination = ".accu_0720_RZC_720"

    elif datatype == "nowpal90_P30":
        acronym = "AZC"
        termination = ".accu_0090_RZC_30"
    elif datatype == "nowpal90_P30_F60":
        acronym = "AZC"
        termination = ".accu_0090_RZC_30_INCA_60"
    elif datatype == "nowpal90_F60":
        acronym = "AZC"
        termination = ".accu_0090_INCA_60"
    elif datatype == "nowpal180_P60":
        acronym = "AZC"
        termination = ".accu_0180_RZC_60"
    elif datatype == "nowpal180_P60_F120":
        acronym = "AZC"
        termination = ".accu_0180_RZC_60_INCA_120"
    elif datatype == "nowpal180_F120":
        acronym = "AZC"
        termination = ".accu_0180_INCA_120"
    elif datatype == "nowpal360_P120":
        acronym = "AZC"
        termination = ".accu_0360_CPC_120"
    elif datatype == "nowpal360_P120_F240":
        acronym = "AZC"
        termination = ".accu_0360_CPC_120_INCA_240"
    elif datatype == "nowpal360_F240":
        acronym = "AZC"
        termination = ".accu_0360_INCA_240"
    elif datatype == "nowpal720_P360":
        acronym = "AZC"
        termination = ".accu_0720_CPC_360"
    elif datatype == "nowpal720_P360_F360":
        acronym = "AZC"
        termination = ".accu_0720_CPC_360_INCA_360"
    elif datatype == "nowpal720_F360":
        acronym = "AZC"
        termination = ".accu_0720_INCA_360"

    elif datatype == "dACC":
        acronym = "ACC"  # Daily precip accumulation using NowPal with CPC
        termination = ".1440"
    elif datatype == "dACCH":
        acronym = "ACC"  # reprocessed after 8 days
        termination = ".1440"
    elif datatype == "dARC":
        acronym = "ARC"  # Daily precip accumulation using NowPal with RZC
        termination = ".1440"

    # rainfall rate products
    elif datatype == "RZC":
        acronym = "RZC"  # rain rate local bias corrected
    elif datatype == "R1F":
        acronym = "R1F"  # RZC including best foreign radars
    elif datatype == "rZC":
        acronym = "rZC"  # local bias not corrected
    elif datatype == "RZF":
        acronym = "RZF"  # RZC including foreign radars
    elif datatype == "dRZC":
        acronym = "RZC"  # Daily maximum rain rate

    # hail products
    elif datatype in ("BZC", "dBZC"):
        acronym = "BZC"  # POH
    elif datatype in ("MZC", "dMZC"):
        acronym = "MZC"  # Maximum expected severe hail size
    elif datatype == "GZC":
        acronym = "GZC"  # Hail probability derived from reflectivity
        termination = ".803"
    elif datatype == "dGZC":
        acronym = "GZC"
        termination = ".824"

    # max echo
    elif datatype in ("CZC", "dCZC"):
        acronym = "CZC"  # max echo
    elif datatype == "HZC":
        acronym = "HZC"  # Max echo height

    # echo tops
    elif datatype in ("EZC15", "dEZC15"):
        acronym = "EZC"  # echo top
        termination = ".815"
    elif datatype == "EZC20":
        acronym = "EZC"  # echo top
        termination = ".820"
    elif datatype in ("EZC45", "dEZC45"):
        acronym = "EZC"
        termination = ".845"
    elif datatype == "EZC50":
        acronym = "EZC"
        termination = ".850"

    # Vertically integrated liquid
    elif datatype in ("LZC", "dLZC"):
        acronym = "LZC"

    # Zh CAPPI
    elif datatype in "OZC01":
        acronym = "OZC"
        termination = ".810"
    elif datatype in "OZC02":
        acronym = "OZC"
        termination = ".820"
    elif datatype in "OZC03":
        acronym = "OZC"
        termination = ".830"
    elif datatype in "OZC04":
        acronym = "OZC"
        termination = ".840"
    elif datatype in "OZC05":
        acronym = "OZC"
        termination = ".850"
    elif datatype in "OZC06":
        acronym = "OZC"
        termination = ".860"
    elif datatype in "OZC07":
        acronym = "OZC"
        termination = ".870"
    elif datatype in "OZC08":
        acronym = "OZC"
        termination = ".880"
    elif datatype in "OZC09":
        acronym = "OZC"
        termination = ".890"
    elif datatype in "OZC10":
        acronym = "OZC"
        termination = ".900"
    elif datatype in "OZC11":
        acronym = "OZC"
        termination = ".910"
    elif datatype in "OZC12":
        acronym = "OZC"
        termination = ".920"
    elif datatype in "OZC13":
        acronym = "OZC"
        termination = ".930"
    elif datatype in "OZC14":
        acronym = "OZC"
        termination = ".940"
    elif datatype in "OZC15":
        acronym = "OZC"
        termination = ".950"
    elif datatype in "OZC16":
        acronym = "OZC"
        termination = ".960"
    elif datatype in "OZC17":
        acronym = "OZC"
        termination = ".970"
    elif datatype in "OZC18":
        acronym = "OZC"
        termination = ".980"

    else:
        raise ValueError("ERROR: Unknown rad4alp product data type " + datatype)

    return acronym, termination


def map_hydro(hydro_data_op):
    """
    maps the operational hydrometeor classification identifiers to the ones
    used by Py-ART

    Parameters
    ----------
    hydro_data_op : numpy array
        The operational hydrometeor classification data

    Returns
    -------
    hydro_data_py : numpy array
        The pyart hydrometeor classification data

    """
    hydro_data_py = deepcopy(hydro_data_op)
    hydro_data_py[hydro_data_op == 25] = 3  # crystals
    hydro_data_py[hydro_data_op == 50] = 2  # aggregate
    hydro_data_py[hydro_data_op == 75] = 4  # light rain
    hydro_data_py[hydro_data_op == 100] = 6  # rain
    hydro_data_py[hydro_data_op == 125] = 5  # graupel
    hydro_data_py[hydro_data_op == 150] = 8  # wet snow
    hydro_data_py[hydro_data_op == 175] = 10  # ice hail
    hydro_data_py[hydro_data_op == 200] = 9  # melting hail

    return hydro_data_py


def mf_sname_to_wmo_number(radar_name):
    """
    returns de WMO radar number when given the short radar name used at MF

    Parameters
    ----------
    radar_name : str
        The radar name

    Returns
    -------
    radar_num : str
        The WMO radar number

    """
    if radar_name == "ABBE":  # Abbeville
        radar_num = "07005"
    elif radar_name == "AJAC":  # Ajaccio
        radar_num = "07760"
    elif radar_name == "ALER":  # Aleria
        radar_num = "07774"
    elif radar_name == "TROY":  # Arcis
        radar_num = "07168"
    elif radar_name == "AVES":  # Avesnes
        radar_num = "07083"
    elif radar_name == "BLAI":  # Blaisy
        radar_num = "07274"
    elif radar_name == "BOLL":  # Bollene
        radar_num = "07569"
    elif radar_name == "BORD":  # Bordeaux
        radar_num = "07510"
    elif radar_name == "BOUR":  # Bourges
        radar_num = "07255"
    elif radar_name == "CHER":  # Cherves
        radar_num = "07336"
    elif radar_name == "COLL":  # Collobrieres
        radar_num = "07671"
    elif radar_name == "CAEN":  # Falaise
        radar_num = "07027"
    elif radar_name == "GREZ":  # Grezes
        radar_num = "07436"
    elif radar_name == "MOMU":  # Momuy
        radar_num = "07606"
    elif radar_name == "MTCY":  # Montancy
        radar_num = "07291"
    elif radar_name == "MCLA":  # Montclar
        radar_num = "07637"
    elif radar_name == "NANC":  # Nancy
        radar_num = "07180"
    elif radar_name == "NIME":  # Nimes
        radar_num = "07645"
    elif radar_name == "OPOU":  # Opoul
        radar_num = "07745"
    elif radar_name == "PLAB":  # Plabennec
        radar_num = "07108"
    elif radar_name == "LEPU":  # Sembadel
        radar_num = "07471"
    elif radar_name == "NIZI":  # StNizier
        radar_num = "07381"
    elif radar_name == "TOUL":  # Toulouse
        radar_num = "07629"
    elif radar_name == "TRAP":  # Trappes
        radar_num = "07145"
    elif radar_name == "TREI":  # Treillieres
        radar_num = "07223"
    elif radar_name == "MAUR":  # Maurel
        radar_num = "07572"
    elif radar_name == "COBI":  # Colombis
        radar_num = "07578"
    elif radar_name == "VARS":  # Vars
        radar_num = "07714"
    elif radar_name == "MOUC":  # Moucherotte
        radar_num = "07468"
    elif radar_name == "REMY":  # StRemy
        radar_num = "07366"
    elif radar_name == "NOYA":  # Noyal
        radar_num = "07122"
    elif radar_name == "LEMO":  # Le Moule
        radar_num = "78891"
    elif radar_name == "DIAM":  # Diamant
        radar_num = "78924"
    elif radar_name == "CORA":  # Colorado
        radar_num = "61979"
    elif radar_name == "PVIL":  # Villers
        radar_num = "61978"
    elif radar_name == "NOUM":  # Noumea
        radar_num = "91592"
    elif radar_name == "LIFO":  # Lifou
        radar_num = "91582"
    elif radar_name == "TIEB":  # Tiebaghi
        radar_num = "91571"
    elif radar_name == "JERS":  # Jersey
        radar_num = "03897"
    elif radar_name == "LDOL":  # Ladole
        radar_num = "06699"
    elif radar_name == "VIAL":  # Vial
        radar_num = "07694"
    else:
        warn("Unable to find radar number for radar name " + radar_name)
        radar_num = ""

    return radar_num


def map_Doppler(Doppler_data_bin, Nyquist_vel):
    """
    maps the binary METRANET Doppler data to actual Doppler velocity

    Parameters
    ----------
    Doppler_data_bin : numpy array
        The binary METRANET data

    Returns
    -------
    Doppler_data : numpy array
        The Doppler veloctiy in [m/s]

    """
    Doppler_data = (Doppler_data_bin - 128.0) / 127.0 * Nyquist_vel

    return Doppler_data


def get_save_dir(
    basepath,
    procname,
    dsname,
    prdname,
    timeinfo=None,
    timeformat="%Y-%m-%d",
    create_dir=True,
):
    """
    obtains the path to a product directory and eventually creates it

    Parameters
    ----------
    basepath : str
        product base path
    procname : str
        name of processing space
    dsname : str
        data set name
    prdname : str
        product name
    timeinfo : datetime
        time info to generate the date directory. If None there is no time
        format in the path
    timeformat : str
        Optional. The time format.
    create_dir : boolean
        If True creates the directory

    Returns
    -------
    savedir : str
        path to product

    """
    if timeinfo is None:
        savedir = "{}{}/{}/{}/".format(basepath, procname, dsname, prdname)
    else:
        savedir = "{}{}/{}/{}/{}/".format(
            basepath, procname, timeinfo.strftime(timeformat), dsname, prdname
        )

    if create_dir is False:
        return savedir

    if not os.path.isdir(savedir):
        os.makedirs(savedir)

    return savedir


def make_filename(
    prdtype,
    dstype,
    dsname,
    ext_list,
    prdcfginfo=None,
    timeinfo=None,
    timeformat="%Y%m%d%H%M%S",
    runinfo=None,
):
    """
    creates a product file name

    Parameters
    ----------
    timeinfo : datetime
        time info to generate the date directory
    prdtype : str
        product type, i.e. 'ppi', etc.
    dstype : str
        data set type, i.e. 'raw', etc.
    dsname : str
        data set name
    ext_list : list of str
        file name extensions, i.e. 'png'
    prdcfginfo : str
        Optional. string to add product configuration information, i.e. 'el0.4'
    timeformat : str
        Optional. The time format
    runinfo : str
        Optional. Additional information about the test (e.g. 'RUN01', 'TS011')

    Returns
    -------
    fname_list : list of str
        list of file names (as many as extensions)

    """
    if timeinfo is None:
        timeinfostr = ""
    else:
        timeinfostr = "{}_".format(timeinfo.strftime(timeformat))

    if prdcfginfo is None:
        cfgstr = ""
    else:
        cfgstr = "_{}".format(prdcfginfo)

    if runinfo is None or runinfo == "":
        runstr = ""
    else:
        runstr = "{}_".format(runinfo)

    fname_list = list()
    for ext in ext_list:
        fname = "{}{}{}_{}_{}{}.{}".format(
            timeinfostr, runstr, prdtype, dstype, dsname, cfgstr, ext
        )
        fname = fname.replace("/", "-")
        fname_list.append(fname)

    return fname_list


def generate_field_name_str(datatype):
    """
    Generates a field name in a nice to read format.

    Parameters
    ----------
    datatype : str
        The data type

    Returns
    -------
    field_str : str
        The field name

    """
    field_name = get_fieldname_pyart(datatype)
    field_dic = get_metadata(field_name)
    field_str = field_dic["standard_name"].replace("_", " ")
    field_str = field_str[0].upper() + field_str[1:]
    field_str += " (" + field_dic["units"] + ")"

    return field_str


def get_field_name(datatype):
    """
    Return long name of datatype.

    Parameters
    ----------
    datatype : str
        The data type

    Returns
    -------
    name : str
        The name

    """
    field_name = get_fieldname_pyart(datatype)
    field_dic = get_metadata(field_name)
    name = field_dic["long_name"].replace("_", " ")
    name = name[0].upper() + name[1:]

    return name


def get_field_unit(datatype):
    """
    Return unit of datatype.

    Parameters
    ----------
    datatype : str
        The data type

    Returns
    -------
    unit : str
        The unit

    """
    field_name = get_fieldname_pyart(datatype)
    field_dic = get_metadata(field_name)

    return field_dic["units"]


def get_datatype_metranet(datatype):
    """
    maps de config file radar data type name into the corresponding metranet
    data type  name and Py-ART field name

    Parameters
    ----------
    datatype : str
        config file radar data type name

    Returns
    -------
    metranet type : dict
        dictionary containing the metranet data type name and its
        corresponding Py-ART field name

    """
    if datatype == "dBZ":
        datatype_metranet = "ZH"
        field_name = "reflectivity"
    elif datatype == "dBZv":
        datatype_metranet = "ZV"
        field_name = "reflectivity_vv"
    elif datatype == "ZDR":
        datatype_metranet = "ZDR"
        field_name = "differential_reflectivity"
    elif datatype == "uRhoHV":
        datatype_metranet = "RHO"
        field_name = "uncorrected_cross_correlation_ratio"
    elif datatype == "uPhiDP":
        datatype_metranet = "PHI"
        field_name = "uncorrected_differential_phase"
    elif datatype == "V":
        datatype_metranet = "VEL"
        field_name = "velocity"
    elif datatype == "W":
        datatype_metranet = "WID"
        field_name = "spectrum_width"
    elif datatype == "CLT":
        datatype_metranet = "CLT"
        field_name = "clutter_exit_code"
    elif datatype == "ST1":
        datatype_metranet = "ST1"
        field_name = "stat_test_lag1"
    elif datatype == "ST2":
        datatype_metranet = "ST2"
        field_name = "stat_test_lag2"
    elif datatype == "WBN":
        datatype_metranet = "WBN"
        field_name = "wide_band_noise"
    elif datatype == "MPH":
        datatype_metranet = "MPH"
        field_name = "mean_phase"
    else:
        raise ValueError("ERROR: Metranet fields do not contain datatype " + datatype)

    return {datatype_metranet: field_name}


def get_datatype_knmi(datatype):
    """
    maps de config file radar data type name into the corresponding KNMI
    data type name and Py-ART field name

    Parameters
    ----------
    datatype : str
        config file radar data type name

    Returns
    -------
    knmi type : dict
        dictionary containing the KNMI data type name and its
        corresponding Py-ART field name

    """
    if datatype == "Raccu":
        datatype_knmi = "ACCUMULATION_[MM]"
        field_name = "rainfall_accumulation"
    elif datatype == "QI":
        datatype_knmi = "QUALITY_[-]"
        field_name = "signal_quality_index"
    elif datatype == "AF":
        datatype_knmi = "ADJUSTMENT_FACTOR_[DB]"
        field_name = "adjustment_factor"
    elif datatype == "RR":
        datatype_knmi = "RAINFALL_RATE_[MM/H]"
        field_name = "radar_estimated_rain_rate"
    else:
        raise ValueError(f"ERROR: KNMI fields do not contain datatype {datatype}")

    return {datatype_knmi: field_name}


def get_datatype_skyecho(datatype):
    """
    maps de config file radar data type name into the corresponding SKYECHO
    data type name and Py-ART field name

    Parameters
    ----------
    datatype : str
        config file radar data type name

    Returns
    -------
    skyecho type : dict
        dictionary containing the skyecho data type name and its
        corresponding Py-ART field name

    """
    if datatype == "dBZ":
        datatype_skyecho = "equivalent_reflectivity_factor_HH"
        field_name = "reflectivity"
    elif datatype == "dBZ_flag":
        datatype_skyecho = "equivalent_reflectivity_factor_HH_flag"
        field_name = "reflectivity_flag"
    elif datatype == "dBZv":
        datatype_skyecho = "equivalent_reflectivity_factor_VV"
        field_name = "reflectivity_vv"
    elif datatype == "dBZhv":
        datatype_skyecho = "equivalent_reflectivity_factor_HV"
        field_name = "reflectivity_hv"
    elif datatype == "dBZvh":
        datatype_skyecho = "equivalent_reflectivity_factor_VH"
        field_name = "reflectivity_vh"

    elif datatype == "V":
        datatype_skyecho = "velocity_HH"
        field_name = "velocity"
    elif datatype == "Vv":
        datatype_skyecho = "velocity_VV"
        field_name = "velocity_vv"
    elif datatype == "Vhv":
        datatype_skyecho = "velocity_HV"
        field_name = "velocity_hv"
    elif datatype == "Vvh":
        datatype_skyecho = "velocity_VH"
        field_name = "velocity_vh"

    elif datatype == "W":
        datatype_skyecho = "spectrum_width_HH"
        field_name = "spectrum_width"
    elif datatype == "Wv":
        datatype_skyecho = "spectrum_width_VV"
        field_name = "spectrum_width_vv"
    elif datatype == "Whv":
        datatype_skyecho = "spectrum_width_HV"
        field_name = "spectrum_width_hv"
    elif datatype == "Wvh":
        datatype_skyecho = "spectrum_width_VH"
        field_name = "spectrum_width_vh"

    elif datatype == "Nh":
        datatype_skyecho = "noise_level_HH"
        field_name = "noisedBZ_hh"
    elif datatype == "Nv":
        datatype_skyecho = "noise_level_VV"
        field_name = "noisedBZ_vv"

    elif datatype == "Nclip_h":
        datatype_skyecho = "noise_clipping_level_HH"
        field_name = "noise_clipping_level_hh_dBZ"
    elif datatype == "Nclip_v":
        datatype_skyecho = "noise_clipping_level_VV"
        field_name = "noise_clipping_level_vv_dBZ"

    elif datatype == "PhiDP":
        datatype_skyecho = "phidp"
        field_name = "differential_phase"
    elif datatype == "PhiDPc":
        datatype_skyecho = "phidp_cor"
        field_name = "corrected_differential_phase"
    elif datatype == "KDP":
        datatype_skyecho = "kdp"
        field_name = "specific_differential_phase"
    elif datatype == "KDPc":
        datatype_skyecho = "kdp_cor"
        field_name = "corrected_specific_differential_phase"
    elif datatype == "ZDR":
        datatype_skyecho = "ZDR"
        field_name = "differential_reflectivity"
    elif datatype == "LDRhv":
        datatype_skyecho = "LDR_HV"
        field_name = "linear_depolarization_ratio_hv"
    elif datatype == "LDRvh":
        datatype_skyecho = "LDR_VH"
        field_name = "linear_depolarization_ratio_vh"
    elif datatype == "RhoHV":
        datatype_skyecho = "copolar_correlation"
        field_name = "cross_correlation_ratio"

    elif datatype == "Ah":
        datatype_skyecho = "attenuationRain"
        field_name = "specific_attenuation"

    elif datatype == "RR":
        datatype_skyecho = "rainfall_rate"
        field_name = "radar_estimated_rain_rate"
    elif datatype == "RR_MP":
        datatype_skyecho = "rainfall_rate_mp1948"
        field_name = "Marshall_Palmer_radar_estimated_rain_rate"
    elif datatype == "RR_Z":
        datatype_skyecho = "rainfall_rate_from_z"
        field_name = "radar_reflectivity_estimated_rain_rate"
    elif datatype == "RR_KDP":
        datatype_skyecho = "rainfall_rate_from_kdp"
        field_name = "radar_kdp_estimated_rain_rate"
    elif datatype == "RR_flag":
        datatype_skyecho = "rainfall_rate_flag"
        field_name = "radar_estimated_rain_rate_flag"

    elif datatype == "EDR":
        datatype_skyecho = "EDR13"
        field_name = "turbulence"

    else:
        raise ValueError("ERROR: SKYECHO fields do not contain datatype " + datatype)

    return {datatype_skyecho: field_name}


def get_datatype_odim(datatype):
    """
    maps the config file radar data type name into the corresponding odim
    data type name and Py-ART field name

    Parameters
    ----------
    datatype : str
        config file radar data type name

    Returns
    -------
    metranet type : dict
        dictionary containing the odim data type name and its
        corresponding Py-ART field name

    """
    if datatype == "dBZ":
        field_name = "reflectivity"
        datatype_odim = "DBZH"
    elif datatype == "sigma_zh":
        field_name = "sigma_zh"
        datatype_odim = "DBZH_DEV"
    elif datatype == "dBuZ":
        field_name = "unfiltered_reflectivity"
        datatype_odim = "TH"
    elif datatype == "dBZc":
        field_name = "corrected_reflectivity"
        datatype_odim = "DBZHC"
    elif datatype == "dBuZc":
        field_name = "corrected_unfiltered_reflectivity"
        datatype_odim = "THC"
    elif datatype == "dBZv":
        field_name = "reflectivity_vv"
        datatype_odim = "DBZV"
    elif datatype == "dBZvc":
        field_name = "corrected_reflectivity_vv"
        datatype_odim = "DBZVC"
    elif datatype == "dBuZv":
        field_name = "unfiltered_reflectivity_vv"
        datatype_odim = "TV"
    elif datatype == "dBuZvc":
        field_name = "corrected_unfiltered_reflectivity_vv"
        datatype_odim = "TVC"
    elif datatype == "dBZ_bias":
        field_name = "reflectivity_bias"
        datatype_odim = "ZBIAS"
    elif datatype == "eta_h":
        field_name = "volumetric_reflectivity"
        datatype_odim = "etah"
    elif datatype == "eta_v":
        field_name = "volumetric_reflectivity_vv"
        datatype_odim = "etav"
    elif datatype == "rcs_h":
        field_name = "radar_cross_section_hh"
        datatype_odim = "RCSH"
    elif datatype == "rcs_v":
        field_name = "radar_cross_section_vv"
        datatype_odim = "RCSV"
    elif datatype == "VPRFEATURES":
        field_name = "vpr_features"
        datatype_odim = "VPRFEATURES"

    elif datatype == "ZDR":
        field_name = "differential_reflectivity"
        datatype_odim = "ZDR"
    elif datatype == "ZDRu":
        field_name = "unfiltered_differential_reflectivity"
        datatype_odim = "ZDRU"
    elif datatype == "ZDRc":
        field_name = "corrected_differential_reflectivity"
        datatype_odim = "ZDRC"
    elif datatype == "ZDRuc":
        field_name = "corrected_unfiltered_differential_reflectivity"
        datatype_odim = "ZDRUC"
    elif datatype == "ZDR_prec":
        field_name = "differential_reflectivity_in_precipitation"
        datatype_odim = "ZDRPREC"
    elif datatype == "ZDR_snow":
        field_name = "differential_reflectivity_in_snow"
        datatype_odim = "ZDRSNOW"

    elif datatype == "dBm":
        field_name = "signal_power_hh"
        datatype_odim = "DBMH"
    elif datatype == "dBmv":
        field_name = "signal_power_vv"
        datatype_odim = "DBMV"
    elif datatype == "Nh":
        field_name = "noisedBZ_hh"
        datatype_odim = "NDBZH"
    elif datatype == "Nv":
        field_name = "noisedBZ_vv"
        datatype_odim = "NDBZV"
    elif datatype == "SNR":
        field_name = "signal_to_noise_ratio"
        datatype_odim = "SNRH"
    elif datatype == "SNRh":
        field_name = "signal_to_noise_ratio_hh"
        datatype_odim = "SNRH"
    elif datatype == "SNRv":
        field_name = "signal_to_noise_ratio_vv"
        datatype_odim = "SNRV"
    elif datatype == "SQI":
        field_name = "normalized_coherent_power"
        datatype_odim = "SQIH"
    elif datatype == "SQIv":
        field_name = "normalized_coherent_power_vv"
        datatype_odim = "SQIV"

    elif datatype == "dBm_sun_hit":
        field_name = "sun_hit_power_h"
        datatype_odim = "DBM_SUNHIT"
    elif datatype == "dBmv_sun_hit":
        field_name = "sun_hit_power_v"
        datatype_odim = "DBMV_SUNHIT"
    elif datatype == "ZDR_sun_hit":
        field_name = "sun_hit_differential_reflectivity"
        datatype_odim = "ZDR_SUNHIT"
    elif datatype == "dBm_sun_est":
        field_name = "sun_est_power_h"
        datatype_odim = "DBM_SUNEST"
    elif datatype == "dBmv_sun_est":
        field_name = "sun_est_power_v"
        datatype_odim = "DBMV_SUNEST"
    elif datatype == "ZDR_sun_est":
        field_name = "sun_est_differential_reflectivity"
        datatype_odim = "ZDR_SUNEST"
    elif datatype == "sun_pos_h":
        field_name = "sun_hit_h"
        datatype_odim = "POSH_SUNHIT"
    elif datatype == "sun_pos_v":
        field_name = "sun_hit_v"
        datatype_odim = "POSV_SUNHIT"
    elif datatype == "sun_pos_zdr":
        field_name = "sun_hit_zdr"
        datatype_odim = "POSZDR_SUNHIT"

    elif datatype == "RhoHV":
        field_name = "cross_correlation_ratio"
        datatype_odim = "RHOHV"
    elif datatype == "uRhoHV":
        field_name = "uncorrected_cross_correlation_ratio"
        datatype_odim = "URHOHV"
    elif datatype == "RhoHVc":
        field_name = "corrected_cross_correlation_ratio"
        datatype_odim = "RHOHVC"
    elif datatype == "RhoHV_rain":
        field_name = "cross_correlation_ratio_in_rain"
        datatype_odim = "RHOHVRAIN"
    elif datatype == "L":
        field_name = "logarithmic_cross_correlation_ratio"
        datatype_odim = "LRHOHV"
    elif datatype == "CDR":
        field_name = "circular_depolarization_ratio"
        datatype_odim = "CDR"
    elif datatype == "LDR":
        field_name = "linear_polarization_ratio"
        datatype_odim = "LDR"

    elif datatype == "PhiDP":
        field_name = "differential_phase"
        datatype_odim = "PHIDP"
    elif datatype == "uPhiDP":
        field_name = "uncorrected_differential_phase"
        datatype_odim = "UPHIDP"
    elif datatype == "PhiDPc":
        field_name = "corrected_differential_phase"
        datatype_odim = "PHIDPC"
    elif datatype == "PhiDP0":
        field_name = "system_differential_phase"
        datatype_odim = "PHIDP0"
    elif datatype == "PhiDP0_bin":
        field_name = "first_gate_differential_phase"
        datatype_odim = "PHIDP0_BIN"
    elif datatype == "KDP":
        field_name = "specific_differential_phase"
        datatype_odim = "KDP"
    elif datatype == "KDPc":
        field_name = "corrected_specific_differential_phase"
        datatype_odim = "KDPC"

    elif datatype == "V":
        field_name = "velocity"
        datatype_odim = "VRADH"
    elif datatype == "Vh":
        field_name = "velocity"
        datatype_odim = "VRADH"
    elif datatype == "dealV":
        field_name = "dealiased_velocity"
        datatype_odim = "VRADDH"
    elif datatype == "Vc":
        field_name = "corrected_velocity"
        datatype_odim = "VRADHC"
    elif datatype == "dealVc":
        field_name = "dealiased_corrected_velocity"
        datatype_odim = "VRADDHC"
    elif datatype == "estV":
        field_name = "retrieved_velocity"
        datatype_odim = "VRADEST"
    elif datatype == "stdV":
        field_name = "retrieved_velocity_std"
        datatype_odim = "sd_vvp"
    elif datatype == "diffV":
        field_name = "velocity_difference"
        datatype_odim = "VDIFF"
    elif datatype == "Vv":
        field_name = "velocity_vv"
        datatype_odim = "VRADV"
    elif datatype == "dealVv":
        field_name = "dealiased_velocity_vv"
        datatype_odim = "VRADDV"
    elif datatype == "W":
        field_name = "spectrum_width"
        datatype_odim = "WRADH"
    elif datatype == "Wc":
        field_name = "corrected_spectrum_width"
        datatype_odim = "WRADHC"
    elif datatype == "Wv":
        field_name = "spectrum_width_vv"
        datatype_odim = "WRADV"
    elif datatype == "wind_vel_h_az":
        field_name = "azimuthal_horizontal_wind_component"
        datatype_odim = "AHWND"
    elif datatype == "wind_vel_v":
        field_name = "vertical_wind_component"
        datatype_odim = "w"
    elif datatype == "wind_vel_h_u":
        field_name = "eastward_wind_component"
        datatype_odim = "UWND"
    elif datatype == "wind_vel_h_v":
        field_name = "northward_wind_component"
        datatype_odim = "VWND"
    elif datatype == "windshear_v":
        field_name = "vertical_wind_shear"
        datatype_odim = "VSHR"
    elif datatype == "WIND_SPEED":
        field_name = "wind_speed"
        datatype_odim = "ff"
    elif datatype == "WIND_DIRECTION":
        field_name = "wind_direction"
        datatype_odim = "dd"

    elif datatype == "Ah":
        field_name = "specific_attenuation"
        datatype_odim = "AH"
    elif datatype == "Ahc":
        field_name = "corrected_specific_attenuation"
        datatype_odim = "AHC"
    elif datatype == "PIA":
        field_name = "path_integrated_attenuation"
        datatype_odim = "PIA"
    elif datatype == "PIAc":
        field_name = "corrected_path_integrated_attenuation"
        datatype_odim = "PIAC"
    elif datatype == "Adp":
        field_name = "specific_differential_attenuation"
        datatype_odim = "ADP"
    elif datatype == "Adpc":
        field_name = "corrected_specific_differential_attenuation"
        datatype_odim = "ADPC"
    elif datatype == "PIDA":
        field_name = "path_integrated_differential_attenuation"
        datatype_odim = "PIDA"
    elif datatype == "PIDAc":
        field_name = "corrected_path_integrated_differential_attenuation"
        datatype_odim = "PIDAC"

    elif datatype == "TEMP":
        field_name = "temperature"
        datatype_odim = "TEMP"
    elif datatype == "ISO0":
        field_name = "iso0"
        datatype_odim = "ISO0"
    elif datatype == "H_ISO0":
        field_name = "height_over_iso0"
        datatype_odim = "HISO0"
    elif datatype == "icon_index":
        field_name = "icon_index"
        datatype_odim = "ICONIND"
    elif datatype == "hzt_index":
        field_name = "hzt_index"
        datatype_odim = "HZTIND"
    elif datatype == "ml":
        field_name = "melting_layer"
        datatype_odim = "ML"

    elif datatype == "VIS":
        field_name = "visibility"
        datatype_odim = "VIS"
    elif datatype == "HGHT":
        field_name = "height"
        datatype_odim = "HGHT"

    elif datatype == "echoID":
        field_name = "radar_echo_id"
        datatype_odim = "ECHOID"
    elif datatype == "CLT":
        field_name = "clutter_exit_code"
        datatype_odim = "CLT"
    elif datatype == "occurrence":
        field_name = "occurrence"
        datatype_odim = "OCC"
    elif datatype == "freq_occu":
        field_name = "frequency_of_occurrence"
        datatype_odim = "OCCFREQ"
    elif datatype == "RR":
        field_name = "radar_estimated_rain_rate"
        datatype_odim = "RATE"
    elif datatype == "Raccu":
        field_name = "rainfall_accumulation"
        datatype_odim = "ACRR"
    elif datatype == "RaccuMF":
        field_name = "rainfall_accumulation"
        datatype_odim = "ACRR_hund_mm"
    elif datatype == "QIMF":
        field_name = "signal_quality_index"
        datatype_odim = "QIND2"

    elif datatype == "hydro":
        field_name = "radar_echo_classification"
        datatype_odim = "CLASS"
    elif datatype == "hydroMF":
        field_name = "radar_echo_classification_MF"
        datatype_odim = "CLASS"
    elif datatype == "entropy":
        field_name = "hydroclass_entropy"
        datatype_odim = "ENTROPY"
    elif datatype == "propAG":
        field_name = "proportion_AG"
        datatype_odim = "propAG"
    elif datatype == "propCR":
        field_name = "proportion_CR"
        datatype_odim = "propCR"
    elif datatype == "propLR":
        field_name = "proportion_LR"
        datatype_odim = "propLR"
    elif datatype == "propRP":
        field_name = "proportion_RP"
        datatype_odim = "propRP"
    elif datatype == "propRN":
        field_name = "proportion_RN"
        datatype_odim = "propRN"
    elif datatype == "propVI":
        field_name = "proportion_VI"
        datatype_odim = "propVI"
    elif datatype == "propWS":
        field_name = "proportion_WS"
        datatype_odim = "propWS"
    elif datatype == "propMH":
        field_name = "proportion_MH"
        datatype_odim = "propMH"
    elif datatype == "propIH":
        field_name = "proportion_IH"
        datatype_odim = "propIH"

    elif datatype == "time_avg_flag":
        field_name = "time_avg_flag"
        datatype_odim = "TAFLAG"
    elif datatype == "colocated_gates":
        field_name = "colocated_gates"
        datatype_odim = "COLGATES"
    elif datatype == "nsamples":
        field_name = "number_of_samples"
        datatype_odim = "ns"
    elif datatype == "bird_density":
        field_name = "bird_density"
        datatype_odim = "dens"
    elif datatype == "std":
        field_name = "standard_deviation"
        datatype_odim = "STD"
    elif datatype == "sum":
        field_name = "sum"
        datatype_odim = "SUM"
    elif datatype == "sum2":
        field_name = "sum_squared"
        datatype_odim = "SUM2"

    # vol2bird field names
    elif datatype == "ff":
        field_name = "wind_speed"
        datatype_odim = "ff"
    elif datatype == "dd":
        field_name = "wind_direction"
        datatype_odim = "dd"
    elif datatype == "u":
        field_name = "eastward_wind_component"
        datatype_odim = "UWND"
    elif datatype == "v":
        field_name = "northward_wind_component"
        datatype_odim = "VWND"
    elif datatype == "w":
        field_name = "vertical_wind_component"
        datatype_odim = "w"
    elif datatype == "width":
        field_name = "height_resolution"
        datatype_odim = "width"
    elif datatype == "gap":
        field_name = "gap"
        datatype_odim = "gap"
    elif datatype == "dbz":
        field_name = "bird_reflectivity"
        datatype_odim = "eta"
    elif datatype == "eta":
        field_name = "volumetric_reflectivity"
        datatype_odim = "etah"
    elif datatype == "dens":
        field_name = "bird_density"
        datatype_odim = "dens"
    elif datatype == "n":
        field_name = "number_of_samples_velocity"
        datatype_odim = "n"
    elif datatype == "n_dbz":
        field_name = "number_of_samples_reflectivity"
        datatype_odim = "n_dbz"
    elif datatype == "sd_vvp":
        field_name = "retrieved_velocity_std"
        datatype_odim = "sd_vvp"
    elif datatype == "DBZH":
        field_name = "reflectivity"
        datatype_odim = "DBZH"
    elif datatype == "n_all":
        field_name = "number_of_samples_velocity_all"
        datatype_odim = "n_all"
    elif datatype == "n_dbz_all":
        field_name = "number_of_samples_reflectivity_all"
        datatype_odim = "n_dbz_all"
    elif datatype == "VOL2BIRD_CLASS":
        field_name = "vol2bird_echo_classification"
        datatype_odim = "CELL"
    elif datatype == "VOL2BIRD_WEATHER":
        field_name = "vol2bird_weather"
        datatype_odim = "WEATHER"
    elif datatype == "VOL2BIRD_BACKGROUND":
        field_name = "vol2bird_background"
        datatype_odim = "BACKGROUND"
    elif datatype == "VOL2BIRD_BIOLOGY":
        field_name = "vol2bird_biology"
        datatype_odim = "BIOLOGY"
    else:
        raise ValueError("ERROR: ODIM fields do not contain datatype " + datatype)

    return {datatype_odim: field_name}


# Function to map datatype to Py-ART field name
def get_fieldname_pyart(datatype):
    """
    Maps the config file radar data type name into the corresponding Py-ART field name
    """
    return pyrad_to_pyart_keys_dict.get(datatype)


# Inverted function to map Py-ART field name back to datatype
def get_datatype_from_pyart(pyart_name):
    """
    Maps the Py-ART field name into the corresponding config file radar data type name
    """
    # Reverse the original dictionary to create a new mapping
    field_to_datatype = {v: k for k, v in pyrad_to_pyart_keys_dict.items()}

    return field_to_datatype.get(pyart_name)


def get_fieldname_icon(field_name):
    """
    maps the Py-ART field name into the corresponding ICON variable name

    Parameters
    ----------
    field_name : str
        Py-ART field name

    Returns
    -------
    icon_name : str
        Py-ART variable name

    """
    if field_name == "temperature":
        icon_name = "T"
    elif field_name == "wind_speed":
        icon_name = "FF"
    elif field_name == "wind_direction":
        icon_name = "DD"
    elif field_name == "vertical_wind_shear":
        icon_name = "WSHEAR"
    else:
        raise ValueError("ERROR: Unknown field name " + field_name)

    return icon_name


def get_scan_files_to_merge(
    basepath,
    scan_list,
    radar_name,
    radar_res,
    voltime,
    dataset_list,
    path_convention="ODIM",
    master_scan_time_tol=0,
    scan_period=0,
):
    """
    gets the list of files to merge into a single radar object

    Parameters
    ----------
    basepath : str
        base path of radar data
    scan_list : list
        list of scans
    radar_name : str
        radar name
    radar_res : str
        radar resolution type
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    dataset_list : list
        list of datasets. Used to get path
    path_convention : str
        The path convention to use. Can be LTE, MCH, ODIM, RADARV, RT
    master_scan_time_tol : int
        time tolerance with respect to the master scan time. Can be
        0 (No tolerance), 1 (time between master scan time and master scan
        time plus scan period), -1 (time between master scan time and master
        scan time minus scan period
    scan_period : int
        scan period (minutes)

    Returns
    -------
    filelist : list of strings
        list of files to merge. If it could not find them
        returns an empty list
    scan_list_aux : list of strings
        list of scans for which a file was found
    """
    fname_list = []
    scan_list_aux = []
    dayinfo = voltime.strftime("%y%j")
    timeinfo = voltime.strftime("%H%M")
    for scan in scan_list:
        file_id = "M"
        if radar_name is not None and radar_res is not None:
            basename = f"{file_id}{radar_res}{radar_name}{dayinfo}"
        if path_convention == "LTE":
            yy = dayinfo[0:2]
            dy = dayinfo[2:]
            subf = f"{file_id}{radar_res}{radar_name}{yy}hdf{dy}"
            datapath = f"{basepath}{subf}/"
            pattern = f"{datapath}{basename}{timeinfo}*{scan}*"
            filename = glob.glob(pattern)
            if not filename:
                file_id = "P"
                basename = f"{file_id}{radar_res}{radar_name}{dayinfo}"
                subf = f"{file_id}{radar_res}{radar_name}{yy}hdf{dy}"
                datapath = f"{basepath}{subf}/"
                pattern = f"{datapath}{basename}{timeinfo}*{scan}*"
                filename = glob.glob(pattern)
        elif path_convention == "MCH":
            datapath = f"{basepath}{dayinfo}/{basename}/"
            pattern = f"{datapath}{basename}{timeinfo}*{scan}*"
            filename = glob.glob(pattern)
            if not filename:
                file_id = "P"
                basename = f"{file_id}{radar_res}{radar_name}{dayinfo}"
                datapath = f"{basepath}{dayinfo}/{basename}/"
                pattern = f"{datapath}{basename}{timeinfo}*{scan}*"
                filename = glob.glob(pattern)
        elif path_convention == "ODIM":
            fpath_strf = dataset_list[0][
                dataset_list[0].find("D") + 2 : dataset_list[0].find("F") - 2
            ]
            fdate_strf = dataset_list[0][dataset_list[0].find("F") + 2 : -1]
            datapath = f"{basepath}{voltime.strftime(fpath_strf)}/"
            pattern = f"{datapath}*{scan}*"
            filenames = glob.glob(pattern)
            filename = []
            for filename_aux in filenames:
                fdatetime = find_date_in_file_name(filename_aux, date_format=fdate_strf)
                if master_scan_time_tol == 0:
                    if fdatetime == voltime:
                        filename = [filename_aux]
                        break
                elif master_scan_time_tol == 1:
                    if (
                        voltime
                        <= fdatetime
                        < voltime + datetime.timedelta(minutes=scan_period)
                    ):
                        filename = [filename_aux]
                        break
                else:
                    if (
                        voltime - datetime.timedelta(minutes=scan_period)
                        < fdatetime
                        <= voltime
                    ):
                        filename = [filename_aux]
                        break
        elif path_convention == "SkyEcho":
            fpath_strf = dataset_list[0][
                dataset_list[0].find("D") + 2 : dataset_list[0].find("F") - 2
            ]
            fdate_strf = dataset_list[0][dataset_list[0].find("F") + 2 : -1]
            datapath = f"{basepath}{voltime.strftime(fpath_strf)}/"
            pattern = f"{datapath}*{scan}*"
            filenames = glob.glob(pattern)
            filename = []
            time_ref = datetime.datetime.strptime(
                voltime.strftime("%Y%m%d000000"), "%Y%m%d%H%M%S"
            )
            for filename_aux in filenames:
                fdatetime = find_date_in_file_name(filename_aux, date_format=fdate_strf)
                if fdatetime >= time_ref:
                    filename = [filename_aux]
        elif path_convention == "RADARV":
            fpath_strf = dataset_list[0][
                dataset_list[0].find("D") + 2 : dataset_list[0].find("F") - 2
            ]
            fdate_strf = dataset_list[0][dataset_list[0].find("F") + 2 : -1]
            pattern = f"{basepath}{scan}{voltime.strftime(fpath_strf)}/*"
            filenames = glob.glob(pattern)
            filename = []
            for filename_aux in filenames:
                fdatetime = find_date_in_file_name(filename_aux, date_format=fdate_strf)
                if master_scan_time_tol == 0:
                    if fdatetime == voltime:
                        filename = [filename_aux]
                        break
                elif master_scan_time_tol == 1:
                    if (
                        voltime
                        <= fdatetime
                        < voltime + datetime.timedelta(minutes=scan_period)
                    ):
                        filename = [filename_aux]
                        break
                else:
                    if (
                        voltime - datetime.timedelta(minutes=scan_period)
                        < fdatetime
                        <= voltime
                    ):
                        filename = [filename_aux]
                        break
        else:
            datapath = f"{basepath}{file_id}{radar_res}{radar_name}/"
            pattern = f"{datapath}{basename}{timeinfo}*{scan}*"
            filename = glob.glob(pattern)
            if not filename:
                file_id = "P"
                basename = f"{file_id}{radar_res}{radar_name}{dayinfo}"
                datapath = f"{basepath}{file_id}{radar_res}{radar_name}/"
                pattern = f"{datapath}{basename}{timeinfo}*{scan}*"
                filename = glob.glob(pattern)

        if not filename:
            warn(f"No file found in {pattern}")
            continue
        fname_list.append(filename[0])
        scan_list_aux.append(scan)

    return fname_list, scan_list_aux


def get_scan_files_to_merge_s3(
    basepath,
    scan_list,
    radar_name,
    radar_res,
    voltime,
    dataset_list,
    cfg,
    path_convention="ODIM",
    master_scan_time_tol=0,
    scan_period=0,
):
    """
    gets the list of files to merge into a single radar object

    Parameters
    ----------
    basepath : str
        base path of radar data
    scan_list : list
        list of scans
    radar_name : str
        radar name
    radar_res : str
        radar resolution type
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    dataset_list : list
        list of datasets. Used to get path
    cfg : dict
        configuration dictionary
    path_convention : str
        The path convention to use. Can be LTE, MCH, ODIM, RADARV, RT
    master_scan_time_tol : int
        time tolerance with respect to the master scan time. Can be
        0 (No tolerance), 1 (time between master scan time and master scan
        time plus scan period), -1 (time between master scan time and master
        scan time minus scan period
    scan_period : int
        scan period (minutes)

    Returns
    -------
    filelist : list of strings
        list of files to merge. If it could not find them
        returns an empty list
    scan_list_aux : list of strings
        list of scans for which a file was found
    """
    fname_list = []
    scan_list_aux = []
    if not _BOTO3_AVAILABLE:
        warn("boto3 not installed")
        return fname_list

    s3_client = boto3.client(
        "s3",
        endpoint_url=cfg["s3EndpointRead"],
        aws_access_key_id=cfg["s3KeyRead"],
        aws_secret_access_key=cfg["s3SecretRead"],
        verify=False,
    )
    response = s3_client.list_objects_v2(Bucket=cfg["s3BucketRead"])

    dayinfo = voltime.strftime("%y%j")
    timeinfo = voltime.strftime("%H%M")
    for scan in scan_list:
        file_id = "M"
        if radar_name is not None and radar_res is not None:
            basename = f"{file_id}{radar_res}{radar_name}{dayinfo}"
        if path_convention == "LTE":
            yy = dayinfo[0:2]
            dy = dayinfo[2:]
            subf = f"{file_id}{radar_res}{radar_name}{yy}hdf{dy}"
            datapath = f"{basepath}{subf}/"
            pattern = f"{datapath}{basename}{timeinfo}*{scan}*"
            found = False
            for content in response["Contents"]:
                if fnmatch.fnmatch(content["Key"], pattern):
                    fname_list.append(content["Key"])
                    scan_list_aux.append(scan)
                    found = True
            if not found:
                file_id = "P"
                basename = f"{file_id}{radar_res}{radar_name}{dayinfo}"
                subf = f"{file_id}{radar_res}{radar_name}{yy}hdf{dy}"
                datapath = f"{basepath}{subf}/"
                pattern = f"{datapath}{basename}{timeinfo}*{scan}*"
                for content in response["Contents"]:
                    if fnmatch.fnmatch(content["Key"], pattern):
                        fname_list.append(content["Key"])
                        scan_list_aux.append(scan)
                        found = True
            if not found:
                warn(f"No file found in {pattern}")
        elif path_convention == "MCH":
            datapath = f"{basepath}{dayinfo}/{basename}/"
            pattern = f"{datapath}{basename}{timeinfo}*{scan}*"
            found = False
            for content in response["Contents"]:
                if fnmatch.fnmatch(content["Key"], pattern):
                    fname_list.append(content["Key"])
                    scan_list_aux.append(scan)
                    found = True
            if not found:
                file_id = "P"
                basename = f"{file_id}{radar_res}{radar_name}{dayinfo}"
                datapath = f"{basepath}{dayinfo}/{basename}/"
                pattern = f"{datapath}{basename}{timeinfo}*{scan}*"
                for content in response["Contents"]:
                    if fnmatch.fnmatch(content["Key"], pattern):
                        fname_list.append(content["Key"])
                        scan_list_aux.append(scan)
                        found = True
            if not found:
                warn(f"No file found in {pattern}")
        elif path_convention == "ODIM":
            fpath_strf = dataset_list[0][
                dataset_list[0].find("D") + 2 : dataset_list[0].find("F") - 2
            ]
            fdate_strf = dataset_list[0][dataset_list[0].find("F") + 2 : -1]
            daydir = voltime.strftime(fpath_strf)
            if daydir == "":
                pattern = f"{basepath}*{scan}*"
            else:
                pattern = f"{basepath}{daydir}/*{scan}*"

            found = False
            for content in response["Contents"]:
                if fnmatch.fnmatch(content["Key"], pattern):
                    fdatetime = find_date_in_file_name(
                        content["Key"], date_format=fdate_strf
                    )
                    if master_scan_time_tol == 0:
                        if fdatetime == voltime:
                            fname_list.append(content["Key"])
                            scan_list_aux.append(scan)
                            found = True
                            break
                    elif master_scan_time_tol == 1:
                        if (
                            voltime
                            <= fdatetime
                            < voltime + datetime.timedelta(minutes=scan_period)
                        ):
                            fname_list.append(content["Key"])
                            scan_list_aux.append(scan)
                            found = True
                            break
                    else:
                        if (
                            voltime - datetime.timedelta(minutes=scan_period)
                            < fdatetime
                            <= voltime
                        ):
                            fname_list.append(content["Key"])
                            scan_list_aux.append(scan)
                            found = True
                            break
            if not found:
                warn(f"No file found in {pattern}")
        elif path_convention == "SkyEcho":
            fpath_strf = dataset_list[0][
                dataset_list[0].find("D") + 2 : dataset_list[0].find("F") - 2
            ]
            fdate_strf = dataset_list[0][dataset_list[0].find("F") + 2 : -1]
            daydir = voltime.strftime(fpath_strf)
            if daydir == "":
                pattern = f"{basepath}*{scan}*"
            else:
                pattern = f"{basepath}{daydir}/*{scan}*"
            time_ref = datetime.datetime.strptime(
                voltime.strftime("%Y%m%d000000"), "%Y%m%d%H%M%S"
            )

            found = False
            fname_list_aux = []
            for content in response["Contents"]:
                if fnmatch.fnmatch(content["Key"], pattern):
                    fdatetime = find_date_in_file_name(
                        content["Key"], date_format=fdate_strf
                    )
                    if fdatetime >= time_ref:
                        fname_list_aux.append(content["Key"])
                        found = True
                        break
            if not found:
                warn(f"No file found in {pattern}")
            else:
                fname_list.append(fname_list_aux[-1])
                scan_list_aux.append(scan)
        elif path_convention == "RADARV":
            fpath_strf = dataset_list[0][
                dataset_list[0].find("D") + 2 : dataset_list[0].find("F") - 2
            ]
            fdate_strf = dataset_list[0][dataset_list[0].find("F") + 2 : -1]
            daydir = voltime.strftime(fpath_strf)
            if daydir == "":
                pattern = f"{basepath}{scan}/*"
            else:
                pattern = f"{basepath}{scan}/{daydir}/*"

            found = False
            for content in response["Contents"]:
                if fnmatch.fnmatch(content["Key"], pattern):
                    fdatetime = find_date_in_file_name(
                        content["Key"], date_format=fdate_strf
                    )
                    if master_scan_time_tol == 0:
                        if fdatetime == voltime:
                            fname_list.append(content["Key"])
                            scan_list_aux.append(scan)
                            found = True
                            break
                    elif master_scan_time_tol == 1:
                        if (
                            voltime
                            <= fdatetime
                            < voltime + datetime.timedelta(minutes=scan_period)
                        ):
                            fname_list.append(content["Key"])
                            scan_list_aux.append(scan)
                            found = True
                            break
                    else:
                        if (
                            voltime - datetime.timedelta(minutes=scan_period)
                            < fdatetime
                            <= voltime
                        ):
                            fname_list.append(content["Key"])
                            scan_list_aux.append(scan)
                            found = True
                            break
            if not found:
                warn(f"No file found in {pattern}")
        else:
            datapath = f"{basepath}{file_id}{radar_res}{radar_name}/"
            pattern = f"{datapath}{basename}{timeinfo}*{scan}*"
            found = False
            for content in response["Contents"]:
                if fnmatch.fnmatch(content["Key"], pattern):
                    fname_list.append(content["Key"])
                    scan_list_aux.append(scan)
                    found = True
            if not found:
                file_id = "P"
                basename = f"{file_id}{radar_res}{radar_name}{dayinfo}"
                datapath = f"{basepath}{file_id}{radar_res}{radar_name}/"
                pattern = f"{datapath}{basename}{timeinfo}*{scan}*"
                for content in response["Contents"]:
                    if fnmatch.fnmatch(content["Key"], pattern):
                        fname_list.append(content["Key"])
                        scan_list_aux.append(scan)
                        found = True
            if not found:
                warn(f"No file found in {pattern}")

    return fname_list, scan_list_aux


def get_file_list(datadescriptor, starttimes, endtimes, cfg, scan=None):
    """
    gets the list of files within a time period

    Parameters
    ----------
    datadescriptor : str
        radar field type. Format : [radar file type]:[datatype]
    startimes : array of datetime objects
        start of time periods
    endtimes : array of datetime object
        end of time periods
    cfg: dictionary of dictionaries
        configuration info to figure out where the data is
    scan : str
        scan name

    Returns
    -------
    filelist : list of strings
        list of files within the time period. If it could not find them
        returns an empty list
    """
    radarnr, datagroup, datatype, dataset, product = get_datatype_fields(datadescriptor)

    ind_rad = int(radarnr[5:8]) - 1
    if datatype in ("Nh", "Nv"):
        datatype = "dBZ"

    filelist = []
    for starttime, endtime in zip(starttimes, endtimes):
        startdate = starttime.replace(hour=0, minute=0, second=0, microsecond=0)
        enddate = endtime.replace(hour=0, minute=0, second=0, microsecond=0)
        ndays = int((enddate - startdate).days) + 1
        t_filelist = []
        pattern = None
        for i in range(ndays):
            if datagroup == "RAINBOW":
                if scan is None:
                    warn("Unknown scan name")
                    return []
                daydir = (starttime + datetime.timedelta(days=i)).strftime("%Y-%m-%d")
                dayinfo = (starttime + datetime.timedelta(days=i)).strftime("%Y%m%d")
                datapath = os.path.join(cfg["datapath"][ind_rad], scan, daydir)
                if not os.path.isdir(datapath):
                    warn(f"WARNING: Unknown datapath '{datapath}'")
                    continue
                pattern = os.path.join(datapath, f"{dayinfo}*00{datatype}.*")
                dayfilelist = glob.glob(pattern)
                for filename in dayfilelist:
                    t_filelist.append(filename)
            elif datagroup == "RAD4ALP":
                if scan is None:
                    warn("Unknown scan name")
                    return []

                datapath, basename = get_rad4alp_dir(
                    cfg["datapath"][ind_rad],
                    starttime + datetime.timedelta(days=i),
                    radar_name=cfg["RadarName"][ind_rad],
                    radar_res=cfg["RadarRes"][ind_rad],
                    scan=scan,
                    path_convention=cfg["path_convention"][ind_rad],
                )

                if not os.path.isdir(datapath):
                    warn(f"WARNING: Unknown datapath '{datapath}'")
                    continue
                pattern = os.path.join(datapath, f"{basename}*.{scan}*")
                dayfilelist = glob.glob(pattern)
                for filename in dayfilelist:
                    t_filelist.append(filename)
            elif datagroup in ("RAD4ALPGRID", "RAD4ALPGIF", "RAD4ALPBIN"):
                acronym, termination = get_rad4alp_prod_fname(datatype)
                dir_day = starttime + datetime.timedelta(days=i)
                dayinfo = dir_day.strftime("%y%j")
                basename = f"{acronym}{dayinfo}"

                datapath = get_rad4alp_grid_dir(
                    cfg["datapath"][ind_rad],
                    dir_day,
                    datatype,
                    acronym,
                    path_convention=cfg["path_convention"][ind_rad],
                )

                if not os.path.isdir(datapath):
                    warn(f"WARNING: Unknown datapath '{datapath}'")
                    continue
                pattern = os.path.join(datapath, f"{basename}*{termination}")
                dayfilelist = glob.glob(pattern)
                for filename in dayfilelist:
                    t_filelist.append(filename)
            elif datagroup == "SATGRID":
                daydir = os.path.join(
                    (starttime + datetime.timedelta(days=i)).strftime("%Y/%m/%d")
                )
                dayinfo = (starttime + datetime.timedelta(days=i)).strftime("%Y%m%d")
                datapath = os.path.join(cfg["satpath"][ind_rad], daydir)
                if not os.path.isdir(datapath):
                    continue
                pattern = os.path.join(datapath, f"MSG?_ccs4_{dayinfo}*_rad_PLAX.nc")
                dayfilelist = glob.glob(pattern)
                for filename in dayfilelist:
                    t_filelist.append(filename)
            elif datagroup == "SKYECHO":
                try:
                    fpath_strf = dataset[dataset.find("D") + 2 : dataset.find("F") - 2]
                except AttributeError:
                    warn(
                        "Unknown directory and/or date "
                        + "convention, check product config file"
                    )
                daydir = (starttime + datetime.timedelta(days=i)).strftime(fpath_strf)
                datapath = os.path.join(cfg["datapath"][ind_rad], daydir)
                pattern = os.path.join(f"{datapath}*", f"*{scan}*")
                dayfilelist = glob.glob(pattern)
                for filename in dayfilelist:
                    t_filelist.append(filename)
            elif datagroup in (
                "ODIM",
                "ODIMBIRDS",
                "CFRADIAL",
                "CFRADIAL2",
                "CF1",
                "NEXRADII",
                "GAMIC",
                "ODIMGRID",
                "KNMIH5GRID",
            ):
                if scan is None:
                    warn("Unknown scan name")
                    return []

                if cfg["path_convention"][ind_rad] == "MCH":
                    dayinfo = (starttime + datetime.timedelta(days=i)).strftime("%y%j")
                    basename = (
                        f'M{cfg["RadarRes"][ind_rad]}'
                        f'{cfg["RadarName"][ind_rad]}{dayinfo}'
                    )
                    datapath = os.path.join(cfg["datapath"][ind_rad], dayinfo, basename)

                    # check that M files exist. if not search P files
                    pattern = os.path.join(datapath, f"{basename}*{scan}*")
                    dayfilelist = glob.glob(pattern)
                    if not dayfilelist:
                        basename = (
                            f'P{cfg["RadarRes"][ind_rad]}'
                            f'{cfg["RadarName"][ind_rad]}{dayinfo}'
                        )
                        datapath = os.path.join(
                            cfg["datapath"][ind_rad], dayinfo, basename
                        )
                        dayfilelist = glob.glob(pattern)

                    if not os.path.isdir(datapath):
                        warn("WARNING: Unknown datapath '%s'" % datapath)
                        continue

                elif cfg["path_convention"][ind_rad] == "ODIM":
                    try:
                        fpath_strf = dataset[
                            dataset.find("D") + 2 : dataset.find("F") - 2
                        ]
                    except AttributeError:
                        warn(
                            "Unknown ODIM directory and/or date "
                            + "convention, check product config file"
                        )
                    daydir = (starttime + datetime.timedelta(days=i)).strftime(
                        fpath_strf
                    )
                    datapath = os.path.join(cfg["datapath"][ind_rad], daydir)
                    pattern = os.path.join(datapath, f"*{scan}*")
                    dayfilelist = glob.glob(pattern)

                elif cfg["path_convention"][ind_rad] == "RADARV":
                    try:
                        fpath_strf = dataset[
                            dataset.find("D") + 2 : dataset.find("F") - 2
                        ]
                    except AttributeError:
                        warn(
                            "Unknown ODIM directory and/or date "
                            + "convention, check product config file"
                        )
                    daydir = (starttime + datetime.timedelta(days=i)).strftime(
                        fpath_strf
                    )
                    datapath = os.path.join(cfg["datapath"][ind_rad], scan)
                    pattern = os.path.join(datapath, daydir, "*")
                    dayfilelist = glob.glob(pattern)

                else:
                    dayinfo = (starttime + datetime.timedelta(days=i)).strftime("%y%j")
                    basename = (
                        "M"
                        + cfg["RadarRes"][ind_rad]
                        + cfg["RadarName"][ind_rad]
                        + dayinfo
                    )
                    datapath = os.path.join(
                        cfg["datapath"][ind_rad],
                        "M" + cfg["RadarRes"][ind_rad] + cfg["RadarName"][ind_rad],
                    )

                    # check that M files exist. if not search P files
                    pattern = os.path.join(datapath, f"{basename}*{scan}*")
                    dayfilelist = glob.glob(pattern)
                    if not dayfilelist:
                        basename = (
                            "P"
                            + cfg["RadarRes"][ind_rad]
                            + cfg["RadarName"][ind_rad]
                            + dayinfo
                        )
                        datapath = os.path.join(
                            cfg["datapath"][ind_rad],
                            "P" + cfg["RadarRes"][ind_rad] + cfg["RadarName"][ind_rad],
                        )
                        dayfilelist = glob.glob(pattern)

                    if not os.path.isdir(datapath):
                        warn("WARNING: Unknown datapath '%s'" % datapath)
                        continue

                for filename in dayfilelist:
                    t_filelist.append(filename)

            elif datagroup in (
                "CFRADIALPYRAD",
                "ODIMPYRAD",
                "PYRADGRID",
                "ODIMPYRADGRID",
                "NETCDFSPECTRA",
                "CSV",
            ):
                termination = ".nc"
                if datagroup in ("ODIMPYRAD", "ODIMPYRADGRID"):
                    termination = ".h*"
                elif datagroup == "CSV":
                    termination = ".csv"

                daydir = (starttime + datetime.timedelta(days=i)).strftime("%Y-%m-%d")
                dayinfo = (starttime + datetime.timedelta(days=i)).strftime("%Y%m%d")
                datapath = os.path.join(
                    cfg["loadbasepath"][ind_rad],
                    cfg["loadname"][ind_rad],
                    daydir,
                    dataset,
                    product,
                )
                if not os.path.isdir(datapath):
                    warn("WARNING: Unknown datapath '%s'" % datapath)
                    continue
                pattern = os.path.join(datapath, f"{dayinfo}*{datatype}{termination}")
                dayfilelist = glob.glob(pattern)
                for filename in dayfilelist:
                    t_filelist.append(filename)

            elif datagroup in ("GECSX"):
                termination = ".nc"
                datapath = os.path.join(
                    cfg["gecsxbasepath"][ind_rad],
                    cfg["gecsxname"][ind_rad],
                    dataset,
                    product,
                )
                if not os.path.isdir(datapath):
                    warn("WARNING: Unknown datapath '%s'" % datapath)
                    continue
                pattern = os.path.join(datapath, "*" + datatype + termination)
                dayfilelist = glob.glob(pattern)
                for filename in dayfilelist:
                    t_filelist.append(filename)
            elif datagroup in (
                "MFCFRADIAL",
                "MFBIN",
                "MFPNG",
                "MFGRIB",
                "MFDAT",
                "MFCF",
            ):
                try:
                    fpath_strf = dataset[dataset.find("D") + 2 : dataset.find("F") - 2]
                except AttributeError:
                    warn(
                        "Unknown directory and/or date "
                        + "convention, check product config file"
                    )
                    continue
                daydir = (starttime + datetime.timedelta(days=i)).strftime(fpath_strf)
                datapath = os.path.join(cfg["datapath"][ind_rad], daydir)
                pattern = os.path.join(datapath, "*" + scan + "*")
                dayfilelist = glob.glob(pattern)

                for filename in dayfilelist:
                    t_filelist.append(filename)
            elif datagroup == "MXPOL":
                if scan is None:
                    warn("Unknown scan name")
                    return []
                if cfg["path_convention"][ind_rad] == "LTE":
                    sub1 = str(starttime.year)
                    sub2 = starttime.strftime("%m")
                    sub3 = starttime.strftime("%d")
                    datapath = os.path.join(cfg["datapath"][ind_rad], sub1, sub2, sub3)
                    basename = (
                        "MXPol-polar-"
                        + starttime.strftime("%Y%m%d")
                        + "-*-"
                        + scan
                        + "*"
                    )
                    pattern = os.path.join(datapath, basename)
                    dayfilelist = glob.glob(pattern)
                else:
                    daydir = (starttime + datetime.timedelta(days=i)).strftime(
                        "%Y-%m-%d"
                    )
                    dayinfo = (starttime + datetime.timedelta(days=i)).strftime(
                        "%Y%m%d"
                    )
                    datapath = os.path.join(cfg["datapath"][ind_rad], scan, daydir)
                    if not os.path.isdir(datapath):
                        warn("WARNING: Unknown datapath '%s'" % datapath)
                        continue
                    pattern = os.path.join(
                        datapath, f"MXPol-polar-{dayinfo}-*-{scan}.nc"
                    )
                    dayfilelist = glob.glob(pattern)
                for filename in dayfilelist:
                    t_filelist.append(filename)

            elif datagroup == "ICONRAW":
                daydir = (starttime + datetime.timedelta(days=i)).strftime("%Y-%m-%d")
                dayinfo = (starttime + datetime.timedelta(days=i)).strftime("%Y%m%d")

                # check that base directory exists
                datapath = os.path.join(cfg["iconpath"][ind_rad], datatype, "raw")
                if not os.path.isdir(datapath):
                    datapath = os.path.join(cfg["iconpath"][ind_rad], datatype, "raw1")
                    if not os.path.isdir(datapath):
                        warn("WARNING: Unknown datapath '%s'" % datapath)
                        continue
                if cfg["path_convention"][ind_rad] == "MCH":
                    datapath = os.path.join(datapath, daydir)

                if not os.path.isdir(datapath):
                    warn("WARNING: Unknown datapath '%s'" % datapath)
                    continue
                pattern = os.path.join(datapath, f"*{dayinfo}*.nc")
                dayfilelist = glob.glob(pattern)
                for filename in dayfilelist:
                    t_filelist.append(filename)

        if datagroup == "GECSX":
            # For GECSX we ignore time, since the visibility is static
            filelist = t_filelist
        elif datagroup == "SKYECHO":
            # Each file contains multiple scans
            for filename in t_filelist:
                _, tend_sweeps, _, _ = get_sweep_time_coverage(filename)
                for tend in tend_sweeps:
                    if starttime <= tend <= endtime:
                        filelist.append(
                            f"{str(filename)}::{tend.strftime('%Y-%m-%dT%H:%M:%S.%f')}"
                        )
        else:
            for filename in t_filelist:
                filenamestr = str(filename)
                fdatetime = get_datetime(filenamestr, datadescriptor)
                if fdatetime is not None:
                    if starttime <= fdatetime <= endtime:
                        if filenamestr not in filelist:
                            filelist.append(filenamestr)

        if not filelist:
            if pattern is not None:
                warn(
                    "WARNING: No file with pattern {:s} could be found between ".format(
                        pattern
                    )
                    + "starttime {:s} and endtime {:s}".format(
                        str(starttime), str(endtime)
                    )
                )
    return sorted(filelist)


def get_file_list_s3(datadescriptor, starttimes, endtimes, cfg, scan=None):
    """
    gets the list of files within a time period from an s3 bucket

    Parameters
    ----------
    datadescriptor : str
        radar field type. Format : [radar file type]:[datatype]
    startimes : array of datetime objects
        start of time periods
    endtimes : array of datetime object
        end of time periods
    cfg: dictionary of dictionaries
        configuration info to figure out where the data is
    scan : str
        scan name

    Returns
    -------
    filelist : list of strings
        list of files within the time period. If it could not find them
        returns an empty list
    """
    filelist = []
    if not _BOTO3_AVAILABLE:
        warn("boto3 not installed")
        return filelist
    radarnr, datagroup, datatype, dataset, product = get_datatype_fields(datadescriptor)

    ind_rad = int(radarnr[5:8]) - 1
    if datatype in ("Nh", "Nv"):
        datatype = "dBZ"

    s3_client = boto3.client(
        "s3",
        endpoint_url=cfg["s3EndpointRead"],
        aws_access_key_id=cfg["s3KeyRead"],
        aws_secret_access_key=cfg["s3SecretRead"],
        verify=False,
    )
    response = s3_client.list_objects_v2(Bucket=cfg["s3BucketRead"])

    for starttime, endtime in zip(starttimes, endtimes):
        startdate = starttime.replace(hour=0, minute=0, second=0, microsecond=0)
        enddate = endtime.replace(hour=0, minute=0, second=0, microsecond=0)
        ndays = int((enddate - startdate).days) + 1
        t_filelist = []
        pattern = None
        for i in range(ndays):
            if datagroup == "SKYECHO":
                try:
                    fpath_strf = dataset[dataset.find("D") + 2 : dataset.find("F") - 2]
                except AttributeError:
                    warn(
                        "Unknown directory and/or date "
                        + "convention, check product config file"
                    )
                daydir = (starttime + datetime.timedelta(days=i)).strftime(fpath_strf)
                if daydir == "":
                    pattern = f'{cfg["s3PathRead"]}*{scan}*'
                else:
                    pattern = f'{cfg["s3PathRead"]}{daydir}/*{scan}*'
                dayfilelist = []
                for content in response["Contents"]:
                    if fnmatch.fnmatch(content["Key"], pattern):
                        dayfilelist.append(content["Key"])

                for filename in dayfilelist:
                    t_filelist.append(filename)
            elif datagroup in (
                "ODIM",
                "ODIMBIRDS",
                "CFRADIAL",
                "CFRADIAL2",
                "CF1",
                "NEXRADII",
                "GAMIC",
                "ODIMGRID",
                "KNMIH5GRID",
            ):
                if scan is None:
                    warn("Unknown scan name")
                    return []
                if cfg["path_convention"][ind_rad] == "MCH":
                    dayinfo = (starttime + datetime.timedelta(days=i)).strftime("%y%j")
                    basename = (
                        f'M{cfg["RadarRes"][ind_rad]}'
                        f'{cfg["RadarName"][ind_rad]}{dayinfo}'
                    )
                    datapath = f'{cfg["s3PathRead"]}{dayinfo}/{basename}/'

                    # check that M files exist. if not search P files
                    pattern = f"{datapath}{basename}*{scan}*"
                    dayfilelist = []
                    for content in response["Contents"]:
                        if fnmatch.fnmatch(content["Key"], pattern):
                            dayfilelist.append(content["Key"])
                    if not dayfilelist:
                        basename = (
                            f'P{cfg["RadarRes"][ind_rad]}'
                            f'{cfg["RadarName"][ind_rad]}{dayinfo}'
                        )
                        datapath = f'{cfg["s3PathRead"]}{dayinfo}/{basename}/'
                        pattern = f"{datapath}{basename}*{scan}*"
                        for content in response["Contents"]:
                            if fnmatch.fnmatch(content["Key"], pattern):
                                dayfilelist.append(content["Key"])
                elif cfg["path_convention"][ind_rad] == "ODIM":
                    try:
                        fpath_strf = dataset[
                            dataset.find("D") + 2 : dataset.find("F") - 2
                        ]
                    except AttributeError:
                        warn(
                            "Unknown ODIM directory and/or date "
                            + "convention, check product config file"
                        )
                    daydir = (starttime + datetime.timedelta(days=i)).strftime(
                        fpath_strf
                    )
                    if daydir == "":
                        pattern = f'{cfg["s3PathRead"]}*{scan}*'
                    else:
                        pattern = f'{cfg["s3PathRead"]}{daydir}/*{scan}*'

                    dayfilelist = []
                    for content in response["Contents"]:
                        if fnmatch.fnmatch(content["Key"], pattern):
                            dayfilelist.append(content["Key"])
                elif cfg["path_convention"][ind_rad] == "RADARV":
                    try:
                        fpath_strf = dataset[
                            dataset.find("D") + 2 : dataset.find("F") - 2
                        ]
                    except AttributeError:
                        warn(
                            "Unknown ODIM directory and/or date "
                            + "convention, check product config file"
                        )
                    daydir = (starttime + datetime.timedelta(days=i)).strftime(
                        fpath_strf
                    )
                    if daydir == "":
                        pattern = f'{cfg["s3PathRead"]}{scan}/*'
                    else:
                        pattern = f'{cfg["s3PathRead"]}{scan}/{daydir}/*'

                    dayfilelist = []
                    for content in response["Contents"]:
                        if fnmatch.fnmatch(content["Key"], pattern):
                            dayfilelist.append(content["Key"])
                else:
                    dayinfo = (starttime + datetime.timedelta(days=i)).strftime("%y%j")
                    basename = (
                        f'M{cfg["RadarRes"][ind_rad]}'
                        f'{cfg["RadarName"][ind_rad]}{dayinfo}'
                    )
                    datapath = (
                        f'{cfg["s3PathRead"]}M{cfg["RadarRes"][ind_rad]}'
                        f'{cfg["RadarName"][ind_rad]}'
                    )

                    # check that M files exist. if not search P files
                    pattern = f"{datapath}{basename}*{scan}*"
                    dayfilelist = []
                    for content in response["Contents"]:
                        if fnmatch.fnmatch(content["Key"], pattern):
                            dayfilelist.append(content["Key"])
                    if not dayfilelist:
                        basename = (
                            f'P{cfg["RadarRes"][ind_rad]}'
                            f'{cfg["RadarName"][ind_rad]}{dayinfo}'
                        )
                        datapath = (
                            f'{cfg["s3PathRead"]}P{cfg["RadarRes"][ind_rad]}'
                            f'{cfg["RadarName"][ind_rad]}'
                        )
                        pattern = f"{datapath}{basename}*{scan}*"
                        for content in response["Contents"]:
                            if fnmatch.fnmatch(content["Key"], pattern):
                                dayfilelist.append(content["Key"])

                for filename in dayfilelist:
                    t_filelist.append(filename)

        if datagroup == "SKYECHO":
            # Each file contains multiple scans
            for filename in t_filelist:
                # we need to download the master file to be able to know the
                # scan coverage
                datapath = f'{cfg["datapath"][ind_rad]}'
                if not os.path.isdir(datapath):
                    os.makedirs(datapath)
                fname_aux = f"{datapath}{os.path.basename(filename)}"
                s3_client.download_file(cfg["s3PathRead"], filename, fname_aux)
                _, tend_sweeps, _, _ = get_sweep_time_coverage(fname_aux)
                for tend in tend_sweeps:
                    if starttime <= tend <= endtime:
                        filelist.append(
                            f"{str(filename)}::{tend.strftime('%Y-%m-%dT%H:%M:%S.%f')}"
                        )
        else:
            for filename in t_filelist:
                filenamestr = str(filename)
                fdatetime = get_datetime(filenamestr, datadescriptor)
                if fdatetime is not None:
                    if starttime <= fdatetime <= endtime:
                        if filenamestr not in filelist:
                            filelist.append(filenamestr)
        if not filelist:
            if pattern is not None:
                warn(
                    "WARNING: No file with pattern {:s} could be found between ".format(
                        pattern
                    )
                    + "starttime {:s} and endtime {:s}".format(
                        str(starttime), str(endtime)
                    )
                )
    return sorted(filelist)


def get_rad4alp_dir(
    basepath, voltime, radar_name="A", radar_res="L", scan="001", path_convention="MCH"
):
    """
    gets the directory where rad4alp data is stored

    Parameters
    ----------
    basepath : str
        base path
    voltime : datetime object
        nominal time
    radar_name : str
        radar name (A, D, L, P, W)
    radar_res : str
        radar resolution (H, L)
    scan : str
        scan
    path_convention : str
        The path convention. Can be 'LTE', 'MCH' or 'RT'

    Returns
    -------
    datapath : str
        The data path
    basename : str
        The base name. ex: PHA17213

    """
    dayinfo = voltime.strftime("%y%j")
    basename = "M" + radar_res + radar_name + dayinfo
    if path_convention == "LTE":
        yy = dayinfo[0:2]
        dy = dayinfo[2:]
        subf = "M" + radar_res + radar_name + yy + "hdf" + dy
        datapath = basepath + subf + "/"

        # check that M files exist. if not search P files
        dayfilelist = glob.glob(datapath + basename + "*." + scan + "*")
        if not dayfilelist:
            subf = "P" + radar_res + radar_name + yy + "hdf" + dy
            datapath = basepath + subf + "/"
            basename = "P" + radar_res + radar_name + dayinfo
    elif path_convention == "MCH":
        datapath = basepath + dayinfo + "/" + basename + "/"

        # check that M files exist. if not search P files
        dayfilelist = glob.glob(datapath + basename + "*." + scan + "*")
        if not dayfilelist:
            basename = "P" + radar_res + radar_name + dayinfo
            datapath = basepath + dayinfo + "/" + basename + "/"
    else:
        datapath = basepath + "M" + radar_res + radar_name + "/"

        # check that M files exist. if not search P files
        dayfilelist = glob.glob(datapath + basename + "*." + scan + "*")
        if not dayfilelist:
            basename = "P" + radar_res + radar_name + dayinfo
            datapath = basepath + "P" + radar_res + radar_name + "/"

    return datapath, basename


def get_rad4alp_grid_dir(basepath, voltime, datatype, acronym, path_convention="MCH"):
    """
    gets the directory where rad4alp grid data is stored

    Parameters
    ----------
    basepath : str
        base path
    voltime : datetime object
        nominal time
    datatype : str
        data type
    acronym : str
        acronym identifying the data type
    path_convention : str
        The path convention. Can be 'LTE', 'MCH' or 'RT'

    Returns
    -------
    datapath : str
        The data path

    """
    nowpal_accu = (
        "nowpal60_P60",
        "nowpal90_P90",
        "nowpal180_P180",
        "nowpal360_P360",
        "nowpal720_P720",
    )
    nowpal = (
        "nowpal90_P30",
        "nowpal90_P30_F60",
        "nowpal90_F60",
        "nowpal180_P60",
        "nowpal180_P60_F120",
        "nowpal180_F120",
        "nowpal360_P120",
        "nowpal360_P120_F240",
        "nowpal360_F240",
        "nowpal720_P360",
        "nowpal720_P360_F360",
        "nowpal720_F360",
    )

    cpch = (
        "CPCH0005",
        "CPCH0060",
        "CPCH0180",
        "CPCH0360",
        "CPCH0720",
        "CPCH1440",
        "CPCH2880",
        "CPCH4320",
    )

    dayinfo = voltime.strftime("%y%j")
    if datatype in nowpal_accu:
        dirbase = "nowpal_accu"
    elif datatype in nowpal:
        dirbase = "nowpal"
    elif datatype.startswith("d") and datatype != "dGZC":
        dirbase = "d" + acronym
        if datatype.endswith("H"):
            dirbase = dirbase + "H"
    elif datatype in cpch:
        dirbase = acronym + "H"
    else:
        dirbase = acronym

    if path_convention == "LTE":
        yy = dayinfo[0:2]
        dy = dayinfo[2:]
        subf = acronym + yy + "hdf" + dy
        datapath = basepath + subf + "/"
    elif path_convention == "MCH":
        datapath = basepath + dayinfo + "/" + dirbase + dayinfo + "/"
    else:
        datapath = basepath + dirbase + "/"

    return datapath


def get_trtfile_list(basepath, starttime, endtime):
    """
    gets the list of TRT files with a time period

    Parameters
    ----------
    datapath : str
        directory where to look for data
    startime : datetime object
        start of time period
    endtime : datetime object
        end of time period

    Returns
    -------
    filelist : list of strings
        list of files within the time period

    """
    startdate = starttime.date()
    enddate = endtime.date()
    ndays = int((enddate - startdate).days) + 1

    t_filelist = []
    for i in range(ndays):
        daydir = (startdate + datetime.timedelta(days=i)).strftime("%y%j")
        datapath = basepath + daydir + "/TRTC" + daydir + "/"
        dayfilelist = glob.glob(datapath + "CZC*0T.trt")
        if not dayfilelist:
            warn("No TRT files in " + datapath)
            continue
        t_filelist.extend(dayfilelist)

    filelist = []
    for filename in t_filelist:
        bfile = os.path.basename(filename)
        datetimestr = bfile[3:12]
        fdatetime = datetime.datetime.strptime(datetimestr, "%y%j%H%M")
        if starttime <= fdatetime <= endtime:
            filelist.append(filename)
        # filelist.append(filename)

    return sorted(filelist)


def get_scan_list(scandescriptor_list):
    """
    determine which is the scan list for each radar

    Parameters
    ----------
    scandescriptor : list of string
        the list of all scans for all radars

    Returns
    -------
    scan_list : list of lists
        the list of scans corresponding to each radar

    """
    descrfields = scandescriptor_list[0].split(":")
    if len(descrfields) == 1:
        # one radar
        return [scandescriptor_list]

    # one or more radars
    # check how many radars are there
    radar_list = set()
    for scandescriptor in scandescriptor_list:
        radar_list.add(scandescriptor.split(":")[0])
    nradar = len(radar_list)

    # create the list of lists
    scan_list = [[] for i in range(nradar)]
    for scandescriptor in scandescriptor_list:
        descrfields = scandescriptor.split(":")
        ind_rad = int(descrfields[0][5:8]) - 1
        scan_list[ind_rad].append(descrfields[1])

    return scan_list


def get_new_rainbow_file_name(master_fname, master_datadescriptor, datatype):
    """
    get the rainbow file name containing datatype from a master file name
    and data type

    Parameters
    ----------
    master_fname : str
        the master file name
    master_datadescriptor : str
        the master data type descriptor
    datatype : str
        the data type of the new file name to be created

    Returns
    -------
    new_fname : str
        the new file name

    """
    _, _, master_datatype, _, _ = get_datatype_fields(master_datadescriptor)
    datapath = os.path.dirname(master_fname)
    voltime = get_datetime(master_fname, master_datatype)
    voltype = os.path.basename(master_fname).split(".")[1]

    return (
        datapath
        + "/"
        + voltime.strftime("%Y%m%d%H%M%S")
        + "00"
        + datatype
        + "."
        + voltype
    )


def get_datatype_fields(datadescriptor):
    """
    splits the data type descriptor and provides each individual member

    Parameters
    ----------
    datadescriptor : str
        radar field type. Format : [radar file type]:[datatype]

    Returns
    -------
    radarnr : str
        radar number, i.e. RADAR1, RADAR2, ...
    datagroup : str
        data type group, i.e. RAINBOW, RAD4ALP, ODIM, CFRADIAL, ICON,
        MXPOL ...
    datatype : str
        data type, i.e. dBZ, ZDR, ISO0, ...
    dataset : str
        dataset type (for saved data only)
    product : str
        product type (for saved data only)

    """
    descrfields = datadescriptor.split(":")
    if len(descrfields) == 1:
        radarnr = "RADAR001"
        datagroup = "RAINBOW"
        datatype = descrfields[0]
        dataset = None
        product = None
    elif descrfields[0].startswith("RADAR"):
        radarnr = descrfields[0]
        if len(descrfields) == 2:
            radarnr = descrfields[0]
            datagroup = "RAINBOW"
            datatype = descrfields[1]
            dataset = None
            product = None
        else:
            datagroup = descrfields[1]
            if datagroup in (
                "CFRADIALPYRAD",
                "ODIMPYRAD",
                "PYRADGRID",
                "ODIMPYRADGRID",
                "NETCDFSPECTRA",
                "CSV",
                "GECSX",
            ):
                descrfields2 = descrfields[2].split(",")
                datatype = descrfields2[0]
                dataset = descrfields2[1]
                product = descrfields2[2]
            elif datagroup == "CFRADIALICON":
                descrfields2 = descrfields[2].split(",")
                datatype = descrfields2[0]
                dataset = descrfields2[1]
                product = None
            elif datagroup == "MXPOL":
                datatype = descrfields[2]
                dataset = None
                product = None
            elif datagroup in (
                "ODIM",
                "ODIMBIRDS",
                "MFCFRADIAL",
                "MFBIN",
                "CFRADIAL2",
                "CF1",
                "NEXRADII",
                "MFPNG",
                "MFGRIB",
                "MFDAT",
                "MFCF",
                "GAMIC",
                "CFRADIAL",
                "ODIMGRID",
                "SKYECHO",
                "KNMIH5GRID",
            ):
                descrfields2 = descrfields[2].split(",")
                datatype = descrfields2[0]
                product = None
                dataset = None
                if np.size(descrfields2) == 2:
                    dataset = descrfields2[1]
            else:
                datatype = descrfields[2]
                dataset = None
                product = None
    else:
        radarnr = "RADAR001"
        datagroup = descrfields[0]
        if datagroup in (
            "CFRADIALPYRAD",
            "ODIMPYRAD",
            "PYRADGRID",
            "ODIMPYRADGRID",
            "NETCDFSPECTRA",
            "CSV",
        ):
            descrfields2 = descrfields[1].split(",")
            datatype = descrfields2[0]
            dataset = descrfields2[1]
            product = descrfields2[2]
        elif datagroup == "CFRADIALICON":
            descrfields2 = descrfields[1].split(",")
            datatype = descrfields2[0]
            dataset = descrfields2[1]
            product = None
        elif datagroup == "MXPOL":
            datatype = descrfields[1]
            dataset = None
            product = None
        elif datagroup in (
            "ODIM",
            "ODIMBIRDS",
            "MFCFRADIAL",
            "MFBIN",
            "NEXRADII",
            "MFPNG",
            "MFGRIB",
            "MFDAT",
            "MFCF",
            "CFRADIAL2",
            "CF1",
            "GAMIC",
            "CFRADIAL",
            "ODIMGRID",
            "SKYECHO",
            "KNMIH5GRID",
        ):
            descrfields2 = descrfields[1].split(",")
            # warn(" descrfields2:  '%s'" % descrfields2[1])
            if len(descrfields2) == 2:
                datatype = descrfields2[0]
                dataset = descrfields2[1]
                product = None
                # warn(" dataset:  '%s'" % dataset)
            else:
                datatype = descrfields[1]
                dataset = None
                product = None
        else:
            datatype = descrfields[1]
            dataset = None
            product = None
    # warn(" dataset:  '%s'" % dataset)
    return radarnr, datagroup, datatype, dataset, product


def get_dataset_fields(datasetdescr):
    """
    splits the dataset type descriptor and provides each individual member

    Parameters
    ----------
    datasetdescr : str
        dataset type. Format : [processing level]:[dataset type]

    Returns
    -------
    proclevel : str
        dataset processing level

    dataset : str
        dataset type, i.e. dBZ, ZDR, ISO0, ...

    """
    descrfields = datasetdescr.split(":")
    if len(descrfields) == 1:
        proclevel = "l00"
        dataset = descrfields[0]
    else:
        proclevel = descrfields[0]
        dataset = descrfields[1]
        if len(proclevel) == 2:
            proclevel = proclevel[0] + "0" + proclevel[1]

    return proclevel, dataset


def get_datetime(fname, datadescriptor):
    """
    Given a data descriptor gets date and time from file name

    Parameters
    ----------
    fname : str
        file name
    datadescriptor : str
        radar field type. Format : [radar file type]:[datatype]

    Returns
    -------
    fdatetime : datetime object
        date and time in file name

    """
    _, datagroup, _, dataset, _ = get_datatype_fields(datadescriptor)

    return _get_datetime(fname, datagroup, ftime_format=dataset)


def find_icon_file(voltime, datatype, cfg, scanid, ind_rad=0):
    """
    Search a ICON file in Rainbow format

    Parameters
    ----------
    voltime : datetime object
        volume scan time
    datatype : str
        type of ICON data to look for
    cfg : dictionary of dictionaries
        configuration info to figure out where the data is
    scanid : str
        name of the scan
    ind_rad : int
        radar index

    Returns
    -------
    fname : str
        Name of ICON file if it exists. None otherwise

    """
    # hour rounded date-time
    fdatetime = voltime.strftime("%Y%m%d%H") + "000000"

    # initial run time to look for
    hvol = int(voltime.strftime("%H"))
    runhour0 = int(hvol / cfg["IconRunFreq"]) * cfg["IconRunFreq"]
    runtime0 = voltime.replace(hour=runhour0, minute=0, second=0)

    # look for icon file
    found = False
    nruns_to_check = int((cfg["IconForecasted"] - 1) / cfg["IconRunFreq"])
    for i in range(nruns_to_check):
        runtime = runtime0 - datetime.timedelta(hours=i * cfg["IconRunFreq"])
        runtimestr = runtime.strftime("%Y%m%d%H") + "000000"

        daydir = runtime.strftime("%Y-%m-%d")
        datapath = cfg["iconpath"][ind_rad] + datatype + "/" + scanid + daydir + "/"

        search_name = (
            datapath + datatype + "_RUN" + runtimestr + "_DX50" + fdatetime + ".*"
        )
        print("Looking for file: " + search_name)
        fname = glob.glob(search_name)
        if fname:
            found = True
            break

    if found:
        return fname[0]

    warn("WARNING: Unable to get ICON " + datatype + " information")
    return None


def find_pyradicon_file(basepath, voltime, datatype, cfg, dataset):
    """
    Search a ICON file in CFRadial or ODIM format

    Parameters
    ----------
    basepath : str
        base path to the ICON file
    voltime : datetime object
        volume scan time
    datatype : str
        type of ICON data to look for
    cfg : dictionary of dictionaries
        configuration info to figure out where the data is
    dataset : str
        name of the folder where the data is stored

    Returns
    -------
    fname : str
        Name of ICON file if it exists. None otherwise

    """
    # hour rounded date-time
    fdatetime = voltime.strftime("%Y%m%d%H") + "0000"

    # initial run time to look for
    hvol = int(voltime.strftime("%H"))
    runhour0 = int(hvol / cfg["IconRunFreq"]) * cfg["IconRunFreq"]
    runtime0 = voltime.replace(hour=runhour0, minute=0, second=0)

    # look for icon file
    found = False
    nruns_to_check = int((cfg["IconForecasted"] - 1) / cfg["IconRunFreq"])
    for i in range(nruns_to_check):
        runtime = runtime0 - datetime.timedelta(hours=i * cfg["IconRunFreq"])
        runtimestr = runtime.strftime("%Y%m%d%H") + "0000"

        daydir = runtime.strftime("%Y-%m-%d")
        datapath = basepath + datatype + "/radar/" + daydir + "/" + dataset + "/"

        search_name = datapath + datatype + "_RUN" + runtimestr + "_" + fdatetime + ".*"
        print("Looking for file: " + search_name)
        fname = glob.glob(search_name)
        if fname:
            found = True
            break

    if found:
        return fname[0]

    warn("WARNING: Unable to get ICON " + datatype + " information")
    return None


def find_raw_icon_file(voltime, datatype, cfg, ind_rad=0):
    """
    Search a ICON file in netcdf format

    Parameters
    ----------
    voltime : datetime object
        volume scan time
    datatype : str
        type of ICON data to look for
    cfg : dictionary of dictionaries
        configuration info to figure out where the data is
    ind_rad : int
        radar index

    Returns
    -------
    fname : str
        Name of ICON file if it exists. None otherwise

    """
    # initial run time to look for
    runhour0 = int(voltime.hour / cfg["IconRunFreq"]) * cfg["IconRunFreq"]
    runtime0 = voltime.replace(hour=runhour0, minute=0, second=0)
    # look for icon file in raw
    found = False
    nruns_to_check = int(cfg["IconForecasted"] / cfg["IconRunFreq"])
    for i in range(nruns_to_check):
        runtime = runtime0 - datetime.timedelta(hours=i * cfg["IconRunFreq"])
        runtimestr = runtime.strftime("%Y_%m_%d_%H")
        daydir = runtime.strftime("%Y-%m-%d")
        datapath = cfg["iconpath"][ind_rad] + datatype + "/raw*/" + daydir + "/"
        for model in ("icon-ch1-eps", "icon-ch2-eps"):
            if datatype == "TEMP":
                search_name = (
                    datapath + runtimestr + "_" + model + "_MDR_3D_m_000" + ".nc"
                )
            elif datatype == "WIND":
                search_name = (
                    datapath + runtimestr + "_" + model + "_MDR_3DWIND_m_000" + ".nc"
                )
            else:
                warn("Unable to get ICON " + datatype + ". Unknown variable")
            print("Looking for file: " + search_name)
            fname = glob.glob(search_name)
            if fname:
                found = True
                break

        if found:
            break

    if found:
        return fname[0]

    warn("WARNING: Unable to get ICON " + datatype + " information")
    return None


def find_hzt_file(voltime, cfg, ind_rad=0):
    """
    Search an ISO-0 degree file in HZT format

    Parameters
    ----------
    voltime : datetime object
        volume scan time
    cfg : dictionary of dictionaries
        configuration info to figure out where the data is
    ind_rad : int
        radar index

    Returns
    -------
    fname : str
        Name of HZT file if it exists. None otherwise

    """
    # initial run time to look for
    runhour0 = int(voltime.hour / cfg["IconRunFreq"]) * cfg["IconRunFreq"]
    runtime0 = voltime.replace(hour=runhour0, minute=0, second=0)

    # look for icon file
    found = False
    nruns_to_check = int((cfg["IconForecasted"] - 1) / cfg["IconRunFreq"])
    for i in range(nruns_to_check):
        runtime = runtime0 - datetime.timedelta(hours=i * cfg["IconRunFreq"])
        target_hour = int((voltime - runtime).total_seconds() / 3600.0)
        runtimestr = runtime.strftime("%y%j%H00")

        daydir = runtime.strftime("%y%j")
        if cfg["path_convention"][ind_rad] == "RT":
            datapath = cfg["iconpath"][ind_rad] + "HZT/"
        else:
            datapath = cfg["iconpath"][ind_rad] + "HZT/" + daydir + "/"
        search_name = (
            datapath + "HZT" + runtimestr + "0L.8" + "{:02d}".format(target_hour)
        )

        print("Looking for file: " + search_name)
        fname = glob.glob(search_name)
        if fname:
            found = True
            break

    if not found:
        warn("WARNING: Unable to find HZT file")
        return None

    return fname[0]


def find_iso0_file(voltime, cfg, ind_rad=0):
    """
    Search an ISO-0 degree file in text format

    Parameters
    ----------
    voltime : datetime object
        volume scan time
    cfg : dictionary of dictionaries
        configuration info to figure out where the data is
    ind_rad : int
        radar index

    Returns
    -------
    fname : str
        Name of iso0 file if it exists. None otherwise

    """
    # initial run time to look for
    runhour0 = int(voltime.hour / cfg["IconRunFreq"]) * cfg["IconRunFreq"]
    runtime0 = voltime.replace(hour=runhour0, minute=0, second=0)

    radar_name = mf_sname_to_wmo_number(cfg["RadarName"][ind_rad])
    # look for file
    found = False
    nruns_to_check = int((cfg["IconForecasted"]) / cfg["IconRunFreq"])
    for i in range(nruns_to_check):
        runtime = runtime0 - datetime.timedelta(hours=i * cfg["IconRunFreq"])
        int((voltime - runtime).total_seconds() / 3600.0)
        runtimestr = runtime.strftime("%Y%m%d%H0000")

        datapath = cfg["iconpath"][ind_rad]
        search_name = datapath + "bdap_iso0_" + radar_name + "_" + runtimestr + ".txt"

        print("Looking for file: " + search_name)
        fname = glob.glob(search_name)
        if fname:
            found = True
            break

    if not found:
        warn("WARNING: Unable to find iso0 file")
        return None

    return fname[0]


def find_iso0_grib_file(voltime, cfg, ind_rad=0):
    """
    Search an ISO-0 degree file in text format

    Parameters
    ----------
    voltime : datetime object
        volume scan time
    cfg : dictionary of dictionaries
        configuration info to figure out where the data is
    ind_rad : int
        radar index

    Returns
    -------
    fname : str
        Name of iso0 file if it exists. None otherwise

    """
    datapath = cfg["iconpath"][ind_rad]

    if cfg["IconRunFreq"] == 0:
        # The date of the NWP file corresponds to the data of the radar
        runtimestr = voltime.strftime("%Y%m%d%H%M")
        search_name = datapath + "ISO_T_PAROME_" + runtimestr + "*.grib"
        print("Looking for file: {}".format(search_name))
        fname = glob.glob(search_name)
        if fname:
            return fname[0]
        warn("WARNING: Unable to find iso0 file")
        return None

    runhour0 = int(voltime.hour / cfg["IconRunFreq"]) * cfg["IconRunFreq"]
    runtime0 = voltime.replace(hour=runhour0, minute=0, second=0)
    nruns_to_check = int((cfg["IconForecasted"]) / cfg["IconRunFreq"])

    # look for file
    found = False
    for i in range(nruns_to_check):
        runtime = runtime0 - datetime.timedelta(hours=i * cfg["IconRunFreq"])
        int((voltime - runtime).total_seconds() / 3600.0)
        runtimestr = runtime.strftime("%Y%m%d%H00")
        search_name = datapath + "ISO_T_PAROME_" + runtimestr + "*.grib"

        print("Looking for file: {}".format(search_name))
        fname = glob.glob(search_name)
        if fname:
            found = True
            break

    if not found:
        warn("WARNING: Unable to find iso0 file")
        return None

    return fname[0]


def find_rad4alpicon_file(voltime, datatype, cfg, scanid, ind_rad=0):
    """
    Search a ICON file

    Parameters
    ----------
    voltime : datetime object
        volume scan time
    datatype : str
        type of ICON data to look for
    cfg: dictionary of dictionaries
        configuration info to figure out where the data is
    ind_rad: int
        radar index

    Returns
    -------
    fname : str
        Name of ICON file if it exists. None otherwise

    scanid: str
        name of the scan

    """
    # hour rounded date-time
    fdatetime = voltime.strftime("%y%j%H") + "00"

    # initial run time to look for
    hvol = int(voltime.strftime("%H"))
    runhour0 = int(hvol / cfg["IconRunFreq"]) * cfg["IconRunFreq"]
    runtime0 = voltime.replace(hour=runhour0, minute=0, second=0)

    # look for icon file
    found = False
    nruns_to_check = int((cfg["IconForecasted"] - 1) / cfg["IconRunFreq"])
    rad_id = "P" + cfg["RadarRes"][ind_rad] + cfg["RadarName"][ind_rad]
    for i in range(nruns_to_check):
        runtime = runtime0 - datetime.timedelta(hours=i * cfg["IconRunFreq"])
        runtimestr = runtime.strftime("%y%j%H") + "00"

        daydir = runtime.strftime("%y%j")
        datapath = (
            cfg["iconpath"][ind_rad] + datatype + "/" + rad_id + "/" + daydir + "/"
        )

        search_name = (
            datapath
            + datatype
            + "_RUN"
            + runtimestr
            + "_"
            + rad_id
            + fdatetime
            + "."
            + scanid
            + ".bin"
        )
        print("Looking for file: " + search_name)
        fname = glob.glob(search_name)
        if fname:
            found = True
            break

    if not found:
        warn("WARNING: Unable to get ICON " + datatype + " information")
        return None

    return fname[0]


def _get_datetime(fname, datagroup, ftime_format=None):
    """
    Given a data group gets date and time from file name

    Parameters
    ----------
    fname : str
        file name
    datadescriptor : str
        radar field type. Format : [radar file type]:[datatype]
    ftime_format : str or None
        if the file is of type ODIM this contain the file time format

    Returns
    -------
    fdatetime : datetime object
        date and time in file name

    """
    bfile = os.path.basename(fname)
    if datagroup in (
        "RAINBOW",
        "CFRADIALPYRAD",
        "ODIMPYRAD",
        "PYRADGRID",
        "ODIMPYRADGRID",
        "NETCDFSPECTRA",
        "CSV",
    ):
        datetimestr = bfile[0:14]
        fdatetime = datetime.datetime.strptime(datetimestr, "%Y%m%d%H%M%S")
    elif datagroup in ("RAD4ALP", "RAD4ALPGRID", "RAD4ALPGIF", "RAD4ALPBIN"):
        datestr = bfile[3:8]
        timestr = bfile[8:12]
        if timestr != "2400":
            fdatetime = datetime.datetime.strptime(datestr + timestr, "%y%j%H%M")
        else:
            fdatetime = datetime.datetime.strptime(
                datestr, "%y%j"
            ) + datetime.timedelta(days=1)
    elif datagroup in (
        "ODIM",
        "ODIMBIRDS",
        "MFCFRADIAL",
        "MFBIN",
        "NEXRADII",
        "MFPNG",
        "MFGRIB",
        "MFDAT",
        "MFCF",
        "CFRADIAL2",
        "CF1",
        "GAMIC",
        "CFRADIAL",
        "ODIMGRID",
        "KNMIH5GRID",
    ):
        if ftime_format is None:
            # we assume is rad4alp format
            datetimestr = bfile[3:12]
            fdatetime = datetime.datetime.strptime(datetimestr, "%y%j%H%M")
        else:
            return find_date_in_file_name(
                bfile, date_format=ftime_format[ftime_format.find("F") + 2 : -1]
            )
    elif datagroup == "MXPOL":
        datetimestr = re.findall(r"([0-9]{8}-[0-9]{6})", bfile)[0]
        fdatetime = datetime.datetime.strptime(datetimestr, "%Y%m%d-%H%M%S")
    elif datagroup == "ICONRAW":
        datetimestr = bfile[-13:-3]
        fdatetime = datetime.datetime.strptime(datetimestr, "%Y%m%d%H")
    elif datagroup == "SATGRID":
        datetimestr = bfile[10:22]
        fdatetime = datetime.datetime.strptime(datetimestr, "%Y%m%d%H%M")
    elif datagroup == "SKYECHO":
        datetimestr = bfile.split("::")[1]
        fdatetime = datetime.datetime.strptime(datetimestr, "%Y-%m-%dT%H:%M:%S.%f")
    else:
        warn("unknown data group")
        return None

    return fdatetime


def find_date_in_file_name(filename, date_format="%Y%m%d%H%M%S"):
    """
    Find a date with date format defined in date_format in a file name.
    If no date is found returns None

    Parameters
    ----------
    filename : str
        file name
    date_format : str
        The time format

    Returns
    -------
    fdatetime : datetime object
        date and time in file name

    """
    today = datetime.datetime.now()
    len_datestr = len(today.strftime(date_format))
    count = 0
    bfile = os.path.basename(filename)
    while True:
        try:
            fdatetime = datetime.datetime.strptime(
                bfile[count : count + len_datestr], date_format
            )
        except ValueError:
            count += 1
            if count + len_datestr > len(bfile):
                warn(
                    f"Unable to find date from string name. Date format "
                    f"{date_format}. File name {bfile}"
                )
                return None
        else:
            # No error, stop the loop
            break

    return fdatetime


def convert_pydda_to_pyart_grid(pydda_grid):
    """
    Converts a PyDDA Dataset back to a Py-ART Grid with the necessary variables.

    Parameters
    ----------
    pydda_grid: xarray.Dataset
        The PyDDA Dataset to convert back to a Py-ART Grid.

    Returns
    -------
    grid: Py-ART Grid
        The Py-ART Grid reconstructed from the PyDDA Dataset.
    """

    # Extract the basic grid properties
    pydda_grid["time"].values
    z = pydda_grid["z"].values
    y = pydda_grid["y"].values
    x = pydda_grid["x"].values
    pydda_grid["point_latitude"].values
    pydda_grid["point_longitude"].values

    # Reconstruct fields
    fields = {}
    for var_name in pydda_grid.data_vars:
        if var_name not in [
            "time",
            "z",
            "y",
            "x",
            "point_latitude",
            "point_longitude",
            "point_altitude",
            "ROI",
            "projection",
            "radar_latitude",
            "radar_longitude",
            "radar_altitude",
            "origin_longitude",
            "origin_latitude",
            "origin_altitude",
            "point_x",
            "point_z",
            "point_y",
            "AZ",
            "EL",
        ]:
            field_data = pydda_grid[var_name].values.squeeze()
            field_metadata = {
                k: v for k, v in pydda_grid[var_name].attrs.items() if k != "data"
            }
            fields[var_name] = {"data": field_data}
            fields[var_name].update(field_metadata)

    # Extract the time information
    pydda_time = pydda_grid["time"].values[0]
    rounded_time = pydda_time.replace(microsecond=0)

    # Convert cftime.DatetimeGregorian to seconds since a reference time
    reference_time_str = pydda_time.strftime("%Y-%m-%dT%H:%M:%SZ")
    time_in_seconds = np.array(
        [(pydda_time - rounded_time).total_seconds()], dtype=np.float32
    )

    time_dict = {
        "units": f"seconds since {reference_time_str}",
        "standard_name": "time",
        "long_name": "Time of grid",
        "calendar": "gregorian",
        "data": time_in_seconds,
    }

    # Reconstructing the Py-ART Grid object
    grid = pyart.core.Grid(
        time=time_dict,
        fields=fields,
        metadata=pydda_grid.attrs,
        origin_latitude={"data": pydda_grid["origin_latitude"].values},
        origin_longitude={"data": pydda_grid["origin_longitude"].values},
        origin_altitude={"data": pydda_grid["origin_altitude"].values},
        x={"data": x},
        y={"data": y},
        z={"data": z},
        projection=pydda_grid["projection"].attrs,
        radar_latitude={"data": pydda_grid["radar_latitude"].values},
        radar_longitude={"data": pydda_grid["radar_longitude"].values},
        radar_altitude={"data": pydda_grid["radar_altitude"].values},
    )
    return grid
