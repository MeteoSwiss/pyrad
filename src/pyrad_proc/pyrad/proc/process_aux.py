"""
pyrad.proc.process_aux
======================

Auxiliary functions. Functions to determine the process type, pass raw data to
the product generation functions, save radar data and extract data at
determined points or regions of interest.

.. autosummary::
    :toctree: generated/

    get_process_func
    process_raw
    process_vol_to_grid
    process_save_radar
    process_fixed_rng
    process_fixed_rng_span
    process_roi
    process_azimuthal_average
    process_moving_azimuthal_average
    process_radar_resampling
    _get_values_antenna_pattern
    _create_target_radar

"""

from time import time

from copy import deepcopy
from warnings import warn
import numpy as np
from scipy.spatial import cKDTree

from pyart.util import cut_radar, compute_directional_stats
from pyart.util import find_neighbour_gates
from pyart.util import compute_azimuthal_average
from pyart.config import get_metadata
from pyart.core import Radar, geographic_to_cartesian_aeqd
from pyart.core import cartesian_to_geographic_aeqd
from pyart.map import grid_from_radars

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..io.read_data_sensor import read_trt_traj_data
from ..io.read_data_other import read_antenna_pattern
from ..io.read_data_cosmo import _put_radar_in_swiss_coord
from ..util.radar_utils import belongs_roi_indices
from ..util.radar_utils import get_fixed_rng_data, get_cercle_coords
from ..util.stat_utils import quantiles_weighted
from ..proc.process_traj import _get_gates_antenna_pattern


def get_process_func(dataset_type, dsname):
    """
    Maps the dataset type into its processing function and data set format
    associated.

    Parameters
    ----------
    dataset_type : str
        The following is a list of data set types ordered by type of output
        dataset with the function they call. For details of what they do check
        the function documentation:
            'VOL' format output:
                'ATTENUATION': process_attenuation
                'AZI_AVG': process_azimuthal_average
                'MOVING_AZI_AVG': process_moving_azimuthal_average
                'BIAS_CORRECTION': process_correct_bias
                'BIRDS_ID': process_birds_id
                'BIRD_DENSITY': process_bird_density
                'CCOR': process_ccor
                'CDF': process_cdf
                'CDR': process_cdr
                'CLT_TO_SAN': process_clt_to_echo_id
                'COSMO': process_cosmo
                'COSMO_LOOKUP': process_cosmo_lookup_table
                'DEM': process_dem
                'DEALIAS_FOURDD': process_dealias_fourdd
                'DEALIAS_REGION': process_dealias_region_based
                'DEALIAS_UNWRAP': process_dealias_unwrap_phase
                'DOPPLER_VELOCITY': process_Doppler_velocity
                'DOPPLER_VELOCITY_IQ': process_Doppler_velocity_iq
                'DOPPLER_WIDTH': process_Doppler_width
                'DOPPLER_WIDTH_IQ': process_Doppler_width_iq
                'ECHO_FILTER': process_echo_filter
                'FIELDS_DIFF': process_fields_diff
                'FIXED_RNG': process_fixed_rng
                'FIXED_RNG_SPAN': process_fixed_rng_span
                'GECSX' : process_gecsx
                'hydroMF_to_hydro': process_hydro_mf_to_hydro
                'hydroMF_to_SAN: process_hydro_mf_to_echo_id
                'HYDROCLASS': process_hydroclass
                'HZT': process_hzt
                'HZT_LOOKUP': process_hzt_lookup_table
                'ISO0_GRIB': process_iso0_grib
                'ISO0_MF': process_iso0_mf
                'KDP_LEASTSQUARE_1W': process_kdp_leastsquare_single_window
                'KDP_LEASTSQUARE_2W': process_kdp_leastsquare_double_window
                'L': process_l
                'MEAN_PHASE_IQ': process_mean_phase_iq
                'NCVOL': process_save_radar
                'NOISE_POWER': process_noise_power
                'OUTLIER_FILTER': process_outlier_filter
                'PhiDP': process_differential_phase
                'PHIDP0_CORRECTION': process_correct_phidp0
                'PHIDP0_ESTIMATE': process_estimate_phidp0
                'PhiDP_IQ': process_differential_phase_iq
                'PHIDP_KDP_KALMAN': process_phidp_kdp_Kalman
                'PHIDP_KDP_LP': process_phidp_kdp_lp
                'PHIDP_KDP_VULPIANI': process_phidp_kdp_Vulpiani
                'PHIDP_SMOOTH_1W': process_smooth_phidp_single_window
                'PHIDP_SMOOTH_2W': process_smooth_phidp_double_window
                'POL_VARIABLES': process_pol_variables
                'POL_VARIABLES_IQ': process_pol_variables_iq
                'PWR': process_signal_power
                'RADAR_RESAMPLING': process_radar_resampling
                'RADIAL_NOISE_HS': process_radial_noise_hs
                'RADIAL_NOISE_IVIC': process_radial_noise_ivic
                'RADIAL_VELOCITY': process_radial_velocity
                'RAINRATE': process_rainrate
                'RAW': process_raw
                'REFLECTIVITY': process_reflectivity
                'REFLECTIVITY_IQ': process_reflectivity_iq
                'RCS': process_rcs
                'RCS_PR': process_rcs_pr
                'RhoHV': process_rhohv
                'RhoHV_IQ': process_rhohv_iq
                'RHOHV_CORRECTION': process_correct_noise_rhohv
                'RHOHV_RAIN': process_rhohv_rain
                'ROI': process_roi
                'SAN': process_echo_id
                'SELFCONSISTENCY_BIAS': process_selfconsistency_bias
                'SELFCONSISTENCY_BIAS2': process_selfconsistency_bias2
                'SELFCONSISTENCY_KDP_PHIDP': process_selfconsistency_kdp_phidp
                'SNR': process_snr
                'SNR_FILTER': process_filter_snr
                'ST1_IQ': process_st1_iq
                'ST2_IQ': process_st2_iq
                'TRAJ_TRT' : process_traj_trt
                'TRAJ_TRT_CONTOUR' : process_traj_trt_contour
                'TURBULENCE': process_turbulence
                'VAD': process_vad
                'VEL_FILTER': process_filter_vel_diff
                'VIS': process_visibility
                'VIS_FILTER': process_filter_visibility
                'VOL_REFL': process_vol_refl
                'VPR': process_vpr
                'WBN': process_wbn_iq
                'WIND_VEL': process_wind_vel
                'WINDSHEAR': process_windshear
                'ZDR': process_differential_reflectivity
                'ZDR_IQ': process_differential_reflectivity_iq
                'ZDR_PREC': process_zdr_precip
                'ZDR_SNOW': process_zdr_snow
            'SPECTRA' format output:
                'FFT': process_fft
                'FILTER_0DOPPLER': process_filter_0Doppler
                'FILTER_SPECTRA_NOISE': process_filter_spectra_noise
                'IFFT': process_ifft
                'RAW_IQ': process_raw_iq
                'RAW_SPECTRA': process_raw_spectra
                'SPECTRA_ANGULAR_AVERAGE': process_spectra_ang_avg
                'SPECTRA_POINT': process_spectra_point
                'SPECTRAL_NOISE': process_spectral_noise
                'SPECTRAL_PHASE': process_spectral_phase
                'SPECTRAL_POWER': process_spectral_power
                'SPECTRAL_REFLECTIVITY': process_spectral_reflectivity
                'sPhiDP': process_spectral_differential_phase
                'sRhoHV': process_spectral_RhoHV
                'SRHOHV_FILTER': process_filter_srhohv
                'sZDR': process_spectral_differential_reflectivity
            'CENTROIDS' format output:
                'CENTROIDS': process_centroids
            'COLOCATED_GATES' format output:
                'COLOCATED_GATES': process_colocated_gates
            'COSMO_COORD' format output:
                'COSMO_COORD': process_cosmo_coord
                'HZT_COORD': process_hzt_coord
            'COSMO2RADAR' format output:
                'COSMO2RADAR': process_cosmo_to_radar
            'GRID' format output:
                'RAW_GRID': process_raw_grid
                'GECSX' : process_gecsx
                'GRID': process_grid
                'GRID_FIELDS_DIFF': process_grid_fields_diff
                'GRID_MASK': process_grid_mask
                'GRID_TEXTURE': process_grid_texture
                'NORMALIZE_LUMINOSITY': process_normalize_luminosity
                'PIXEL_FILTER': process_pixel_filter
                'VOL2GRID': process_vol_to_grid
            'GRID_TIMEAVG' format output:
                'GRID_TIME_STATS': process_grid_time_stats
                'GRID_TIME_STATS2': process_grid_time_stats2
            'INTERCOMP' format output:
                'INTERCOMP': process_intercomp
                'INTERCOMP_FIELDS': process_intercomp_fields
                'INTERCOMP_TIME_AVG': process_intercomp_time_avg
            'ML' format output:
                'ML_DETECTION': process_melting_layer
            'MONITORING' format output:
                'GC_MONITORING': process_gc_monitoring
                'MONITORING': process_monitoring
            'OCCURRENCE' format output:
                'OCCURRENCE': process_occurrence
                'OCCURRENCE_PERIOD': process_occurrence_period
                'TIMEAVG_STD': process_time_avg_std
            'QVP' format output:
                'EVP': process_evp
                'QVP': process_qvp
                'rQVP': process_rqvp
                'SVP': process_svp
                'TIME_HEIGHT': process_time_height
                'TIME_ALONG_COORD': process_ts_along_coord
            'SPARSE_GRID' format output:
                'ZDR_COLUMN': process_zdr_column
            'SUN_HITS' format output:
                'SUN_HITS': process_sun_hits
            'SUNSCAN' format output:
                'SUNSCAN': process_sunscan
            'TIMEAVG' format output:
                'FLAG_TIME_AVG': process_time_avg_flag
                'TIME_AVG': process_time_avg
                'WEIGHTED_TIME_AVG': process_weighted_time_avg
                'TIME_STATS': process_time_stats
                'TIME_STATS2': process_time_stats2
                'RAIN_ACCU': process_rainfall_accumulation
            'TIMESERIES' format output:
                'GRID_POINT_MEASUREMENT': process_grid_point
                'POINT_MEASUREMENT': 'process_point_measurement'
                'TRAJ_ANTENNA_PATTERN': process_traj_antenna_pattern
                'TRAJ_ATPLANE': process_traj_atplane
                'TRAJ_LIGHTNING': process_traj_lightning
            'TRAJ_ONLY' format output:
                'TRAJ': process_trajectory
    dsname : str
        Name of dataset

    Returns
    -------
    func_name : str or processing function
        pyrad function used to process the data set type
    dsformat : str
        data set format, i.e.: 'VOL', etc.

    """

    dsformat = 'VOL'
    if dataset_type == 'RAW':
        func_name = process_raw
    elif dataset_type == 'AZI_AVG':
        func_name = process_azimuthal_average
    elif dataset_type == 'MOVING_AZI_AVG':
        func_name = process_moving_azimuthal_average
    elif dataset_type == 'RADAR_RESAMPLING':
        func_name = 'process_radar_resampling'
    elif dataset_type == 'CCOR':
        func_name = 'process_ccor'
    elif dataset_type == 'GECSX':
        func_name = 'process_gecsx'
        dsformat = ['GRID', 'VOL']
    elif dataset_type == 'GRID':
        func_name = 'process_grid'
        dsformat = 'GRID'
    elif dataset_type == 'RAW_GRID':
        func_name = 'process_raw_grid'
        dsformat = 'GRID'
    elif dataset_type == 'VOL2GRID':
        func_name = 'process_vol_to_grid'
        dsformat = 'GRID'
    elif dataset_type == 'GRID_FIELDS_DIFF':
        func_name = 'process_grid_fields_diff'
        dsformat = 'GRID'
    elif dataset_type == 'GRID_MASK':
        func_name = 'process_grid_mask'
        dsformat = 'GRID'
    elif dataset_type == 'GRID_TEXTURE':
        func_name = 'process_grid_texture'
        dsformat = 'GRID'
    elif dataset_type == 'NORMALIZE_LUMINOSITY':
        func_name = 'process_normalize_luminosity'
        dsformat = 'GRID'
    elif dataset_type == 'PIXEL_FILTER':
        func_name = 'process_pixel_filter'
        dsformat = 'GRID'
    elif dataset_type == 'RAW_SPECTRA':
        func_name = 'process_raw_spectra'
        dsformat = 'SPECTRA'
    elif dataset_type == 'SPECTRA_POINT':
        func_name = 'process_spectra_point'
        dsformat = 'SPECTRA'
    elif dataset_type == 'IFFT':
        func_name = 'process_ifft'
        dsformat = 'SPECTRA'
    elif dataset_type == 'SPECTRAL_POWER':
        func_name = 'process_spectral_power'
        dsformat = 'SPECTRA'
    elif dataset_type == 'SPECTRAL_NOISE':
        func_name = 'process_spectral_noise'
        dsformat = 'SPECTRA'
    elif dataset_type == 'SPECTRAL_PHASE':
        func_name = 'process_spectral_phase'
        dsformat = 'SPECTRA'
    elif dataset_type == 'SPECTRAL_REFLECTIVITY':
        func_name = 'process_spectral_reflectivity'
        dsformat = 'SPECTRA'
    elif dataset_type == 'sZDR':
        func_name = 'process_spectral_differential_reflectivity'
        dsformat = 'SPECTRA'
    elif dataset_type == 'sPhiDP':
        func_name = 'process_spectral_differential_phase'
        dsformat = 'SPECTRA'
    elif dataset_type == 'sRhoHV':
        func_name = 'process_spectral_rhohv'
        dsformat = 'SPECTRA'
    elif dataset_type == 'FILTER_SPECTRA_NOISE':
        func_name = 'process_filter_spectra_noise'
        dsformat = 'SPECTRA'
    elif dataset_type == 'FILTER_0DOPPLER':
        func_name = 'process_filter_0Doppler'
        dsformat = 'SPECTRA'
    elif dataset_type == 'SRHOHV_FILTER':
        func_name = 'process_filter_srhohv'
        dsformat = 'SPECTRA'
    elif dataset_type == 'SPECTRA_ANGULAR_AVERAGE':
        func_name = 'process_spectra_ang_avg'
        dsformat = 'SPECTRA'
    elif dataset_type == 'FFT':
        func_name = 'process_fft'
        dsformat = 'SPECTRA'
    elif dataset_type == 'RAW_IQ':
        func_name = 'process_raw_iq'
        dsformat = 'SPECTRA'
    elif dataset_type == 'QVP':
        func_name = 'process_qvp'
        dsformat = 'QVP'
    elif dataset_type == 'rQVP':
        func_name = 'process_rqvp'
        dsformat = 'QVP'
    elif dataset_type == 'SVP':
        func_name = 'process_svp'
        dsformat = 'QVP'
    elif dataset_type == 'EVP':
        func_name = 'process_evp'
        dsformat = 'QVP'
    elif dataset_type == 'TIME_HEIGHT':
        func_name = 'process_time_height'
        dsformat = 'QVP'
    elif dataset_type == 'TIME_ALONG_COORD':
        func_name = 'process_ts_along_coord'
        dsformat = 'QVP'
    elif dataset_type == 'CDF':
        func_name = 'process_cdf'
    elif dataset_type == 'NCVOL':
        func_name = process_save_radar
    elif dataset_type == 'PWR':
        func_name = 'process_signal_power'
    elif dataset_type == 'RCS_PR':
        func_name = 'process_rcs_pr'
    elif dataset_type == 'RCS':
        func_name = 'process_rcs'
    elif dataset_type == 'SNR':
        func_name = 'process_snr'
    elif dataset_type == 'RADIAL_NOISE_HS':
        func_name = 'process_radial_noise_hs'
    elif dataset_type == 'RADIAL_NOISE_IVIC':
        func_name = 'process_radial_noise_ivic'
    elif dataset_type == 'VOL_REFL':
        func_name = 'process_vol_refl'
    elif dataset_type == 'VPR':
        func_name = 'process_vpr'
    elif dataset_type == 'BIRD_DENSITY':
        func_name = 'process_bird_density'
    elif dataset_type == 'RHOHV_CORRECTION':
        func_name = 'process_correct_noise_rhohv'
    elif dataset_type == 'BIAS_CORRECTION':
        func_name = 'process_correct_bias'
    elif dataset_type == 'L':
        func_name = 'process_l'
    elif dataset_type == 'CDR':
        func_name = 'process_cdr'
    elif dataset_type == 'SAN':
        func_name = 'process_echo_id'
    elif dataset_type == 'BIRDS_ID':
        func_name = 'process_birds_id'
    elif dataset_type == 'CLT_TO_SAN':
        func_name = 'process_clt_to_echo_id'
    elif dataset_type == 'hydroMF_to_hydro':
        func_name = 'process_hydro_mf_to_hydro'
    elif dataset_type == 'hydroMF_to_SAN':
        func_name = 'process_hydro_mf_to_echo_id'
    elif dataset_type == 'ECHO_FILTER':
        func_name = 'process_echo_filter'
    elif dataset_type == 'ZDR_COLUMN':
        func_name = 'process_zdr_column'
        dsformat = 'SPARSE_GRID'
    elif dataset_type == 'SNR_FILTER':
        func_name = 'process_filter_snr'
    elif dataset_type == 'VEL_FILTER':
        func_name = 'process_filter_vel_diff'
    elif dataset_type == 'VIS_FILTER':
        func_name = 'process_filter_visibility'
    elif dataset_type == 'VIS':
        func_name = 'process_visibility'
    elif dataset_type == 'OUTLIER_FILTER':
        func_name = 'process_outlier_filter'
    elif dataset_type == 'PHIDP0_CORRECTION':
        func_name = 'process_correct_phidp0'
    elif dataset_type == 'PHIDP_SMOOTH_1W':
        func_name = 'process_smooth_phidp_single_window'
    elif dataset_type == 'PHIDP_SMOOTH_2W':
        func_name = 'process_smooth_phidp_double_window'
    elif dataset_type == 'PHIDP_KDP_VULPIANI':
        func_name = 'process_phidp_kdp_Vulpiani'
    elif dataset_type == 'PHIDP_KDP_KALMAN':
        func_name = 'process_phidp_kdp_Kalman'
    elif dataset_type == 'PHIDP_KDP_MAESAKA':
        func_name = 'process_phidp_kdp_Maesaka'
    elif dataset_type == 'PHIDP_KDP_LP':
        func_name = 'process_phidp_kdp_lp'
    elif dataset_type == 'KDP_LEASTSQUARE_1W':
        func_name = 'process_kdp_leastsquare_single_window'
    elif dataset_type == 'KDP_LEASTSQUARE_2W':
        func_name = 'process_kdp_leastsquare_double_window'
    elif dataset_type == 'ATTENUATION':
        func_name = 'process_attenuation'
    elif dataset_type == 'RAINRATE':
        func_name = 'process_rainrate'
    elif dataset_type == 'RAIN_ACCU':
        func_name = 'process_rainfall_accumulation'
        dsformat = 'TIMEAVG'
    elif dataset_type == 'TURBULENCE':
        func_name = 'process_turbulence'
    elif dataset_type == 'DEALIAS_FOURDD':
        func_name = 'process_dealias_fourdd'
    elif dataset_type == 'DEALIAS_REGION':
        func_name = 'process_dealias_region_based'
    elif dataset_type == 'DEALIAS_UNWRAP':
        func_name = 'process_dealias_unwrap_phase'
    elif dataset_type == 'RADIAL_VELOCITY':
        func_name = 'process_radial_velocity'
    elif dataset_type == 'WIND_VEL':
        func_name = 'process_wind_vel'
    elif dataset_type == 'VAD':
        func_name = 'process_vad'
    elif dataset_type == 'WINDSHEAR':
        func_name = 'process_windshear'
    elif dataset_type == 'HYDROCLASS':
        func_name = 'process_hydroclass'
    elif dataset_type == 'CENTROIDS':
        func_name = 'process_centroids'
        dsformat = 'CENTROIDS'
    elif dataset_type == 'ML_DETECTION':
        func_name = 'process_melting_layer'
        dsformat = 'ML'
    elif dataset_type == 'PHIDP0_ESTIMATE':
        func_name = 'process_estimate_phidp0'
    elif dataset_type == 'RHOHV_RAIN':
        func_name = 'process_rhohv_rain'
    elif dataset_type == 'ZDR_PREC':
        func_name = 'process_zdr_precip'
    elif dataset_type == 'ZDR_SNOW':
        func_name = 'process_zdr_snow'
    elif dataset_type == 'POL_VARIABLES':
        func_name = 'process_pol_variables'
    elif dataset_type == 'NOISE_POWER':
        func_name = 'process_noise_power'
    elif dataset_type == 'REFLECTIVITY':
        func_name = 'process_reflectivity'
    elif dataset_type == 'ZDR':
        func_name = 'process_differential_reflectivity'
    elif dataset_type == 'PhiDP':
        func_name = 'process_differential_phase'
    elif dataset_type == 'RhoHV':
        func_name = 'process_rhohv'
    elif dataset_type == 'DOPPLER_VELOCITY':
        func_name = 'process_Doppler_velocity'
    elif dataset_type == 'DOPPLER_WIDTH':
        func_name = 'process_Doppler_width'
    elif dataset_type == 'POL_VARIABLES_IQ':
        func_name = 'process_pol_variables_iq'
    elif dataset_type == 'REFLECTIVITY_IQ':
        func_name = 'process_reflectivity_iq'
    elif dataset_type == 'ZDR_IQ':
        func_name = 'process_differential_reflectivity_iq'
    elif dataset_type == 'PhiDP_IQ':
        func_name = 'process_differential_phase_iq'
    elif dataset_type == 'RhoHV_IQ':
        func_name = 'process_rhohv_iq'
    elif dataset_type == 'DOPPLER_VELOCITY_IQ':
        func_name = 'process_Doppler_velocity_iq'
    elif dataset_type == 'DOPPLER_WIDTH_IQ':
        func_name = 'process_Doppler_width_iq'
    elif dataset_type == 'MEAN_PHASE_IQ':
        func_name = 'process_mean_phase_iq'
    elif dataset_type == 'ST1_IQ':
        func_name = 'process_st1_iq'
    elif dataset_type == 'ST2_IQ':
        func_name = 'process_st2_iq'
    elif dataset_type == 'WBN_IQ':
        func_name = 'process_wbn_iq'
    elif dataset_type == 'SELFCONSISTENCY_KDP_PHIDP':
        func_name = 'process_selfconsistency_kdp_phidp'
    elif dataset_type == 'SELFCONSISTENCY_BIAS':
        func_name = 'process_selfconsistency_bias'
    elif dataset_type == 'SELFCONSISTENCY_BIAS2':
        func_name = 'process_selfconsistency_bias2'
    elif dataset_type == 'COSMO':
        func_name = 'process_cosmo'
    elif dataset_type == 'COSMO_LOOKUP':
        func_name = 'process_cosmo_lookup_table'
    elif dataset_type == 'COSMO_COORD':
        func_name = 'process_cosmo_coord'
        dsformat = 'COSMO_COORD'
    elif dataset_type == 'HZT_COORD':
        func_name = 'process_hzt_coord'
        dsformat = 'COSMO_COORD'
    elif dataset_type == 'COSMO2RADAR':
        func_name = 'process_cosmo_to_radar'
        dsformat = 'COSMO2RADAR'
    elif dataset_type == 'HZT':
        func_name = 'process_hzt'
    elif dataset_type == 'ISO0_MF':
        func_name = 'process_iso0_mf'
    elif dataset_type == 'ISO0_GRIB':
        func_name = 'process_iso0_grib'
    elif dataset_type == 'HZT_LOOKUP':
        func_name = 'process_hzt_lookup_table'
    elif dataset_type == 'DEM':
        func_name = 'process_dem'
    elif dataset_type == 'TIME_AVG':
        func_name = 'process_time_avg'
        dsformat = 'TIMEAVG'
    elif dataset_type == 'WEIGHTED_TIME_AVG':
        func_name = 'process_weighted_time_avg'
        dsformat = 'TIMEAVG'
    elif dataset_type == 'FLAG_TIME_AVG':
        func_name = 'process_time_avg_flag'
        dsformat = 'TIMEAVG'
    elif dataset_type == 'TIME_STATS':
        func_name = 'process_time_stats'
        dsformat = 'TIMEAVG'
    elif dataset_type == 'TIME_STATS2':
        func_name = 'process_time_stats2'
        dsformat = 'TIMEAVG'
    elif dataset_type == 'GRID_TIME_STATS':
        func_name = 'process_grid_time_stats'
        dsformat = 'GRID_TIMEAVG'
    elif dataset_type == 'GRID_TIME_STATS2':
        func_name = 'process_grid_time_stats2'
        dsformat = 'GRID_TIMEAVG'
    elif dataset_type == 'COLOCATED_GATES':
        func_name = 'process_colocated_gates'
        dsformat = 'COLOCATED_GATES'
    elif dataset_type == 'INTERCOMP':
        func_name = 'process_intercomp'
        dsformat = 'INTERCOMP'
    elif dataset_type == 'INTERCOMP_FIELDS':
        func_name = 'process_intercomp_fields'
        dsformat = 'INTERCOMP'
    elif dataset_type == 'INTERCOMP_TIME_AVG':
        func_name = 'process_intercomp_time_avg'
        dsformat = 'INTERCOMP'
    elif dataset_type == 'FIELDS_DIFF':
        func_name = 'process_fields_diff'
    elif dataset_type == 'MONITORING':
        func_name = 'process_monitoring'
        dsformat = 'MONITORING'
    elif dataset_type == 'GC_MONITORING':
        func_name = 'process_gc_monitoring'
        dsformat = 'MONITORING'
    elif dataset_type == 'OCCURRENCE':
        func_name = 'process_occurrence'
        dsformat = 'OCCURRENCE'
    elif dataset_type == 'TIMEAVG_STD':
        func_name = 'process_time_avg_std'
        dsformat = 'OCCURRENCE'
    elif dataset_type == 'OCCURRENCE_PERIOD':
        func_name = 'process_occurrence_period'
        dsformat = 'OCCURRENCE'
    elif dataset_type == 'SUN_HITS':
        func_name = 'process_sun_hits'
        dsformat = 'SUN_HITS'
    elif dataset_type == 'SUNSCAN':
        func_name = 'process_sunscan'
        dsformat = 'SUN_HITS'
    elif dataset_type == 'POINT_MEASUREMENT':
        func_name = 'process_point_measurement'
        dsformat = 'TIMESERIES'
    elif dataset_type == 'GRID_POINT_MEASUREMENT':
        func_name = 'process_grid_point'
        dsformat = 'TIMESERIES'
    elif dataset_type == 'ROI':
        func_name = process_roi
    elif dataset_type == 'TRAJ':
        func_name = 'process_trajectory'
        dsformat = 'TRAJ_ONLY'
    elif dataset_type == 'TRAJ_ATPLANE':
        func_name = 'process_traj_atplane'
        dsformat = 'TIMESERIES'
    elif dataset_type == 'TRAJ_ANTENNA_PATTERN':
        func_name = 'process_traj_antenna_pattern'
        dsformat = 'TIMESERIES'
    elif dataset_type == 'TRAJ_LIGHTNING':
        func_name = 'process_traj_lightning'
        dsformat = 'TIMESERIES'
    elif dataset_type == 'TRAJ_TRT':
        func_name = 'process_traj_trt'
    elif dataset_type == 'TRAJ_TRT_CONTOUR':
        func_name = 'process_traj_trt_contour'
    elif dataset_type == 'FIXED_RNG':
        func_name = process_fixed_rng
    elif dataset_type == 'FIXED_RNG_SPAN':
        func_name = process_fixed_rng_span
    else:
        raise ValueError("ERROR: Unknown dataset type '%s' of dataset '%s'"
                         % (dataset_type, dsname))

    return func_name, dsformat


def process_raw(procstatus, dscfg, radar_list=None):
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


def process_vol_to_grid(procstatus, dscfg, radar_list=None):
    """
    Function to convert polar data into a Cartesian grid

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration
    radar_list : list of Radar objects
        Optional. list of radar objects
    xmin, xmax, ymin, ymax : float
        Horizontal limits of the grid [m from origin]. Default +-20000.
    zmin, zmax : float
        vertical limits of the grid [masl]. Default 1000.
    hres, vres : float
        horizontal and vertical resolution [m]. Default 1000.
    lat0, lon0 : float
        Grid origin [deg]. The default will be the radar position
    alt0 : float
        Grid origin altitude [masl]. Default is 0
    wfunc : str
        Weighting function. Default NEAREST


    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    if radar_list is None:
        warn('ERROR: No valid radar found')
        return None, None

    # Process
    field_names_aux = []
    ind_rads_aux = []
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_names_aux.append(get_fieldname_pyart(datatype))
        ind_rads_aux.append(int(radarnr[5:8])-1)
    field_names_aux = np.array(field_names_aux)
    field_names_aux = np.unique(field_names_aux)
    ind_rads_aux = np.array(ind_rads_aux)
    ind_rads_aux = np.unique(ind_rads_aux)

    radar_list_aux = []
    ind_rads = []
    for ind_rad in ind_rads_aux:
        if radar_list[ind_rad] is None:
            warn('ERROR: Radar {} not valid'.format(ind_rad))
            continue
        radar_list_aux.append(radar_list[ind_rad])
        ind_rads.append(ind_rad)

    if not radar_list_aux:
        warn('ERROR: No valid radar found')
        return None, None

    # keep only fields present in radar object
    field_names = []
    nrad = len(radar_list_aux)
    for field_name in field_names_aux:
        nfields = 0
        for radar in radar_list_aux:
            if field_name not in radar.fields:
                warn('Field name '+field_name+' not available in radar object')
                continue
            nfields += 1
        if nfields == nrad:
            field_names.append(field_name)

    if not field_names:
        warn("Fields not available in radar data")
        return None, None

    xmin = dscfg.get('xmin', -20000.)
    xmax = dscfg.get('xmax', 20000.)
    ymin = dscfg.get('ymin', -20000.)
    ymax = dscfg.get('ymax', 20000.)
    zmin = dscfg.get('zmin', 1000.)
    zmax = dscfg.get('zmax', 1000.)
    hres = dscfg.get('hres', 1000.)
    vres = dscfg.get('vres', 500.)
    lat = dscfg.get('lat0', float(radar_list_aux[0].latitude['data']))
    lon = dscfg.get('lon0', float(radar_list_aux[0].latitude['data']))
    alt = dscfg.get('alt0', 0.)

    wfunc = dscfg.get('wfunc', 'NEAREST')

    # number of grid points
    ny = int((ymax-ymin)/hres)+1
    nx = int((xmax-xmin)/hres)+1
    nz = int((zmax-zmin)/vres)+1

    # parameters to determine the gates to use for each grid point
    if (radar_list_aux[0].instrument_parameters is not None and
            'radar_beam_width_h' in radar_list_aux[0].instrument_parameters):
        beamwidth = radar_list_aux[0].instrument_parameters[
            'radar_beam_width_h']['data'][0]
    else:
        warn('beamwidth not defined. Assumed 1 deg')
        beamwidth = 1.

    if radar_list_aux[0].ray_angle_res is not None:
        beam_spacing = radar_list_aux[0].ray_angle_res['data'][0]
    else:
        warn('beam spacing not defined. Assumed 1 deg')
        beam_spacing = 1.

    # cartesian mapping
    grid = grid_from_radars(
        radar_list_aux, gridding_algo='map_to_grid', weighting_function=wfunc,
        roi_func='dist_beam', h_factor=1.0, nb=beamwidth, bsp=beam_spacing,
        min_radius=hres/2.,
        grid_shape=(nz, ny, nx),
        grid_limits=((zmin, zmax), (ymin, ymax), (xmin, xmax)),
        grid_origin=(lat, lon), grid_origin_alt=alt,
        fields=field_names)

    new_dataset = {'radar_out': grid}

    return new_dataset, ind_rad


def process_save_radar(procstatus, dscfg, radar_list=None):
    """
    Dummy function that allows to save the entire radar object

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


def process_fixed_rng(procstatus, dscfg, radar_list=None):
    """
    Obtains radar data at a fixed range

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of strings. Dataset keyword
            The fields we want to extract
        rng : float. Dataset keyword
            The fixed range [m]
        RngTol : float. Dataset keyword
            The tolerance between the nominal range and the radar range
        ele_min, ele_max, azi_min, azi_max : floats. Dataset keyword
            The azimuth and elevation limits of the data [deg]

    radar_list : list of Radar objects
          Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the data and metadata at the point of interest
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    field_names = []
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(
            datatypedescr)
        field_names.append(get_fieldname_pyart(datatype))
    ind_rad = int(radarnr[5:8])-1

    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    # user defined parameters
    rng_tol = dscfg.get('RngTol', 50.)
    ele_min = dscfg.get('ele_min', None)
    ele_max = dscfg.get('ele_max', None)
    azi_min = dscfg.get('azi_min', None)
    azi_max = dscfg.get('azi_max', None)

    radar_aux = get_fixed_rng_data(
        radar, field_names, dscfg['rng'], rng_tol=rng_tol, ele_min=ele_min,
        ele_max=ele_max, azi_min=azi_min, azi_max=azi_max)

    if radar_aux is None:
        return None

    new_dataset = {'radar_out': radar_aux}

    return new_dataset, ind_rad


def process_fixed_rng_span(procstatus, dscfg, radar_list=None):
    """
    For each azimuth-elevation gets the data within a fixed range span
    and computes a user-defined statistic: mean, min, max, mode, median

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of strings. Dataset keyword
            The fields we want to extract
        rmin, rmax : float. Dataset keyword
            The range limits [m]
        ele_min, ele_max, azi_min, azi_max : floats. Dataset keyword
            The azimuth and elevation limits of the data [deg]

    radar_list : list of Radar objects
          Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the data and metadata at the point of interest
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    field_names = []
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(
            datatypedescr)
        field_names.append(get_fieldname_pyart(datatype))
    ind_rad = int(radarnr[5:8])-1

    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    # user defined parameters
    rmin = dscfg.get('rmin', None)
    rmax = dscfg.get('rmax', None)
    ele_min = dscfg.get('ele_min', None)
    ele_max = dscfg.get('ele_max', None)
    azi_min = dscfg.get('azi_min', None)
    azi_max = dscfg.get('azi_max', None)

    radar_aux = cut_radar(
        radar, field_names, rng_min=rmin, rng_max=rmax, ele_min=ele_min,
        ele_max=ele_max, azi_min=azi_min, azi_max=azi_max)

    if radar_aux is None:
        return None

    new_dataset = {'radar_out': radar_aux}

    return new_dataset, ind_rad


def process_roi(procstatus, dscfg, radar_list=None):
    """
    Obtains the radar data at a region of interest defined by a TRT file or
    by the user.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement
        trtfile : str. Dataset keyword
            TRT file from which to extract the region of interest
        time_tol : float. Dataset keyword
            Time tolerance between the TRT file date and the nominal radar
            volume time
        lon_roi, lat_roi : float array. Dataset keyword
            latitude and longitude positions defining a region of interest
        alt_min, alt_max : float. Dataset keyword
            Minimum and maximum altitude of the region of interest. Can be
            None
        cercle : boolean. Dataset keyword
            If True the region of interest is going to be defined as a cercle
            centered at a particular point. Default False
        lon_centre, lat_centre : Float. Dataset keyword
            The position of the centre of the cercle
        rad_cercle : Float. Dataset keyword
            The radius of the cercle in m. Default 1000.
        res_cercle : int. Dataset keyword
            Number of points used to define a quarter of cercle. Default 16
        use_latlon : Bool. Dataset keyword
            If True the coordinates used to find the radar gates within the
            ROI will be lat/lon. If false it will use Cartesian Coordinates
            with origin the radar position. Default True

    radar_list : list of Radar objects
          Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the data and metadata at the point of interest
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

    # Choose origin of ROI definition
    cercle = dscfg.get('cercle', False)
    use_latlon = dscfg.get('use_latlon', True)
    if 'trtfile' in dscfg:
        (_, yyyymmddHHMM, lon, lat, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
         _, _, _, _, _, _, _, _, _, cell_contour) = read_trt_traj_data(
             dscfg['trtfile'])

        time_tol = dscfg.get('TimeTol', 100.)
        dt = np.empty(yyyymmddHHMM.size, dtype=float)
        for i, time_traj in enumerate(yyyymmddHHMM):
            dt[i] = np.abs((dscfg['timeinfo'] - time_traj).total_seconds())
        if dt.min() > time_tol:
            warn('No TRT data for radar volume time')
            return None, None

        ind = np.argmin(dt)
        lon_roi = cell_contour[ind]['lon']
        lat_roi = cell_contour[ind]['lat']

        x_roi, y_roi = geographic_to_cartesian_aeqd(
            lon_roi, lat_roi, radar.longitude['data'][0],
            radar.longitude['data'][0])
    elif 'cercle':
        lon_centre = dscfg.get('lon_centre', None)
        lat_centre = dscfg.get('lat_centre', None)
        rad_cercle = dscfg.get('rad_cercle', 1000.)  # m
        res_cercle = dscfg.get('res_cercle', 16)

        if lon_centre is None or lat_centre is None:
            warn('Undefined ROI')
            return None, None

        # put data as x,y coordinates
        x_centre, y_centre = geographic_to_cartesian_aeqd(
            lon_centre, lat_centre, radar.longitude['data'][0],
            radar.latitude['data'][0])
        x_roi, y_roi = get_cercle_coords(
            x_centre, y_centre, radius=rad_cercle, resolution=res_cercle)
        if x_roi is None:
            return None, None
        lon_roi, lat_roi = cartesian_to_geographic_aeqd(
            x_roi, y_roi, radar.longitude['data'][0],
            radar.latitude['data'][0])
    else:
        lon_roi = dscfg.get('lon_roi', None)
        lat_roi = dscfg.get('lat_roi', None)

        if lon_roi is None or lat_roi is None:
            warn('Undefined ROI')
            return None, None
        x_roi, y_roi = geographic_to_cartesian_aeqd(
            lon_roi, lat_roi, radar.longitude['data'][0],
            radar.latitude['data'][0])

    alt_min = dscfg.get('alt_min', None)
    alt_max = dscfg.get('alt_max', None)

    roi_dict = {
        'lon': lon_roi,
        'lat': lat_roi,
        'x': x_roi,
        'y': y_roi,
        'alt_min': alt_min,
        'alt_max': alt_max}

    # extract the data within the ROI boundaries
    inds_ray, inds_rng = np.indices(np.shape(radar.gate_longitude['data']))

    if use_latlon:
        mask = np.logical_and(
            np.logical_and(
                radar.gate_latitude['data'] >= roi_dict['lat'].min(),
                radar.gate_latitude['data'] <= roi_dict['lat'].max()),
            np.logical_and(
                radar.gate_longitude['data'] >= roi_dict['lon'].min(),
                radar.gate_longitude['data'] <= roi_dict['lon'].max()))
    else:
        mask = np.logical_and(
            np.logical_and(
                radar.gate_y['data'] >= roi_dict['y'].min(),
                radar.gate_y['data'] <= roi_dict['y'].max()),
            np.logical_and(
                radar.gate_x['data'] >= roi_dict['x'].min(),
                radar.gate_x['data'] <= roi_dict['x'].max()))

    if alt_min is not None:
        mask[radar.gate_altitude['data'] < alt_min] = 0
    if alt_max is not None:
        mask[radar.gate_altitude['data'] > alt_max] = 0

    if np.all(mask == 0):
        warn('No values within ROI')
        return None, None

    inds_ray = inds_ray[mask]
    inds_rng = inds_rng[mask]

    # extract the data inside the ROI
    lat = radar.gate_latitude['data'][mask]
    lon = radar.gate_longitude['data'][mask]
    if use_latlon:
        inds, is_roi = belongs_roi_indices(
            lat, lon, roi_dict, use_latlon=use_latlon)
    else:
        y = radar.gate_y['data'][mask]
        x = radar.gate_x['data'][mask]
        inds, is_roi = belongs_roi_indices(
            y, x, roi_dict, use_latlon=use_latlon)

    if is_roi == 'None':
        warn('No values within ROI')
        return None, None

    inds_ray = inds_ray[inds]
    inds_rng = inds_rng[inds]

    lat = lat[inds][np.newaxis, :]
    lon = lon[inds][np.newaxis, :]
    alt = radar.gate_altitude['data'][inds_ray, inds_rng][np.newaxis, :]

    # prepare new radar object output
    new_dataset = {'radar_out': deepcopy(radar)}

    new_dataset['radar_out'].range['data'] = radar.range['data'][inds_rng]
    new_dataset['radar_out'].ngates = inds_rng.size
    new_dataset['radar_out'].time['data'] = np.asarray(
        [new_dataset['radar_out'].time['data'][0]])
    new_dataset['radar_out'].scan_type = 'roi'
    new_dataset['radar_out'].sweep_mode['data'] = np.array(['roi'])
    new_dataset['radar_out'].sweep_start_ray_index['data'] = np.array(
        [0], dtype='int32')
    new_dataset['radar_out'].fixed_angle['data'] = np.array(
        [], dtype='float64')
    new_dataset['radar_out'].sweep_number['data'] = np.array(
        [0], dtype='int32')
    new_dataset['radar_out'].nsweeps = 1

    if radar.rays_are_indexed is not None:
        new_dataset['radar_out'].rays_are_indexed['data'] = np.array(
            [radar.rays_are_indexed['data'][0]])
    if radar.ray_angle_res is not None:
        new_dataset['radar_out'].ray_angle_res['data'] = np.array(
            [radar.ray_angle_res['data'][0]])

    new_dataset['radar_out'].sweep_end_ray_index['data'] = np.array(
        [1], dtype='int32')
    new_dataset['radar_out'].rays_per_sweep = np.array([1], dtype='int32')
    new_dataset['radar_out'].azimuth['data'] = np.array([], dtype='float64')
    new_dataset['radar_out'].elevation['data'] = np.array([], dtype='float64')
    new_dataset['radar_out'].nrays = 1

    new_dataset['radar_out'].gate_longitude['data'] = lon
    new_dataset['radar_out'].gate_latitude['data'] = lat
    new_dataset['radar_out'].gate_altitude['data'] = alt

    new_dataset['radar_out'].gate_x['data'] = (
        radar.gate_x['data'][inds_ray, inds_rng][np.newaxis, :])
    new_dataset['radar_out'].gate_y['data'] = (
        radar.gate_y['data'][inds_ray, inds_rng][np.newaxis, :])
    new_dataset['radar_out'].gate_z['data'] = (
        radar.gate_z['data'][inds_ray, inds_rng][np.newaxis, :])

    new_dataset['radar_out'].fields = dict()
    for field_name in field_names:
        field_dict = deepcopy(radar.fields[field_name])
        field_dict['data'] = (
            radar.fields[field_name]['data'][inds_ray, inds_rng][
                np.newaxis, :])
        new_dataset['radar_out'].add_field(field_name, field_dict)

    return new_dataset, ind_rad


def process_azimuthal_average(procstatus, dscfg, radar_list=None):
    """
    Averages radar data in azimuth obtaining and RHI as a result

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement
        angle : float or None. Dataset keyword
            The center angle to average. If not set or set to -1 all
            available azimuth angles will be used
        delta_azi : float. Dataset keyword
            The angle span to average. If not set or set to -1 all the
            available azimuth angles will be used
        avg_type : str. Dataset keyword
            Average type. Can be mean or median. Default mean
        nvalid_min : int. Dataset keyword
            the minimum number of valid points to consdier the average valid.
            Default 1


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
    angle = dscfg.get('angle', None)
    delta_azi = dscfg.get('delta_azi', None)
    avg_type = dscfg.get('avg_type', 'mean')
    nvalid_min = dscfg.get('nvalid_min', 1)
    if avg_type not in ('mean', 'median'):
        warn('Unsuported statistics '+avg_type)
        return None, None

    if delta_azi == -1:
        delta_azi = None
    if angle == -1:
        angle = None

    radar_rhi = compute_azimuthal_average(
        radar, field_names, angle=angle, delta_azi=delta_azi,
        avg_type=avg_type, nvalid_min=nvalid_min)

    # prepare for exit
    new_dataset = {'radar_out': radar_rhi}

    return new_dataset, ind_rad


def process_moving_azimuthal_average(procstatus, dscfg, radar_list=None):
    """
    Applies a moving azimuthal average to the radar data

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement
        delta_azi : float. Dataset keyword
            The angle span to average. Default 20
        avg_type : str. Dataset keyword
            Average type. Can be mean or median. Default mean
        nvalid_min : int. Dataset keyword
            the minimum number of valid points to consdier the average valid.
            Default 1


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
    delta_azi = dscfg.get('delta_azi', 20.)
    avg_type = dscfg.get('avg_type', 'mean')
    nvalid_min = dscfg.get('nvalid_min', 1)
    if avg_type not in ('mean', 'median'):
        warn('Unsuported statistics '+avg_type)
        return None, None

    radar_out = deepcopy(radar)
    # At the moment only PPI scans are supported
    if radar_out.scan_type != 'ppi':
        warn('Error: unsupported scan type.')
        return None, None

    # keep only fields of interest
    radar_out.fields = dict()
    for field_name in field_names:
        radar_out.add_field(field_name, radar.fields[field_name])

    # average radar data
    for sweep in range(radar_out.nsweeps):
        radar_aux = deepcopy(radar_out)
        radar_aux = radar_aux.extract_sweeps([sweep])

        for ind_ray_sweep, angle in enumerate(radar_aux.azimuth['data']):
            ind_ray = (
                ind_ray_sweep+radar_out.sweep_start_ray_index['data'][sweep])
            # find neighbouring gates to be selected
            inds_ray, inds_rng = find_neighbour_gates(
                radar_aux, angle, None, delta_azi=delta_azi, delta_rng=None)

            # keep only data we are interested in
            for field_name in field_names:
                field_aux = radar_aux.fields[field_name]['data'][:, inds_rng]
                field_aux = field_aux[inds_ray, :]

                vals, _ = compute_directional_stats(
                    field_aux, avg_type=avg_type, nvalid_min=nvalid_min,
                    axis=0)

                radar_out.fields[field_name]['data'][ind_ray, :] = vals

    # prepare for exit
    new_dataset = {'radar_out': radar_out}

    return new_dataset, ind_rad


def process_radar_resampling(procstatus, dscfg, radar_list=None):
    """
    Resamples the radar data to mimic another radar with different geometry
    and antenna pattern

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
        par_azimuth_antenna : dict. Global keyword
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
        change_antenna_pattern : Bool. Dataset keyword
            If true the target radar has a different antenna pattern than the
            observations radar. Default True
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
        alt_tol : float. Dataset keyword
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
            Dictionary specifying whether the nans have to be used in the
            computation of the statistics for each field. Default False
        nan_value : dict. Dataset keyword
            Dictionary with the value to use to substitute the NaN values when
            computing the statistics of each field. Default 0
        moving_angle_min, moving_angle_max: float. Dataset keyword
            The minimum and maximum azimuth angle (deg) of the target radar.
            Default 0, 360.
        ray_res: float
            Ray resolution (deg). Default 1 deg.
        rng_min, rng_max:
            The minimum and maximum range of the target radar (m).
            Default 0, 100000
        rng_res : float
            The target radar range resolution (m). Default 100.
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the new radar
    ind_rad : int
        radar index
    """

    if procstatus != 1:
        return None, None

    # Process
    field_names_aux = []
    datatypes = []
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_names_aux.append(get_fieldname_pyart(datatype))
        datatypes.append(datatype)

    ind_rad = int(radarnr[5:8])-1
    if ((radar_list is None) or (radar_list[ind_rad] is None)):
        warn('ERROR: No valid radar found')
        return None, None
    radar = deepcopy(radar_list[ind_rad])

    field_names = []
    for field_name in field_names_aux:
        if field_name not in radar.fields:
            warn('Field '+field_name+' not in observations radar object')
            continue
        field_names.append(field_name)

    if not dscfg['initialized']:
        if 'antennaType' not in dscfg:
            raise Exception("ERROR: Undefined 'antennaType' for dataset '%s'"
                            % dscfg['dsname'])
        if 'configpath' not in dscfg:
            raise Exception("ERROR: Undefined 'configpath' for dataset '%s'"
                            % dscfg['dsname'])
        if 'target_radar_pos' not in dscfg:
            radar_antenna_atsameplace = True
            warn('No target radar position specified. ' +
                 'The radars are assumed co-located')
        else:
            radar_antenna_atsameplace = False

        if dscfg['antennaType'] == 'AZIMUTH':
            is_azimuth_antenna = True
            info = 'parAzAnt'
            if 'par_azimuth_antenna' not in dscfg:
                raise Exception("ERROR: Undefined 'par_azimuth_antenna' for"
                                " dataset '%s'" % dscfg['dsname'])

            patternfile = dscfg['configpath'] + 'antenna/' \
                + dscfg['par_azimuth_antenna']['elPatternFile']
            fixed_angle_val = dscfg['par_azimuth_antenna']['fixed_angle']

        elif dscfg['antennaType'] == 'ELEVATION':
            is_azimuth_antenna = False
            info = 'parElAnt'
            if 'par_elevation_antenna' not in dscfg:
                raise Exception("ERROR: Undefined 'par_elevation_antenna' for"
                                " dataset '%s'" % dscfg['dsname'])

            patternfile = dscfg['configpath'] + 'antenna/' \
                + dscfg['par_elevation_antenna']['azPatternFile']
            fixed_angle_val = dscfg['par_elevation_antenna']['fixed_angle']

        elif dscfg['antennaType'] == 'LOWBEAM':
            is_azimuth_antenna = True
            info = 'asrLowBeamAnt'
            if 'asr_lowbeam_antenna' not in dscfg:
                raise Exception("ERROR: Undefined 'asr_lowbeam_antenna' for"
                                " dataset '%s'" % dscfg['dsname'])

            patternfile = dscfg['configpath'] + 'antenna/' \
                + dscfg['asr_lowbeam_antenna']['elPatternFile']
            fixed_angle_val = dscfg['asr_lowbeam_antenna']['fixed_angle']

        elif dscfg['antennaType'] == 'HIGHBEAM':
            is_azimuth_antenna = True
            info = 'asrHighBeamAnt'
            if 'asr_highbeam_antenna' not in dscfg:
                raise Exception("ERROR: Undefined 'asr_highbeam_antenna' for"
                                " dataset '%s'" % dscfg['dsname'])

            patternfile = dscfg['configpath'] + 'antenna/' \
                + dscfg['asr_highbeam_antenna']['elPatternFile']
            patternfile_low = dscfg['configpath'] + 'antenna/' \
                + dscfg['asr_lowbeam_antenna']['elPatternFile']
            fixed_angle_val = dscfg['asr_highbeam_antenna']['fixed_angle']
        else:
            raise Exception("ERROR: Unexpected antenna type '%s' for dataset"
                            " '%s'" % (dscfg['antennaType'], dscfg['dsname']))

        if isinstance(fixed_angle_val, float):
            fixed_angle_val = [fixed_angle_val]

        change_antenna_pattern = dscfg.get('change_antenna_pattern', True)

        # Read dataset config parameters:
        weight_threshold = dscfg.get('pattern_thres', 0.)

        # Config parameters for processing when the weather radar and the
        # antenna are not at the same place:
        rhi_resolution = dscfg.get('rhi_resolution', 0.5)  # [deg]
        max_altitude = dscfg.get('max_altitude', 12000.)  # [m]
        latlon_tol = dscfg.get('latlon_tol', 0.04)  # [deg]
        alt_tol = dscfg.get('alt_tol', 1000.)  # [m]
        distance_upper_bound = dscfg.get('distance_upper_bound', 1000.)
        use_cKDTree = dscfg.get('use_cKDTree', True)
        quants = np.array(dscfg.get(
            'quants', [0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95]))

        target_radar = _create_target_radar(
            radar, dscfg, fixed_angle_val, info, field_names,
            change_antenna_pattern=change_antenna_pattern,
            quantiles=100*quants)

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

        pattern_angles = antpattern['angle'] + fixed_angle_val[0]
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

        data_is_log = dict()
        use_nans = dict()
        nan_value = dict()
        for datatype, field_name in zip(datatypes, field_names):
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

        # Persistent data structure
        trdict = dict({
            'target_radar': target_radar,
            'is_azimuth_antenna': is_azimuth_antenna,
            'info': info,
            'scan_angles': scan_angles,
            'radar_antenna_atsameplace': radar_antenna_atsameplace,
            'weightvec': weightvec,
            'quantiles': quants,
            'use_nans': use_nans,
            'nan_value': nan_value,
            'weight_threshold': weight_threshold,
            'max_altitude': max_altitude,
            'latlon_tol': latlon_tol,
            'alt_tol': alt_tol,
            'distance_upper_bound': distance_upper_bound,
            'use_cKDTree': use_cKDTree,
            'data_is_log': data_is_log,
            'change_antenna_pattern': change_antenna_pattern})

        dscfg['global_data'] = trdict
        dscfg['initialized'] = True
        # end init
    else:
        # init already done
        trdict = dscfg['global_data']

        # update time of target radar
        trdict['target_radar'].time = deepcopy(radar.time)
        time_data = np.sort(trdict['target_radar'].time['data'])
        time_res = time_data[1]-time_data[0]
        trdict['target_radar'].time['data'] = np.arange(
            trdict['target_radar'].nrays)*time_res

        # reset field values
        for field_name in trdict['target_radar'].fields.keys():
            if 'npoints' in field_name:
                trdict['target_radar'].fields[field_name]['data'] = (
                    np.ma.zeros(
                        (trdict['target_radar'].nrays,
                         trdict['target_radar'].ngates), dtype=np.int32))
                continue
            trdict['target_radar'].fields[field_name]['data'] = (
                np.ma.masked_all(
                    (trdict['target_radar'].nrays,
                     trdict['target_radar'].ngates)))

    target_radar = _get_values_antenna_pattern(radar, trdict, field_names)

    if target_radar is None:
        return None, None

    new_dataset = {'radar_out': target_radar}

    return new_dataset, ind_rad


def _get_values_antenna_pattern(radar, tadict, field_names):
    """
    Get the values of a synthetic radar

    Parameters
    ----------
    radar : radar object
        The radar volume with the data
    tadict : dict
        A dictionary containing parameters useful for radar re-sampling
    field_names : list of str
        list of names of the radar field

    Returns
    -------
    target_radar : radar object
        The synthetic radar

    """
    is_azimuth_antenna = tadict['is_azimuth_antenna']
    scan_angles = tadict['scan_angles']
    radar_antenna_atsameplace = tadict['radar_antenna_atsameplace']
    nan_value = tadict['nan_value']
    use_nans = tadict['use_nans']
    weight_threshold = tadict['weight_threshold']
    target_radar = tadict['target_radar']
    max_altitude = tadict['max_altitude']
    latlon_tol = tadict['latlon_tol']
    alt_tol = tadict['alt_tol']
    distance_upper_bound = tadict['distance_upper_bound']
    use_cKDTree = tadict['use_cKDTree']
    data_is_log = tadict['data_is_log']
    change_antenna_pattern = tadict['change_antenna_pattern']

    # find closest radar gate to target
    x_radar, y_radar, z_radar = _put_radar_in_swiss_coord(radar)
    x_target, y_target, z_target = _put_radar_in_swiss_coord(target_radar)

    tree = cKDTree(
        np.transpose(
            (x_radar.flatten(), y_radar.flatten(), z_radar.flatten())),
        compact_nodes=False, balanced_tree=False)
    _, ind_vec = tree.query(np.transpose(
        (x_target.flatten(), y_target.flatten(), z_target.flatten())), k=1)

    # temporary solution to get right time:
    target_radar.time['data'][:] = radar.time['data'][0]

    if not change_antenna_pattern:
        for field_name in field_names:
            if field_name not in radar.fields:
                warn('Field '+field_name+' not in observations radar object')
                continue

            values = radar.fields[field_name]['data'].flatten()
            target_radar.fields[field_name]['data'][:] = (
                values[ind_vec].reshape(
                    target_radar.nrays, target_radar.ngates))

        return target_radar

    # Find closest azimuth and elevation ray to target radar
    rad_ind_rays, rad_ind_rngs = np.unravel_index(
        ind_vec, (radar.nrays, radar.ngates))
    for sample, (rad_ind_ray, rad_ind_rng) in enumerate(
            zip(rad_ind_rays, rad_ind_rngs)):

        # measure time
        tstart = time()

        trad_ind_ray, trad_ind_rng = np.unravel_index(
            sample, (target_radar.nrays, target_radar.ngates))

        if radar_antenna_atsameplace:
            # ==============================================================
            # Radar and scanning antenna are at the SAME place
            # ==============================================================

            # ==============================================================
            # Get sample at bin
            if is_azimuth_antenna:
                angles = radar.azimuth['data']
                angles_scan = radar.elevation['data']
                ray_angle = radar.azimuth['data'][rad_ind_ray]
            else:
                angles = radar.elevation['data']
                angles_scan = radar.azimuth['data']
                ray_angle = radar.elevation['data'][rad_ind_ray]

            d_angle = np.abs(angles - ray_angle)
            ray_inds = np.where(d_angle < 0.09)[0]
            angles_sortind = np.argsort(angles_scan[ray_inds])

            ray_inds = ray_inds[angles_sortind]
            angles_sorted = angles_scan[ray_inds]

            # Set default values
            if ((scan_angles.size != angles_sorted.size) or
                    (np.max(np.abs(scan_angles - angles_sorted)) > 0.1)):
                warn("Scan angle mismatch!")
                continue

            w_vec = tadict['weightvec']
            for field_name in field_names:
                if field_name not in radar.fields:
                    warn("Datatype '%s' not available in radar data" %
                         field_name)
                    continue
                values = radar.fields[field_name]['data'][
                    ray_inds, rad_ind_rng]
                if use_nans[field_name]:
                    values_ma = np.ma.getmaskarray(values)
                    values[values_ma] = nan_value[field_name]

                try:
                    (avg, qvals, nvals_valid) = quantiles_weighted(
                        values,
                        weight_vector=deepcopy(w_vec),
                        quantiles=tadict['quantiles'],
                        weight_threshold=weight_threshold,
                        data_is_log=data_is_log[field_name])
                except Exception as ee:
                    warn(str(ee))
                    continue

                if avg is None:
                    continue

                # average field
                target_radar.fields['avg_'+field_name]['data'][
                    trad_ind_ray, trad_ind_rng] = avg

                # npoints field
                target_radar.fields['npoints_'+field_name]['data'][
                    trad_ind_ray, trad_ind_rng] = nvals_valid

                # quantile fields
                for quant, val in zip(tadict['quantiles'], qvals):
                    if val is None:
                        continue
                    quant_field = (
                        'quant'+'{:02d}'.format(int(100*quant))+'_' +
                        field_name)
                    target_radar.fields[quant_field]['data'][
                        trad_ind_ray, trad_ind_rng] = val

        else:
            # ================================================================
            # Radar and scanning antenna are NOT at the same place
            # ================================================================
            ray_inds, rng_inds, w_inds = _get_gates_antenna_pattern(
                radar, target_radar,
                target_radar.azimuth['data'][trad_ind_ray],
                target_radar.range['data'][trad_ind_rng],
                target_radar.time['data'][trad_ind_ray], scan_angles,
                alt_tol=alt_tol, latlon_tol=latlon_tol,
                max_altitude=max_altitude,
                distance_upper_bound=distance_upper_bound,
                use_cKDTree=use_cKDTree)

            w_vec = tadict['weightvec'][w_inds]
            for field_name in field_names:
                if field_name not in radar.fields:
                    warn("Datatype '%s' not available in radar data" %
                         field_name)
                    continue
                values = radar.fields[field_name]['data'][ray_inds, rng_inds]
                if use_nans[field_name]:
                    values_ma = np.ma.getmaskarray(values)
                    values[values_ma] = nan_value[field_name]

                try:
                    (avg, qvals, nvals_valid) = quantiles_weighted(
                        values,
                        weight_vector=deepcopy(w_vec),
                        quantiles=tadict['quantiles'],
                        weight_threshold=weight_threshold,
                        data_is_log=data_is_log[field_name])
                except Exception as ee:
                    warn(str(ee))
                    continue

                if avg is None:
                    continue

                # average field
                target_radar.fields['avg_'+field_name]['data'][
                    trad_ind_ray, trad_ind_rng] = avg

                # npoints field
                target_radar.fields['npoints_'+field_name]['data'][
                    trad_ind_ray, trad_ind_rng] = nvals_valid

                # quantile fields
                for quant, val in zip(tadict['quantiles'], qvals):
                    if val is None:
                        continue
                    quant_field = (
                        'quant'+'{:02d}'.format(int(100*quant))+'_' +
                        field_name)
                    target_radar.fields[quant_field]['data'][
                        trad_ind_ray, trad_ind_rng] = val

            tend = time()

            print(
                'original radar indices (azi, rng): '+str(rad_ind_ray)+', ' +
                str(rad_ind_rng) +
                ' target radar indices (azi, rng): '+str(trad_ind_ray)+', ' +
                str(trad_ind_rng) +
                ' Samples done: '+str(sample)+'/'+str(rad_ind_rngs.size) +
                ' Time used: '+str(tend-tstart),
                end="\r", flush=True)

    return target_radar


def _create_target_radar(radar, dscfg, fixed_angle_val, info, field_names,
                         change_antenna_pattern=False, quantiles=[50]):
    """
    Creates the target radar

    Parameters
    ----------
    radar : radar object
        the radar object containing the observed data
    dscfg : dict
        dict with the configuration
    fixed_angle_val : array of floats
        array containing the fixed angles
    info : str
        String with info on the type of antenna
    field_names : list of str
        the list of field names that the target radar will contain
    change_antenna_pattern : bool
        Whether the antenna pattern of the target radar is different from the
        observations radar
    quantiles : list of floats
        the quantiles to be computed if the target radar has a different
        antenna pattern

    Returns
    -------
    target_radar : radar object
        The target radar

    """
    # Parameters to create the new radar
    moving_angle_min = dscfg.get('moving_angle_min', 0.)
    moving_angle_max = dscfg.get('moving_angle_max', 359.)
    ray_res = dscfg.get('ray_res', 1.)
    rng_min = dscfg.get('rng_min', 0.)
    rng_max = dscfg.get('rng_max', 100000.)
    rng_res = dscfg.get('rng_res', 100.)

    # metadata needed
    _time = get_metadata('time')
    _range = get_metadata('range')
    sweep_number = get_metadata('sweep_number')
    sweep_mode = get_metadata('sweep_mode')
    fixed_angle = get_metadata('fixed_angle')
    sweep_start_ray_index = get_metadata('sweep_start_ray_index')
    sweep_end_ray_index = get_metadata('sweep_end_ray_index')
    azimuth = get_metadata('azimuth')
    elevation = get_metadata('elevation')
    metadata = dict()

    latitude = deepcopy(radar.latitude)
    longitude = deepcopy(radar.longitude)
    altitude = deepcopy(radar.altitude)
    if 'target_radar_pos' in dscfg:
        latitude['data'] = np.array(
            [dscfg['target_radar_pos']['latitude']], dtype=np.float)
        longitude['data'] = np.array(
            [dscfg['target_radar_pos']['longitude']], dtype=np.float)
        altitude['data'] = np.array(
            [dscfg['target_radar_pos']['altitude']], dtype=np.float)

    _range['data'] = np.arange(rng_min, rng_max+rng_res, rng_res)
    ngates = _range['data'].size

    fixed_angle['data'] = np.array(fixed_angle_val)
    nsweeps = fixed_angle['data'].size
    if info in ('parElAnt', 'asrLowBeamAnt', 'asrHighBeamAnt'):
        scan_type = 'ppi'
        sweep_mode['data'] = np.array(nsweeps*['azimuth_surveillance'])
    else:
        scan_type = 'rhi'
        sweep_mode['data'] = np.array(nsweeps*['elevation_surveillance'])

    moving_angle = np.arange(
        moving_angle_min, moving_angle_max+ray_res, ray_res)

    nrays = moving_angle.size*nsweeps
    sweep_start_ray_index['data'] = np.empty((nsweeps), dtype=np.int32)
    sweep_end_ray_index['data'] = np.empty((nsweeps), dtype=np.int32)
    sweep_number['data'] = np.arange(nsweeps)
    for sweep in range(nsweeps):
        sweep_start_ray_index['data'][sweep] = sweep*nrays
        sweep_end_ray_index['data'][sweep] = (sweep+1)*(nrays-1)

    elevation['data'] = np.empty((nrays), dtype=float)
    azimuth['data'] = np.empty((nrays), dtype=float)
    if scan_type == 'ppi':
        for sweep, (start_ray, end_ray) in enumerate(zip(
                sweep_start_ray_index['data'],
                sweep_end_ray_index['data'])):
            azimuth['data'][start_ray:end_ray+1] = moving_angle
            elevation['data'][start_ray:end_ray+1] = (
                fixed_angle['data'][sweep])
    else:
        for sweep, (start_ray, end_ray) in enumerate(zip(
                sweep_start_ray_index['data'],
                sweep_end_ray_index['data'])):
            elevation['data'][start_ray:end_ray+1] = moving_angle
            azimuth['data'][start_ray:end_ray+1] = fixed_angle['data'][sweep]

    _time = deepcopy(radar.time)
    time_data = np.sort(_time['data'])
    time_res = time_data[1]-time_data[0]
    _time['data'] = np.arange(nrays)*time_res

    fields = dict()
    for field_name in field_names:
        if not change_antenna_pattern:
            fields.update({field_name: get_metadata(field_name)})
            fields[field_name]['data'] = np.ma.masked_all((nrays, ngates))
            continue

        # average field
        field_name_aux = 'avg_'+field_name
        fields.update({field_name_aux: get_metadata(field_name_aux)})
        fields[field_name_aux]['data'] = np.ma.masked_all((nrays, ngates))

        # npoints field
        field_name_aux = 'npoints_'+field_name
        fields.update({field_name_aux: get_metadata(field_name_aux)})
        fields[field_name_aux]['data'] = np.ma.zeros(
            (nrays, ngates), dtype=np.int32)

        # quantile fields
        for quant in quantiles:
            field_name_aux = (
                'quant'+'{:02d}'.format(int(quant))+'_'+field_name)
            fields.update({field_name_aux: get_metadata(field_name_aux)})
            fields[field_name_aux]['data'] = np.ma.masked_all((nrays, ngates))

    target_radar = Radar(
        _time, _range, fields, metadata, scan_type, latitude, longitude,
        altitude, sweep_number, sweep_mode, fixed_angle,
        sweep_start_ray_index, sweep_end_ray_index, azimuth, elevation)

    return target_radar
