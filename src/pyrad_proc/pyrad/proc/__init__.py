"""
======================================================
Dataset processing (:mod:`pyrad.proc`)
======================================================

.. currentmodule:: pyrad.proc

Initiate the dataset processing.

Auxiliary functions
====================

.. autosummary::
    :toctree: generated/

    get_process_func
    process_raw
    process_save_radar
    process_fixed_rng
    process_fixed_rng_span
    process_roi
    process_keep_roi
    process_roi2
    process_azimuthal_average
    process_moving_azimuthal_average
    process_radar_resampling
    process_vol_to_grid

Gridded data functions
======================

.. autosummary::
    :toctree: generated/

    process_raw_grid
    process_grid
    process_grid_stats
    process_grid_point
    process_grid_multiple_points
    process_grid_time_stats
    process_grid_time_stats2
    process_grid_rainfall_accumulation
    process_grid_texture
    process_grid_fields_diff
    process_grid_mask
    process_normalize_luminosity
    process_pixel_filter

Spectral data functions
=======================

.. autosummary::
    :toctree: generated/

    process_raw_spectra
    process_spectra_point
    process_filter_0Doppler
    process_filter_spectra_noise
    process_filter_srhohv
    process_spectra_ang_avg
    process_spectral_power
    process_spectral_noise
    process_spectral_phase
    process_spectral_reflectivity
    process_spectral_differential_reflectivity
    process_spectral_differential_phase
    process_spectral_rhohv
    process_sunscan
    process_pol_variables
    process_noise_power
    process_reflectivity
    process_differential_reflectivity
    process_differential_phase
    process_rhohv
    process_Doppler_velocity
    process_Doppler_width
    process_ifft

IQ data functions
=======================

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

Echo classification and filtering
=================================

.. autosummary::
    :toctree: generated/

    process_echo_id
    process_birds_id
    process_clt_to_echo_id
    process_vstatus_to_echo_id
    process_hydro_mf_to_echo_id
    process_hydro_mf_to_hydro
    process_echo_filter
    process_cdf
    process_filter_snr
    process_filter_visibility
    process_gatefilter
    process_outlier_filter
    process_filter_vol2bird
    process_gate_filter_vol2bird
    process_hydroclass
    process_centroids
    process_melting_layer
    process_filter_vel_diff
    process_zdr_column

Phase processing and attenuation correction
===========================================

.. autosummary::
    :toctree: generated/

    process_correct_phidp0
    process_smooth_phidp_single_window
    process_smooth_phidp_double_window
    process_kdp_leastsquare_single_window
    process_kdp_leastsquare_double_window
    process_phidp_kdp_Vulpiani
    process_phidp_kdp_Kalman
    process_phidp_kdp_Maesaka
    process_phidp_kdp_lp
    process_attenuation

Monitoring, calibration and noise correction
============================================

.. autosummary::
    :toctree: generated/

    process_correct_bias
    process_correct_noise_rhohv
    process_rhohv_rain
    process_zdr_precip
    process_zdr_snow
    process_estimate_phidp0
    process_sun_hits
    process_selfconsistency_kdp_phidp
    process_selfconsistency_bias
    process_selfconsistency_bias2
    process_time_avg_std
    process_occurrence
    process_occurrence_period
    process_monitoring
    process_gc_monitoring
    process_time_avg
    process_weighted_time_avg
    process_time_avg_flag
    process_time_stats
    process_time_stats2
    process_colocated_gates
    process_intercomp
    process_intercomp_time_avg
    process_fields_diff
    process_intercomp_fields

Retrievals
==========

.. autosummary::
    :toctree: generated/

    process_ccor
    process_signal_power
    process_rcs
    process_rcs_pr
    process_radial_noise_hs
    process_radial_noise_ivic
    process_snr
    process_l
    process_cdr
    process_refl_from_zdr
    process_rainrate
    process_rainfall_accumulation
    process_vol_refl
    process_bird_density
    process_vpr

Doppler processing
==================

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

Time series functions
====================

.. autosummary::
    :toctree: generated/

    process_point_measurement
    process_multiple_points
    process_qvp
    process_rqvp
    process_svp
    process_evp
    process_time_height
    process_ts_along_coord

Trajectory functions
====================

.. autosummary::
    :toctree: generated/

    process_trajectory
    process_traj_atplane
    process_traj_antenna_pattern
    process_traj_lightning
    process_traj_trt
    process_traj_trt_contour

icon data
==========

.. autosummary::
    :toctree: generated/

    process_icon
    process_icon_lookup_table
    process_icon_coord
    process_hzt
    process_hzt_lookup_table
    process_hzt_coord
    process_iso0_mf
    process_iso0_grib
    process_icon_to_radar


DEM data
==========

.. autosummary::
    :toctree: generated/

    process_dem
    process_visibility
    process_gecsx
"""

from .process_aux import get_process_func, process_raw, process_save_radar  # noqa
from .process_aux import process_roi, process_roi2, process_azimuthal_average  # noqa
from .process_aux import process_keep_roi  # noqa
from .process_aux import process_fixed_rng, process_fixed_rng_span  # noqa
from .process_aux import process_radar_resampling, process_vol_to_grid  # noqa
from .process_aux import process_moving_azimuthal_average  # noqa

from .process_grid import process_grid, process_raw_grid, process_grid_point  # noqa
from .process_grid import (  # noqa
    process_grid_time_stats,
    process_grid_time_stats2,
    process_grid_stats,
)
from .process_grid import process_grid_fields_diff, process_grid_mask  # noqa
from .process_grid import process_normalize_luminosity, process_pixel_filter  # noqa
from .process_grid import process_grid_texture, process_grid_multiple_points  # noqa
from .process_grid import process_grid_rainfall_accumulation  # noqa

from .process_spectra import process_raw_spectra, process_spectral_power  # noqa
from .process_spectra import process_spectra_point, process_spectral_phase  # noqa
from .process_spectra import process_spectral_reflectivity  # noqa
from .process_spectra import process_spectral_differential_reflectivity  # noqa
from .process_spectra import process_spectral_differential_phase  # noqa
from .process_spectra import process_spectral_rhohv, process_filter_0Doppler  # noqa
from .process_spectra import process_filter_spectra_noise  # noqa
from .process_spectra import process_filter_srhohv, process_ifft  # noqa
from .process_spectra import process_pol_variables, process_reflectivity  # noqa
from .process_spectra import process_differential_reflectivity  # noqa
from .process_spectra import process_differential_phase  # noqa
from .process_spectra import process_rhohv, process_Doppler_velocity  # noqa
from .process_spectra import process_Doppler_width, process_spectra_ang_avg  # noqa
from .process_spectra import process_spectral_noise, process_noise_power  # noqa

from .process_iq import process_raw_iq, process_reflectivity_iq  # noqa
from .process_iq import process_differential_reflectivity_iq  # noqa
from .process_iq import process_rhohv_iq, process_differential_phase_iq  # noqa
from .process_iq import process_Doppler_velocity_iq, process_Doppler_width_iq  # noqa
from .process_iq import process_pol_variables_iq, process_fft  # noqa
from .process_iq import process_mean_phase_iq, process_st1_iq, process_st2_iq  # noqa
from .process_iq import process_wbn_iq  # noqa

from .process_timeseries import process_point_measurement, process_qvp  # noqa
from .process_timeseries import process_rqvp, process_evp, process_svp  # noqa
from .process_timeseries import process_time_height, process_ts_along_coord  # noqa
from .process_timeseries import process_multiple_points  # noqa

from .process_traj import process_trajectory, process_traj_atplane  # noqa
from .process_traj import process_traj_antenna_pattern, process_traj_lightning  # noqa
from .process_traj import process_traj_trt, process_traj_trt_contour  # noqa

from .process_echoclass import process_echo_id, process_birds_id  # noqa
from .process_echoclass import process_echo_filter, process_clt_to_echo_id  # noqa
from .process_echoclass import process_filter_snr, process_filter_visibility  # noqa
from .process_echoclass import process_outlier_filter, process_hydroclass  # noqa
from .process_echoclass import process_centroids, process_hydro_mf_to_hydro  # noqa
from .process_echoclass import process_cdf, process_melting_layer  # noqa
from .process_echoclass import process_filter_vel_diff, process_zdr_column  # noqa
from .process_echoclass import process_hydro_mf_to_echo_id  # noqa
from .process_echoclass import process_filter_vol2bird  # noqa
from .process_echoclass import process_gate_filter_vol2bird  # noqa
from .process_echoclass import process_gatefilter, process_vstatus_to_echo_id  # noqa

from .process_phase import process_correct_phidp0  # noqa
from .process_phase import process_smooth_phidp_single_window  # noqa
from .process_phase import process_smooth_phidp_double_window  # noqa
from .process_phase import process_kdp_leastsquare_single_window  # noqa
from .process_phase import process_kdp_leastsquare_double_window  # noqa
from .process_phase import process_phidp_kdp_lp, process_phidp_kdp_Maesaka  # noqa
from .process_phase import process_phidp_kdp_Vulpiani  # noqa
from .process_phase import process_phidp_kdp_Kalman  # noqa
from .process_phase import process_attenuation  # noqa

from .process_calib import process_correct_bias, process_correct_noise_rhohv  # noqa
from .process_calib import process_occurrence, process_occurrence_period  # noqa
from .process_calib import process_gc_monitoring, process_sun_hits  # noqa
from .process_calib import process_time_avg_std  # noqa
from .process_calib import process_sunscan  # noqa

from .process_intercomp import process_time_avg, process_weighted_time_avg  # noqa
from .process_intercomp import process_time_avg_flag, process_time_stats  # noqa
from .process_intercomp import process_colocated_gates, process_time_stats2  # noqa
from .process_intercomp import process_intercomp, process_intercomp_time_avg  # noqa
from .process_intercomp import process_fields_diff, process_intercomp_fields  # noqa

from .process_monitoring import process_zdr_snow, process_zdr_precip  # noqa
from .process_monitoring import process_estimate_phidp0, process_rhohv_rain  # noqa
from .process_monitoring import process_selfconsistency_kdp_phidp  # noqa
from .process_monitoring import process_selfconsistency_bias  # noqa
from .process_monitoring import process_selfconsistency_bias2  # noqa
from .process_monitoring import process_monitoring  # noqa

from .process_retrieve import process_signal_power, process_snr, process_ccor  # noqa
from .process_retrieve import (  # noqa
    process_l,
    process_cdr,
    process_bird_density,
    process_refl_from_zdr,
)
from .process_retrieve import process_rainrate, process_vol_refl, process_rcs  # noqa
from .process_retrieve import process_rcs_pr, process_rainfall_accumulation  # noqa
from .process_retrieve import process_radial_noise_hs, process_vpr  # noqa
from .process_retrieve import process_radial_noise_ivic  # noqa

from .process_Doppler import process_wind_vel, process_windshear  # noqa
from .process_Doppler import process_radial_velocity, process_dealias_fourdd  # noqa
from .process_Doppler import process_dealias_region_based  # noqa
from .process_Doppler import process_dealias_unwrap_phase  # noqa
from .process_Doppler import process_vad, process_turbulence, process_dda  # noqa
from .process_Doppler import process_windshear_lidar  # noqa

from .process_icon import process_icon, process_icon_lookup_table  # noqa
from .process_icon import process_icon_coord, process_hzt, process_iso0_mf  # noqa
from .process_icon import process_hzt_lookup_table, process_hzt_coord  # noqa
from .process_icon import process_icon_to_radar, process_iso0_grib  # noqa

from .process_dem import process_dem, process_visibility, process_gecsx  # noqa

__all__ = [s for s in dir() if not s.startswith("_")]
