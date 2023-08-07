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
    process_hydro_mf_to_echo_id
    process_hydro_mf_to_hydro
    process_echo_filter
    process_cdf
    process_filter_snr
    process_filter_visibility
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
    process_vad

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

COSMO data
==========

.. autosummary::
    :toctree: generated/

    process_cosmo
    process_cosmo_lookup_table
    process_cosmo_coord
    process_hzt
    process_hzt_lookup_table
    process_hzt_coord
    process_iso0_mf
    process_iso0_grib
    process_cosmo_to_radar


DEM data
==========

.. autosummary::
    :toctree: generated/

    process_dem
    process_visibility
    process_gecsx
"""
















__all__ = [s for s in dir() if not s.startswith('_')]
