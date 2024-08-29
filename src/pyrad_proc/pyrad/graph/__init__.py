"""
==================================
Plotting (:mod:`pyrad.graph`)
==================================

.. currentmodule:: pyrad.graph

Functions to plot graphics.

Plots
=====

.. autosummary::
    :toctree: generated/

    plot_ray
    plot_surface_raw
    plot_surface
    plot_latitude_slice
    plot_longitude_slice
    plot_cross_section
    plot_dda_map
    plot_dda_slice
    plot_ppi
    plot_ppi_contour
    plot_ppi_map
    plot_rhi
    plot_rhi_contour
    plot_bscope
    plot_fixed_rng
    plot_fixed_rng_span
    plot_fixed_rng_sun
    plot_time_range
    plot_rhi_profile
    plot_along_coord
    plot_field_coverage
    plot_density
    plot_scatter
    plot_centroids
    plot_cappi
    plot_traj
    plot_pos
    plot_pos_map
    plot_quantiles
    plot_histogram
    plot_histogram2
    plot_antenna_pattern
    plot_selfconsistency
    plot_selfconsistency_instrument
    plot_timeseries
    plot_timeseries_comp
    plot_monitoring_ts
    plot_scatter_comp
    plot_intercomp_scores_ts
    plot_ml_ts
    plot_sun_hits
    plot_sun_retrieval_ts
    plot_Doppler
    plot_complex_Doppler
    plot_amp_phase_Doppler
    plot_range_Doppler
    plot_complex_range_Doppler
    plot_amp_phase_range_Doppler
    plot_angle_Doppler
    plot_complex_angle_Doppler
    plot_amp_phase_angle_Doppler
    plot_time_Doppler
    plot_complex_time_Doppler
    plot_amp_phase_time_Doppler
    plot_roi_contour
    get_colobar_label
    get_field_name
    _plot_time_range

"""

from .plots import plot_histogram, plot_histogram2, plot_density, plot_scatter  # noqa
from .plots import plot_pos, plot_pos_map, plot_sun_hits, plot_antenna_pattern  # noqa
from .plots import plot_scatter_comp, plot_quantiles, plot_selfconsistency  # noqa
from .plots import plot_selfconsistency_instrument, plot_centroids  # noqa
from .plots import _plot_time_range  # noqa

from .plots_vol import plot_ppi, plot_ppi_map, plot_rhi, plot_bscope  # noqa
from .plots_vol import plot_time_range, plot_cappi, plot_rhi_profile  # noqa
from .plots_vol import plot_along_coord, plot_field_coverage, plot_traj  # noqa
from .plots_vol import plot_rhi_contour, plot_ppi_contour, plot_ray  # noqa
from .plots_vol import plot_fixed_rng, plot_fixed_rng_span, plot_roi_contour  # noqa
from .plots_vol import plot_fixed_rng_sun  # noqa

from .plots_grid import plot_surface_raw, plot_surface, plot_latitude_slice  # noqa
from .plots_grid import plot_longitude_slice, plot_cross_section  # noqa
from .plots_grid import plot_dda_slice, plot_dda_map  # noqa

from .plots_spectra import plot_range_Doppler, plot_complex_range_Doppler  # noqa
from .plots_spectra import plot_amp_phase_range_Doppler, plot_Doppler  # noqa
from .plots_spectra import plot_complex_Doppler, plot_amp_phase_Doppler  # noqa
from .plots_spectra import plot_time_Doppler, plot_complex_time_Doppler  # noqa
from .plots_spectra import plot_amp_phase_time_Doppler, plot_angle_Doppler  # noqa
from .plots_spectra import plot_complex_angle_Doppler  # noqa
from .plots_spectra import plot_amp_phase_angle_Doppler  # noqa

from .plots_timeseries import plot_timeseries, plot_timeseries_comp  # noqa
from .plots_timeseries import plot_monitoring_ts, plot_intercomp_scores_ts  # noqa
from .plots_timeseries import plot_sun_retrieval_ts, plot_ml_ts  # noqa

from .plots_aux import get_colobar_label, get_field_name  # noqa


__all__ = [s for s in dir() if not s.startswith("_")]
