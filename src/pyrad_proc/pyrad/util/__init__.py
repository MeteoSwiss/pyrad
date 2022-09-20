"""
==================================
Utilities (:mod:`pyrad.util`)
==================================

.. currentmodule:: pyrad.util

Functions to read and write data and configuration files.

Radar Utilities
===============

.. autosummary::
    :toctree: generated/

    get_data_along_rng
    get_data_along_azi
    get_data_along_ele
    get_ROI
    rainfall_accumulation
    time_series_statistics
    find_contiguous_times
    join_time_series
    get_range_bins_to_avg
    find_ray_index
    find_rng_index
    find_nearest_gate
    find_colocated_indexes
    get_fixed_rng_data
    time_avg_range
    get_closest_solar_flux
    create_sun_hits_field
    create_sun_retrieval_field
    compute_quantiles
    compute_quantiles_from_hist
    compute_quantiles_sweep
    compute_2d_hist
    compute_1d_stats
    compute_2d_stats
    compute_histogram
    compute_histogram_sweep
    belongs_roi_indices
    get_cercle_coords
    compute_profile_stats
    project_to_vertical

    quantiles_weighted
    ratio_bootstrapping
"""

from .radar_utils import time_avg_range, get_closest_solar_flux
from .radar_utils import create_sun_hits_field, create_sun_retrieval_field
from .radar_utils import compute_histogram, compute_histogram_sweep
from .radar_utils import compute_quantiles, compute_quantiles_sweep
from .radar_utils import compute_quantiles_from_hist, get_range_bins_to_avg
from .radar_utils import find_ray_index, find_rng_index, find_nearest_gate
from .radar_utils import find_colocated_indexes, find_contiguous_times
from .radar_utils import compute_2d_hist, compute_1d_stats, compute_2d_stats
from .radar_utils import time_series_statistics, join_time_series
from .radar_utils import rainfall_accumulation, get_ROI, belongs_roi_indices
from .radar_utils import project_to_vertical, get_data_along_rng
from .radar_utils import get_data_along_azi, get_data_along_ele
from .radar_utils import get_fixed_rng_data, get_cercle_coords
from .radar_utils import compute_profile_stats

from .stat_utils import quantiles_weighted, ratio_bootstrapping

__all__ = [s for s in dir() if not s.startswith('_')]
