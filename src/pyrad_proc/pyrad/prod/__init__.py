"""
======================================================
Products generation (:mod:`pyrad.prod`)
======================================================

.. currentmodule:: pyrad.prod

Initiate the products generation.

Auxiliary functions
===================

.. autosummary::
    :toctree: generated/

    get_dsformat_func

Product generation
==================

.. autosummary::
    :toctree: generated/

    generate_occurrence_products
    generate_cosmo_coord_products
    generate_cosmo_to_radar_products
    generate_sun_hits_products
    generate_intercomp_products
    generate_colocated_gates_products
    generate_time_avg_products
    generate_qvp_products
    generate_vol_products
    generate_timeseries_products
    generate_monitoring_products
    generate_spectra_products
    generate_grid_products
    generate_grid_time_avg_products
    generate_traj_product
    generate_ml_products
    generate_centroids_products

"""

from .product_aux import get_prodgen_func

from .process_product import generate_sun_hits_products
from .process_product import generate_cosmo_coord_products
from .process_product import generate_qvp_products
from .process_product import generate_occurrence_products
from .process_product import generate_ml_products
from .process_product import generate_cosmo_to_radar_products
from .process_product import generate_centroids_products

from .process_vol_products import generate_vol_products
from .process_grid_products import generate_grid_products
from .process_grid_products import generate_grid_time_avg_products
from .process_spectra_products import generate_spectra_products
from .process_timeseries_products import generate_timeseries_products
from .process_traj_products import generate_traj_product
from .process_monitoring_products import generate_monitoring_products
from .process_intercomp_products import generate_intercomp_products
from .process_intercomp_products import generate_colocated_gates_products
from .process_intercomp_products import generate_time_avg_products

__all__ = [s for s in dir() if not s.startswith('_')]
