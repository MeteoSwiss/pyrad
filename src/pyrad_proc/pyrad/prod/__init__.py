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




__all__ = [s for s in dir() if not s.startswith('_')]
