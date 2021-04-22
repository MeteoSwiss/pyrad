"""
pyrad.prod.product_aux
======================

Auxiliary functions to generate products

.. autosummary::
    :toctree: generated/

    get_prodgen_func

"""
from .process_product import generate_sun_hits_products
from .process_product import generate_cosmo_coord_products
from .process_product import generate_occurrence_products
from .process_product import generate_qvp_products
from .process_product import generate_ml_products
from .process_product import generate_cosmo_to_radar_products
from .process_product import generate_centroids_products

from .process_vol_products import generate_vol_products
from .process_grid_products import generate_grid_products
from .process_grid_products import generate_sparse_grid_products
from .process_grid_products import generate_grid_time_avg_products
from .process_spectra_products import generate_spectra_products
from .process_timeseries_products import generate_timeseries_products
from .process_traj_products import generate_traj_product
from .process_monitoring_products import generate_monitoring_products
from .process_intercomp_products import generate_colocated_gates_products
from .process_intercomp_products import generate_intercomp_products
from .process_intercomp_products import generate_time_avg_products


def get_prodgen_func(dsformat, dsname, dstype):
    """
    maps the dataset format into its processing function

    Parameters
    ----------
    dsformat : str
        dataset group. The following is a list of dataset groups with the
        function that is called to generate their products. For details about
        what the functions do check the function documentation:
            'VOL': generate_vol_products
            'CENTROIDS': generate_centroids_products
            'COLOCATED_GATES': generate_colocated_gates_products
            'COSMO_COORD': generate_cosmo_coord_products
            'COSMO2RADAR': generate_cosmo_to_radar_products
            'GRID': generate_grid_products
            'SPECTRA': generate_spectra_products
            'GRID_TIMEAVG': generate_grid_time_avg_products
            'INTERCOMP': generate_intercomp_products
            'ML': generate_ml_products
            'MONITORING': generate_monitoring_products
            'OCCURRENCE': generate_occurrence_products
            'QVP': generate_qvp_products
            'SPARSE_GRID': generate_sparse_grid_products
            'SUN_HITS': generate_sun_hits_products
            'TIMEAVG': generate_time_avg_products
            'TIMESERIES': generate_timeseries_products
            'TRAJ_ONLY': generate_traj_product

    Returns
    -------
    func : function
        pyrad function used to generate the products

    """

    if dsformat == 'VOL':
        func = generate_vol_products
    elif dsformat == 'TIMESERIES':
        func = generate_timeseries_products
    elif dsformat == 'GRID':
        func = generate_grid_products
    elif dsformat == 'SPECTRA':
        func = generate_spectra_products
    elif dsformat == 'SPARSE_GRID':
        func = generate_sparse_grid_products
    elif dsformat == 'TIMEAVG':
        func = generate_time_avg_products
    elif dsformat == 'GRID_TIMEAVG':
        func = generate_grid_time_avg_products
    elif dsformat == 'SUN_HITS':
        func = generate_sun_hits_products
    elif dsformat == 'MONITORING':
        func = generate_monitoring_products
    elif dsformat == 'OCCURRENCE':
        func = generate_occurrence_products
    elif dsformat == 'INTERCOMP':
        func = generate_intercomp_products
    elif dsformat == 'TRAJ_ONLY':
        func = generate_traj_product
    elif dsformat == 'COLOCATED_GATES':
        func = generate_colocated_gates_products
    elif dsformat == 'COSMO_COORD':
        func = generate_cosmo_coord_products
    elif dsformat == 'COSMO2RADAR':
        func = generate_cosmo_to_radar_products
    elif dsformat == 'QVP':
        func = generate_qvp_products
    elif dsformat == 'ML':
        func = generate_ml_products
    elif dsformat == 'CENTROIDS':
        func = generate_centroids_products
    else:
        raise ValueError("ERROR: Unknown dataset format '%s' of dataset '%s'"
                         "(dataset type '%s')" % (dsformat, dsname, dstype))

    return func
