"""
pyrad.prod.process_intercomp_products
=====================================

Functions for obtaining Pyrad products from datasets used in the
intercomparison process

.. autosummary::
    :toctree: generated/

    generate_time_avg_products

"""

from .process_vol_products import generate_vol_products


def generate_time_avg_products(dataset, prdcfg):
    """
    generates time average products. Accepted product types:
        All the products of the 'VOL' dataset group

    Parameters
    ----------
    dataset : tuple
        radar objects and colocated gates dictionary

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    filename : str
        the name of the file created. None otherwise

    """
    prdcfg["timeinfo"] = dataset["timeinfo"]

    return generate_vol_products(dataset, prdcfg)
