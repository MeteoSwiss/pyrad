"""
pyrad.util.stat_utils
======================

Miscellaneous functions dealing with statistics

.. autosummary::
    :toctree: generated/

    quantiles_weighted
    ratio_bootstrapping

"""

import numpy as np


def quantiles_weighted(values, weight_vector=None, quantiles=np.array([0.5]),
                       weight_threshold=None, data_is_log=False,
                       nvalid_min=3):
    """
    Given a set of values and weights, compute the weighted quantile(s) and
    average.

    Parameters
    ----------
    values : array of floats
        Array containing the values. Can be 2-dimensional
    weight_vector : array of floats or None
        array containing the weights to apply. If None it will be an array
        of ones (uniform weight). If values is a 2D array it will be repeated
        for the second dimension
    quantiles : array of floats
        The quantiles to be computed
    weight_threshold : float or None
        If weight_threshold is set quantiles will be computed only if the
        total weight (sum of the weights of valid data) exceeds this threshold
    data_is_log : Bool
        If true the values will be considered to be in logarithmic scale and
        transformed into linear scale before computing the quantiles and
        average
    nvalid_min : int
        Minimum number of valid points to consider the computation valid

    Returns
    -------
    avg : float
        the weighted average
    quants : array of floats
        an array containing the weighted quantiles in the same order as the
        quantiles vector
    nvalid : int
        Number of valid points in the computation of the statistics

    """

    if weight_vector is not None:
        if weight_vector.size != values.shape[0]:
            raise Exception(
                "ERROR: Unexpected size of weight vector "
                "(%d instead of %d)" % (weight_vector.size, values.shape[0]))
    else:
        weight_vector = np.ones(values.shape[0], dtype=float)

    if len(values.shape) > 1:
        # repeat weight vec
        weight_vector = np.repeat(weight_vector, values.shape[1]) \
            .reshape(weight_vector.size, values.shape[1])

        values = values.reshape(-1)
        weight_vector = weight_vector.reshape(-1)

    # there must be more than 3 valid values
    mask = np.ma.getmaskarray(values)
    nvalid = np.count_nonzero(np.logical_not(mask))
    if nvalid < nvalid_min:
        return (None, np.array([None] * quantiles.size), None)

    # mask weights in non-valid data
    weight_vector[mask] = np.ma.masked

    total_weight = np.ma.sum(weight_vector)

    if data_is_log:
        # Convert log to lin
        values = 10.**(values/10.)

    # Average
    avg = np.ma.sum(values*weight_vector) / total_weight

    if weight_threshold is not None:
        if total_weight < weight_threshold:
            if data_is_log:
                # Convert lin to log
                avg = 10. * np.log10(avg)
            return (avg, np.array([None] * quantiles.size), nvalid)

    # sort the valid data
    values = values[~mask]
    weight_vector = weight_vector[~mask]

    sorter = np.argsort(values, axis=None)
    values = values[sorter]
    weight_vector = weight_vector[sorter]

    weighted_quantiles = np.cumsum(weight_vector) - 0.5 * weight_vector

    weighted_quantiles /= total_weight

    # As done by np.percentile():
    # weighted_quantiles -= weighted_quantiles[0]
    # weighted_quantiles /= weighted_quantiles[-1]

    # Note: Does not extrapolate
    quants = np.interp(quantiles, weighted_quantiles, values)

    if data_is_log:
        # Convert lin to log
        avg = 10. * np.log10(avg)
        quants = 10. * np.log10(quants)

    return (avg, quants, nvalid)


def ratio_bootstrapping(nominator, denominator, nsamples=1000):
    """
    Computes a set of samples obtained as sum(nominator)/sum(denominator)
    where the nominator and the denominator are randomly sampled with
    replacement.

    Parameters
    ----------
    nominator, denominator : 1D array
        The data points in the nominator and the denominator. Nominator
        and denominator are not independent, i.e. data point i in the
        nominator is linked to data point i in the denominator
    nsamples : int
        Number of iteration, i.e. number of samples desired

    Returns
    -------
    samples : 1D array
        the resultant samples

    """
    ind_values = np.arange(nominator.size)

    samples = np.ma.masked_all(nsamples)
    for i in range(nsamples):
        # this is for version of numpy from 1.7 to 1.15. for higher versions
        # np.random.Generator.choice should be used
        ind_sample = np.random.choice(
            ind_values, size=ind_values.size, replace=True, p=None)
        samples[i] = (np.ma.sum(nominator[ind_sample]) /
                      np.ma.sum(denominator[ind_sample]))
    return samples
