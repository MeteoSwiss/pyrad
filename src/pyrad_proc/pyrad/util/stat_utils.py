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
import ast
import operator as op
from collections import OrderedDict
from scipy.stats import energy_distance


def quantiles_weighted(
    values,
    weight_vector=None,
    quantiles=np.array([0.5]),
    weight_threshold=None,
    data_is_log=False,
    nvalid_min=3,
):
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

    if isinstance(quantiles, list):
        quantiles = np.array(quantiles)

    if weight_vector is not None:
        if weight_vector.size != values.shape[0]:
            raise Exception(
                "ERROR: Unexpected size of weight vector "
                "(%d instead of %d)" % (weight_vector.size, values.shape[0])
            )
    else:
        weight_vector = np.ones(values.shape[0], dtype=float)

    if len(values.shape) > 1:
        # repeat weight vec
        weight_vector = np.repeat(weight_vector, values.shape[1]).reshape(
            weight_vector.size, values.shape[1]
        )

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
        values = 10.0 ** (values / 10.0)

    # Average
    avg = np.ma.sum(values * weight_vector) / total_weight

    if weight_threshold is not None:
        if total_weight < weight_threshold:
            if data_is_log:
                # Convert lin to log
                avg = 10.0 * np.log10(avg)
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
        avg = 10.0 * np.log10(avg)
        quants = 10.0 * np.log10(quants)

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
            ind_values, size=ind_values.size, replace=True, p=None
        )
        samples[i] = np.ma.sum(nominator[ind_sample]) / np.ma.sum(
            denominator[ind_sample]
        )
    return samples


# Supported operators mapped to NumPy operations
operators = {
    ast.Add: op.add,
    ast.Sub: op.sub,
    ast.Mult: op.mul,
    ast.Div: op.truediv,
    ast.Pow: op.pow,
    ast.USub: op.neg,
}

# Supported NumPy mathematical functions
functions = {
    "sin": np.sin,
    "cos": np.cos,
    "tan": np.tan,
    "exp": np.exp,
    "sqrt": np.sqrt,
    "log": np.log,
    "log10": np.log10,
    "log2": np.log2,
    "abs": np.abs,
}


def parse_math_expression(expr):
    """
    Parses a mathematical expression into a vectorized NumPy-compatible function.

    Supports element-wise operations on both scalars and NumPy arrays.

    Supported operations:
        - Basic arithmetic: +, -, *, /, **
        - Unary negation: -x
        - Mathematical functions: sin, cos, tan, exp, sqrt, log, abs

    Parameters
    ----------
    expr : str
        A string representing the mathematical expression (e.g., `"sin(x) + x**2"`).

    Returns
    -------
    function
        A callable function `f(x)` that evaluates the expression using NumPy.
        Accepts scalars or NumPy arrays for vectorized operations.
    """

    def eval_(node, x):
        if isinstance(node, ast.Constant):  # Numeric literals
            return node.value
        elif isinstance(node, ast.BinOp):  # Binary operations (e.g., +, -, *, **)
            return operators[type(node.op)](eval_(node.left, x), eval_(node.right, x))
        elif isinstance(node, ast.UnaryOp):  # Unary operations (e.g., -x)
            return operators[type(node.op)](eval_(node.operand, x))
        elif isinstance(node, ast.Call):  # Function calls (e.g., sin(x))
            func = functions.get(node.func.id)
            if func:
                args = [eval_(arg, x) for arg in node.args]
                return func(*args)
            else:
                raise ValueError(f"Unsupported function: {node.func.id}")
        elif isinstance(node, ast.Name):  # Variable (e.g., x)
            if node.id == "x":
                return x
            else:
                raise ValueError(f"Unknown variable: {node.id}")
        else:
            raise TypeError(f"Unsupported operation: {node}")

    def parsed_function(x):
        """
        Evaluates the parsed expression using scalar or array inputs.

        Parameters
        ----------
        x : float, int, or numpy.ndarray
            Input value(s) to substitute for 'x' in the expression.

        Returns
        -------
        float or numpy.ndarray
            Result of evaluating the expression at input `x`.
        """
        parsed_expr = ast.parse(expr, mode="eval").body
        return eval_(parsed_expr, x)

    return parsed_function


def perfscores(
    est_data, ref_data, bounds=None, array=False, doublecond_thresh=0.1, linearize=False
):
    """
    Computes a set of precipitation performance scores, on different data ranges.
    The scores are
        - scatter: 0.5 * (Qw84(x) - Qw16(x)), where Qw is a quantile weighted
          by ref_data / sum(ref_data) and x is est_data / ref_data in dB scale
        - RMSE: root mean  square error (linear error)
        - bias:  (ME/mean(ref_data) + 1) in dB
        - ED: the energy distance which is a measure of the distance between
          two distributions (https://en.wikipedia.org/wiki/Energy_distance)

    Parameters
    ----------
    est_data : ndarray
        array of estimates (ex. precip from QPE)
    ref_data : ndarray
        array of reference (ex. precip from gauge)
    bounds : list (optional)
        list of bounds on ref_data for which to compute the error metrics,
        by default all data will be used (unbounded), note that even if you
        prescribe bounds the scores for the overall data will always be
        added in the output
    array: boolean (optional)
        Whether or not to convert the output dict to a numpy array
    doublecond_thresh: float
        By default scores are computed with a double-conditional threshold.
        Set this value to 0 to have unconditional scores.
    linearize: boolean (optional)
        Whether the input data should be linearized before computing the scores.
        Default False.
    Returns
    -------
    all_metrics : dict or ndarray
        a dictionary containing all the scores, organized in the following way
        all_metrics[bound][score]
    """
    all_metrics = OrderedDict()

    valid = np.logical_and(est_data >= 0, ref_data >= 0)
    est_data = est_data[valid > 0]
    ref_data = ref_data[valid > 0]

    est = est_data
    ref = ref_data

    all_metrics["all"] = _perfscores(est, ref)

    if bounds is not None:
        for i in range(len(bounds) - 1):
            bound_str = "{:2.1f}-{:2.1f}".format(bounds[i], bounds[i + 1])
            cond = np.logical_and(ref_data < bounds[i + 1], ref_data >= bounds[i])
            if np.sum(cond) > 0:
                est = est_data[cond]
                ref = ref_data[cond]

                all_metrics[bound_str] = _perfscores(est, ref, linearize)

    if array:
        arr = []
        for k in all_metrics:
            arr.append(list(all_metrics[k].values()))
        arr = np.array(arr)
        all_metrics = np.array(arr)

    return all_metrics


def _perfscores(est_data, ref_data, doublecond_thresh=-np.inf, linearize=False):
    """
    Robust performance scores.
    If any computation fails, the corresponding metric is set to np.nan.
    If something catastrophic happens, all metrics are np.nan.

    Notes
    -----
    - if linearize is True, the input data will first be linearized before computing the scores
    """
    # Ensure arrays
    est_data = np.asarray(est_data)
    ref_data = np.asarray(ref_data)

    # Template output with NaNs
    metrics = {
        "RMSE": np.nan,
        "scatter": np.nan,
        "logBias": np.nan,
        "ED": np.nan,
        "corr": np.nan,
        "NP": 0,
        "NP_all": int(ref_data.size),
        "est_mean": np.nan,
        "ref_mean": np.nan,
        "est_std": np.nan,
        "ref_std": np.nan,
    }

    try:
        # Optional linearization
        if linearize:
            est_data_lin = 10.0 ** (0.1 * est_data.astype(float))
            ref_data_lin = 10.0 ** (0.1 * ref_data.astype(float))
        else:
            est_data_lin = est_data.astype(float)
            ref_data_lin = ref_data.astype(float)

        # Basic sanity
        if est_data.shape != ref_data.shape:
            return metrics

        # Build condition mask
        doublecond = (ref_data > doublecond_thresh) & (est_data > doublecond_thresh)
        metrics["NP"] = int(np.count_nonzero(doublecond))

        if metrics["NP"] == 0:
            return metrics

        est_dc = est_data_lin[doublecond]
        ref_dc = ref_data_lin[doublecond]

        # RMSE
        try:
            metrics["RMSE"] = float(np.sqrt(np.nanmean((est_dc - ref_dc) ** 2)))
        except Exception:
            metrics["RMSE"] = np.nan

        # Means/stds
        try:
            metrics["est_mean"] = float(np.nanmean(est_dc))
        except Exception:
            metrics["est_mean"] = np.nan

        try:
            metrics["ref_mean"] = float(np.nanmean(ref_dc))
        except Exception:
            metrics["ref_mean"] = np.nan

        try:
            metrics["est_std"] = float(np.nanstd(est_dc))
        except Exception:
            metrics["est_std"] = np.nan

        try:
            metrics["ref_std"] = float(np.nanstd(ref_dc))
        except Exception:
            metrics["ref_std"] = np.nan

        # Pearson correlation
        try:
            good_corr = np.isfinite(est_dc) & np.isfinite(ref_dc)
            if np.count_nonzero(good_corr) >= 2:
                x = est_dc[good_corr]
                y = ref_dc[good_corr]
                if np.nanstd(x) > 0 and np.nanstd(y) > 0:
                    metrics["corr"] = float(np.corrcoef(x, y)[0, 1])
        except Exception:
            metrics["corr"] = np.nan

        # Scatter from weighted quantiles of dB error
        try:
            ratio = est_dc / ref_dc
            good_ratio = (
                np.isfinite(ratio) & (ratio > 0) & np.isfinite(ref_dc) & (ref_dc > 0)
            )
            if np.any(good_ratio):
                db_err = 10.0 * np.log10(ratio[good_ratio])
                w = ref_dc[good_ratio]
                wsum = np.nansum(w)
                if np.isfinite(wsum) and wsum > 0:
                    weights = w / wsum
                    _, q, _ = quantiles_weighted(db_err, weights, [0.16, 0.84])
                    metrics["scatter"] = float(0.5 * (q[1] - q[0]))
        except Exception:
            metrics["scatter"] = np.nan

        # logBias
        try:
            s_est = np.nansum(est_dc)
            s_ref = np.nansum(ref_dc)
            if np.isfinite(s_est) and np.isfinite(s_ref) and s_est > 0 and s_ref > 0:
                metrics["logBias"] = float(10.0 * np.log10(s_est / s_ref))
        except Exception:
            metrics["logBias"] = np.nan

        # Energy distance on original variables
        try:
            finite = np.isfinite(est_data) & np.isfinite(ref_data)
            e_fin = est_data[finite]
            r_fin = ref_data[finite]
            if e_fin.size > 0:
                metrics["ED"] = float(energy_distance(e_fin, r_fin))
        except Exception:
            metrics["ED"] = np.nan

        return metrics

    except Exception:
        return metrics
