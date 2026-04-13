"""
pyrad.graph.plot_timeseries
===========================

Functions to plot Pyrad datasets

.. autosummary::
    :toctree: generated/

    plot_timeseries
    plot_timeseries_comp
    plot_monitoring_ts
    plot_intercomp_scores_ts
    plot_ml_ts
    plot_sun_retrieval_ts

"""

import pyart
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import re

import numpy as np
import datetime

import matplotlib as mpl

from ..util import warn

mpl.use("Agg")

# Increase a bit font size
mpl.rcParams.update({"font.size": 16})
mpl.rcParams.update({"font.family": "sans-serif"})


def plot_timeseries(
    tvec,
    data_list,
    fname_list,
    labelx="Time [UTC]",
    labely="Value",
    labels=None,  # <- changed default
    title="Time Series",
    period=0,
    timeformat=None,
    colors=None,
    linestyles=None,
    markers=None,
    ymin=None,
    ymax=None,
    dpi=72,
):
    """
    Plot one or multiple time series.

    Behavior
    --------
    - Legacy: if `data_list` is a list of 1D arrays (one per series), plot them.
    - Multi-timeseries: if `data_list` is a list of lists/tuples, each inner list
      represents a "group" (e.g., multiple radars) and all series in the group
      are plotted as lines. If labels is None, defaults to RADAR001, RADAR002, ...

      Examples:
        data_list = [series1, series2]                 -> 2 lines (legacy)
        data_list = [[r1, r2, r3]]                     -> 3 lines
        data_list = [[r1, r2], [r1b, r2b]]             -> 4 lines (all plotted)

    Parameters
    ----------
    tvec : array-like of datetime
        Time vector.
    data_list : list
        List of series (legacy) OR list of lists of series (multi-timeseries).
    fname_list : list of str
        Output filenames.
    labels : list of str or None
        Legend labels (one per plotted line). If None and multi-timeseries is
        detected, defaults to RADAR001, RADAR002, ... up to number of lines.
        If None in legacy mode, no legend is shown.
    period : float
        Measurement period in seconds used to compute accumulation. If 0 no
        accumulation is computed.
    timeformat : str, optional
        Datetime formatter string for x-axis.
    colors, linestyles, markers : list, optional
        Styling per plotted line (length must match number of plotted lines),
        or None to use matplotlib defaults.
    ymin, ymax : float, optional
        Y-axis limits.
    dpi : int
        Figure DPI.

    Returns
    -------
    fname_list : list of str
        Output filenames.
    """
    tvec = np.array(tvec)  # convert from pandas if needed
    if tvec.size == 0:
        raise ValueError("tvec is empty")

    # ---- normalize data_list into a flat list of series ----
    def _is_sequence_of_series(obj):
        # True if obj looks like a list/tuple of arrays (multi-timeseries group)
        if not isinstance(obj, (list, tuple)):
            return False
        if len(obj) == 0:
            return False
        return isinstance(obj[0], (list, tuple, np.ndarray))

    multi_mode = _is_sequence_of_series(data_list) and isinstance(
        data_list[0], (list, tuple)
    )
    if multi_mode:
        # Flatten groups: [[a,b],[c]] -> [a,b,c]
        series = []
        for grp in data_list:
            series.extend(list(grp))
        data_series = series
    else:
        # Legacy: [a,b,c] -> [a,b,c]
        data_series = list(data_list)

    n_lines = len(data_series)
    if n_lines == 0:
        raise ValueError("data_list contains no series")

    # ---- accumulation if requested ----
    if period > 0:
        for i, data in enumerate(data_series):
            d = np.asanyarray(data)
            d = d * (period / 3600.0)
            data_series[i] = np.ma.cumsum(d)

    # ---- default labels ----
    if labels is None and multi_mode:
        labels = [f"RADAR{idx+1:03d}" for idx in range(n_lines)]

    # ---- validate style arrays ----
    def _get_style(arr, i, default):
        if arr is None:
            return default
        if i >= len(arr):
            return default
        return arr[i]

    fig, ax = plt.subplots(figsize=[10, 6], dpi=dpi)

    any_label = False
    for i, data in enumerate(data_series):
        lab = labels[i] if labels is not None and i < len(labels) else None
        if lab is not None:
            any_label = True

        col = _get_style(colors, i, None)
        lstyle = _get_style(linestyles, i, "--")
        marker = _get_style(markers, i, "o")

        ax.plot(tvec, data, label=lab, color=col, linestyle=lstyle, marker=marker)

    ax.set_title(title)
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_ylim(bottom=ymin, top=ymax)
    ax.set_xlim([tvec[0], tvec[-1]])
    ax.grid(True)

    if timeformat is not None:
        ax.xaxis.set_major_formatter(mdates.DateFormatter(timeformat))

    fig.autofmt_xdate()
    fig.tight_layout()

    # Show legend only if we actually have labels
    if any_label:
        ax.legend(loc="best")

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def plot_timeseries_comp(
    dates,
    values,
    fname_list,
    labelx="Time [UTC]",
    labely="Value",
    labels=None,
    titl="Time Series Comparison",
    periods=None,
    ymin=None,
    ymax=None,
    dpi=72,
    colors=None,
    linestyles=None,
    markers=None,
):
    """
    Plot an arbitrary number of time series on the same axes.

    Parameters
    ----------
    dates : sequence of array-like of datetime
        One time vector per series.
    values : sequence of array-like
        One value vector per series (same length as corresponding dates entry).
    fname_list : list of str
        Output filenames.
    labelx, labely : str
        Axis labels.
    labels : sequence of str or None, optional
        Legend labels (one per series). If None, labels are "Series 1", ...
    titl : str
        Figure title.
    periods : sequence of float or None, optional
        Measurement period (seconds) per series used to compute accumulation.
        If provided and periods[i] > 0, series i is converted to accumulation:
            cumsum(values[i] * periods[i] / 3600)
        If None, no accumulation is applied.
    ymin, ymax : float, optional
        Y-axis limits.
    dpi : int
        Figure DPI.
    colors, linestyles, markers : sequence or None, optional
        Style per series. If None, matplotlib defaults are used.

    Returns
    -------
    fname_list : list of str
        Output filenames.
    """
    if len(dates) != len(values):
        raise ValueError("dates and values must have the same length")

    n = len(dates)
    if n == 0:
        raise ValueError("No time series provided")

    if labels is None:
        labels = [f"Series {i+1}" for i in range(n)]
    if periods is None:
        periods = [0] * n

    def _style(arr, i, default):
        if arr is None or i >= len(arr):
            return default
        return arr[i]

    # Apply accumulation per-series if requested
    values_proc = []
    for i, v in enumerate(values):
        v = np.asanyarray(v)
        p = periods[i] if periods is not None else 0
        if p and p > 0:
            v = v * (p / 3600.0)
            v = np.ma.cumsum(v)
        values_proc.append(v)

    fig, ax = plt.subplots(figsize=[10, 6.5], dpi=dpi)

    # Plot all series
    for i in range(n):
        col = _style(colors, i, None)
        ls = _style(linestyles, i, "--")
        mk = _style(markers, i, "o")
        ax.plot(
            dates[i],
            values_proc[i],
            label=labels[i],
            color=col,
            linestyle=ls,
            marker=mk,
        )

    ax.legend(loc="best")
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_title(titl)
    ax.grid(True)

    ax.set_ylim(bottom=ymin, top=ymax)

    # Set x-limits to the min/max over all series (robust)
    x0 = min(np.min(d) for d in dates if len(d) > 0)
    x1 = max(np.max(d) for d in dates if len(d) > 0)
    ax.set_xlim([x0, x1])

    fig.autofmt_xdate()
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def plot_monitoring_ts(
    date,
    np_t,
    val_mean_and_quant,
    quantiles,
    field_name,
    fname_list,
    ref_value=None,
    vmin=None,
    vmax=None,
    np_min=0,
    labelx="Time [UTC]",
    labely="Value",
    titl="Time Series",
    dpi=72,
    plot_until_year_end=False,
):
    """
    Plots a time series of monitoring data.

    Parameters
    ----------
    date : datetime object
        Time of the time series.
    np_t : int array
        Number of points.
    val_mean_and_quant : float array
        Values of the mean in dB, mean in dB computed from linear units,
        followed by all quantiles.
    quantiles : 1D ndarray
        The computed quantiles.
    field_name : str
        Name of the field.
    fname_list : list of str
        List of filenames where to store the plot.
    ref_value : float, optional
        The reference value.
    vmin, vmax : float, optional
        Limits of the y-axis.
    np_min : int, optional
        Minimum number of points to consider the sample plottable.
    labelx : str, optional
        Label for the X-axis.
    labely : str, optional
        Label for the Y-axis.
    titl : str, optional
        The figure title.
    dpi : int, optional
        Dots per inch.
    plot_until_year_end : bool, optional
        If True, sets the maximum x-tick to the end of the year.

    Returns
    -------
    fname_list : list of str
        List of filenames of the created plots.
    """

    # Get default Py-ART field limits
    vmin_pyart, vmax_pyart = pyart.config.get_field_limits(field_name)
    if vmin is None:
        vmin = vmin_pyart
    if vmax is None:
        vmax = vmax_pyart

    date = np.array(date)

    # Separate mean values and quantiles from val_mean_and_quant
    geometric_mean_dB = val_mean_and_quant[:, 0]
    linear_mean_dB = val_mean_and_quant[:, 1]
    quantile_values = val_mean_and_quant[:, 2:]

    # Define colors for plotting
    quantile_colors = plt.cm.viridis_r(np.linspace(0, 1, len(quantiles)))

    # Determine valid data points
    isvalid = np.logical_not(np.ma.getmaskarray(geometric_mean_dB))
    if np_min > 0:
        has_np = np_t > np_min
        isvalid = np.logical_and(isvalid, has_np)

    if not np.any(isvalid):
        warn(
            f"No valid row found in dataset (no row with num points > minimum ({np_min}))"
        )
        return

    date_plt = date[isvalid]
    mean_dB_plt = geometric_mean_dB[isvalid]
    mean_dB_from_lin_plt = linear_mean_dB[isvalid]
    quantile_values_plt = quantile_values[isvalid]

    fig = plt.figure(figsize=[15, 13], dpi=dpi)

    ax = fig.add_subplot(2, 1, 1)

    # Plot means
    ax.plot(date_plt, mean_dB_plt, "rx--", label="Geometric mean")
    ax.plot(date_plt, mean_dB_from_lin_plt, "bx--", label="Linear mean")

    # Plot all quantiles dynamically
    for i, q in enumerate(quantiles):
        ax.plot(
            date_plt,
            quantile_values_plt[:, i],
            "x-",
            color=quantile_colors[i],
            label=f"Q{q}",
        )

    # Plot reference value if provided
    if ref_value is not None:
        ax.axhline(ref_value, color="k", linestyle="--", label="Reference")

    ax.set_ylabel(labely)
    ax.set_title(titl)
    ax.set_ylim([vmin, vmax])
    ax.grid(True)
    ax.legend(loc="best")

    # Adjust x-axis limits
    if plot_until_year_end:
        t0 = date_plt[0]
        tend = datetime.datetime(year=t0.year, month=12, day=31, hour=23, minute=59)
        ax.set_xlim([t0, tend])
    else:
        ax.autoscale(enable=True, axis="x", tight=True)
        ax.set_xlim([date_plt[0], date_plt[-1]])

    # Second subplot for NP values
    ax = fig.add_subplot(2, 1, 2)
    ax.plot(date, np_t, "x-", label="Number of Samples")

    if np_min is not None:
        ax.axhline(np_min, color="k", linestyle="--", label=f"Min NP ({np_min})")

    if plot_until_year_end:
        ax.set_xlim([t0, tend])
    else:
        ax.autoscale(enable=True, axis="x", tight=True)
        ax.set_xlim([date_plt[0], date_plt[-1]])

    ax.set_ylabel("Number of Samples")
    ax.set_xlabel(labelx)
    ax.legend(loc="best")
    fig.autofmt_xdate()

    plt.tight_layout()

    # Save figure
    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def plot_sensor_scores_timeseries(
    df,
    fname_list,
    bounds,
    double_conditional_threshold,
    np_min=0,
    labelx="Time [UTC]",
    labely="Value",
    titl=None,
    dpi=72,
    plot_until_year_end=False,
    labels=None,
):
    """
    Multi-radar capable time series plot of sensor scores in compact stacked subplots.

    Behavior
    --------
    - Single-radar: metric columns are plain (e.g. "RMSE", "N", ...). One line per subplot.
    - Multi-radar: metric columns have suffixes like "_radar001", "_radar002", ... (case-insensitive).
      For each base metric (e.g. "RMSE"), all radar-specific columns are plotted as separate lines
      on the same subplot, with a legend.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain column "date" (datetime64 or convertible). Other columns are metrics.
        Multi-radar format example: "RMSE_radar001", "RMSE_radar002", ...
    fname_list : list of str
        Output filenames.
    bounds : str
        Bounds metadata (used only for title/annotation if desired).
    double_conditional_threshold : float
        Double conditional threshold metadata (used only for title/annotation if desired).
    np_min : int, optional
        Minimum number of points required to plot (based on finite "date" count).
    labelx, labely : str
        Axis labels (x label used; y label is per-subplot metric name).
    titl : str, optional
        Figure title.
    dpi : int, optional
        Figure DPI.
    plot_until_year_end : bool, optional
        If True, set xlim to end of the year (requires df["date"] non-empty).
    labels : dict or list or None, optional
        Legend labels to use for radars.
        - If dict: maps radar suffix (e.g. "radar001") -> label string.
        - If list/tuple: labels in the order of detected radar suffixes (sorted).
        - If None: legend labels default to the detected suffixes (e.g. "RADAR001").

    Returns
    -------
    fname_list : list of str
        Output filenames (same as input).
    """
    if "date" not in df.columns:
        raise ValueError("Expected a 'date' column in df.")

    # Ensure datetime
    df = df.copy()
    df["date"] = np.array(df["date"])
    # Drop rows with invalid dates
    ok_date = np.isfinite(df["date"].astype("datetime64[ns]").astype("int64"))
    df = df.loc[ok_date].sort_values("date")

    if df.empty or len(df) < max(1, np_min):
        # Nothing to plot
        return fname_list

    # --- detect multi-radar metric columns ---
    # Pattern: <metric>_<suffix>, where suffix like radar001 (case-insensitive)
    suffix_re = re.compile(r"^(?P<base>.+?)_(?P<suf>radar\d{3})$", re.IGNORECASE)

    metric_cols = [c for c in df.columns if c != "date"]

    # Group columns by base metric:
    # base -> list of (suffix, colname)
    groups = {}
    for c in metric_cols:
        m = suffix_re.match(c)
        if m:
            base = m.group("base")
            suf = m.group("suf").lower()  # radar001
        else:
            base = c
            suf = None
        groups.setdefault(base, []).append((suf, c))

    # If every group only has one column and all suffixes are None -> single-radar mode
    # But mixed is allowed; we just plot what's there.
    base_metrics = list(groups.keys())
    if not base_metrics:
        raise ValueError("No metric columns found (expected columns besides 'date').")

    # Determine radar suffix ordering for legend
    detected_suffixes = sorted(
        {suf for lst in groups.values() for (suf, _) in lst if suf is not None}
    )
    is_multiradar = len(detected_suffixes) > 0

    def _legend_label_for_suffix(suf, idx):
        # suf is like "radar001"
        if labels is None:
            return suf.upper()  # RADAR001
        if isinstance(labels, dict):
            return labels.get(suf, suf.upper())
        if isinstance(labels, (list, tuple)):
            if idx < len(labels):
                return labels[idx]
            return suf.upper()
        return suf.upper()

    # --- figure setup ---
    n = len(base_metrics)
    fig_h = max(3.0, 2.0 * n)
    fig, axes = plt.subplots(
        nrows=n, ncols=1, sharex=True, figsize=(10, fig_h), constrained_layout=True
    )
    fig.get_layout_engine().set(h_pad=2 / 72, w_pad=2 / 72, hspace=0.0, wspace=0.0)

    if n == 1:
        axes = [axes]

    # x-axis limits
    t0 = df["date"].iloc[0]
    t1 = df["date"].iloc[-1]
    if plot_until_year_end:
        # end of year based on first date
        year = t0.year
        tend = t0.__class__(year, 12, 31, 23, 59, 59)
        xlim = (t0, tend)
    else:
        xlim = (t0, t1)

    # --- plotting ---
    for ax, base in zip(axes, base_metrics):
        series = groups[base]

        # Plot per radar line if multi, otherwise single line
        if is_multiradar:
            # Keep a stable order: by detected_suffixes first, then None-suffix at end
            # Build map suf->col
            suf_to_col = {suf: col for (suf, col) in series if suf is not None}
            # plot known suffixes
            for j, suf in enumerate(detected_suffixes):
                col = suf_to_col.get(suf, None)
                if col is None:
                    continue
                ax.plot(
                    df["date"],
                    df[col],
                    marker="o",
                    linewidth=1,
                    label=_legend_label_for_suffix(suf, j),
                )
            # plot any non-suffixed column if present (rare/mixed)
            for suf, col in series:
                if suf is None:
                    ax.plot(df["date"], df[col], marker="o", linewidth=1, label=col)
        else:
            # single radar: just plot the column(s) in this base group
            # (usually exactly one)
            for _, col in series:
                ax.plot(df["date"], df[col], marker="o", linewidth=1)

        ax.set_ylabel(base, rotation=0, ha="right", va="center")
        ax.grid(True, alpha=0.3)
        ax.tick_params(axis="both", which="major", labelsize=9)
        ax.set_xlim(xlim)

        if is_multiradar:
            ax.legend(loc="best", fontsize=9)

    axes[-1].set_xlabel(labelx)

    # Title: if not provided, build something informative
    if titl is None:
        titl = f"Sensor scores (bounds={bounds}, dct={double_conditional_threshold})"
    if titl:
        fig.suptitle(titl)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def plot_intercomp_scores_ts(
    date_vec,
    np_vec,
    meanbias_vec,
    medianbias_vec,
    quant25bias_vec,
    quant75bias_vec,
    modebias_vec,
    corr_vec,
    slope_vec,
    intercep_vec,
    intercep_slope1_vec,
    fname_list,
    ref_value=0.0,
    np_min=0,
    corr_min=0.0,
    labelx="Time UTC",
    titl="RADAR001-RADAR002 intercomparison",
    dpi=72,
):
    """
    plots a time series of radar intercomparison scores

    Parameters
    ----------
    date_vec : datetime object
        time of the time series
    np_vec : int array
        number of points
    meanbias_vec, medianbias_vec, modebias_vec : float array
        mean, median and mode bias
    quant25bias_vec, quant75bias_vec: 25th and 75th percentile of the bias
    corr_vec : float array
        correlation
    slope_vec, intercep_vec : float array
        slope and intercep of a linear regression
    intercep_slope1_vec : float
        the intercep point of a inear regression of slope 1
    ref_value : float
        the reference value
    np_min : int
        The minimum number of points to consider the result valid
    corr_min : float
        The minimum correlation to consider the results valid
    labelx : str
        The label of the X axis
    titl : str
        The figure title

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    # plot only valid data (but keep first and last date)
    date2 = np.array(date_vec)
    isvalid = np.logical_not(np.ma.getmaskarray(meanbias_vec))
    isvalid_corr = np.logical_not(np.ma.getmaskarray(corr_vec))
    if np_min > 0:
        has_np = np_vec > np_min
        isvalid = np.logical_and(isvalid, has_np)
    if corr_min > 0:
        has_corr_min = corr_vec > corr_min
        isvalid = np.logical_and(isvalid, has_corr_min)

    meanbias_plt = meanbias_vec[isvalid]
    medianbias_plt = medianbias_vec[isvalid]
    quant25bias_plt = quant25bias_vec[isvalid]
    quant75bias_plt = quant75bias_vec[isvalid]
    modebias_plt = modebias_vec[isvalid]
    intercep_plt = intercep_slope1_vec[isvalid]
    corr_plt = corr_vec[isvalid_corr]
    date_corr = date2[isvalid_corr]
    date_plt = date2[isvalid]
    if not isvalid[0]:
        meanbias_plt = np.ma.append(np.ma.masked, meanbias_plt)
        medianbias_plt = np.ma.append(np.ma.masked, medianbias_plt)
        quant25bias_plt = np.ma.append(np.ma.masked, quant25bias_plt)
        quant75bias_plt = np.ma.append(np.ma.masked, quant75bias_plt)
        modebias_plt = np.ma.append(np.ma.masked, modebias_plt)
        intercep_plt = np.ma.append(np.ma.masked, intercep_plt)
        date_plt = np.ma.append(date2[0], date_plt)
    if not isvalid[-1]:
        meanbias_plt = np.ma.append(meanbias_plt, np.ma.masked)
        medianbias_plt = np.ma.append(medianbias_plt, np.ma.masked)
        quant25bias_plt = np.ma.append(quant25bias_plt, np.ma.masked)
        quant75bias_plt = np.ma.append(quant75bias_plt, np.ma.masked)
        modebias_plt = np.ma.append(modebias_plt, np.ma.masked)
        intercep_plt = np.ma.append(intercep_plt, np.ma.masked)
        date_plt = np.ma.append(date_plt, date2[-1])

    if not isvalid_corr[0]:
        corr_plt = np.ma.append(np.ma.masked, corr_plt)
        date_corr = np.ma.append(date2[0], date_corr)
    if not isvalid_corr[-1]:
        corr_plt = np.ma.append(corr_plt, np.ma.masked)
        date_corr = np.ma.append(date_corr, date2[-1])

    fig = plt.figure(figsize=[10, 20], dpi=dpi)

    ax = fig.add_subplot(4, 1, 1)
    ax.plot(date_plt, medianbias_plt, "bx-", label="median")
    ax.plot(date_plt, meanbias_plt, "rx-", label="mean")
    ax.plot(date_plt, modebias_plt, "gx-", label="mode")
    ax.plot(date_plt, intercep_plt, "yx-", label="intercep of slope 1 LR")
    if ref_value is not None:
        ax.plot(date_plt, np.zeros(len(date_plt)) + ref_value, "k--")
    # plt.legend(loc='best')
    ax.set_ylabel("bias [dB]")
    ax.set_title(titl)
    ax.set_ylim([-5.0, 5.0])

    # tight x axis
    ax.autoscale(enable=True, axis="x", tight=True)
    ax.grid(True)

    ax = fig.add_subplot(4, 1, 2)
    ax.plot(date_plt, medianbias_plt, "bx-", label="median")
    ax.plot(date_plt, quant25bias_plt, "rx-", label="25-percentile")
    ax.plot(date_plt, quant75bias_plt, "rx-", label="75-percentile")
    if ref_value is not None:
        ax.plot(date_plt, np.zeros(len(date_plt)) + ref_value, "k--")
    # plt.legend(loc='best')
    ax.set_ylabel("bias [dB]")
    ax.set_ylim([-5.0, 5.0])

    # tight x axis
    ax.autoscale(enable=True, axis="x", tight=True)
    ax.grid(True)

    ax = fig.add_subplot(4, 1, 3)
    ax.plot(date_corr, corr_plt, "bx-")

    if corr_min > 0:
        ax.plot(date_corr, np.zeros(len(date_corr)) + corr_min, "k--")

    ax.set_ylabel("correlation")
    ax.set_ylim([0.0, 1.0])

    # tight x axis
    ax.autoscale(enable=True, axis="x", tight=True)
    ax.grid(True)

    ax = fig.add_subplot(4, 1, 4)
    ax.plot(date2, np_vec, "bx-")

    if np_min > 0:
        ax.plot(date2, np.zeros(len(date2)) + np_min, "k--")

    ax.set_ylabel("Number of Samples")
    ax.set_xlabel(labelx)

    # tight x axis
    ax.autoscale(enable=True, axis="x", tight=True)

    # rotates and right aligns the x labels, and moves the bottom of the
    # axes up to make room for them
    fig.autofmt_xdate()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return fname_list


def plot_ml_ts(
    dt_ml_arr,
    ml_top_avg_arr,
    ml_top_std_arr,
    thick_avg_arr,
    thick_std_arr,
    nrays_valid_arr,
    nrays_total_arr,
    fname_list,
    labelx="Time UTC",
    titl="Melting layer time series",
    dpi=72,
):
    """
    plots a time series of melting layer data

    Parameters
    ----------
    dt_ml_arr : datetime object
        time of the time series
    np_vec : int array
        number of points
    meanbias_vec, medianbias_vec, modebias_vec : float array
        mean, median and mode bias
    quant25bias_vec, quant75bias_vec: 25th and 75th percentile of the bias
    corr_vec : float array
        correlation
    slope_vec, intercep_vec : float array
        slope and intercep of a linear regression
    intercep_slope1_vec : float
        the intercep point of a inear regression of slope 1
    ref_value : float
        the reference value
    np_min : int
        The minimum number of points to consider the result valid
    corr_min : float
        The minimum correlation to consider the results valid
    labelx : str
        The label of the X axis
    titl : str
        The figure title

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    fig = plt.figure(figsize=[10, 15], dpi=dpi)

    ax = fig.add_subplot(3, 1, 1)
    ax.plot(dt_ml_arr, ml_top_avg_arr, "bx-", label="avg")
    ax.plot(dt_ml_arr, ml_top_avg_arr + ml_top_std_arr, "rx-", label="avg+std")
    ax.plot(dt_ml_arr, ml_top_avg_arr - ml_top_std_arr, "rx-", label="avg-std")
    # plt.legend(loc='best')
    ax.set_ylabel("Top height [m MSL]")
    ax.set_title(titl)
    ax.set_ylim([0.0, 6000.0])
    ax.set_xlim([dt_ml_arr[0], dt_ml_arr[-1]])

    # tight x axis
    ax.autoscale(enable=True, axis="x", tight=True)
    ax.grid(True)

    ax = fig.add_subplot(3, 1, 2)
    ax.plot(dt_ml_arr, thick_avg_arr, "bx-", label="avg")
    ax.plot(dt_ml_arr, thick_avg_arr + thick_std_arr, "rx-", label="avg+std")
    ax.plot(dt_ml_arr, thick_avg_arr - thick_std_arr, "rx-", label="avg-std")
    # plt.legend(loc='best')
    ax.set_ylabel("Thickness [m]")
    ax.set_ylim([0.0, 3000.0])
    ax.set_xlim([dt_ml_arr[0], dt_ml_arr[-1]])

    # tight x axis
    ax.autoscale(enable=True, axis="x", tight=True)
    ax.grid(True)

    ax = fig.add_subplot(3, 1, 3)
    ax.plot(dt_ml_arr, nrays_valid_arr, "bx-", label="N valid rays")
    ax.plot(dt_ml_arr, nrays_total_arr, "rx-", label="rays total")
    # plt.legend(loc='best')
    ax.set_ylabel("Rays")
    ax.set_xlabel(labelx)
    ax.set_ylim([0, np.max(nrays_total_arr) + 5])
    ax.set_xlim([dt_ml_arr[0], dt_ml_arr[-1]])

    # tight x axis
    ax.autoscale(enable=True, axis="x", tight=True)
    ax.grid(True)

    # rotates and right aligns the x labels, and moves the bottom of the
    # axes up to make room for them
    fig.autofmt_xdate()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def plot_sun_retrieval_ts(
    sun_retrieval,
    data_type,
    fname_list,
    labelx="Date",
    titl="Sun retrieval Time Series",
    dpi=72,
):
    """
    plots sun retrieval time series series

    Parameters
    ----------
    sun_retrieval : tuple
        tuple containing the retrieved parameters
    data_type : str
        parameter to be plotted
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        the x label
    titl : str
        the title of the plot
    dpi : int
        dots per inch

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    value_std = None
    ref = None
    date = sun_retrieval[1]
    if data_type == "nhits_h":
        value = sun_retrieval[2]
        labely = "Number of sun hits H channel"
        vmin = 0
        vmax = np.max(sun_retrieval[2]) + 1
    elif data_type == "el_width_h":
        value = sun_retrieval[3]
        labely = "Elevation beamwidth H channel (Deg)"
        vmin = 0.0
        vmax = 4.0
    elif data_type == "az_width_h":
        value = sun_retrieval[4]
        labely = "Azimuth beamwidth H channel (Deg)"
        vmin = 0.0
        vmax = 4.0
    elif data_type == "el_bias_h":
        value = sun_retrieval[5]
        ref = np.zeros(len(value))
        labely = "Elevation pointing bias H channel (Deg)"
        vmin = -2.0
        vmax = 2.0
    elif data_type == "az_bias_h":
        value = sun_retrieval[6]
        ref = np.zeros(len(value))
        labely = "Azimuth pointing bias H channel (Deg)"
        vmin = -2.0
        vmax = 2.0
    elif data_type == "dBm_sun_est":
        value = sun_retrieval[7]
        value_std = sun_retrieval[8]
        labely = "Sun Power H channel (dBm)"
        vmin = -110.0
        vmax = -90.0
    elif data_type == "rx_bias_h":
        value = 10.0 * np.ma.log10(sun_retrieval[9]) - 10.0 * np.ma.log10(
            sun_retrieval[21]
        )
        value_std = sun_retrieval[8]
        ref = np.zeros(len(value))
        labely = "Receiver bias H channel (dB)"
        vmin = -5.0
        vmax = 5.0
    elif data_type == "sf_h":
        value = 10.0 * np.ma.log10(sun_retrieval[9])
        # value_std = sun_retrieval[8]
        ref = 10.0 * np.ma.log10(sun_retrieval[21])
        labely = "Observed solar flux H channel (dB(sfu))"
        vmin = 15.0
        vmax = 30.0
    elif data_type == "nhits_v":
        value = sun_retrieval[10]
        labely = "Number of sun hits V channel"
        vmin = 0
        vmax = np.max(sun_retrieval[10]) + 1
    elif data_type == "el_width_v":
        value = sun_retrieval[11]
        labely = "Elevation beamwidth V channel (Deg)"
        vmin = 0.0
        vmax = 4.0
    elif data_type == "az_width_v":
        value = sun_retrieval[12]
        labely = "Azimuth beamwidth V channel (Deg)"
        vmin = 0.0
        vmax = 4.0
    elif data_type == "el_bias_v":
        value = sun_retrieval[13]
        ref = np.zeros(len(value))
        labely = "Elevation pointing bias V channel (Deg)"
        vmin = -2.0
        vmax = 2.0
    elif data_type == "az_bias_v":
        value = sun_retrieval[14]
        ref = np.zeros(len(value))
        labely = "Azimuth pointing bias V channel (Deg)"
        vmin = -2.0
        vmax = 2.0
    elif data_type == "dBmv_sun_est":
        value = sun_retrieval[15]
        value_std = sun_retrieval[16]
        labely = "Sun Power V channel (dBm)"
        vmin = -110.0
        vmax = -90.0
    elif data_type == "rx_bias_v":
        value = 10.0 * np.ma.log10(sun_retrieval[17]) - 10.0 * np.ma.log10(
            sun_retrieval[21]
        )
        value_std = sun_retrieval[16]
        ref = np.zeros(len(value))
        labely = "Receiver bias V channel (dB)"
        vmin = -5.0
        vmax = 5.0
    elif data_type == "sf_v":
        value = 10.0 * np.ma.log10(sun_retrieval[17])
        # value_std = sun_retrieval[16]
        ref = 10.0 * np.ma.log10(sun_retrieval[21])
        labely = "Observed solar flux V channel (dB(sfu))"
        vmin = 15.0
        vmax = 30.0
    elif data_type == "nhits_zdr":
        value = sun_retrieval[18]
        labely = "Number of sun hits ZDR"
        vmin = 0
        vmax = np.max(sun_retrieval[18]) + 1
    elif data_type == "ZDR_sun_est":
        value = sun_retrieval[19]
        value_std = sun_retrieval[20]
        ref = np.zeros(len(value))
        labely = "Sun ZDR (dB)"
        vmin = -2.0
        vmax = 2.0

    mask = np.ma.getmaskarray(value)
    if mask.all():
        warn("Unable to create figure " + " ".join(fname_list) + ". No valid data")
        return None

    # plot only valid data (but keep first and last date)
    isvalid = np.logical_not(mask)
    date2 = np.array(date)

    value_plt = value[isvalid]
    date_plt = date2[isvalid]
    if not isvalid[0]:
        value_plt = np.ma.append(np.ma.masked, value_plt)
        date_plt = np.ma.append(date2[0], date_plt)
    if not isvalid[-1]:
        value_plt = np.ma.append(value_plt, np.ma.masked)
        date_plt = np.ma.append(date_plt, date2[-1])

    fig, ax = plt.subplots(figsize=[10, 6], dpi=dpi)
    ax.plot(date_plt, value_plt, "x-")
    if value_std is not None:
        value_std_plt = value_std[isvalid]
        if not isvalid[0]:
            value_std_plt = np.ma.append(np.ma.masked, value_std_plt)
        if not isvalid[-1]:
            value_std_plt = np.ma.append(value_std_plt, np.ma.masked)

        ax.plot(date_plt, value_plt + value_std_plt, "rx-")
        ax.plot(date_plt, value_plt - value_std_plt, "rx-")
    if ref is not None:
        ref_plt = ref[isvalid]
        if not isvalid[0]:
            ref_plt = np.ma.append(ref[0], ref_plt)
        if not isvalid[-1]:
            ref_plt = np.ma.append(ref_plt, ref[-1])
        ax.plot(date_plt, ref_plt, "k--")
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_title(titl)
    ax.set_ylim([vmin, vmax])
    ax.set_xlim([date_plt[0], date_plt[-1]])
    # tight x axis
    ax.autoscale(enable=True, axis="x", tight=True)
    ax.grid(True)

    # rotates and right aligns the x labels, and moves the bottom of the
    # axes up to make room for them
    fig.autofmt_xdate()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list
