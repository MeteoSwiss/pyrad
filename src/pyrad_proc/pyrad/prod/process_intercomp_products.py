"""
pyrad.prod.process_intercomp_products
=====================================

Functions for obtaining Pyrad products from datasets used in the
intercomparison process

.. autosummary::
    :toctree: generated/

    generate_intercomp_products
    generate_colocated_gates_products
    generate_time_avg_products

"""

from ..util import warn

import numpy as np

from .process_vol_products import generate_vol_products

from ..io.io_aux import get_fieldname_pyart
from ..io.io_aux import get_save_dir, make_filename

from ..io.read_data_other import read_intercomp_scores_ts

from ..io.write_data import write_colocated_gates, write_colocated_data
from ..io.write_data import write_intercomp_scores_ts

from ..graph.plots import plot_scatter, plot_histogram
from ..graph.plots_timeseries import plot_intercomp_scores_ts
from ..graph import plot_colocated_gates
from ..graph import plot_bias_over_range

from ..util.radar_utils import compute_2d_stats
from ..util.stat_utils import parse_math_expression


def generate_intercomp_products(dataset, prdcfg):
    """
    Generates radar intercomparison products. Accepted product types:
        'PLOT_AND_WRITE_INTERCOMP_TS': Writes statistics of radar
            intercomparison in a file and plots the time series of the
            statistics. Note that in order to use this product you need
            to also define the product WRITE_INTERCOMP for the dataset in question.
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                add_date_in_fname: Bool
                    If true adds the year in the csv file containing the
                    statistics. Default False
                sort_by_date: Bool
                    If true sorts the statistics by date when reading the
                    csv file containing the statistics. Default False
                step: float
                    The quantization step of the data. If none it will be
                    computed using the Py-ART config file. Default None
                rewrite: Bool
                    If true rewrites the csv file containing the statistics.
                    Default False
                npoints_min: int
                    The minimum number of points to consider the statistics
                    valid and therefore use the data point in the plotting.
                    Default 0
                transform: str
                    A transform to apply to the data before computing the statistics
                    and generating the plots. Any mathematical function with argument "x"
                    is accepted for example "sin(x) + x**2" or "10 * log10(x)"
                    Default is to use no transform
                range_bins : list of float, optional
                    Range bin edges in meters. Separate files and plots are created for every set of range
                    bins.
                    Default is [0, np.inf].
                corr_min: float
                    The minimum correlation to consider the statistics
                    valid and therefore use the data point in the plotting.
                    Default 0.
        'INTERCOMP_HISTOGRAMS': Plots one-dimensional histograms of colocated intercomparison data from two
            radars. The distributions of radar 1 and radar 2 are shown on the same plot,
            optionally for different range intervals.
            User defined parameters:
                histogram_type : str, optional
                    Type of histogram product. If set to "cumulative", the plot is only
                    generated when dataset["final"] is True. Default is "cumulative".
                range_bins : list of float, optional
                    Range bin edges in meters. One histogram is produced for each range
                    interval.
                    Default is [0, np.inf].
                transform : str, optional
                    Mathematical expression applied to both radar values before plotting.
                    Default is "x".
                plot_diff : bool, optional
                    If True, plots the histogram of the difference between radar 1 and radar 2
                    values instead of the histograms of the two radars. Default is False.
                step : float, optional
                    Histogram bin width. If None, bin edges are determined automatically.
                vmin, vmax : float, optional
                    Minimum and maximum histogram limits. If not specified, limits are inferred
                    from the selected data.
                cap_limits : bool, optional
                    If True, values outside [vmin, vmax] are discarded. If False, values are
                    clipped to the nearest limit. Default is False.
                double_conditional_threshold : float, optional
                    If set, only samples where both radar values are larger than this threshold
                    are used.
                binwidth_equal : bool, optional
                    If True, histogram bars are shown with equal visual width regardless of the
                    actual bin size. Default is False.
        'PLOT_SCATTER_INTERCOMP': Plots a density plot (2D histogram) with the points of
            radar 1 versus the points of radar 2. Note that in order to use this product you need
            to also define the product WRITE_INTERCOMP_TIME_AVG for the dataset in question.
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                step: float
                    The quantization step of the data. If none it will be
                    computed using the Py-ART config file. Default None
                vmin: float
                    Minimum value of the density plot. If vmin or vmax are not defined the range limits
                    of the field as defined in the Py-ART config file are going to be used.
                vmax: float
                    Maximum value of the density plot. If vmin or vmax are not defined the range limits
                    of the field as defined in the Py-ART config file are going to be used.
                min_cnt: int
                    Minimum number of samples in a 2D histogram bin for the bin to be displayed.
                    Default is 0 so all histogram bins with a last one sample are displayed.
                bins_transform: str
                    Either "linear" or "log". Specifies whether the counts in the colorbar of the scatter-plots
                    are in linear or in log scale. Default is "linear"
                cap_limits: bool
                    If set to 1, all values larger than vmax or smaller than vmin will be discarded from the
                    scatter plot. If set to 0, they will be kept and assigned to the smallest/largest bin
                    of the histogram. Default is 0.
                double_conditional_threshold: float
                    If set a double conditional threshold will be applied in the computation of the scores,
                    e.g. only samples where both radars have measured a value > threshold will be used.
                scatter_type: str
                    Type of scatter plot. Can be a plot for each radar volume
                    (instant) or at the end of the processing period
                    (cumulative). Default is cumulative
                marginals: bool
                    If set to 1, will also plot the two marginal distributions on the sides of the scatter-plot
                range_bins: list of floats
                    Range bins in meters to radar 1 to consider, a different plot will be generated for
                    all range bins. This can be used to differentiate the analysis in different
                    radar ranges.
                    Default is to generate only one figure, regardless of the distance to radar 1.
                transform: str
                    A transform to apply to the data before computing the statistics
                    and generating the plots. Any mathematical function with argument "x"
                    is accepted for example "sin(x) + x**2" or "10 * log10(x)"
                    Default is to use no transform
        'WRITE_INTERCOMP': Writes the instantaneously intercompared data (gate positions, values, etc.) in a csv file.
            Beware that if the file already exists for a given processing day, the data will be appended to the existing file
            and you might get duplicates. To avoid this delete the existing file first.
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset

    Parameters
    ----------
    dataset : tuple
        values of colocated gates dictionary

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    filename : str
        the name of the file created. None otherwise

    """

    dssavedir = prdcfg["dsname"]
    if "dssavename" in prdcfg:
        dssavedir = prdcfg["dssavename"]

    if prdcfg["type"] == "WRITE_INTERCOMP":
        if dataset["final"]:
            return None
        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=dataset["timeinfo"],
        )

        fname = make_filename(
            "colocated_data",
            prdcfg["dstype"],
            prdcfg["voltype"],
            ["csv"],
            timeinfo=dataset["timeinfo"],
            timeformat="%Y%m%d",
        )

        fname = savedir + fname[0]

        write_colocated_data(dataset["intercomp_dict"], fname)
        print("saved colocated data file: " + fname)

        return fname
    if prdcfg["type"] == "INTERCOMP_HISTOGRAMS":
        histogram_type = prdcfg.get("histogram_type", "cumulative")
        if histogram_type == "cumulative" and not dataset["final"]:
            return None

        timeformat = "%Y%m%d"
        colname = "val"
        if f"rad1_{colname}" not in dataset["intercomp_dict"]:
            return None

        field_name = get_fieldname_pyart(prdcfg["voltype"])
        transform_str = prdcfg.get("transform", "x")
        transform = parse_math_expression(transform_str)

        range_bins = prdcfg.get("range_bins", [0, np.inf])
        step = prdcfg.get("step", None)
        vmin = prdcfg.get("vmin", None)
        vmax = prdcfg.get("vmax", None)
        cap_limits = prdcfg.get("cap_limits", False)
        binwidth_equal = prdcfg.get("binwidth_equal", False)
        double_conditional_threshold = prdcfg.get("double_conditional_threshold", None)
        plot_diff = prdcfg.get("plot_diff", False)

        fname_list = []

        for i in range(len(range_bins) - 1):
            if len(range_bins) > 2:
                rangebin_info = (
                    f"range_bin_{range_bins[i]:.0f}_{range_bins[i + 1]:.0f}m"
                )
                rangebin_info_title = (
                    f"range_bin {range_bins[i]:.0f}-{range_bins[i + 1]:.0f}m"
                )
            else:
                rangebin_info = None
                rangebin_info_title = ""

            savedir = get_save_dir(
                prdcfg["basepath"],
                prdcfg["procname"],
                dssavedir,
                prdcfg["prdname"],
                subdir=rangebin_info,
                timeinfo=dataset["timeinfo"],
            )

            f_list = make_filename(
                "histogram",
                prdcfg["dstype"],
                prdcfg["voltype"],
                prdcfg["imgformat"],
                prdcfginfo=rangebin_info,
                timeinfo=dataset["timeinfo"],
                timeformat=timeformat,
            )

            for j, fname in enumerate(f_list):
                f_list[j] = savedir + fname

            rad1 = dataset["intercomp_dict"]["rad1_name"]
            rad2 = dataset["intercomp_dict"]["rad2_name"]

            titl = (
                f"colocated radar gates {rad1}-{rad2}\n"
                f"{rangebin_info_title} " + dataset["timeinfo"].strftime(timeformat)
            )
            if transform_str != "x":
                titl += f"\n f(x)={transform_str}"

            selection = np.logical_and(
                dataset["intercomp_dict"]["rad1_rng"] >= range_bins[i],
                dataset["intercomp_dict"]["rad1_rng"] < range_bins[i + 1],
            )
            if not np.sum(selection):
                warn(
                    "No data in range bin "
                    + rangebin_info_title
                    + ", skipping histogram"
                )

            rad1_values = np.ma.asarray(
                dataset["intercomp_dict"][f"rad1_{colname}"][selection]
            )
            rad2_values = np.ma.asarray(
                dataset["intercomp_dict"][f"rad2_{colname}"][selection]
            )

            if transform is not None:
                rad1_values = transform(rad1_values)
                rad2_values = transform(rad2_values)

            if double_conditional_threshold is not None:
                valid = np.logical_and(
                    rad1_values > double_conditional_threshold,
                    rad2_values > double_conditional_threshold,
                )
                rad1_values = rad1_values[valid]
                rad2_values = rad2_values[valid]

            if rad1_values.size == 0 or rad2_values.size == 0:
                continue

            if vmin is None:
                vmin_aux = min(
                    np.ma.min(rad1_values),
                    np.ma.min(rad2_values),
                )
            else:
                vmin_aux = vmin

            if vmax is None:
                vmax_aux = max(
                    np.ma.max(rad1_values),
                    np.ma.max(rad2_values),
                )
            else:
                vmax_aux = vmax

            if cap_limits:
                valid = np.logical_and.reduce(
                    (
                        rad1_values >= vmin_aux,
                        rad1_values <= vmax_aux,
                        rad2_values >= vmin_aux,
                        rad2_values <= vmax_aux,
                    )
                )
                rad1_values = rad1_values[valid]
                rad2_values = rad2_values[valid]
            else:
                rad1_values = np.ma.clip(rad1_values, vmin_aux, vmax_aux)
                rad2_values = np.ma.clip(rad2_values, vmin_aux, vmax_aux)

            if plot_diff:
                diff = rad2_values.compressed() - rad1_values.compressed()
                bin_edges = np.histogram_bin_edges(
                    diff,
                    bins="auto",
                )
            else:
                if step is None:
                    bin_edges = np.histogram_bin_edges(
                        np.ma.concatenate(
                            [rad1_values.compressed(), rad2_values.compressed()]
                        ),
                        bins="auto",
                        range=(vmin_aux, vmax_aux),
                    )
                else:
                    bin_edges = np.arange(vmin_aux, vmax_aux + step, step)

            labelx = field_name
            if transform_str != "x":
                labelx = f"f({field_name})"

            if plot_diff:
                data_histogram = rad2_values.compressed() - rad1_values.compressed()
                labels = [f"{rad2} - {rad1}"]
                labelx += " differences"
                vert_line = 0
            else:
                data_histogram = [rad1_values.compressed(), rad2_values.compressed()]
                labels = [rad1, rad2]
                vert_line = None

            plot_histogram(
                bin_edges,
                data_histogram,
                f_list,
                labelx=labelx,
                labely="Number of Samples",
                titl=titl,
                binwidth_equal=binwidth_equal,
                dpi=prdcfg.get("dpi", 72),
                labels=labels,
                vert_line=vert_line,
            )

            fname_list.extend(f_list)

            print("----- save to " + " ".join(f_list))

        return fname_list
    if prdcfg["type"] == "PLOT_SCATTER_INTERCOMP":
        scatter_type = prdcfg.get("scatter_type", "cumulative")
        if scatter_type == "cumulative" and not dataset["final"]:
            return None
        timeformat = "%Y%m%d"
        colname = "val"
        if f"rad1_{colname}" not in dataset["intercomp_dict"]:
            return

        field_name = get_fieldname_pyart(prdcfg["voltype"])
        transform_str = prdcfg.get("transform", "x")
        transform = parse_math_expression(transform_str)
        range_bins = prdcfg.get("range_bins", [0, np.inf])
        step = prdcfg.get("step", None)
        vmin = prdcfg.get("vmin", None)
        vmax = prdcfg.get("vmax", None)
        marginals = prdcfg.get("marginals", False)
        min_cnt = prdcfg.get("min_cnt", 0)
        bins_transform = prdcfg.get("bins_transform", 0)
        cap_limits = prdcfg.get("cap_limits", 0)
        double_conditional_threshold = prdcfg.get("double_conditional_threshold", None)

        fname_list = []
        for i in range(len(range_bins) - 1):  # loop on range bins
            if len(range_bins) > 2:
                rangebin_info = f"range_bin_{range_bins[i]:.0f}_{range_bins[i+1]:.0f}m"
                rangebin_info_title = (
                    f"range_bin {range_bins[i]:.0f}-{range_bins[i+1]:.0f}m"
                )
            else:
                rangebin_info = None  # only one range bin, leave empty
                rangebin_info_title = ""

            savedir = get_save_dir(
                prdcfg["basepath"],
                prdcfg["procname"],
                dssavedir,
                prdcfg["prdname"],
                subdir=rangebin_info,
                timeinfo=dataset["timeinfo"],
            )

            f_list = make_filename(
                "scatter",
                prdcfg["dstype"],
                prdcfg["voltype"],
                prdcfg["imgformat"],
                prdcfginfo=rangebin_info,
                timeinfo=dataset["timeinfo"],
                timeformat=timeformat,
            )

            for j, fname in enumerate(f_list):
                f_list[j] = savedir + fname

            rad1 = dataset["intercomp_dict"]["rad1_name"]
            rad2 = dataset["intercomp_dict"]["rad2_name"]
            titl = (
                f"colocated radar gates {rad2}-{rad1} \n {rangebin_info_title} "
                + dataset["timeinfo"].strftime(timeformat)
            )
            if transform_str != "x":
                titl += f"\n f(x)={transform_str}"

            selection = np.logical_and(
                dataset["intercomp_dict"]["rad1_rng"] >= range_bins[i],
                dataset["intercomp_dict"]["rad1_rng"] < range_bins[i + 1],
            )
            if not np.sum(selection):  # skip if selection empty
                warn(
                    "No data in range bin "
                    + rangebin_info_title
                    + ", skipping scatter plot"
                )

            hist_2d, bin_edges1, bin_edges2, stats = compute_2d_stats(
                np.ma.asarray(dataset["intercomp_dict"][f"rad1_{colname}"][selection]),
                np.ma.asarray(dataset["intercomp_dict"][f"rad2_{colname}"][selection]),
                field_name,
                field_name,
                step1=step,
                step2=step,
                transform=transform,
                vmin=vmin,
                vmax=vmax,
                cap_limits=cap_limits,
                double_conditional_threshold=double_conditional_threshold,
            )
            if hist_2d is None:
                return None

            metadata = (
                "npoints: "
                + str(stats["npoints"])
                + "\n"
                + "mode bias: "
                + "{:.2f}".format(float(stats["modebias"]))
                + "\n"
                + "median bias: "
                + "{:.2f}".format(float(stats["medianbias"]))
                + "\n"
                + "mean bias: "
                + "{:.2f}".format(float(stats["meanbias"]))
                + "\n"
                + "intercep slope 1: "
                + "{:.2f}".format(float(stats["intercep_slope_1"]))
                + "\n"
                + "corr: "
                + "{:.2f}".format(float(stats["corr"]))
                + "\n"
                + "slope: "
                + "{:.2f}".format(float(stats["slope"]))
                + "\n"
                + "intercep: "
                + "{:.2f}".format(float(stats["intercep"]))
                + "\n"
            )
            if transform_str != "x":
                field_name = f"f({field_name})"
            plot_scatter(
                bin_edges1,
                bin_edges2,
                np.ma.asarray(hist_2d),
                field_name,
                field_name,
                f_list,
                prdcfg,
                titl=titl,
                metadata=metadata,
                lin_regr=[stats["slope"], stats["intercep"]],
                lin_regr_slope1=stats["intercep_slope_1"],
                rad1_name=dataset["intercomp_dict"]["rad1_name"],
                rad2_name=dataset["intercomp_dict"]["rad2_name"],
                vmin=vmin,
                vmax=vmax,
                min_cnt=min_cnt,
                marginals=marginals,
                bins_transform=bins_transform,
            )

            fname_list.extend(f_list)

            print("----- save to " + " ".join(f_list))

        return fname_list

    if prdcfg["type"] == "BIAS_OVER_RANGE":
        """
        Plots radar intercomparison bias as a function of range.

        The product uses only data already available in:
            dataset["intercomp_dict"][f"rad1_{colname}"]
            dataset["intercomp_dict"][f"rad2_{colname}"]

        User defined parameters
        -----------------------
        voltype : str
            Pyrad variable name, used for labels and filenames.
        colname : str, optional
            Suffix of the intercomparison value columns.
            Default: "val", i.e. rad1_val and rad2_val.
        reference : str, optional
            Bias convention. Accepted values are:
                "rad2-rad1"
                "rad1-rad2"
            Default: "rad2-rad1".
        bias_type : str, optional
            Type of bias-over-range product. Accepted values are:
                "cumulative":
                    The plot is only generated when dataset["final"] is True.
                "instant":
                    The plot is generated for the current dataset.
            Default: "cumulative".
        range_min : float, optional
            Minimum range bin edge in meters.
            Used together with range_max and range_step.
            Default: 0.
        range_max : float, optional
            Maximum range bin edge in meters.
            Used together with range_min and range_step.
            The value is interpreted as the last bin edge.
            Example:
                range_min = 1300
                range_max = 40300
                range_step = 1000
            gives bins:
                1300, 2300, ..., 39300, 40300
        range_step : float, optional
            Range bin width in meters.
            Used together with range_min and range_max.
        range_bins : list of float, optional
            Explicit range bin edges in meters.
            If provided, this takes precedence over range_min, range_max
            and range_step.
            This option is kept for backward compatibility.
        percentiles : list of float, optional
            Percentiles of the point-wise offset to plot.
            The offset is computed according to the selected reference.
            Default: [25, 50, 75].
        mean_method : str, optional
            Method used to compute the mean bias. Accepted values are:
                "linear_ratio":
                    Computes mean bias as:
                        10*log10(mean(10**(a/10)) / mean(10**(b/10)))
                    This is appropriate for dBZ-like quantities.

                "offset":
                    Computes the arithmetic mean of the point-wise offsets:
                        mean(a-b)

            Default: "linear_ratio".
        min_samples_per_range : int, optional
            Minimum number of valid samples required in a range bin.
            Range bins with fewer samples are set to NaN.
            Default: 0.
        transform : str, optional
            Mathematical expression applied to both radar values before
            computing statistics.
            The expression must be a valid function of x.
            Examples:
                "x"
                "10 * log10(x)"
            Default: "x".
        double_conditional_threshold : float, optional
            If set, only samples where both radar values are larger than this
            threshold are used.
        range_marker : float, optional
            If set, plots a vertical marker at this range in meters.
        plot_npoints : bool, optional
            If True, plots the number of samples on a secondary y-axis.
            Default: True.
        write_csv : bool, optional
            If True, writes the bias-over-range statistics to a CSV file.
            Default: True.
        bias_vmin, bias_vmax : float, optional
            y-axis limits for the bias axis.
        figsize : tuple, optional
            Figure size passed to the plotting function.
            Default: (12, 7).
        dpi : int, optional
            Figure resolution.
            Default: 72.

        Returns
        -------
        fname_list : list of str
            List of generated CSV and figure filenames.
            None if no product is generated.
        """

        bias_type = prdcfg.get("bias_type", "cumulative")
        if bias_type == "cumulative" and not dataset["final"]:
            return None

        timeformat = "%Y%m%d"

        intercomp_dict = dataset["intercomp_dict"]

        colname = prdcfg.get("colname", "val")
        rad1_key = f"rad1_{colname}"
        rad2_key = f"rad2_{colname}"

        if rad1_key not in intercomp_dict:
            warn(f"Unable to generate BIAS_OVER_RANGE. Missing key: {rad1_key}")
            return None

        if rad2_key not in intercomp_dict:
            warn(f"Unable to generate BIAS_OVER_RANGE. Missing key: {rad2_key}")
            return None

        if "rad1_rng" not in intercomp_dict:
            warn("Unable to generate BIAS_OVER_RANGE. Missing key: rad1_rng")
            return None

        field_name = get_fieldname_pyart(prdcfg["voltype"])

        transform_str = prdcfg.get("transform", "x")
        transform = parse_math_expression(transform_str)

        if "range_bins" in prdcfg:
            # Backward-compatible explicit bin edges
            range_bins = np.asarray(prdcfg["range_bins"], dtype=float)
        else:
            range_min = prdcfg.get("range_min", 0.0)
            range_max = prdcfg.get("range_max", np.inf)
            range_step = prdcfg.get("range_step", None)

            if range_step is None:
                if np.isinf(range_max):
                    range_bins = np.asarray([range_min, np.inf], dtype=float)
                else:
                    raise ValueError(
                        "BIAS_OVER_RANGE requires either range_bins or "
                        "range_min, range_max and range_step"
                    )
            else:
                range_min = float(range_min)
                range_max = float(range_max)
                range_step = float(range_step)

                if range_step <= 0.0:
                    raise ValueError("range_step must be > 0")

                if range_max <= range_min:
                    raise ValueError("range_max must be > range_min")

                # range_max is interpreted as the last bin edge.
                # Example: 1300, 40300, 1000 -> 1300, 2300, ..., 40300
                range_bins = np.arange(
                    range_min,
                    range_max + 0.5 * range_step,
                    range_step,
                    dtype=float,
                )

                # Avoid overshooting too much in case of floating-point steps
                range_bins = range_bins[range_bins <= range_max + 1e-6]

                if range_bins[-1] < range_max:
                    range_bins = np.append(range_bins, range_max)

        if range_bins.size < 2:
            raise ValueError("At least two range bin edges are required")

        range_centers = range_bins[:-1] + np.diff(range_bins) / 2.0

        reference = prdcfg.get("reference", "rad2-rad1")
        if reference == "rad2-rad1":
            reference_label = (
                f"{intercomp_dict['rad2_name']} - {intercomp_dict['rad1_name']}"
            )
            first_key = rad2_key
            second_key = rad1_key
        elif reference == "rad1-rad2":
            reference_label = (
                f"{intercomp_dict['rad1_name']} - {intercomp_dict['rad2_name']}"
            )
            first_key = rad1_key
            second_key = rad2_key
        else:
            raise ValueError("reference must be either 'rad2-rad1' or 'rad1-rad2'")

        percentiles = np.asarray(
            prdcfg.get("percentiles", [25.0, 50.0, 75.0]),
            dtype=float,
        )

        if np.any(percentiles < 0.0) or np.any(percentiles > 100.0):
            raise ValueError("percentiles must be between 0 and 100")

        mean_method = prdcfg.get("mean_method", "linear_ratio")
        if mean_method not in ("linear_ratio", "offset"):
            raise ValueError("mean_method must be 'linear_ratio' or 'offset'")

        min_samples_per_range = prdcfg.get("min_samples_per_range", 0)
        double_conditional_threshold = prdcfg.get("double_conditional_threshold", None)

        rad1_rng = np.ma.asarray(intercomp_dict["rad1_rng"])

        first_values_all = np.ma.asarray(intercomp_dict[first_key])
        second_values_all = np.ma.asarray(intercomp_dict[second_key])

        if transform is not None:
            first_values_all = transform(first_values_all)
            second_values_all = transform(second_values_all)

        first_values_all = np.ma.masked_invalid(first_values_all)
        second_values_all = np.ma.masked_invalid(second_values_all)

        meanbias_vec = []
        nsamples_vec = []
        percentile_vecs = {perc: [] for perc in percentiles}

        for i in range(len(range_bins) - 1):
            selection = np.logical_and(
                rad1_rng >= range_bins[i],
                rad1_rng < range_bins[i + 1],
            )

            first_values = first_values_all[selection]
            second_values = second_values_all[selection]

            valid = np.logical_not(
                np.logical_or(
                    np.ma.getmaskarray(first_values),
                    np.ma.getmaskarray(second_values),
                )
            )

            first_values = np.asarray(first_values.filled(np.nan), dtype=float)
            second_values = np.asarray(second_values.filled(np.nan), dtype=float)

            valid = np.logical_and.reduce(
                (
                    valid,
                    np.isfinite(first_values),
                    np.isfinite(second_values),
                )
            )

            if double_conditional_threshold is not None:
                valid = np.logical_and.reduce(
                    (
                        valid,
                        first_values > double_conditional_threshold,
                        second_values > double_conditional_threshold,
                    )
                )

            first_values = first_values[valid]
            second_values = second_values[valid]

            npoints = first_values.size
            nsamples_vec.append(npoints)

            if npoints < min_samples_per_range:
                meanbias_vec.append(np.nan)
                for perc in percentiles:
                    percentile_vecs[perc].append(np.nan)
                continue

            offset = first_values - second_values

            if mean_method == "linear_ratio":
                mean_first = np.nanmean(10.0 ** (first_values / 10.0))
                mean_second = np.nanmean(10.0 ** (second_values / 10.0))

                if (
                    np.isfinite(mean_first)
                    and np.isfinite(mean_second)
                    and mean_first > 0.0
                    and mean_second > 0.0
                ):
                    meanbias = 10.0 * np.log10(mean_first / mean_second)
                else:
                    meanbias = np.nan
            else:
                meanbias = np.nanmean(offset)

            meanbias_vec.append(meanbias)

            for perc in percentiles:
                percentile_vecs[perc].append(np.nanpercentile(offset, perc))

        meanbias_vec = np.asarray(meanbias_vec)
        nsamples_vec = np.asarray(nsamples_vec)

        for perc in percentiles:
            percentile_vecs[perc] = np.asarray(percentile_vecs[perc])

        rad1_name = intercomp_dict["rad1_name"]
        rad2_name = intercomp_dict["rad2_name"]

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=dataset["timeinfo"],
        )

        prdcfginfo = f"{rad1_name}-{rad2_name}_{reference}"

        # ------------------------------------------------------------------
        # Write CSV
        # ------------------------------------------------------------------
        fname_list = []

        if prdcfg.get("write_csv", True):
            csvfname = make_filename(
                "bias_over_range",
                prdcfg["dstype"],
                prdcfg["voltype"],
                ["csv"],
                prdcfginfo=prdcfginfo,
                timeinfo=dataset["timeinfo"],
                timeformat=timeformat,
            )[0]

            csvfname = savedir + csvfname

            csv_header = [
                "range_start_m",
                "range_end_m",
                "range_center_m",
                "npoints",
                "meanbias",
            ]

            csv_data = [
                range_bins[:-1],
                range_bins[1:],
                range_centers,
                nsamples_vec,
                meanbias_vec,
            ]

            for perc in percentiles:
                csv_header.append(f"percentile_{perc:g}")
                csv_data.append(percentile_vecs[perc])

            csv_data = np.column_stack(csv_data)

            np.savetxt(
                csvfname,
                csv_data,
                delimiter=",",
                header=",".join(csv_header),
                comments="",
            )

            print("saved CSV file: " + csvfname)
            fname_list.append(csvfname)

        # ------------------------------------------------------------------
        # Plot
        # ------------------------------------------------------------------
        figfname_list = make_filename(
            "bias_over_range",
            prdcfg["dstype"],
            prdcfg["voltype"],
            prdcfg.get("imgformat", ["png"]),
            prdcfginfo=prdcfginfo,
            timeinfo=dataset["timeinfo"],
            timeformat=timeformat,
        )

        for j, figfname in enumerate(figfname_list):
            figfname_list[j] = savedir + figfname

        plot_field_name = field_name
        if transform_str != "x":
            plot_field_name = f"f({field_name})"

        titl = (
            f"{rad1_name}-{rad2_name} {plot_field_name} bias over range\n"
            f"{reference_label}, "
            f"{dataset['timeinfo'].strftime(timeformat)}"
        )

        if transform_str != "x":
            titl += f"\n f(x)={transform_str}"

        percentile_dict = {perc: percentile_vecs[perc] for perc in percentiles}

        plot_bias_over_range(
            range_centers,
            meanbias_vec,
            percentile_dict,
            figfname_list,
            nsamples=nsamples_vec,
            reference_label=reference_label,
            field_name=plot_field_name,
            titl=titl,
            mean_label=f"mean bias ({mean_method})",
            range_marker=prdcfg.get("range_marker", None),
            plot_npoints=prdcfg.get("plot_npoints", True),
            bias_vmin=prdcfg.get("bias_vmin", None),
            bias_vmax=prdcfg.get("bias_vmax", None),
            figsize=prdcfg.get("figsize", (12, 7)),
            dpi=prdcfg.get("dpi", 72),
        )

        print("----- save to " + " ".join(figfname_list))
        fname_list.extend(figfname_list)

        return fname_list

    if prdcfg["type"] == "PLOT_AND_WRITE_INTERCOMP_TS":
        if not dataset["final"]:
            return None

        field_name = get_fieldname_pyart(prdcfg["voltype"])
        step = prdcfg.get("step", None)
        transform_str = prdcfg.get("transform", "x")
        transform = parse_math_expression(transform_str)

        range_bins = prdcfg.get("range_bins", [0, np.inf])
        vmin = prdcfg.get("vmin", None)
        vmax = prdcfg.get("vmax", None)
        cap_limits = prdcfg.get("cap_limits", False)
        double_conditional_threshold = prdcfg.get("double_conditional_threshold", None)

        rad1_name = dataset["intercomp_dict"]["rad1_name"]
        rad2_name = dataset["intercomp_dict"]["rad2_name"]

        # put time info in file path and name
        csvtimeinfo_file = None
        timeformat = None
        sort_by_date = prdcfg.get("sort_by_date", False)
        rewrite = prdcfg.get("rewrite", False)

        if prdcfg.get("add_date_in_fname", False):
            csvtimeinfo_file = dataset["timeinfo"]
            timeformat = "%Y"

        fname_list = []

        for i in range(len(range_bins) - 1):
            if len(range_bins) > 2:
                rangebin_info = (
                    f"range_bin_{range_bins[i]:.0f}_{range_bins[i + 1]:.0f}m"
                )
                rangebin_info_title = (
                    f" range {range_bins[i]:.0f}-{range_bins[i + 1]:.0f} m"
                )
            else:
                rangebin_info = None
                rangebin_info_title = ""

            savedir = get_save_dir(
                prdcfg["basepath"],
                prdcfg["procname"],
                dssavedir,
                prdcfg["prdname"],
                subdir=rangebin_info,
                timeinfo=None,
            )

            selection = np.logical_and(
                dataset["intercomp_dict"]["rad1_rng"] >= range_bins[i],
                dataset["intercomp_dict"]["rad1_rng"] < range_bins[i + 1],
            )

            if not np.sum(selection):
                warn(
                    "No data in range bin "
                    + rangebin_info_title
                    + ", skipping statistics for timeseries plot"
                )
                continue

            hist_2d, bin_edges1, bin_edges2, stats = compute_2d_stats(
                np.ma.asarray(dataset["intercomp_dict"]["rad1_val"][selection]),
                np.ma.asarray(dataset["intercomp_dict"]["rad2_val"][selection]),
                field_name,
                field_name,
                step1=step,
                step2=step,
                transform=transform,
                vmin=vmin,
                vmax=vmax,
                cap_limits=cap_limits,
                double_conditional_threshold=double_conditional_threshold,
            )

            if hist_2d is None:
                continue

            prdcfginfo = rad1_name + "-" + rad2_name
            if rangebin_info is not None:
                prdcfginfo += "_" + rangebin_info

            csvfname = make_filename(
                "ts",
                prdcfg["dstype"],
                prdcfg["voltype"],
                ["csv"],
                prdcfginfo=prdcfginfo,
                timeinfo=csvtimeinfo_file,
                timeformat=timeformat,
            )[0]

            csvfname = savedir + csvfname
            fname_list.append(csvfname)
            write_intercomp_scores_ts(
                dataset["timeinfo"],
                stats,
                field_name,
                csvfname,
                rad1_name=rad1_name,
                rad2_name=rad2_name,
            )

            print("saved CSV file: " + csvfname)

            (
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
            ) = read_intercomp_scores_ts(csvfname, sort_by_date=sort_by_date)

            if date_vec is None:
                warn("Unable to plot time series. No valid data")
                continue

            if len(date_vec) < 2:
                warn("Unable to plot time series. Not enough points")
                continue

            if rewrite:
                stats_rewrite = {
                    "npoints": np_vec,
                    "meanbias": meanbias_vec,
                    "medianbias": medianbias_vec,
                    "quant25bias": quant25bias_vec,
                    "quant75bias": quant75bias_vec,
                    "modebias": modebias_vec,
                    "corr": corr_vec,
                    "slope": slope_vec,
                    "intercep": intercep_vec,
                    "intercep_slope_1": intercep_slope1_vec,
                }

                write_intercomp_scores_ts(
                    date_vec,
                    stats_rewrite,
                    field_name,
                    csvfname,
                    rad1_name=rad1_name,
                    rad2_name=rad2_name,
                )

            figtimeinfo = None
            fig_timeformat = None
            titldate = (
                date_vec.tolist()[0].strftime("%Y%m%d")
                + "-"
                + date_vec.tolist()[-1].strftime("%Y%m%d")
            )

            if prdcfg.get("add_date_in_fname", False):
                figtimeinfo = date_vec[0]
                fig_timeformat = "%Y"

            figfname_list = make_filename(
                "ts",
                prdcfg["dstype"],
                prdcfg["voltype"],
                prdcfg["imgformat"],
                prdcfginfo=prdcfginfo,
                timeinfo=figtimeinfo,
                timeformat=fig_timeformat,
            )

            for j, figfname in enumerate(figfname_list):
                figfname_list[j] = savedir + figfname

            np_min = prdcfg.get("npoints_min", 0)
            corr_min = prdcfg.get("corr_min", 0.0)

            plot_field_name = field_name
            if transform_str != "x":
                plot_field_name = f"f({field_name})"

            titl = (
                rad1_name
                + "-"
                + rad2_name
                + " "
                + plot_field_name
                + " intercomparison "
                + titldate
                + rangebin_info_title
            )

            plot_intercomp_scores_ts(
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
                figfname_list,
                ref_value=0.0,
                np_min=np_min,
                corr_min=corr_min,
                labelx="Time UTC",
                titl=titl,
            )

            print("----- save to " + " ".join(figfname_list))

            fname_list.extend(figfname_list)

        return fname_list

    else:
        warn(" Unsupported product type: " + prdcfg["type"])
    return None


def generate_colocated_gates_products(dataset, prdcfg):
    """
    Generates colocated gates products. Accepted product types:
        'WRITE_COLOCATED_GATES': Writes the position of the co-located gates
            in a csv file.
        'PLOT_COLOCATED_GATES': Plots the colocated gates over a geo-referenced map,
            as well as the coordinates of both radars. The plot performs a gridding and
            plots the number of cologates gates within every grid-cell as well as their
            average height. All other plotting parameters are obtained from the
            gridMapImageConfig structure of the loc file.
            User defined parameters:
                grid_size: float
                    Size of the grid in degrees (lat/lon). If not provided,
                    the grid size defined in gridMapImageConfig with
                    lat_min, lon_min, lat_max, lon_max and lonstep, latstep
                    will be used.
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
    if prdcfg["type"] == "WRITE_COLOCATED_GATES":
        if prdcfg["radar"] not in dataset:
            return None
        if "coloc_dict" not in dataset[prdcfg["radar"]]:
            return None

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            "colocated_gates",
            prdcfg["prdname"],
            timeinfo=None,
        )

        fname = make_filename(
            "info", prdcfg["dstype"], prdcfg["prdname"], ["csv"], timeinfo=None
        )

        fname = savedir + fname[0]

        write_colocated_gates(dataset[prdcfg["radar"]]["coloc_dict"], fname)

        print("saved colocated gates file: " + fname)

        return fname

    if prdcfg["type"] == "PLOT_COLOCATED_GATES":
        dssavedir = prdcfg["dsname"]
        if "dssavename" in prdcfg:
            dssavedir = prdcfg["dssavename"]
        prdsavedir = prdcfg["prdname"]

        grid_size = prdcfg.get("grid_size", None)

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdsavedir,
            timeinfo=None,
        )
        fname_list = make_filename(
            "map",
            prdcfg["dstype"],
            "counts",
            prdcfg["imgformat"],
        )

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir + fname

        plot_colocated_gates(dataset, prdcfg, fname_list, grid_size=grid_size)

        print("----- save to " + " ".join(fname_list))

        return fname_list

    if prdcfg["radar"] not in dataset:
        return None
    if "radar_out" not in dataset[prdcfg["radar"]]:
        return None

    prdcfg["timeinfo"] = None
    return generate_vol_products(dataset[prdcfg["radar"]], prdcfg)


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
