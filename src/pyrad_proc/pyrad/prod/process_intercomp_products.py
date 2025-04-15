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

from warnings import warn

import numpy as np

from .process_vol_products import generate_vol_products

from ..io.io_aux import get_fieldname_pyart
from ..io.io_aux import get_save_dir, make_filename

from ..io.read_data_other import read_intercomp_scores_ts

from ..io.write_data import write_colocated_gates, write_colocated_data
from ..io.write_data import write_colocated_data_time_avg
from ..io.write_data import write_intercomp_scores_ts

from ..graph.plots import plot_scatter
from ..graph.plots_timeseries import plot_intercomp_scores_ts

from ..util.radar_utils import compute_2d_stats
from ..util.stat_utils import parse_math_expression


def generate_intercomp_products(dataset, prdcfg):
    """
    Generates radar intercomparison products. Accepted product types:
        'PLOT_AND_WRITE_INTERCOMP_TS': Writes statistics of radar
            intercomparison in a file and plots the time series of the
            statistics.
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
                corr_min: float
                    The minimum correlation to consider the statistics
                    valid and therefore use the data point in the plotting.
                    Default 0.
        'PLOT_SCATTER_INTERCOMP': Plots a density plot with the points of
            radar 1 versus the points of radar 2
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                step: float
                    The quantization step of the data. If none it will be
                    computed using the Py-ART config file. Default None
                scatter_type: str
                    Type of scatter plot. Can be a plot for each radar volume
                    (instant) or at the end of the processing period
                    (cumulative). Default is cumulative
                range_bins: list of floats
                    Range bins to radar 1 to consider, a different plot will be generated for
                    all range bins.This can be used to differentiate the analysis in different
                    radar ranges.
                    Default is to generate only one figure, regardless of the distance to radar 1.
                transform: str
                    A transform to apply to the data before computing the statistics
                    and generating the plots. Any mathematical function with argument "x"
                    is accepted for example "sin(x) + x**2" or "10 * log10(x)"
                    Default is to use no transform
        'WRITE_INTERCOMP': Writes the instantaneously intercompared data
            (gate positions, values, etc.) in a csv file.
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
        'WRITE_INTERCOMP_TIME_AVG': Writes the time-averaged intercompared
            data (gate positions, values, etc.) in a csv file.
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

    if prdcfg["type"] == "WRITE_INTERCOMP_TIME_AVG":
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
        write_colocated_data_time_avg(dataset["intercomp_dict"], fname)
        print("saved colocated time averaged data file: " + fname)

        return fname

    if prdcfg["type"] == "PLOT_SCATTER_INTERCOMP":
        scatter_type = prdcfg.get("scatter_type", "cumulative")
        if scatter_type == "cumulative" and not dataset["final"]:
            return None

        timeformat = "%Y%m%d"
        if scatter_type == "instant":
            timeformat = "%Y%m%d%H%M%S"

        field_name = get_fieldname_pyart(prdcfg["voltype"])
        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=dataset["timeinfo"],
        )

        transform_str = prdcfg.get("transform", "x")
        transform = parse_math_expression(transform_str)
        range_bins = prdcfg.get("range_bins", [0, np.inf])
        step = prdcfg.get("step", None)

        fname_list = []
        for i in range(len(range_bins) - 1):  # loop on range bins
            if len(range_bins) > 2:
                rangebin_info = f"range_bin_{range_bins[i]:.0f}_{range_bins[i+1]:.0f}m"
                rangebin_info_title = (
                    f"range_bin {range_bins[i]:.0f}-{range_bins[i+1]:.0f}m"
                )
            else:
                rangebin_info = ""  # only one range bin, leave empty
                rangebin_info_title = rangebin_info

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

            titl = f"colocated radar gates \n {rangebin_info_title} " + dataset[
                "timeinfo"
            ].strftime(timeformat)

            selection = np.logical_and(
                dataset["intercomp_dict"]["rad1_rng"] >= range_bins[i],
                dataset["intercomp_dict"]["rad1_rng"] < range_bins[i + 1],
            )
            if not np.sum(selection):  # skip if selection empty
                continue

            hist_2d, bin_edges1, bin_edges2, stats = compute_2d_stats(
                np.ma.asarray(dataset["intercomp_dict"]["rad1_val"][selection]),
                np.ma.asarray(dataset["intercomp_dict"]["rad2_val"][selection]),
                field_name,
                field_name,
                step1=step,
                step2=step,
                transform=transform,
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
                field_name = f"{transform_str} of {field_name}"
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
            )

            fname_list.extend(f_list)

            print("----- save to " + " ".join(f_list))

        return fname_list

    if prdcfg["type"] == "PLOT_AND_WRITE_INTERCOMP_TS":
        if not dataset["final"]:
            return None

        field_name = get_fieldname_pyart(prdcfg["voltype"])
        step = prdcfg.get("step", None)
        transform = parse_math_expression(prdcfg.get("transform", "x"))
        rad1_name = dataset["intercomp_dict"]["rad1_name"]
        rad2_name = dataset["intercomp_dict"]["rad2_name"]

        hist_2d, bin_edges1, bin_edges2, stats = compute_2d_stats(
            transform(np.ma.asarray(dataset["intercomp_dict"]["rad1_val"])),
            transform(np.ma.asarray(dataset["intercomp_dict"]["rad2_val"])),
            field_name,
            field_name,
            step1=step,
            step2=step,
        )

        # put time info in file path and name
        csvtimeinfo_file = None
        timeformat = None
        sort_by_date = False
        rewrite = False
        if "add_date_in_fname" in prdcfg:
            if prdcfg["add_date_in_fname"]:
                csvtimeinfo_file = dataset["timeinfo"]
                timeformat = "%Y"
        if "sort_by_date" in prdcfg:
            sort_by_date = prdcfg["sort_by_date"]
        if "rewrite" in prdcfg:
            rewrite = prdcfg["rewrite"]

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=None,
        )

        csvfname = make_filename(
            "ts",
            prdcfg["dstype"],
            prdcfg["voltype"],
            ["csv"],
            prdcfginfo=rad1_name + "-" + rad2_name,
            timeinfo=csvtimeinfo_file,
            timeformat=timeformat,
        )[0]

        csvfname = savedir + csvfname

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
            return None

        if len(date_vec) < 2:
            warn("Unable to plot time series. Not enough points")
            return None

        if rewrite:
            stats = {
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
                stats,
                field_name,
                csvfname,
                rad1_name=rad1_name,
                rad2_name=rad2_name,
            )

        figtimeinfo = None
        titldate = (
            date_vec[0].strftime("%Y%m%d") + "-" + date_vec[-1].strftime("%Y%m%d")
        )
        if "add_date_in_fname" in prdcfg:
            if prdcfg["add_date_in_fname"]:
                figtimeinfo = date_vec[0]
                timeformat = "%Y"

        figfname_list = make_filename(
            "ts",
            prdcfg["dstype"],
            prdcfg["voltype"],
            prdcfg["imgformat"],
            prdcfginfo=rad1_name + "-" + rad2_name,
            timeinfo=figtimeinfo,
            timeformat=timeformat,
        )

        for i, figfname in enumerate(figfname_list):
            figfname_list[i] = savedir + figfname

        np_min = prdcfg.get("npoints_min", 0)
        corr_min = prdcfg.get("corr_min", 0.0)

        titl = (
            rad1_name
            + "-"
            + rad2_name
            + " "
            + field_name
            + " intercomparison "
            + titldate
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

        return figfname_list

    warn(" Unsupported product type: " + prdcfg["type"])
    return None


def generate_colocated_gates_products(dataset, prdcfg):
    """
    Generates colocated gates products. Accepted product types:
        WRITE_COLOCATED_GATES: Writes the position of the co-located gates
            in a csv file
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
