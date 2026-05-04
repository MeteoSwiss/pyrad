"""
pyrad.prod.process_timeseries_products
======================================

Functions for obtaining Pyrad products from a time series datasets

.. autosummary::
    :toctree: generated/

    generate_timeseries_products

"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from netCDF4 import num2date
import pyart

from ..util import warn

from ..io.io_aux import get_save_dir, make_filename, get_fieldname_pyart
from ..io.io_aux import generate_field_name_str

from ..io.read_data_sensor import get_sensor_data
from ..io.read_data_other import read_sensor_scores_ts

from ..io.write_data import write_ts_polar_data
from ..io.write_data import write_ts_grid_data, write_multiple_points
from ..io.write_data import write_multiple_points_grid, write_sensor_scores

from ..graph.plots_timeseries import (
    plot_timeseries,
    plot_timeseries_comp,
    plot_sensor_scores_timeseries,
)
from ..graph.plots_vol import plot_cappi, plot_traj
from ..graph.plots import plot_scatter_comp

from ..util.radar_utils import rainfall_accumulation
from ..util.stat_utils import perfscores

# Matplotlib color cycler
mpl_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]


def generate_timeseries_products(dataset, prdcfg):
    """
    Generates time series products. Accepted product types:
        'POINT_TS': Writes and plots point time series for one or multiple radars.

            It generates a CSV file containing the time series values at a single
            point or radar gate. If multiple radar entries are present in the
            dataset, all radar series are written to the same CSV file and plotted
            together in one figure. If the timestamps of the different radars don't agree
            an outer join will be performed and missing timestamps will be filled with
            NaNs

            User defined parameters:
                set_time_info: bool
                    If True, date information is added to the output directory and
                    filenames. Default is True.
                do_only_final: bool
                    If True, the product is only generated when dataset['final'] is
                    True. Default is True.
                RadarName: list of str
                    Names of the radars to use as plot labels. The order should
                    match the RADAR entries in the dataset.
                basepath: str
                    Base path where the product output directory is created.
                procname: str
                    Processing name used to build the output directory.
                dstype: str
                    Dataset type used in the output filename.
                imgformat: list of str
                    List of image formats to generate.
                dpi: int
                    Figure resolution in dots per inch. Default is 72.
                vmin: float or None
                    Minimum y-axis value used in the plot. If None, the limit is
                    determined automatically.
                vmax: float or None
                    Maximum y-axis value used in the plot. If None, the limit is
                    determined automatically.
                label: str
                    Optional label added to the plot title.

    'COMPARE_POINT_SENSOR': Plots and/or compares radar and co-located gauge (sensor)
        data at a point of interest. The product always writes a combined CSV file
        containing both radar and sensor data. Optionally, the radar and sensor
        time series can be time-averaged or time-accumulated to a user-defined
        interval (if no interval is provided, native sampling is used). The plot
        can be either a dual time series comparison or a scatter plot of radar vs
        sensor.
        User defined parameters:
            set_time_info: bool
                If True, date information is added to the output directory and
                filenames. Default is True.
            do_only_final: bool
                If True, the product is only generated when dataset['final'] is
                True. Default is True.
            dpi: int
                The pixel density of the plot. Default 72
            vmin, vmax: float
                The limits of the Y-axis. If none they will be obtained
                from the Py-ART config file.
            plot_type: str
                The type of plot to generate. Can be 'timeseries' (default) or
                'scatter'. Scatter-plots also show the RMSE, corr and bias.
            double_conditional_threshold: float
                When using plot_type == "scatter", the data can be filtered with a double-conditional threshold.
                Default is -infinity (no thresholding)
            time_interval: int or None
                Time interval in seconds used to resample both radar and sensor
                series before comparison. If None, the native timestamps are used.
            agg_mode: str
                The aggregation method used when time_interval is set. Can be
                'mean' (time-averaging) or 'sum' (time-accumulation). Default
                'mean'.
            align: str
                How to align radar and sensor timestamps when producing paired
                values (used for CSV pairing and scatter plots). Can be 'nearest'
                (default) or 'inner'. 'inner' requires that the timestamps of the
                radar and gauge match perfectly.
            nearest_tol: int or None
                Maximum allowed time difference in seconds when align='nearest'.
                If None, no tolerance is applied (nearest_tol = infinity), which is not recommended when there are holes in the
                radar timeseries, as it will lead to erroneous interpolation.
            sensor: str
                The sensor type. Can be 'rgage', 'rgage_knmi' or 'disdro'
            sensorid: str
                The sensor ID.
            varname: str
                Name of the variable to use in the sensor files. It corresponds to the column you want to
                compare with your radar measurements.
            sensor_scaling: float
                Scaling factor for the sensor data. Use 6 for example to convert from mm for 10-min measurements
                to mm/hr (which is the radar reference). Default is 1 (no scaling)
            dir_template: str
                The directory template where to find the file. dir_template should be a directory template string
                (relative to rgpath/disdropath or absolute) that can contain {sensorid} and any strftime date format codes
                (e.g. %Y, %m, %d, %Y%m%d, optionally written as {%Y%m%d}), which will be expanded using the provided sensorid
                and date before performing the recursive search, for example "{%Y-%m-%d}/{sensorid}/"
                If not provided a potentially deep recursive search will be performed which can take some time.
    'COMPARE_SCORES_SENSOR': Computes performance scores using sensor data as a reference and radar
        data as the estimate.  It will generate a csv file containing the daily scores for one or multiple radars (if present)
        as well as update a plot of the time series of performance scores over time.
        User defined parameters:
            scores: list of str
                Which scores to compute, currently the following are supported:
                - corr (pearson correlation)
                - RMSE (root mean square error)
                - logBias (bias in log scale)
                - scatter (see https://doi.org/10.1256/qj.05.190)
                Default is RMSE and corr
            bounds: list of float
                If set will compute the scores in the specified bounds based on the sensor value.
                e.g. 0,1,10 will compute the scores in the [0,1] and [1,10] will be computed separately
                and will be written to the csv file.
                If not provided no bounds will be used
            double_conditional_threshold: float
                Optionally a double conditional threshold can be applied in the computation of the scores
                e.g. only points where both the estimate and the reference exceed this threshold will be used
                in the score computation
            time_interval: int or None
                Time interval in seconds used to resample both radar and sensor
                series before comparison. If None, the native timestamps are used.
            agg_mode: str
                The aggregation method used when time_interval is set. Can be
                'mean' (time-averaging) or 'sum' (time-accumulation). Default
                'mean'.
            align: str
                How to align radar and sensor timestamps when producing paired
                values (used for CSV pairing and scatter plots). Can be 'nearest'
                (default) or 'inner'. 'inner' requires that the timestamps of the
                radar and gauge match perfectly.
            nearest_tol: int or None
                Maximum allowed time difference in seconds when align='nearest'.
                If None, no tolerance is applied (nearest_tol = infinity), which is not recommended when there are holes in the
                radar timeseries, as it will lead to erroneous interpolation.
            sensor: str
                The sensor type. Can be 'rgage', 'rgage_knmi' or 'disdro'
            sensorid: str
                The sensor ID.
            varname: str
                Name of the variable to use in the sensor files. It corresponds to the column you want to
                compare with your radar measurements.
            sensor_scaling: float
                Scaling factor for the sensor data. Use 6 for example to convert from mm for 10-min measurements
                to mm/hr (which is the radar reference). Default is 1 (no scaling)
            dir_template: str
                The directory template where to find the file. dir_template should be a directory template string
                (relative to rgpath/disdropath or absolute) that can contain {sensorid} and any strftime date format codes
                (e.g. %Y, %m, %d, %Y%m%d, optionally written as {%Y%m%d}), which will be expanded using the provided sensorid
                and date before performing the recursive search, for example "{%Y-%m-%d}/{sensorid}/"
                If not provided a potentially deep recursive search will be performed which can take some time.
        'PLOT_AND_WRITE_TRAJ': Writes and plots a trajectory time series.
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                ymin, ymax: float
                    The minimum and maximum value of the Y-axis. If none it
                    will be obtained from the Py-ART config file.
        'PLOT_HIST_TRAJ': plots and writes a histogram of all the data gathered
            during the trajectory processing
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                step: float or None
                    The quantization step of the data. If None it will be
                    obtained from the Py-ART config file
        'CAPPI_IMAGE_TRAJ': Creates a CAPPI image with the trajectory position
            overplot on it.
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                color_ref: str
                    The meaning of the color code with which the trajectory is
                    plotted. Can be 'None', 'altitude' (the absolute
                    altitude), 'rel_altitude' (altitude relative to the CAPPI
                    altitude), 'time' (trajectory time respect of the start of
                    the radar scan leading to the CAPPI)
                altitude: float
                    The CAPPI altitude [m]
                wfunc: str
                    Function used in the gridding of the radar data. The
                    function types are defined in pyart.map.grid_from_radars.
                    Default 'NEAREST'
                res: float
                    The CAPPI resolution [m]. Default 500.
        'WRITE_MULTIPLE_POINTS': Writes a csv file containing data from
            multiple points obtained from a radar object
        'WRITE_MULTIPLE_POINTS_GRID': Writes a csv file containing data from
            multiple points obtained from a grid object

    Parameters
    ----------
    dataset : dictionary
        radar object

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    no return

    """
    dssavedir = prdcfg["dsname"]
    if "dssavename" in prdcfg:
        dssavedir = prdcfg["dssavename"]

    prdsavedir = prdcfg["prdname"]
    if "prdsavedir" in prdcfg:
        prdsavedir = prdcfg["prdsavedir"]

    if prdcfg["type"] == "WRITE_MULTIPLE_POINTS":
        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdsavedir,
            timeinfo=dataset["ref_time"],
        )

        csvfname = make_filename(
            "POI",
            prdcfg["dstype"],
            dataset["datatype"],
            ["csv"],
            timeinfo=dataset["ref_time"],
        )[0]

        csvfname = savedir + csvfname
        write_multiple_points(dataset, csvfname)
        print("saved CSV file: " + csvfname)

        return csvfname

    if prdcfg["type"] == "WRITE_MULTIPLE_POINTS_GRID":
        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdsavedir,
            timeinfo=dataset["ref_time"],
        )

        csvfname = make_filename(
            "POI",
            prdcfg["dstype"],
            dataset["datatype"],
            ["csv"],
            timeinfo=dataset["ref_time"],
        )[0]

        csvfname = savedir + csvfname
        write_multiple_points_grid(dataset, csvfname)
        print("saved CSV file: " + csvfname)

        return csvfname

    if prdcfg["type"] == "POINT_TS":
        # ---- configuration ----
        set_time_info = prdcfg.get("set_time_info", True)
        do_only_final = prdcfg.get("do_only_final", True)
        radar_names = prdcfg["RadarName"]
        basepath = prdcfg["basepath"]
        procname = prdcfg["procname"]
        dstype = prdcfg["dstype"]
        imgformat = prdcfg["imgformat"]
        label = prdcfg.get("label", "")
        dpi = prdcfg.get("dpi", 72)
        vmin = prdcfg.get("vmin", None)
        vmax = prdcfg.get("vmax", None)

        if do_only_final and not dataset["final"]:
            return None

        # ---- radar entries ----
        radar_keys = sorted(
            key
            for key in dataset.keys()
            if isinstance(key, str) and key.startswith("RADAR")
        )

        if not radar_keys:
            warn("Unable to plot time series. No radar entries in dataset")
            return None

        # Reference radar for metadata, filenames and time vector
        ds0 = dataset[radar_keys[0]]
        datatype = ds0["datatype"]

        # ---- gate / point information ----
        is_single_radar = len(radar_names) == 1
        has_antenna_coordinates = "antenna_coordinates_az_el_r" in ds0

        if has_antenna_coordinates and is_single_radar:
            az, el, rng = ds0["antenna_coordinates_az_el_r"]
            prdcfginfo = f"az{az:.1f}r{rng:.1f}el{el:.1f}"
        else:
            lon, lat, alt = ds0["point_coordinates_WGS84_lon_lat_alt"]
            prdcfginfo = f"lon{lon:.3f}lat{lat:.3f}alt{alt:.1f}"

        # ---- time information ----
        timeformat = "%Y%m%d" if set_time_info else None
        timeinfo = ds0.get("ref_time", None) if set_time_info else None

        date = ds0.get("time", None)
        if date is None or len(date) == 0:
            warn("Unable to plot time series. No valid time vector")
            return None

        timeinfo_fig = date[0] if set_time_info else None

        # ---- output paths ----
        savedir = get_save_dir(
            basepath,
            procname,
            dssavedir,
            prdsavedir,
            timeinfo=timeinfo,
        )

        csvfname = make_filename(
            "ts",
            dstype,
            datatype,
            ["csv"],
            prdcfginfo=prdcfginfo,
            timeinfo=timeinfo,
            timeformat=timeformat,
        )[0]
        csvfname = savedir + csvfname

        figfname_list = make_filename(
            "ts",
            dstype,
            datatype,
            imgformat,
            prdcfginfo=prdcfginfo,
            timeinfo=timeinfo_fig,
            timeformat=timeformat,
        )
        figfname_list = [savedir + fname for fname in figfname_list]

        # ---- write CSV ----
        if has_antenna_coordinates:
            write_ts_polar_data(dataset, csvfname)
        else:
            write_ts_grid_data(dataset, csvfname)

        print("saved CSV file: " + csvfname)

        # ---- extract radar time series ----
        series_list = []
        dates_list = []
        for radar_key in radar_keys:
            values = dataset[radar_key].get("value", None)
            times = dataset[radar_key].get("time", None)
            if values is None:
                series_list.append(np.full(len(date), np.nan))
                continue

            if np.ma.isMaskedArray(values):
                values = values.filled(np.nan)
            else:
                try:
                    values = np.array(
                        [
                            np.nan if value is np.ma.masked else value
                            for value in values
                        ],
                        dtype=float,
                    )
                except Exception:
                    values = np.asarray(values)

            series_list.append(values)
            dates_list.append(times)

        # ---- plot ----
        title = f"Time Series {label} {date[0].strftime('%Y-%m-%d')}"
        labely = generate_field_name_str(datatype)

        plot_timeseries(
            dates_list,
            series_list,
            figfname_list,
            labelx="Time UTC",
            labely=labely,
            labels=radar_names,
            title=title,
            dpi=dpi,
            ymin=vmin,
            ymax=vmax,
        )

        print("----- save to " + " ".join(figfname_list))

        figfname_list.append(csvfname)
        return figfname_list

    if prdcfg["type"] == "COMPARE_POINT_SENSOR":
        # -------------------------
        # 0) common options / guards
        # -------------------------
        dpi = prdcfg.get("dpi", 72)
        vmin = prdcfg.get("vmin", None)
        vmax = prdcfg.get("vmax", None)
        do_only_final = prdcfg.get("do_only_final", True)
        set_time_info = prdcfg.get("set_time_info", True)

        plot_type = prdcfg.get("plot_type", "timeseries")  # "timeseries" or "scatter"
        time_interval = prdcfg.get("time_interval", None)  # seconds or None
        base_time = prdcfg.get("base_time", 0)  # seconds
        agg_mode = prdcfg.get("agg_mode", "mean")  # "mean" or "sum"
        sensor_scaling = prdcfg.get("sensor_scaling", 1)
        align = prdcfg.get("align", "nearest")  # "nearest" or "inner"
        nearest_tol = prdcfg.get("nearest_tol", None)  # seconds or None
        sensor_id = prdcfg.get("sensorid", "")  # sensor ID
        doublecond = prdcfg.get("double_conditional_threshold", -np.inf)

        radar_keys = sorted(
            [k for k in dataset.keys() if isinstance(k, str) and k.startswith("RADAR")]
        )
        if not radar_keys:
            warn("COMPARE_POINT_SENSOR: no RADAR### keys found in dataset")
            return None

        ds0 = dataset[radar_keys[0]]
        if do_only_final and not dataset["final"]:
            return None

        # -------------------------
        # 1) point / gateinfo (use reference radar)
        # -------------------------
        is_single_radar = len(prdcfg["RadarName"]) == 1
        has_antenna = "antenna_coordinates_az_el_r" in ds0

        if has_antenna and is_single_radar:
            az, el, r = ds0["antenna_coordinates_az_el_r"]
            gateinfo = f"az{az:.1f}r{r:.1f}el{el:.1f}"
            label_radar_base = f"Radar (az, el, r): ({az:.1f}, {el:.1f}, {r:.1f})"
        else:
            lon, lat, alt = ds0["point_coordinates_WGS84_lon_lat_alt"]
            gateinfo = f"lon{lon:.3f}lat{lat:.3f}alt{alt:.1f}"
            prefix = "Radar" if has_antenna else "Grid"
            label_radar_base = (
                f"{prefix} (lon, lat, alt): ({lon:.3f}, {lat:.3f}, {alt:.1f})"
            )

        timeformat = None
        timeinfo = None
        if set_time_info:
            timeformat = "%Y%m%d"
            timeinfo = ds0.get("ref_time", None)

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdsavedir,
            timeinfo=timeinfo,
        )

        # -------------------------
        # 2) build sensor dataframe once
        # -------------------------
        # Use last available radar time as query anchor (from reference radar)
        rad0_date = ds0.get("time", None)
        if rad0_date is None or len(rad0_date) == 0:
            warn("COMPARE_POINT_SENSOR: reference radar has no date vector")
            return None

        sensordate, sensorvalue, sensortype, _ = get_sensor_data(
            rad0_date[0], rad0_date[-1], ds0["datatype"], prdcfg
        )

        if sensordate is None:
            warn("Unable to compare at POI. No valid sensor data")
            return None

        label_sensor = f"{sensortype} {sensor_id}"
        labely = generate_field_name_str(ds0["datatype"])

        df_s = (
            pd.DataFrame(
                {
                    "date": pd.to_datetime(sensordate, errors="coerce", format="mixed"),
                    "sensor": np.asarray(sensorvalue, dtype=float)
                    * float(sensor_scaling),
                }
            )
            .dropna(subset=["date"])
            .sort_values("date")
        )

        if df_s.empty:
            warn("Sensor dataframe empty after cleaning")
            return None

        # optional resampling/accumulation on sensor
        if time_interval is not None:
            t_out_vec, val_out_vec, np_vec = rainfall_accumulation(
                df_s["date"].to_numpy(),
                df_s["sensor"].to_numpy(),
                time_interval,
                base_time,
            )
            df_s = pd.DataFrame(
                {
                    "date": pd.to_datetime(t_out_vec, errors="coerce", format="mixed"),
                    "sensor": val_out_vec,
                    "n_agg_pts_sensor": np_vec,
                }
            )

        df_s = df_s.dropna(subset=["date"]).sort_values("date").reset_index(drop=True)

        # -------------------------
        # 3) build combined dataframe (sensor axis) + collect series for plotting
        # -------------------------
        df_out = df_s[["date", "sensor"]].copy()

        # for plotting (keep original/resampled series per radar + sensor)
        plot_dates = [
            pd.to_datetime(df_s["date"], errors="coerce", format="mixed").to_list()
        ]
        plot_values = [df_s["sensor"].to_numpy()]
        plot_labels = [label_sensor]
        # merge_asof tolerance
        tol = None
        if nearest_tol is not None:
            tol = pd.Timedelta(seconds=float(nearest_tol))

        for radarnr in radar_keys:
            ds = dataset[radarnr]
            radarvalue = ds.get("value", None)
            radardate = ds.get("time", None)  # as requested
            if radarvalue is None or radardate is None:
                warn(f"{radarnr}: missing date/value -> column will be NaN")
                df_out[f"radar_{radarnr.lower()}"] = np.nan
                continue

            # masked -> nan
            if np.ma.isMaskedArray(radarvalue):
                radarvalue = radarvalue.filled(np.nan)
            else:
                try:
                    radarvalue = [
                        np.nan if x is np.ma.masked else x for x in radarvalue
                    ]
                except Exception:
                    radarvalue = np.asarray(radarvalue)

            df_r = (
                pd.DataFrame(
                    {
                        "date": pd.to_datetime(
                            radardate, errors="coerce", format="mixed"
                        ),
                        "radar": np.asarray(radarvalue, dtype=float),
                    }
                )
                .dropna(subset=["date"])
                .sort_values("date")
            )

            if df_r.empty:
                warn(f"{radarnr}: radar dataframe empty -> column will be NaN")
                df_out[f"radar_{radarnr.lower()}"] = np.nan
                continue

            # optional resampling/accumulation on radar
            if time_interval is not None:
                t_out_vec, val_out_vec, np_vec = rainfall_accumulation(
                    df_r["date"].to_numpy(),
                    df_r["radar"].to_numpy(),
                    time_interval,
                    base_time,
                )
                df_r = pd.DataFrame(
                    {
                        "date": pd.to_datetime(
                            t_out_vec, errors="coerce", format="mixed"
                        ),
                        "radar": val_out_vec,
                        "n_agg_pts": np_vec,
                    }
                )

            df_r = (
                df_r.dropna(subset=["date"]).sort_values("date").reset_index(drop=True)
            )

            # align radar to sensor dates (sensor-axis table)
            if align == "inner":
                # exact timestamp match; keep sensor axis and place NaN where no match
                df_tmp = pd.merge(df_out[["date"]], df_r, on="date", how="left")
            else:
                df_tmp = pd.merge_asof(
                    df_out[["date"]].sort_values("date"),
                    df_r.sort_values("date"),
                    on="date",
                    direction="nearest",
                    tolerance=tol,
                )

            df_out[f"radar_{radarnr.lower()}"] = df_tmp["radar"].to_numpy()

            # Collect for plotting
            plot_dates.append(
                pd.to_datetime(
                    df_out["date"], errors="coerce", format="mixed"
                ).to_list()
            )
            plot_values.append(df_tmp["radar"].to_numpy())

        # Replace NaN
        df_out.fillna(pyart.config.get_fillvalue(), inplace=True)
        # Update plot labels
        plot_labels.extend(prdcfg["RadarName"])

        # -------------------------
        # 4) write one combined CSV
        # -------------------------
        csv_tag = "ts_comp"
        if time_interval is not None:
            csv_tag = f"{int(time_interval)}s_{agg_mode}_ts_comp"

        timeinfo_csv = rad0_date[-1] if set_time_info else None

        comp_csvfname = make_filename(
            csv_tag,
            prdcfg["dstype"],
            ds0["datatype"],
            ["csv"],
            prdcfginfo=gateinfo,
            timeinfo=timeinfo_csv,
            timeformat=timeformat,
        )[0]
        comp_csvfname = savedir + comp_csvfname

        header = (
            "# Weather radar/sensor comparison timeseries data file\n"
            '# Comment lines are preceded by "#"\n'
            "# Description: \n"
            "# Time series of radar and sensor data over a fixed location.\n"
            f"# Location: {label_radar_base}\n"
            f"# {sensortype}: {sensor_id}\n"
            f"# Time aggregation period (s): {time_interval}\n"
            f"# Time aggregation method: {agg_mode}\n"
            f"# Sensor/Radar timestamp alignment method: {align}\n"
            f"# Sensor/Radar timestamp alignment tolerance (s): {nearest_tol}\n"
            f"# Data: {generate_field_name_str(ds0.get('datatype'))}\n"
            f"# Fill Value: {pyart.config.get_fillvalue()}\n"
            f"# Start: {rad0_date[0].strftime('%Y-%m-%d %H:%M:%S UTC')}\n"
            "#\n"
        )
        with open(comp_csvfname, "w", encoding="utf-8", newline="") as f:
            f.write(header)
            df_out.to_csv(f, index=False)

        print("saved radar+sensor CSV file (all radars): " + comp_csvfname)

        # -------------------------
        # 5) plot one figure (sensor + all radars)
        # -------------------------
        timeinfo_fig = timeinfo_csv

        fig_tag = "ts_comp" if plot_type == "timeseries" else "scatter_comp"
        if time_interval is not None:
            fig_tag = f"{int(time_interval)}s_{agg_mode}_" + fig_tag

        figfname_list = make_filename(
            fig_tag,
            prdcfg["dstype"],
            ds0["datatype"],
            prdcfg["imgformat"],
            prdcfginfo=gateinfo,
            timeinfo=timeinfo_fig,
            timeformat=timeformat,
        )
        figfname_list = [savedir + f for f in figfname_list]

        if plot_type == "scatter":
            x = plot_values[0]
            ys = plot_values[1:]

            if len(ys) == 0:
                warn("No radar series available to plot scatter")
                return None

            # Create list of x,y pairs to plot
            pairs = []
            for y in ys:
                # Keep only valid x,y
                ok = np.isfinite(x)
                x_pair = x[ok]
                y_pair = y[ok]

                # Use double conditional threshold
                ok = (x_pair >= doublecond) & (y_pair >= doublecond)
                x_pair = x_pair[ok]
                y_pair = y_pair[ok]

                if np.any(np.isfinite(x_pair) & np.isfinite(y_pair)):
                    pairs.append((x_pair, y_pair))

            if not len(pairs):
                warn("No matched radar/sensor pairs to plot scatter")
                return None

            labelx = f"{label_sensor} ({labely})" if labely else label_sensor
            labely2 = f"{label_radar_base}\n({labely})" if labely else label_radar_base

            titl = "Radar vs Sensor"
            if time_interval is not None:
                titl += f" ({int(time_interval)}s {agg_mode})"
            titl += " " + pd.to_datetime(
                df_out["date"].iloc[0], errors="coerce", format="mixed"
            ).strftime("%Y-%m-%d")

            plot_scatter_comp(
                pairs,
                figfname_list,
                labelx=labelx,
                labely=labely2,
                labels=plot_labels[1:],
                titl=titl,
                axis="equal",
                dpi=dpi,
            )

        else:
            # timeseries comparison
            # Use the existing compare plot helper (expects separate arrays)
            titl = f"Time Series comparison with {sensortype} sensor {sensor_id}"
            if time_interval is not None:
                titl += f" ({int(time_interval)}s {agg_mode})"
            titl += " " + pd.to_datetime(
                df_s["date"].iloc[0], format="mixed", errors="coerce"
            ).strftime("%Y-%m-%d")

            colors = ["#000000"] + [mpl_colors[i] for i in range(len(plot_dates))]
            linestyles = ["--"] + ["-"] * len(plot_dates)
            plot_timeseries_comp(
                plot_dates,
                plot_values,
                figfname_list,
                labelx="Time UTC",
                labely=labely,
                labels=plot_labels,  # sensor + RADAR001.. lines
                titl=titl,
                ymin=vmin,
                ymax=vmax,
                dpi=dpi,
                colors=colors,
                linestyles=linestyles,
            )

        print("----- save to " + " ".join(figfname_list))

        # Return outputs
        figfname_list.append(comp_csvfname)
        return figfname_list

    if prdcfg["type"] == "COMPARE_SCORES_SENSOR":
        dpi = prdcfg.get("dpi", 72)
        bounds = prdcfg.get("bounds", None)
        time_interval = prdcfg.get("time_interval", None)  # seconds or None
        base_time = prdcfg.get("base_time", 0)  # seconds
        agg_mode = prdcfg.get("agg_mode", "mean")  # "mean" or "sum"
        sensor_scaling = prdcfg.get("sensor_scaling", 1)
        align = prdcfg.get("align", "nearest")  # "nearest" or "inner"
        nearest_tol = prdcfg.get("nearest_tol", None)  # seconds or None
        doublecond = prdcfg.get("double_conditional_threshold", -np.inf)
        selected_scores = prdcfg.get("scores", ["RMSE", "corr"])
        linearize = prdcfg.get("linearize", False)

        if not dataset["final"]:
            return None

        radar_keys = sorted(
            [k for k in dataset.keys() if isinstance(k, str) and k.startswith("RADAR")]
        )
        # Reference radar for metadata/filenames
        ds0 = dataset[radar_keys[0]]

        # -------------------------
        # 1) point / gateinfo (use reference radar)
        # -------------------------
        if (
            "antenna_coordinates_az_el_r" in ds0 and len(prdcfg["RadarName"]) == 1
        ):  # Single radar case
            az = "{:.1f}".format(ds0["antenna_coordinates_az_el_r"][0])
            el = "{:.1f}".format(ds0["antenna_coordinates_az_el_r"][1])
            r = "{:.1f}".format(ds0["antenna_coordinates_az_el_r"][2])
            gateinfo = "az" + az + "r" + r + "el" + el
        else:
            lon = "{:.3f}".format(ds0["point_coordinates_WGS84_lon_lat_alt"][0])
            lat = "{:.3f}".format(ds0["point_coordinates_WGS84_lon_lat_alt"][1])
            alt = "{:.1f}".format(ds0["point_coordinates_WGS84_lon_lat_alt"][2])
            gateinfo = "lon" + lon + "lat" + lat + "alt" + alt

        # -------------------------
        # 2) load sensor data ONCE (anchor on latest time from reference radar)
        # -------------------------
        rad0_date = ds0.get("time", None)  # or ds0.get("time", None)
        if rad0_date is None or len(rad0_date) == 0:
            warn("COMPARE_SCORES_SENSOR: reference radar has no date vector")
            return None

        sensordate, sensorvalue, sensortype, _ = get_sensor_data(
            rad0_date[0],
            rad0_date[-1],
            ds0["datatype"],
            prdcfg,
        )
        if sensordate is None:
            warn("Unable to compare at POI. No valid sensor data")
            return None

        label_sensor = f"{sensortype} {prdcfg.get('sensorid','')}"

        # Sensor dataframe
        df_s = (
            pd.DataFrame(
                {
                    "date": pd.to_datetime(sensordate, errors="coerce", format="mixed"),
                    "sensor": np.asarray(sensorvalue, dtype=float)
                    * float(sensor_scaling),
                }
            )
            .dropna(subset=["date"])
            .sort_values("date")
        )

        if df_s.empty:
            warn("Sensor dataframe empty after cleaning")
            return None

        # Optional aggregation on sensor
        if time_interval is not None:
            t_out_vec, val_out_vec, np_vec = rainfall_accumulation(
                df_s["date"].to_numpy(),
                df_s["sensor"].to_numpy(),
                time_interval,
                base_time,
            )
            df_s = pd.DataFrame(
                {
                    "date": pd.to_datetime(t_out_vec, errors="coerce", format="mixed"),
                    "sensor": val_out_vec,
                    "n_agg_pts": np_vec,
                }
            ).sort_values("date")

        # -------------------------
        # 3) per-radar: build radar df, optionally aggregate, align, compute scores
        # -------------------------
        tol = None
        if nearest_tol is not None:
            tol = pd.Timedelta(seconds=float(nearest_tol))

        # This is what we will feed to write_sensor_scores in multi-radar mode
        # scores_out[radarnr][bounds][metric] = scalar OR length-1 array
        scores_out = {}

        # Use last radar time across radars for file naming (most recent)
        last_times = []

        for radarnr in radar_keys:
            ds = dataset[radarnr]

            radarvalue = ds.get("value", None)
            radardate = ds.get("time", None)

            if radarvalue is None or radardate is None or len(radardate) == 0:
                warn(f"{radarnr}: missing date/value -> skipping")
                continue

            last_times.append(radardate[-1])

            # masked -> nan
            if np.ma.isMaskedArray(radarvalue):
                radarvalue = radarvalue.filled(np.nan)
            else:
                try:
                    radarvalue = [
                        np.nan if x is np.ma.masked else x for x in radarvalue
                    ]
                except Exception:
                    radarvalue = np.asarray(radarvalue)

            df_r = (
                pd.DataFrame(
                    {
                        "date": pd.to_datetime(
                            radardate, errors="coerce", format="mixed"
                        ),
                        "radar": np.asarray(radarvalue, dtype=float),
                    }
                )
                .dropna(subset=["date"])
                .sort_values("date")
            )

            if df_r.empty:
                warn(f"{radarnr}: radar dataframe empty after cleaning -> skipping")
                continue

            # Optional aggregation on radar
            if time_interval is not None:
                t_out_vec, val_out_vec, np_vec = rainfall_accumulation(
                    df_r["date"].to_numpy(),
                    df_r["radar"].to_numpy(),
                    time_interval,
                    base_time,
                )
                df_r = pd.DataFrame(
                    {
                        "date": pd.to_datetime(
                            t_out_vec, errors="coerce", format="mixed"
                        ),
                        "radar": val_out_vec,
                        "n_agg_pts": np_vec,
                    }
                ).sort_values("date")

            # Align
            if align == "inner":
                df = pd.merge(df_r, df_s, on="date", how="inner")
            else:
                df = pd.merge_asof(
                    df_r.sort_values("date"),
                    df_s.sort_values("date"),
                    on="date",
                    direction="nearest",
                    tolerance=tol,
                )

            df = df.sort_values("date").reset_index(drop=True)

            # Compute scores
            all_scores = perfscores(
                df["radar"].to_numpy(),
                df["sensor"].to_numpy(),
                doublecond_thresh=doublecond,
                bounds=bounds,
                linearize=linearize,
            )

            # keep user-selected + NP (as in your original)
            scores_out[radarnr] = {}
            for b in all_scores:
                keep = selected_scores + ["NP"]
                scores_out[radarnr][b] = {
                    k: all_scores[b][k] for k in keep if k in all_scores[b]
                }

        if not scores_out:
            warn("COMPARE_SCORES_SENSOR: no radar produced valid scores")
            return None

        # -------------------------
        # 4) write single scores CSV (multi-radar)
        # -------------------------
        csvtimeinfo_path = None
        csvtimeinfo_file = None
        timeformat = None
        if prdcfg.get("add_date_in_fname", False):
            csvtimeinfo_file = ds0.get("timeinfo", None)
            timeformat = "%Y"

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=csvtimeinfo_path,
        )

        filenames = make_filename(
            "ts_scores",
            prdcfg["dstype"],
            ds0["datatype"],
            ext_list=["csv"] + prdcfg["imgformat"],
            timeinfo=csvtimeinfo_file,
            timeformat=timeformat,
            runinfo=prdcfg["runinfo"],
        )
        csvfname = savedir + filenames[0]
        figfname_list = [savedir + f for f in filenames[1:]]

        # choose one start_time for the score entry (use reference radar last time by default)
        # if you want "now", use dataset["timeinfo"] or similar
        start_time_for_row = max(last_times) if last_times else rad0_date[-1]

        write_sensor_scores(
            start_time_for_row,
            scores_out,
            doublecond,
            ds0["datatype"],
            csvfname,
            rewrite=False,
        )
        print("saved CSV file: " + csvfname)

        # -------------------------
        # 5) reread scores timeseries + plot
        # -------------------------
        scores_ts, meta = read_sensor_scores_ts(csvfname)

        titldate = (
            scores_ts["date"].iloc[0].strftime("%Y%m%d")
            + "-"
            + scores_ts["date"].iloc[-1].strftime("%Y%m%d")
        )
        titl = (
            f"{prdcfg['runinfo']} Sensor scores {titldate}\n"
            f"Bounds = {meta['bounds']}, Double conditional threshold = {meta['doublecond_thresh']}"
        )

        # Multi-radar plot: this function now plots all *_radar### metrics as separate lines
        plot_sensor_scores_timeseries(
            scores_ts,
            figfname_list,
            meta["bounds"],
            meta["doublecond_thresh"],
            titl=titl,
            dpi=dpi,
            labels=prdcfg["RadarName"],
        )
        print("----- save to " + " ".join(figfname_list))

        return [csvfname] + figfname_list

    # ================================================================
    if prdcfg["type"] == "PLOT_AND_WRITE_TRAJ":
        if not dataset["final"]:
            return None

        field_name = get_fieldname_pyart(prdcfg["voltype"])
        if field_name not in dataset["ts_dict"]:
            warn(
                " Field type "
                + field_name
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        ts = dataset["ts_dict"][field_name]
        timeinfo = ts.time_vector[0]

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdsavedir,
            timeinfo=timeinfo,
        )

        dstype_str = prdcfg["dstype"].lower().replace("_", "")
        csvfname = make_filename(
            "ts",
            dstype_str,
            ts.datatype,
            ["csv"],
            prdcfginfo=None,
            timeinfo=timeinfo,
            timeformat="%Y%m%d%H%M%S",
            runinfo=prdcfg["runinfo"],
        )[0]
        csvfname = savedir + csvfname

        ts.write(csvfname)

        figfname_list = make_filename(
            "ts",
            dstype_str,
            ts.datatype,
            prdcfg["imgformat"],
            prdcfginfo=None,
            timeinfo=timeinfo,
            timeformat="%Y%m%d%H%M%S",
            runinfo=prdcfg["runinfo"],
        )
        for i, figfname in enumerate(figfname_list):
            figfname_list[i] = savedir + figfname

        ymin = prdcfg.get("ymin", None)
        ymax = prdcfg.get("ymax", None)

        ts.plot(figfname_list[0], ymin=ymin, ymax=ymax)

        figfname_list.append(csvfname)
        return figfname_list

    if prdcfg["type"] == "PLOT_HIST_TRAJ":
        if not dataset["final"]:
            return None

        field_name = get_fieldname_pyart(prdcfg["voltype"])
        if field_name not in dataset["ts_dict"]:
            warn(
                " Field type "
                + field_name
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        ts = dataset["ts_dict"][field_name]
        timeinfo = ts.time_vector[0]

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdsavedir,
            timeinfo=timeinfo,
        )

        dstype_str = prdcfg["dstype"].lower().replace("_", "")

        figfname_list = make_filename(
            "hist",
            dstype_str,
            ts.datatype,
            prdcfg["imgformat"],
            prdcfginfo=None,
            timeinfo=timeinfo,
            timeformat="%Y%m%d%H%M%S",
            runinfo=prdcfg["runinfo"],
        )
        for i, figfname in enumerate(figfname_list):
            figfname_list[i] = savedir + figfname

        step = prdcfg.get("step", None)

        ts.plot_hist(figfname_list[0], step=step)

        return figfname_list

    if prdcfg["type"] == "CAPPI_IMAGE_TRAJ":
        if dataset["final"]:
            return None

        field_name = get_fieldname_pyart(prdcfg["voltype"])
        if field_name not in dataset["radar"].fields:
            warn(
                " Field type "
                + field_name
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdsavedir,
            timeinfo=prdcfg["timeinfo"],
        )

        color_ref = prdcfg.get("color_ref", "None")
        if color_ref == "altitude":
            prdtype = "cappi_alt"
        elif color_ref == "rel_altitude":
            prdtype = "cappi_rel_alt"
        elif color_ref == "time":
            prdtype = "cappi_time"
        else:
            prdtype = "cappi"

        figfname_list = make_filename(
            prdtype,
            prdcfg["dstype"],
            prdcfg["voltype"],
            prdcfg["imgformat"],
            prdcfginfo="alt" + "{:.1f}".format(prdcfg["altitude"]),
            timeinfo=prdcfg["timeinfo"],
            runinfo=prdcfg["runinfo"],
        )

        for i, fname in enumerate(figfname_list):
            figfname_list[i] = savedir + fname

        fig, ax = plot_cappi(
            dataset["radar"],
            field_name,
            prdcfg["altitude"],
            prdcfg,
            figfname_list,
            save_fig=False,
        )

        # start time of radar object
        t_start = dataset["radar"].time["data"].min()
        dt_start = num2date(
            t_start, dataset["radar"].time["units"], dataset["radar"].time["calendar"]
        )

        figfname_list = plot_traj(
            dataset["rng_traj"],
            dataset["azi_traj"],
            dataset["ele_traj"],
            dataset["time_traj"],
            prdcfg,
            figfname_list,
            rad_alt=dataset["radar"].altitude["data"],
            rad_tstart=dt_start,
            ax=ax,
            fig=fig,
            save_fig=True,
        )

        print("----- saved to " + " ".join(figfname_list))

        return figfname_list

    else:
        warn(" Unsupported product type: " + prdcfg["type"])
    return None
