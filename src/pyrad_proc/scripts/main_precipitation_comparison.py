#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_precipitation_comparison
================================================

This program compares radar data with a point measurement sensor.

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import glob
import os
import argparse

import numpy as np

from pyrad.io import read_ts_cum
from pyrad.graph import plot_scatter_comp
from pyrad.util import compute_1d_stats

print(__doc__)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Compare radar data with a point measurement sensor."
    )
    parser.add_argument(
        "directories",
        help="Directories where the precipitation products are stored, multiple directories (corresponding typically to multiple radars) can be specified (comma-separated)"
        + "If multiple directories are specified, both an intercomparison between directories will be performed, as well as a "
        + "separate analysis of cumulated precipitation for every directory",
    )
    parser.add_argument(
        "--parameters", default=["RR"], nargs="+", help="List of parameters to check"
    )

    return parser.parse_args()


def find_max_min_dates(file_paths):
    dates = []

    for file_path in file_paths:
        # Extract the date component from the filename
        filename = os.path.basename(file_path)
        date_str = filename.split("_")[0]

        try:
            # Convert the date string to a datetime object
            date = datetime.datetime.strptime(date_str, "%Y%m%d")
            dates.append(date)
        except ValueError:
            # Handle files that don't have a valid date format
            print(f"Warning: '{file_path}' does not contain a valid date format.")

    if dates:
        # Find the minimum and maximum dates
        min_date = min(dates)
        max_date = max(dates)
        return min_date, max_date
    else:
        return None, None


def main():
    """Main function."""
    args = parse_args()

    fpaths = args.directories
    fpaths = fpaths.split(",")
    params = args.parameters

    np_radar_min = 6
    np_sensor_min = 6
    min_val = 0.2

    img_ext = "png"
    avg_time = 3600

    print(
        "====== precipitation comparison started: %s"
        % datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
    )

    for fpath in fpaths:
        print(f"Processing directory {fpath}")

        # Get name of pyrad output
        pyrad_name = [k for k in fpath.split("/") if len(k)][-1]

        # Find all products within file structure
        all_RR_files = glob.glob(
            os.path.join(fpath, "**", "*acc_ts_comp_POINT_MEASUREMENT_RR*.csv"),
            recursive=True,
        )

        print(f"Found {len(all_RR_files)} files with acc. precip.")

        startdate, enddate = find_max_min_dates(all_RR_files)
        print(f"Start date: {startdate}")
        print(f"End date: {enddate}")

        # Sort files by parameter
        list_params = sorted(params)[::-1]
        RR_files_by_param = {}
        for file in all_RR_files:
            for param in list_params:
                if param in file:
                    if param not in RR_files_by_param:
                        RR_files_by_param[param] = []
                    RR_files_by_param[param].append(file)

        for param in list_params:
            ts_vec = np.array([])
            val_radar = np.ma.array([])
            np_radar = np.array([])
            val_sensor = np.ma.array([])
            np_sensor = np.array([])
            for file in RR_files_by_param[param]:
                (
                    ts_aux,
                    np_radar_aux,
                    radar_value_aux,
                    np_sensor_aux,
                    sensor_value_aux,
                ) = read_ts_cum(file)
                ts_vec = np.append(ts_vec, ts_aux)
                val_radar = np.ma.append(val_radar, radar_value_aux)
                np_radar = np.append(np_radar, np_radar_aux)
                val_sensor = np.ma.append(val_sensor, sensor_value_aux)
                np_sensor = np.append(np_sensor, np_sensor_aux)
            import pdb

            pdb.set_trace()
            # filter out undesired data
            ind = np.where(
                np.logical_and(
                    np.logical_and(
                        np_radar >= np_radar_min, np_sensor >= np_sensor_min
                    ),
                    np.logical_and(val_sensor >= min_val, val_radar >= min_val),
                )
            )[0]

            val_sensor = val_sensor[ind]
            val_radar = val_radar[ind]

            # compute statistics
            stats = compute_1d_stats(val_sensor, val_radar)

            # create output image
            fpath = fpath + "RR/"
            if os.path.isdir(fpath):
                pass
            else:
                os.makedirs(fpath)
            print(f"Saving outputs to {fpath}")
            figfname = [
                startdate.strftime("%Y%m%d")
                + "-"
                + enddate.strftime("%Y%m%d")
                + "_"
                + str(avg_time)
                + "s_acc_ts_comp_"
                + param
                + "."
                + img_ext
            ]

            for i in range(len(figfname)):
                figfname[i] = fpath + figfname[i]

            labelx = "RG (mm)"
            labely = "Radar (mm)"
            titl = (
                str(avg_time)
                + " s Acc. Comp. "
                + startdate.strftime("%Y%m%d")
                + "-"
                + enddate.strftime("%Y%m%d")
            )

            metadata = (
                "npoints: "
                + str(stats["npoints"])
                + "\n"
                + "NB: "
                + "{:.2f}".format(float(stats["NB"]))
                + "\n"
                + "corr: "
                + "{:.2f}".format(float(stats["corr"]))
                + "\n"
                + "RMS: "
                + "{:.2f}".format(float(stats["RMS"]))
                + "\n"
                + "Nash: "
                + "{:.2f}".format(float(stats["Nash"]))
                + "\n"
            )

            plot_scatter_comp(
                val_sensor,
                val_radar,
                figfname,
                labelx=labelx,
                labely=labely,
                titl=titl,
                axis="equal",
                metadata=metadata,
                dpi=300,
            )


def _print_end_msg(text):
    """Prints end message."""
    print(
        text
        + datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
    )


# ---------------------------------------------------------
# Start main:
# ---------------------------------------------------------
if __name__ == "__main__":
    main()
