"""
pyrad.proc.process_timeavg
=============================

Functions to obtain time averages of radar data mostly for intercomparison purposes

.. autosummary::
    :toctree: generated/

    process_time_stats
    process_time_stats2
    process_time_avg
    process_weighted_time_avg
    process_time_avg_flag

"""


import datetime
from ..util import warn
import numpy as np
import scipy
import pyart
from copy import deepcopy

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..util import time_avg_range
from ..io.read_data_radar import interpol_field


def process_time_stats(procstatus, dscfg, radar_list=None):
    """
    computes the temporal statistics of a field

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            Arbitrary data type supported by pyrad and contained in the radar data
        period : float. Dataset keyword
            the period to average [s]. If -1 the statistics are going to be
            performed over the entire data. Default 3600.
        start_average : float. Dataset keyword
            when to start the average [s from midnight UTC]. Default 0.
        lin_trans: int. Dataset keyword
            If 1 apply linear transformation before averaging
        use_nan : bool. Dataset keyword
            If true non valid data will be used
        nan_value : float. Dataset keyword
            The value of the non valid data. Default 0
        stat: string. Dataset keyword
            Statistic to compute: Can be mean, std, cov, min, max. Default
            mean
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the statistic computed on the input field, as well as
        "nsamples", as well as
        "sum2" (sum-squared) if stat in (cov, std), as well as

    ind_rad : int
        radar index

    """
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_name = get_fieldname_pyart(datatype)
        break
    ind_rad = int(radarnr[5:8]) - 1

    start_average = dscfg.get("start_average", 0.0)
    period = dscfg.get("period", 3600.0)
    lin_trans = dscfg.get("lin_trans", 0)
    use_nan = dscfg.get("use_nan", 0)
    nan_value = dscfg.get("nan_value", 0.0)
    stat = dscfg.get("stat", "mean")

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn("No valid radar")
            return None, None
        radar = radar_list[ind_rad]

        if field_name not in radar.fields:
            warn(field_name + " not available.")
            return None, None

        # Prepare auxiliary radar
        field = deepcopy(radar.fields[field_name])
        if stat in ("mean", "std", "cov"):
            if lin_trans:
                field["data"] = np.ma.power(10.0, 0.1 * field["data"])

            if use_nan:
                field["data"] = np.ma.asarray(field["data"].filled(nan_value))

            if stat in ("std", "cov"):
                sum2_dict = pyart.config.get_metadata("sum_squared")
                sum2_dict["data"] = field["data"] * field["data"]
        else:
            if use_nan:
                field["data"] = np.ma.asarray(field["data"].filled(nan_value))

        npoints_dict = pyart.config.get_metadata("number_of_samples")
        npoints_dict["data"] = np.ma.asarray(
            np.logical_not(np.ma.getmaskarray(field["data"])), dtype=int
        )

        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()
        radar_aux.add_field(field_name, field)
        radar_aux.add_field("number_of_samples", npoints_dict)

        if stat in ("std", "cov"):
            radar_aux.add_field("sum_squared", sum2_dict)

        # first volume: initialize start and end time of averaging
        if dscfg["initialized"] == 0:
            avg_par = dict()
            if period != -1:
                date_00 = dscfg["timeinfo"].replace(
                    hour=0, minute=0, second=0, microsecond=0
                )

                avg_par.update(
                    {"starttime": date_00 + datetime.timedelta(seconds=start_average)}
                )
                avg_par.update(
                    {
                        "endtime": avg_par["starttime"]
                        + datetime.timedelta(seconds=period)
                    }
                )
            else:
                avg_par.update({"starttime": dscfg["timeinfo"]})
                avg_par.update({"endtime": dscfg["timeinfo"]})

            avg_par.update({"timeinfo": dscfg["timeinfo"]})
            dscfg["global_data"] = avg_par
            dscfg["initialized"] = 1

        if dscfg["initialized"] == 0:
            return None, None

        dscfg["global_data"]["timeinfo"] = dscfg["timeinfo"]
        # no radar object in global data: create it
        if "radar_out" not in dscfg["global_data"]:
            if period != -1:
                # get start and stop times of new radar object
                (
                    dscfg["global_data"]["starttime"],
                    dscfg["global_data"]["endtime"],
                ) = time_avg_range(
                    dscfg["timeinfo"],
                    dscfg["global_data"]["starttime"],
                    dscfg["global_data"]["endtime"],
                    period,
                )

                # check if volume time older than starttime
                if dscfg["timeinfo"] > dscfg["global_data"]["starttime"]:
                    dscfg["global_data"].update({"radar_out": radar_aux})
            else:
                dscfg["global_data"].update({"radar_out": radar_aux})

            return None, None

        # still accumulating: add field to global field
        if period == -1 or dscfg["timeinfo"] < dscfg["global_data"]["endtime"]:
            if period == -1:
                dscfg["global_data"]["endtime"] = dscfg["timeinfo"]

            field_interp = interpol_field(
                dscfg["global_data"]["radar_out"], radar_aux, field_name
            )
            npoints_interp = interpol_field(
                dscfg["global_data"]["radar_out"], radar_aux, "number_of_samples"
            )

            if use_nan:
                field_interp["data"] = np.ma.asarray(
                    field_interp["data"].filled(nan_value)
                )
                dscfg["global_data"]["radar_out"].fields["number_of_samples"][
                    "data"
                ] += np.ma.asarray(
                    npoints_interp["data"].filled(fill_value=1), dtype=int
                )
            else:
                dscfg["global_data"]["radar_out"].fields["number_of_samples"][
                    "data"
                ] += np.ma.asarray(
                    npoints_interp["data"].filled(fill_value=0), dtype=int
                )

            if stat in ("mean", "std", "cov"):
                masked_sum = np.ma.getmaskarray(
                    dscfg["global_data"]["radar_out"].fields[field_name]["data"]
                )
                valid_sum = np.logical_and(
                    np.logical_not(masked_sum),
                    np.logical_not(np.ma.getmaskarray(field_interp["data"])),
                )

                dscfg["global_data"]["radar_out"].fields[field_name]["data"][
                    masked_sum
                ] = field_interp["data"][masked_sum]

                dscfg["global_data"]["radar_out"].fields[field_name]["data"][
                    valid_sum
                ] += field_interp["data"][valid_sum]

                if stat in ("cov", "std"):
                    dscfg["global_data"]["radar_out"].fields["sum_squared"]["data"][
                        masked_sum
                    ] = (
                        field_interp["data"][masked_sum]
                        * field_interp["data"][masked_sum]
                    )

                    dscfg["global_data"]["radar_out"].fields["sum_squared"]["data"][
                        valid_sum
                    ] += (
                        field_interp["data"][valid_sum]
                        * field_interp["data"][valid_sum]
                    )

            elif stat == "max":
                dscfg["global_data"]["radar_out"].fields[field_name][
                    "data"
                ] = np.maximum(
                    dscfg["global_data"]["radar_out"]
                    .fields[field_name]["data"]
                    .filled(fill_value=-1.0e300),
                    field_interp["data"].filled(fill_value=-1.0e300),
                )

                dscfg["global_data"]["radar_out"].fields[field_name][
                    "data"
                ] = np.ma.masked_values(
                    dscfg["global_data"]["radar_out"].fields[field_name]["data"],
                    -1.0e300,
                )
            elif stat == "min":
                dscfg["global_data"]["radar_out"].fields[field_name][
                    "data"
                ] = np.minimum(
                    dscfg["global_data"]["radar_out"]
                    .fields[field_name]["data"]
                    .filled(fill_value=1.0e300),
                    field_interp["data"].filled(fill_value=1.0e300),
                )

                dscfg["global_data"]["radar_out"].fields[field_name][
                    "data"
                ] = np.ma.masked_values(
                    dscfg["global_data"]["radar_out"].fields[field_name]["data"],
                    1.0e300,
                )

            return None, None

        # we have reached the end of the accumulation period: do the averaging
        # and start a new object (only reachable if period != -1)
        if stat in ("mean", "std", "cov"):
            field_mean = (
                dscfg["global_data"]["radar_out"].fields[field_name]["data"]
                / dscfg["global_data"]["radar_out"].fields["number_of_samples"]["data"]
            )

            if stat == "mean":
                if lin_trans:
                    dscfg["global_data"]["radar_out"].fields[field_name][
                        "data"
                    ] = 10.0 * np.ma.log10(field_mean)
                else:
                    dscfg["global_data"]["radar_out"].fields[field_name][
                        "data"
                    ] = field_mean
            elif stat in ("std", "cov"):
                field_std = np.ma.sqrt(
                    dscfg["global_data"]["radar_out"].fields["sum_squared"]["data"]
                    / dscfg["global_data"]["radar_out"].fields["number_of_samples"][
                        "data"
                    ]
                    - field_mean * field_mean
                )

                if stat == "std":
                    if lin_trans:
                        dscfg["global_data"]["radar_out"].fields[field_name][
                            "data"
                        ] = 10.0 * np.ma.log10(field_std)
                    else:
                        dscfg["global_data"]["radar_out"].fields[field_name][
                            "data"
                        ] = field_std
                else:
                    if lin_trans:
                        dscfg["global_data"]["radar_out"].fields[field_name][
                            "data"
                        ] = 10.0 * np.ma.log10(field_std / field_mean)
                    else:
                        dscfg["global_data"]["radar_out"].fields[field_name]["data"] = (
                            field_std / field_mean
                        )

        new_dataset = {
            "radar_out": deepcopy(dscfg["global_data"]["radar_out"]),
            "timeinfo": dscfg["global_data"]["endtime"],
        }

        dscfg["global_data"]["starttime"] += datetime.timedelta(seconds=period)
        dscfg["global_data"]["endtime"] += datetime.timedelta(seconds=period)

        # remove old radar object from global_data dictionary
        dscfg["global_data"].pop("radar_out", None)

        # get start and stop times of new radar object
        (
            dscfg["global_data"]["starttime"],
            dscfg["global_data"]["endtime"],
        ) = time_avg_range(
            dscfg["timeinfo"],
            dscfg["global_data"]["starttime"],
            dscfg["global_data"]["endtime"],
            period,
        )

        # check if volume time older than starttime
        if dscfg["timeinfo"] > dscfg["global_data"]["starttime"]:
            dscfg["global_data"].update({"radar_out": radar_aux})

        return new_dataset, ind_rad

    # no more files to process if there is global data pack it up
    if procstatus == 2:
        if dscfg["initialized"] == 0:
            return None, None
        if "radar_out" not in dscfg["global_data"]:
            return None, None

        if stat in ("mean", "std", "cov"):
            field_mean = (
                dscfg["global_data"]["radar_out"].fields[field_name]["data"]
                / dscfg["global_data"]["radar_out"].fields["number_of_samples"]["data"]
            )

            if stat == "mean":
                if lin_trans:
                    dscfg["global_data"]["radar_out"].fields[field_name][
                        "data"
                    ] = 10.0 * np.ma.log10(field_mean)
                else:
                    dscfg["global_data"]["radar_out"].fields[field_name][
                        "data"
                    ] = field_mean

            elif stat in ("std", "cov"):
                field_std = np.ma.sqrt(
                    dscfg["global_data"]["radar_out"].fields["sum_squared"]["data"]
                    / dscfg["global_data"]["radar_out"].fields["number_of_samples"][
                        "data"
                    ]
                    - field_mean * field_mean
                )
                if stat == "std":
                    if lin_trans:
                        dscfg["global_data"]["radar_out"].fields[field_name][
                            "data"
                        ] = 10.0 * np.ma.log10(field_std)
                    else:
                        dscfg["global_data"]["radar_out"].fields[field_name][
                            "data"
                        ] = field_std
                else:
                    dscfg["global_data"]["radar_out"].fields[field_name]["data"] = (
                        field_std / field_mean
                    )

        new_dataset = {
            "radar_out": deepcopy(dscfg["global_data"]["radar_out"]),
            "timeinfo": dscfg["global_data"]["endtime"],
        }

        return new_dataset, ind_rad


def process_time_stats2(procstatus, dscfg, radar_list=None):
    """
    computes the temporal mean of a field

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            Arbitrary data type supported by pyrad and contianed in the radar data
        period : float. Dataset keyword
            the period to average [s]. If -1 the statistics are going to be
            performed over the entire data. Default 3600.
        start_average : float. Dataset keyword
            when to start the average [s from midnight UTC]. Default 0.
        stat: string. Dataset keyword
            Statistic to compute: Can be median, mode, percentileXX
        use_nan : bool. Dataset keyword
            If true non valid data will be used
        nan_value : float. Dataset keyword
            The value of the non valid data. Default 0
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the statistic computed on the input field, as well as
        "nsamples"

    ind_rad : int
        radar index

    """
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_name = get_fieldname_pyart(datatype)
        break
    ind_rad = int(radarnr[5:8]) - 1

    start_average = dscfg.get("start_average", 0.0)
    period = dscfg.get("period", 3600.0)
    use_nan = dscfg.get("use_nan", 0)
    nan_value = dscfg.get("nan_value", 0.0)
    stat = dscfg.get("stat", "median")
    if "percentile" in stat:
        percentile = float(stat.replace("percentile", ""))

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn("No valid radar")
            return None, None
        radar = radar_list[ind_rad]

        if field_name not in radar.fields:
            warn(field_name + " not available.")
            return None, None

        # prepare auxiliary radar
        field = deepcopy(radar.fields[field_name])
        if use_nan:
            field["data"] = np.ma.asarray(field["data"].filled(nan_value))
        npoints_dict = pyart.config.get_metadata("number_of_samples")
        npoints_dict["data"] = np.ma.asarray(
            np.logical_not(np.ma.getmaskarray(field["data"])), dtype=int
        )

        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()
        radar_aux.add_field(field_name, field)
        radar_aux.add_field("number_of_samples", npoints_dict)

        # first volume: initialize start and end time of averaging
        if dscfg["initialized"] == 0:
            avg_par = dict()
            if period != -1:
                date_00 = dscfg["timeinfo"].replace(
                    hour=0, minute=0, second=0, microsecond=0
                )

                avg_par.update(
                    {"starttime": date_00 + datetime.timedelta(seconds=start_average)}
                )
                avg_par.update(
                    {
                        "endtime": avg_par["starttime"]
                        + datetime.timedelta(seconds=period)
                    }
                )
            else:
                avg_par.update({"starttime": dscfg["timeinfo"]})
                avg_par.update({"endtime": dscfg["timeinfo"]})

            avg_par.update({"timeinfo": dscfg["timeinfo"]})
            dscfg["global_data"] = avg_par
            dscfg["initialized"] = 1

        if dscfg["initialized"] == 0:
            return None, None

        dscfg["global_data"]["timeinfo"] = dscfg["timeinfo"]
        # no radar object in global data: create it
        if "radar_out" not in dscfg["global_data"]:
            if period != -1:
                # get start and stop times of new radar object
                (
                    dscfg["global_data"]["starttime"],
                    dscfg["global_data"]["endtime"],
                ) = time_avg_range(
                    dscfg["timeinfo"],
                    dscfg["global_data"]["starttime"],
                    dscfg["global_data"]["endtime"],
                    period,
                )

                # check if volume time older than starttime
                if dscfg["timeinfo"] > dscfg["global_data"]["starttime"]:
                    dscfg["global_data"].update({"radar_out": radar_aux})
                    dscfg["global_data"].update(
                        {
                            "field_data": np.atleast_3d(
                                radar_aux.fields[field_name]["data"]
                            )
                        }
                    )
            else:
                dscfg["global_data"].update({"radar_out": radar_aux})
                dscfg["global_data"].update(
                    {"field_data": np.atleast_3d(radar_aux.fields[field_name]["data"])}
                )

            return None, None

        # still accumulating: add field to global field
        if period == -1 or dscfg["timeinfo"] < dscfg["global_data"]["endtime"]:
            if period == -1:
                dscfg["global_data"]["endtime"] = dscfg["timeinfo"]

            field_interp = interpol_field(
                dscfg["global_data"]["radar_out"], radar_aux, field_name
            )
            npoints_interp = interpol_field(
                dscfg["global_data"]["radar_out"], radar_aux, "number_of_samples"
            )

            if use_nan:
                field_interp["data"] = np.ma.asarray(
                    field_interp["data"].filled(nan_value)
                )
                dscfg["global_data"]["radar_out"].fields["number_of_samples"][
                    "data"
                ] += np.ma.asarray(
                    npoints_interp["data"].filled(fill_value=1), dtype=int
                )
            else:
                dscfg["global_data"]["radar_out"].fields["number_of_samples"][
                    "data"
                ] += np.ma.asarray(
                    npoints_interp["data"].filled(fill_value=0), dtype=int
                )

            dscfg["global_data"]["field_data"] = np.ma.append(
                dscfg["global_data"]["field_data"],
                np.atleast_3d(field_interp["data"]),
                axis=2,
            )

            return None, None

        # we have reached the end of the accumulation period: do the averaging
        # and start a new object (only reachable if period != -1)
        if stat == "median":
            dscfg["global_data"]["radar_out"].fields[field_name]["data"] = np.ma.median(
                dscfg["global_data"]["field_data"], axis=2
            )
        elif stat == "mode":
            mode_data, _ = scipy.stats.mode(
                dscfg["global_data"]["field_data"].filled(fill_value=np.nan),
                axis=2,
                nan_policy="omit",
            )
            dscfg["global_data"]["radar_out"].fields[field_name][
                "data"
            ] = np.ma.masked_invalid(np.squeeze(mode_data, axis=2))
        elif "percentile" in stat:
            percent_data = np.nanpercentile(
                dscfg["global_data"]["field_data"].filled(fill_value=np.nan),
                percentile,
                axis=2,
            )
            dscfg["global_data"]["radar_out"].fields[field_name][
                "data"
            ] = np.ma.masked_invalid(percent_data)

        new_dataset = {
            "radar_out": deepcopy(dscfg["global_data"]["radar_out"]),
            "timeinfo": dscfg["global_data"]["endtime"],
        }

        dscfg["global_data"]["starttime"] += datetime.timedelta(seconds=period)
        dscfg["global_data"]["endtime"] += datetime.timedelta(seconds=period)

        # remove old radar object from global_data dictionary
        dscfg["global_data"].pop("radar_out", None)

        # get start and stop times of new radar object
        (
            dscfg["global_data"]["starttime"],
            dscfg["global_data"]["endtime"],
        ) = time_avg_range(
            dscfg["timeinfo"],
            dscfg["global_data"]["starttime"],
            dscfg["global_data"]["endtime"],
            period,
        )

        # check if volume time older than starttime
        if dscfg["timeinfo"] > dscfg["global_data"]["starttime"]:
            dscfg["global_data"].update({"radar_out": radar_aux})

        return new_dataset, ind_rad

    # no more files to process if there is global data pack it up
    if procstatus == 2:
        if dscfg["initialized"] == 0:
            return None, None
        if "radar_out" not in dscfg["global_data"]:
            return None, None

        if stat == "median":
            dscfg["global_data"]["radar_out"].fields[field_name]["data"] = np.ma.median(
                dscfg["global_data"]["field_data"], axis=2
            )
        elif stat == "mode":
            mode_data, _ = scipy.stats.mode(
                dscfg["global_data"]["field_data"].filled(fill_value=np.nan),
                axis=2,
                nan_policy="omit",
            )
            dscfg["global_data"]["radar_out"].fields[field_name][
                "data"
            ] = np.ma.masked_invalid(np.squeeze(mode_data, axis=2))
        elif "percentile" in stat:
            percent_data = np.nanpercentile(
                dscfg["global_data"]["field_data"].filled(fill_value=np.nan),
                percentile,
                axis=2,
            )
            dscfg["global_data"]["radar_out"].fields[field_name][
                "data"
            ] = np.ma.masked_invalid(percent_data)

        new_dataset = {
            "radar_out": deepcopy(dscfg["global_data"]["radar_out"]),
            "timeinfo": dscfg["global_data"]["endtime"],
        }

        return new_dataset, ind_rad


def process_time_avg(procstatus, dscfg, radar_list=None):
    """
    computes the temporal mean of a field

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            Arbitrary data type supported by pyrad and contained in the radar data
        period : float. Dataset keyword
            the period to average [s]. Default 3600.
        start_average : float. Dataset keyword
            when to start the average [s from midnight UTC]. Default 0.
        lin_trans: int. Dataset keyword
            If 1 apply linear transformation before averaging
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the statistic computed on the input field, as well as
        "nsamples"
    ind_rad : int
        radar index

    """
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_name = get_fieldname_pyart(datatype)
        break
    ind_rad = int(radarnr[5:8]) - 1

    lin_trans = dscfg.get("lin_trans", 0)

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn("No valid radar")
            return None, None
        radar = radar_list[ind_rad]

        if field_name not in radar.fields:
            warn(field_name + " not available.")
            return None, None

        period = dscfg.get("period", 3600.0)

        field = deepcopy(radar.fields[field_name])
        if lin_trans:
            field["data"] = np.ma.power(10.0, 0.1 * field["data"])

        field["data"] = field["data"].filled(fill_value=0.0)
        field["data"] = np.ma.asarray(field["data"])

        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()
        radar_aux.add_field(field_name, field)
        npoints_dict = pyart.config.get_metadata("number_of_samples")
        npoints_dict["data"] = np.ma.ones((radar.nrays, radar.ngates), dtype=int)
        radar_aux.add_field("number_of_samples", npoints_dict)

        # first volume: initialize start and end time of averaging
        if dscfg["initialized"] == 0:
            start_average = dscfg.get("start_average", 0.0)

            date_00 = dscfg["timeinfo"].replace(
                hour=0, minute=0, second=0, microsecond=0
            )

            avg_par = dict()
            avg_par.update(
                {"starttime": date_00 + datetime.timedelta(seconds=start_average)}
            )
            avg_par.update(
                {"endtime": avg_par["starttime"] + datetime.timedelta(seconds=period)}
            )
            avg_par.update({"timeinfo": dscfg["timeinfo"]})
            dscfg["global_data"] = avg_par
            dscfg["initialized"] = 1

        if dscfg["initialized"] == 0:
            return None, None

        dscfg["global_data"]["timeinfo"] = dscfg["timeinfo"]
        # no radar object in global data: create it
        if "radar_out" not in dscfg["global_data"]:
            # get start and stop times of new radar object
            (
                dscfg["global_data"]["starttime"],
                dscfg["global_data"]["endtime"],
            ) = time_avg_range(
                dscfg["timeinfo"],
                dscfg["global_data"]["starttime"],
                dscfg["global_data"]["endtime"],
                period,
            )

            # check if volume time older than starttime
            if dscfg["timeinfo"] > dscfg["global_data"]["starttime"]:
                dscfg["global_data"].update({"radar_out": radar_aux})

            return None, None

        # still accumulating: add field to global field
        if dscfg["timeinfo"] < dscfg["global_data"]["endtime"]:
            field_interp = interpol_field(
                dscfg["global_data"]["radar_out"], radar_aux, field_name
            )
            npoints_interp = interpol_field(
                dscfg["global_data"]["radar_out"], radar_aux, "number_of_samples"
            )
            dscfg["global_data"]["radar_out"].fields[field_name][
                "data"
            ] += field_interp["data"].filled(fill_value=0)
            dscfg["global_data"]["radar_out"].fields["number_of_samples"]["data"] += (
                npoints_interp["data"].filled(fill_value=0)
            ).astype("int")

            return None, None

        # we have reached the end of the accumulation period: do the averaging
        # and start a new object
        dscfg["global_data"]["radar_out"].fields[field_name]["data"] /= dscfg[
            "global_data"
        ]["radar_out"].fields["number_of_samples"]["data"]
        if lin_trans:
            dscfg["global_data"]["radar_out"].fields[field_name][
                "data"
            ] = 10.0 * np.ma.log10(
                dscfg["global_data"]["radar_out"].fields[field_name]["data"]
            )

        new_dataset = {
            "radar_out": deepcopy(dscfg["global_data"]["radar_out"]),
            "timeinfo": dscfg["global_data"]["endtime"],
        }

        dscfg["global_data"]["starttime"] += datetime.timedelta(seconds=period)
        dscfg["global_data"]["endtime"] += datetime.timedelta(seconds=period)

        # remove old radar object from global_data dictionary
        dscfg["global_data"].pop("radar_out", None)

        # get start and stop times of new radar object
        (
            dscfg["global_data"]["starttime"],
            dscfg["global_data"]["endtime"],
        ) = time_avg_range(
            dscfg["timeinfo"],
            dscfg["global_data"]["starttime"],
            dscfg["global_data"]["endtime"],
            period,
        )

        # check if volume time older than starttime
        if dscfg["timeinfo"] > dscfg["global_data"]["starttime"]:
            dscfg["global_data"].update({"radar_out": radar_aux})

        return new_dataset, ind_rad

    # no more files to process if there is global data pack it up
    if procstatus == 2:
        if dscfg["initialized"] == 0:
            return None, None
        if "radar_out" not in dscfg["global_data"]:
            return None, None

        (dscfg["global_data"]["radar_out"].fields[field_name]["data"]) /= dscfg[
            "global_data"
        ]["radar_out"].fields["number_of_samples"]["data"]
        if lin_trans:
            dscfg["global_data"]["radar_out"].fields[field_name][
                "data"
            ] = 10.0 * np.ma.log10(
                dscfg["global_data"]["radar_out"].fields[field_name]["data"]
            )

        new_dataset = {
            "radar_out": deepcopy(dscfg["global_data"]["radar_out"]),
            "timeinfo": dscfg["global_data"]["endtime"],
        }

        return new_dataset, ind_rad


def process_weighted_time_avg(procstatus, dscfg, radar_list=None):
    """
    computes the temporal mean of a field weighted by the reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            Arbitrary data type supported by pyrad and contained in the radar data, as well as
            "dBZ" or "dBZc" or "dBuZ" or "dBZv" or "dBZvc" or "dBuZv" (refl. weighting)
        period : float. Dataset keyword
            the period to average [s]. Default 3600.
        start_average : float. Dataset keyword
            when to start the average [s from midnight UTC]. Default 0.
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the statistic computed on the input field
    ind_rad : int
        radar index

    """
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("dBZ", "dBZc", "dBuZ", "dBZv", "dBZvc", "dBuZv"):
            refl_name = get_fieldname_pyart(datatype)
        else:
            field_name = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8]) - 1

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn("No valid radar")
            return None, None
        radar = radar_list[ind_rad]
        if field_name not in radar.fields or refl_name not in radar.fields:
            warn("Unable to compute weighted average. Missing data")
            return None, None

        period = dscfg.get("period", 3600.0)

        field = deepcopy(radar.fields[field_name])
        field["data"] = field["data"].filled(fill_value=0.0)
        field["data"] = np.ma.asarray(field["data"])

        refl_field = deepcopy(radar.fields[refl_name])
        refl_field["data"] = np.ma.power(10.0, 0.1 * refl_field["data"])
        refl_field["data"] = refl_field["data"].filled(fill_value=0.0)
        refl_field["data"] = np.ma.asarray(refl_field["data"])

        field["data"] *= refl_field["data"]

        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()
        radar_aux.add_field(field_name, field)
        radar_aux.add_field(refl_name, refl_field)

        # first volume: initialize start and end time of averaging
        if dscfg["initialized"] == 0:
            start_average = dscfg.get("start_average", 0.0)

            date_00 = dscfg["timeinfo"].replace(
                hour=0, minute=0, second=0, microsecond=0
            )

            avg_par = dict()
            avg_par.update(
                {"starttime": date_00 + datetime.timedelta(seconds=start_average)}
            )
            avg_par.update(
                {"endtime": avg_par["starttime"] + datetime.timedelta(seconds=period)}
            )
            avg_par.update({"timeinfo": dscfg["timeinfo"]})
            dscfg["global_data"] = avg_par
            dscfg["initialized"] = 1

        if dscfg["initialized"] == 0:
            return None, None

        dscfg["global_data"]["timeinfo"] = dscfg["timeinfo"]
        # no radar object in global data: create it
        if "radar_out" not in dscfg["global_data"]:
            # get start and stop times of new radar object
            (
                dscfg["global_data"]["starttime"],
                dscfg["global_data"]["endtime"],
            ) = time_avg_range(
                dscfg["timeinfo"],
                dscfg["global_data"]["starttime"],
                dscfg["global_data"]["endtime"],
                period,
            )

            # check if volume time older than starttime
            if dscfg["timeinfo"] > dscfg["global_data"]["starttime"]:
                dscfg["global_data"].update({"radar_out": radar_aux})

            return None, None

        # still accumulating: add field to global field
        if dscfg["timeinfo"] < dscfg["global_data"]["endtime"]:
            field_interp = interpol_field(
                dscfg["global_data"]["radar_out"], radar_aux, field_name
            )
            dscfg["global_data"]["radar_out"].fields[field_name][
                "data"
            ] += field_interp["data"].filled(fill_value=0)

            refl_interp = interpol_field(
                dscfg["global_data"]["radar_out"], radar_aux, refl_name
            )
            dscfg["global_data"]["radar_out"].fields[refl_name]["data"] += refl_interp[
                "data"
            ].filled(fill_value=0)

            return None, None

        # we have reached the end of the accumulation period: do the averaging
        # and start a new object
        dscfg["global_data"]["radar_out"].fields[field_name]["data"] /= dscfg[
            "global_data"
        ]["radar_out"].fields[refl_name]["data"]

        new_dataset = {
            "radar_out": deepcopy(dscfg["global_data"]["radar_out"]),
            "timeinfo": dscfg["global_data"]["endtime"],
        }

        dscfg["global_data"]["starttime"] += datetime.timedelta(seconds=period)
        dscfg["global_data"]["endtime"] += datetime.timedelta(seconds=period)

        # remove old radar object from global_data dictionary
        dscfg["global_data"].pop("radar_out", None)

        # get start and stop times of new radar object
        (
            dscfg["global_data"]["starttime"],
            dscfg["global_data"]["endtime"],
        ) = time_avg_range(
            dscfg["timeinfo"],
            dscfg["global_data"]["starttime"],
            dscfg["global_data"]["endtime"],
            period,
        )

        # check if volume time older than starttime
        if dscfg["timeinfo"] > dscfg["global_data"]["starttime"]:
            dscfg["global_data"].update({"radar_out": radar_aux})

        return new_dataset, ind_rad

    # no more files to process if there is global data pack it up
    if procstatus == 2:
        if dscfg["initialized"] == 0:
            return None, None
        if "radar_out" not in dscfg["global_data"]:
            return None, None

        dscfg["global_data"]["radar_out"].fields[field_name]["data"] /= dscfg[
            "global_data"
        ]["radar_out"].fields[refl_name]["data"]

        new_dataset = {
            "radar_out": deepcopy(dscfg["global_data"]["radar_out"]),
            "timeinfo": dscfg["global_data"]["endtime"],
        }

        return new_dataset, ind_rad


def process_time_avg_flag(procstatus, dscfg, radar_list=None):
    """
    computes a flag field describing the conditions of the data used while
    averaging. The flag is an integer that tracks up to 999 occurrences of the number of samples as well as
    of three conditions during data accumulation: 𝜙𝑑𝑝 exceeding a threshold (𝜙𝑑𝑝ₘₐₓ), clutter,
    and non-rain precipitation. It is a packed representation of 4 numbers into one large integer.
    Flags are encoded by adding +1 for nsamples, +1E3 for 𝜙𝑑𝑝ₘₐₓ exceedance, + 1E6 for clutter,
    and +1E9 for non-rain.
    Inputs include the 𝜙𝑑𝑝 field, a dynamic clutter map, and either a hydrometeor
    classification or temperature field to identify precipitation phase.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types, must be
            "PhiDP" or "PhiDPc" (Optional, for PhiDP flagging), and,
            "echoID" (Optional, for echoID flagging), and,
            "hydro" (Optional, for no rain flagging), and,
            "TEMP" (Optional, for solid precip flagging), and,
            "H_ISO0" (Optional, also for solid precip flagging)
        period : float. Dataset keyword
            the period to average [s]. Default 3600.
        start_average : float. Dataset keyword
            when to start the average [s from midnight UTC]. Default 0.
        phidpmax: float. Dataset keyword
            maximum PhiDP
        beamwidth : float. Dataset keyword
            the antenna beamwidth [deg]. If None that of the keys
            radar_beam_width_h or radar_beam_width_v in attribute
            instrument_parameters of the radar object will be used. If the key
            or the attribute are not present the beamwidth will be set to None
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the field "time_avg_flag"
    ind_rad : int
        radar index

    """
    temp_name = None
    hydro_name = None
    iso0_name = None
    echo_name = None
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ("PhiDP", "PhiDPc"):
            phidp_name = get_fieldname_pyart(datatype)
        elif datatype == "echoID":
            echo_name = get_fieldname_pyart(datatype)
        elif datatype == "hydro":
            hydro_name = get_fieldname_pyart(datatype)
        elif datatype == "TEMP":
            temp_name = get_fieldname_pyart(datatype)
        elif datatype == "H_ISO0":
            iso0_name = "height_over_iso0"

    ind_rad = int(radarnr[5:8]) - 1

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn("No valid radar")
            return None, None
        radar = radar_list[ind_rad]

        phidpmax = dscfg.get("phidpmax", 60.0)
        period = dscfg.get("period", 3600.0)

        time_avg_flag = pyart.config.get_metadata("time_avg_flag")
        time_avg_flag["data"] = np.ma.zeros((radar.nrays, radar.ngates), dtype=int)

        time_avg_flag["data"] += 1  # number of samples

        if phidp_name not in radar.fields:
            warn("Missing PhiDP data")
            time_avg_flag["data"] += int(1e3)
        else:
            phidp_field = radar.fields[phidp_name]
            time_avg_flag["data"][phidp_field["data"] > phidpmax] += int(1e3)

        if echo_name is not None:
            if echo_name not in radar.fields:
                warn("Missing echo ID data")
                time_avg_flag["data"] += int(1e6)
            else:
                echo_field = radar.fields[echo_name]
                time_avg_flag["data"][echo_field["data"] == 2] += int(1e6)

        if hydro_name is not None and echo_name is not None:
            if (hydro_name not in radar.fields) or (echo_name not in radar.fields):
                warn("Missing hydrometeor classification data")
                time_avg_flag["data"] += int(1e9)
            else:
                hydro_field = radar.fields[hydro_name]
                # check where is no rain
                is_not_rain = np.logical_and(
                    hydro_field["data"] != 4, hydro_field["data"] != 6
                )
                # where is no rain should be precip
                is_not_rain = np.logical_and(is_not_rain, echo_field["data"] == 3)
                time_avg_flag["data"][is_not_rain] += int(1e9)
        elif temp_name is not None:
            if temp_name not in radar.fields:
                warn("Missing temperature data")
                time_avg_flag["data"] += int(1e9)
            else:
                beamwidth = dscfg.get("beamwidth", None)
                if beamwidth is None:
                    if radar.instrument_parameters is not None:
                        if "radar_beam_width_h" in radar.instrument_parameters:
                            beamwidth = radar.instrument_parameters[
                                "radar_beam_width_h"
                            ]["data"][0]
                        elif "radar_beam_width_v" in radar.instrument_parameters:
                            beamwidth = radar.instrument_parameters[
                                "radar_beam_width_v"
                            ]["data"][0]
                if beamwidth is None:
                    warn("Antenna beam width unknown.")

                mask_fzl, _ = pyart.correct.get_mask_fzl(
                    radar,
                    fzl=None,
                    doc=None,
                    min_temp=0.0,
                    max_h_iso0=0.0,
                    thickness=700.0,
                    beamwidth=beamwidth,
                    temp_field=temp_name,
                    iso0_field=iso0_name,
                    temp_ref="temperature",
                )
                time_avg_flag["data"][mask_fzl] += int(1e6)
        elif iso0_name is not None:
            if iso0_name not in radar.fields:
                warn("Missing height relative to iso0 data")
                time_avg_flag["data"] += int(1e6)
            else:
                beamwidth = dscfg.get("beamwidth", None)
                if beamwidth is None:
                    if radar.instrument_parameters is not None:
                        if "radar_beam_width_h" in radar.instrument_parameters:
                            beamwidth = radar.instrument_parameters[
                                "radar_beam_width_h"
                            ]["data"][0]
                        elif "radar_beam_width_v" in radar.instrument_parameters:
                            beamwidth = radar.instrument_parameters[
                                "radar_beam_width_v"
                            ]["data"][0]
                if beamwidth is None:
                    warn("Antenna beam width unknown.")

                mask_fzl, _ = pyart.correct.get_mask_fzl(
                    radar,
                    fzl=None,
                    doc=None,
                    min_temp=0.0,
                    max_h_iso0=0.0,
                    thickness=700.0,
                    beamwidth=beamwidth,
                    temp_field=temp_name,
                    iso0_field=iso0_name,
                    temp_ref="height_over_iso0",
                )
                time_avg_flag["data"][mask_fzl] += int(1e9)

        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()
        radar_aux.add_field("time_avg_flag", time_avg_flag)

        # first volume: initialize start and end time of averaging
        if dscfg["initialized"] == 0:
            start_average = dscfg.get("start_average", 0.0)

            date_00 = dscfg["timeinfo"].replace(
                hour=0, minute=0, second=0, microsecond=0
            )

            avg_par = dict()
            avg_par.update(
                {"starttime": date_00 + datetime.timedelta(seconds=start_average)}
            )
            avg_par.update(
                {"endtime": avg_par["starttime"] + datetime.timedelta(seconds=period)}
            )
            avg_par.update({"timeinfo": dscfg["timeinfo"]})
            dscfg["global_data"] = avg_par
            dscfg["initialized"] = 1

        if dscfg["initialized"] == 0:
            return None, None

        dscfg["global_data"]["timeinfo"] = dscfg["timeinfo"]
        # no radar object in global data: create it
        if "radar_out" not in dscfg["global_data"]:
            # get start and stop times of new radar object
            (
                dscfg["global_data"]["starttime"],
                dscfg["global_data"]["endtime"],
            ) = time_avg_range(
                dscfg["timeinfo"],
                dscfg["global_data"]["starttime"],
                dscfg["global_data"]["endtime"],
                period,
            )

            # check if volume time older than starttime
            if dscfg["timeinfo"] > dscfg["global_data"]["starttime"]:
                dscfg["global_data"].update({"radar_out": radar_aux})

            return None, None

        # still accumulating: add field to global field
        if dscfg["timeinfo"] < dscfg["global_data"]["endtime"]:
            flag_interp = interpol_field(
                dscfg["global_data"]["radar_out"], radar_aux, "time_avg_flag"
            )
            dscfg["global_data"]["radar_out"].fields["time_avg_flag"]["data"] += (
                flag_interp["data"].filled(fill_value=0)
            ).astype(int)

            return None, None

        # we have reached the end of the accumulation: start a new object
        new_dataset = {
            "radar_out": deepcopy(dscfg["global_data"]["radar_out"]),
            "timeinfo": dscfg["global_data"]["endtime"],
        }

        dscfg["global_data"]["starttime"] += datetime.timedelta(seconds=period)
        dscfg["global_data"]["endtime"] += datetime.timedelta(seconds=period)

        # remove old radar object from global_data dictionary
        dscfg["global_data"].pop("radar_out", None)

        # get start and stop times of new radar object
        (
            dscfg["global_data"]["starttime"],
            dscfg["global_data"]["endtime"],
        ) = time_avg_range(
            dscfg["timeinfo"],
            dscfg["global_data"]["starttime"],
            dscfg["global_data"]["endtime"],
            period,
        )

        # check if volume time older than starttime
        if dscfg["timeinfo"] > dscfg["global_data"]["starttime"]:
            dscfg["global_data"].update({"radar_out": radar_aux})

        return new_dataset, ind_rad

    # no more files to process if there is global data pack it up
    if procstatus == 2:
        if dscfg["initialized"] == 0:
            return None, None
        if "radar_out" not in dscfg["global_data"]:
            return None, None

        new_dataset = {
            "radar_out": deepcopy(dscfg["global_data"]["radar_out"]),
            "timeinfo": dscfg["global_data"]["endtime"],
        }

        return new_dataset, ind_rad
