"""
pyrad.io.read_data_other
========================

Functions for reading auxiliary data

.. autosummary::
    :toctree: generated/

    read_centroids_npz
    read_centroids
    read_proc_periods
    read_profile_ts
    read_histogram_ts
    read_quantiles_ts
    read_vpr_theo_parameters
    read_rhi_profile
    read_last_state
    read_status
    read_mch_xml_vad
    read_rad4alp_icon
    read_rad4alp_vis
    read_mf_vis
    read_histogram
    read_quantiles
    read_excess_gates
    read_colocated_gates
    read_colocated_data
    read_colocated_data_time_avg
    read_timeseries
    read_ts_cum
    read_ml_ts
    read_monitoring_ts
    read_monitoring_ts_old
    read_intercomp_scores_ts
    read_intercomp_scores_ts_old
    read_intercomp_scores_ts_old_v0
    read_selfconsistency
    read_antenna_pattern
    read_xls

"""

import os
import glob
import datetime

import pandas as pd
import csv
import xml.etree.ElementTree as et
from warnings import warn

try:
    import fcntl

    FCNTL_AVAIL = True
except ModuleNotFoundError:
    import msvcrt  # For windows

    FCNTL_AVAIL = False
import time
import errno

import numpy as np

import pyart
from pyart.config import get_fillvalue, get_metadata

from .io_aux import get_fieldname_pyart, _get_datetime
from .flock_utils import lock_file, unlock_file


def read_centroids_npz(fname):
    """
    Reads an npz file containing data used in centroid computation

    Parameters
    ----------
    fname : str
        name of the file to read

    Returns
    -------
    centroids_data : dict
        Dictionary with the data

    """
    try:
        with np.load(fname, allow_pickle=True) as npzfile:
            return npzfile["centroids_data"].item()
    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None


def read_centroids(fname):
    """
    Reads a file containing the centroids for the semi-supervised
    classification

    Parameters
    ----------
    fname : str
        name of the file to read

    Returns
    -------
    mass_centers : 2D-array
        Matrix with the centroids as (nclasses, nvariables)

    """
    try:
        with open(fname, "r", newline="") as csvfile:
            # first count the lines
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            nhydro = sum(1 for row in reader)
            nvars = len(reader.fieldnames) - 1
            mass_centers = np.empty((nhydro, nvars))
            hydro_names = []
            var_names = reader.fieldnames[1:]

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            for i, row in enumerate(reader):
                hydro_names.append(row[reader.fieldnames[0]])
                for j, col in enumerate(var_names):
                    mass_centers[i, j] = float(row[col])

            return mass_centers, hydro_names, var_names
    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None, None, None


def read_proc_periods(fname):
    """
    Reads a file containing the start and stop times of periods to process

    Parameters
    ----------
    fname : str
        name of the file to read

    Returns
    -------
    starttimes, endtimes : array of datetime objects or None
        The start and end times of the periods to process if the reading has
        been successful

    """
    try:
        with open(fname, "r", newline="") as csvfile:
            # first count the lines
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            nrows = sum(1 for row in reader)
            startimes = np.empty(nrows, dtype=datetime.datetime)
            endtimes = np.empty(nrows, dtype=datetime.datetime)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            for i, row in enumerate(reader):
                startimes[i] = datetime.datetime.strptime(
                    row["start_time"], "%Y%m%d%H%M%S"
                )
                endtimes[i] = datetime.datetime.strptime(
                    row["end_time"], "%Y%m%d%H%M%S"
                )

            return startimes, endtimes
    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None, None


def read_profile_ts(fname_list, labels, hres=None, label_nr=0, t_res=300.0):
    """
    Reads a colection of profile data file and creates a time series

    Parameters
    ----------
    fname_list : str
        list of files to read
    labels : list of str
        The data labels
    hres : float
        Height resolution
    label_nr : int
        the label nr of the data that will be used in the time series
    t_res : float
        time resolution [s]. If None the time resolution is taken as the
        median

    Returns
    -------
    tbin_edges, hbin_edges, np_ma, data_ma, datetime_arr : tupple
        The read data. None otherwise

    """
    data_ma = []
    np_ma = []
    datetime_arr = np.ma.array([], dtype=datetime.datetime)
    for fname in fname_list:
        datetime_arr = np.append(datetime_arr, _get_datetime(fname, "RAINBOW"))
        height, np_t, vals = read_rhi_profile(fname, labels)
        if hres is None:
            hres = np.mean(height[1:] - height[:-1])
        hbin_edges = np.append(height - hres / 2, height[-1] + hres / 2)
        val = vals[:, label_nr]
        data_ma.append(val)
        np_ma.append(np_t)
    data_ma = np.ma.asarray(data_ma)
    np_ma = np.asarray(np_ma, dtype=int)

    # sort data as a function of time
    ind = np.argsort(datetime_arr)
    datetime_arr = datetime_arr[ind]
    data_ma = data_ma[ind, :]
    np_ma = np_ma[ind, :]

    # put date time array as seconds from start of TRT cell
    dt_s = np.empty(datetime_arr.size, dtype=float)
    for j, dt in enumerate(datetime_arr):
        dt_s[j] = (dt - datetime_arr[0]).total_seconds()

    if t_res is None:
        if dt_s.size > 1:
            t_res = np.median(dt_s[1:] - dt_s[:-1])
        else:
            t_res = 300.0
    tbin_edges = np.append(dt_s - t_res, dt_s[-1])

    return tbin_edges, hbin_edges, np_ma, data_ma, datetime_arr


def read_histogram_ts(fname_list, datatype, t_res=300.0):
    """
    Reads a colection of histogram data file and creates a time series

    Parameters
    ----------
    fname_list : str
        list of files to read
    datatype : str
        The data type (dBZ, ZDR, etc.)
    t_res : float
        time resolution [s]. If None the time resolution is taken as the
        median

    Returns
    -------
    tbin_edges, bin_edges, data_ma, datetime_arr : tupple
        The read data. None otherwise

    """
    data_ma = []
    datetime_arr = np.ma.array([], dtype=datetime.datetime)
    for fname in fname_list:
        datetime_arr = np.append(datetime_arr, _get_datetime(fname, "RAINBOW"))
        hist, bin_edges = read_histogram(fname)
        if datatype == "hydro":
            hist[0] = 0
        elif datatype == "RhoHVc":
            hist = hist[bin_edges[0:-1] > 0.95]
            bin_edges = bin_edges[bin_edges > 0.95]
        npoints = np.sum(hist)
        if npoints > 0:
            data_ma.append(hist / np.sum(hist) * 100.0)
        else:
            data_ma.append(hist)
    data_ma = np.ma.asarray(data_ma)

    # sort data as a function of time
    ind = np.argsort(datetime_arr)
    datetime_arr = datetime_arr[ind]
    data_ma = data_ma[ind, :]

    # put date time array as seconds from start of TRT cell
    dt_s = np.empty(datetime_arr.size, dtype=float)
    for j, dt in enumerate(datetime_arr):
        dt_s[j] = (dt - datetime_arr[0]).total_seconds()

    if t_res is None:
        if dt_s.size > 1:
            t_res = np.median(dt_s[1:] - dt_s[:-1])
        else:
            t_res = 300.0
    tbin_edges = np.append(dt_s - t_res, dt_s[-1])

    return tbin_edges, bin_edges, data_ma, datetime_arr


def read_quantiles_ts(fname_list, step=5.0, qmin=0.0, qmax=100.0, t_res=300.0):
    """
    Reads a colection of quantiles data file and creates a time series

    Parameters
    ----------
    fname_list : str
        list of files to read
    step, qmin, qmax : float
        The minimum, maximum and step quantiles
    t_res : float
        time resolution [s]. If None the time resolution is taken as the
        median

    Returns
    -------
    tbin_edges, qbin_edges, data_ma, datetime_arr : tupple
        The read data. None otherwise

    """
    data_ma = []
    datetime_arr = np.ma.array([], dtype=datetime.datetime)
    qbin_edges = np.arange(
        qmin - step / 2, qmax + step / 2 + step / 2, step=step, dtype=float
    )
    qbin_centers = np.arange(qmin, qmax + step / 2, step=step)
    for fname in fname_list:
        values_aux = np.ma.masked_all(qbin_edges.size - 1, dtype=float)
        datetime_arr = np.append(datetime_arr, _get_datetime(fname, "RAINBOW"))

        quantiles, values = read_quantiles(fname)
        if quantiles is not None:
            for i, qbin_center in enumerate(qbin_centers):
                val_aux = values[quantiles == qbin_center]
                if val_aux.size > 0:
                    values_aux[i] = val_aux
        data_ma.append(values_aux)
    data_ma = np.ma.asarray(data_ma)

    # sort data as a function of time
    ind = np.argsort(datetime_arr)
    datetime_arr = datetime_arr[ind]
    data_ma = data_ma[ind, :]

    # put date time array as seconds from start of TRT cell
    dt_s = np.empty(datetime_arr.size, dtype=float)
    for j, dt in enumerate(datetime_arr):
        dt_s[j] = (dt - datetime_arr[0]).total_seconds()

    if t_res is None:
        if dt_s.size > 1:
            t_res = np.median(dt_s[1:] - dt_s[:-1])
        else:
            t_res = 300.0
    tbin_edges = np.append(dt_s - t_res, dt_s[-1])

    return tbin_edges, qbin_edges, data_ma, datetime_arr


def read_vpr_theo_parameters(fname):
    """
    Reads the parameters defining a theoretical VPR profile from a csv file

    Parameters
    ----------
    fname : str
        file name where to store the data

    Returns
    -------
    vpr_theo_dict : dict
        Dictionary containing the parameters

    """
    try:
        with open(fname, "r", newline="") as csvfile:
            reader = csv.DictReader(csvfile)
            vpr_theo_dict = next(reader)
            for k, v in vpr_theo_dict.items():
                vpr_theo_dict[k] = float(v)
            return vpr_theo_dict

    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None


def read_rhi_profile(
    fname, labels=["50.0-percentile", "25.0-percentile", "75.0-percentile"]
):
    """
    Reads a monitoring time series contained in a csv file

    Parameters
    ----------
    fname : str
        path of time series file
    labels : list of str
        The data labels

    Returns
    -------
    height, np_t, vals : tupple
        The read data. None otherwise

    """
    try:
        with open(fname, "r", newline="") as csvfile:
            nfields = len(labels)
            # first count the lines
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            nrows = sum(1 for row in reader)

            height = np.empty(nrows, dtype=float)
            np_t = np.zeros(nrows, dtype=int)
            vals = np.ma.empty((nrows, nfields), dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            for i, row in enumerate(reader):
                height[i] = float(row["Altitude [m MSL]"])
                if "N valid" in row:
                    np_t[i] = int(row["N valid"])
                for j, label in enumerate(labels):
                    vals[i, j] = float(row[label])

            vals = np.ma.masked_values(vals, get_fillvalue())

            return height, np_t, vals
    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None, None, None


def read_last_state(fname):
    """
    Reads a file containing the date of acquisition of the last volume
    processed

    Parameters
    ----------
    fname : str
        name of the file to read

    Returns
    -------
    last_state : datetime object
        the date

    """
    try:
        with open(fname, "r", newline="") as txtfile:
            line = txtfile.readline()
            txtfile.close()
            if not line:
                warn("File " + fname + " is empty.")
                return None
            try:
                last_state = datetime.datetime.strptime(line, "%Y%m%d%H%M%S")
                return last_state
            except ValueError:
                warn("File " + fname + " does not contain a valid date.")
                return None
    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None


def read_status(voltime, cfg, ind_rad=0):
    """
    Reads rad4alp xml status file.

    Parameters
    ----------
    voltime : datetime object
        volume scan time

    cfg: dictionary of dictionaries
        configuration info to figure out where the data is
    ind_rad: int
        radar index

    Returns
    -------
    root : root element object
        The information contained in the status file

    """
    if cfg["RadarName"] is None:
        raise ValueError(
            "ERROR: Radar Name not specified in config file. \
            Unable to read status data"
        )

    dayinfo = voltime.strftime("%y%j")
    timeinfo = voltime.strftime("%H%M")
    basename = "ST" + cfg["RadarName"][ind_rad] + dayinfo
    if cfg["path_convention"][ind_rad] == "RT":
        datapath = cfg["datapath"][ind_rad] + "ST" + cfg["RadarName"][ind_rad] + "/"
    else:
        datapath = cfg["datapath"][ind_rad] + dayinfo + "/" + basename + "/"
    filename = glob.glob(datapath + basename + timeinfo + "*.xml")
    if not filename:
        warn(
            "rad4alp status file " + datapath + basename + timeinfo + "*.xml not found"
        )
        return None
    try:
        root = et.parse(filename[0]).getroot()
        return root
    except et.ParseError as ee:
        warn(str(ee))
        warn("Unable to read file " + filename[0])

        return None


def read_mch_xml_vad(voltime, cfg, ind_rad=0):
    """
    Reads a MeteoSwiss VAD file (VAx xml files)

    Parameters
    ----------
    voltime : datetime object
        volume scan time

    cfg: dictionary of dictionaries
        configuration info to figure out where the data is
    ind_rad: int
        radar index

    Returns
    -------
    pyart.core.HorizontalWindProfile
        pyart Wind profile object containing the VAD data

    """

    if cfg["RadarName"] is None:
        raise ValueError(
            "ERROR: Radar Name not specified in config file. \
            Unable to read status data"
        )

    dayinfo = voltime.strftime("%y%j")
    timeinfo = voltime.strftime("%H%M")
    basename = "VA" + cfg["RadarName"][ind_rad] + dayinfo
    if cfg["path_convention"][ind_rad] == "RT":
        datapath = cfg["datapath"][ind_rad] + "VA" + cfg["RadarName"][ind_rad] + "/"
    else:
        datapath = cfg["datapath"][ind_rad] + dayinfo + "/" + basename + "/"
    filename = glob.glob(datapath + basename + timeinfo + "*.xml")

    return pyart.aux_io.read_mch_vad(filename)


def read_rad4alp_icon(fname, datatype, ngates=0):
    """
    Reads rad4alp ICON data binary file.

    Parameters
    ----------
    fname : str
        name of the file to read
    datatype : str
        name of the data type
    ngates : int
        maximum number of range gates per ray. If larger than 0 the radar
        field will be cut accordingly.

    Returns
    -------
    field : dictionary
        The data field

    """
    try:
        bindata = np.fromfile(fname, dtype="uint8", count=-1)
        nbins = bindata.size
        naz = 360
        nr = int(nbins / naz)
        bindata = np.reshape(bindata, (naz, nr))
        mask = bindata == 0

        if datatype == "TEMP":
            field_name = get_fieldname_pyart(datatype)
            field_data = np.ma.masked_where(
                mask, (bindata - 1).astype(float) * 0.5 - 87.0
            )

            field = get_metadata(field_name)
            field["data"] = field_data
            if ngates > 0:
                field["data"] = field["data"][:, :ngates]
            return field
        if datatype == "ISO0":
            field_name = get_fieldname_pyart(datatype)
            field_data = np.ma.masked_where(mask, (bindata - 1).astype(float))

            field = get_metadata(field_name)
            field["data"] = field_data
            if ngates > 0:
                field["data"] = field["data"][:, :ngates]
            return field

        warn("Unknown COSMO data type " + datatype)
        return None
    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None


def read_rad4alp_vis(fname, datatype):
    """
    Reads rad4alp visibility data binary file.

    Parameters
    ----------
    fname : str
        name of the file to read

    datatype : str
        name of the data type

    Returns
    -------
    field_list : list of dictionaries
        A data field. Each element of the list corresponds to one elevation

    """
    if datatype != "VIS":
        warn("Unknown DEM data type " + datatype)
        return None

    header_size = 64
    naz = 360
    nrngs = [
        492,
        420,
        492,
        324,
        366,
        324,
        292,
        324,
        280,
        242,
        222,
        200,
        174,
        150,
        124,
        100,
        82,
        68,
        60,
        54,
    ]

    try:
        with open(fname, "rb") as visfile:
            field_list = list()
            # skip header
            visfile.seek(header_size, os.SEEK_SET)

            bindata = np.fromfile(visfile, dtype="uint8", count=-1)

            sweep_start_index = 0
            for nrng in nrngs:
                sweep_end_index = sweep_start_index + naz * nrng - 1
                field_data = np.reshape(
                    bindata[sweep_start_index : sweep_end_index + 1], (naz, nrng)
                ).astype(float)
                sweep_start_index = sweep_end_index + 1

                field_name = get_fieldname_pyart(datatype)
                field = get_metadata(field_name)
                field["data"] = field_data
                field_list.append(field)

            visfile.close()

            return field_list

    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None


def read_mf_vis(fname, datatype):
    """
    Reads Météo-France visibility data binary file.

    Parameters
    ----------
    fname : str
        name of the file to read

    datatype : str
        name of the data type

    Returns
    -------
    field : dict
        A data field

    """
    datatype = datatype.strip()
    if datatype != "VIS":
        warn("Unknown DEM data type " + datatype)
        return None

    try:
        with open(fname, "rb") as visfile:
            if os.path.getsize(fname) == 262144:
                bindata = np.fromfile(visfile, dtype="uint8", count=-1)
                bindata = np.reshape(bindata, (512, 512)).astype(float)
                field_data = np.ma.masked_values(bindata, 255)
                if "coefmasqZ" not in fname:
                    field_data = 100.0 * np.ma.power(field_data / 100.0, 1.6)
                field_data = 100.0 - np.broadcast_to(field_data[::-1, :], (1, 512, 512))
            else:
                bindata = 100.0 - np.fromfile(visfile, dtype="float32", count=-1)
                field_data = np.reshape(bindata, (720, 270))

            field_name = get_fieldname_pyart(datatype)
            field = get_metadata(field_name)
            field["data"] = field_data

            visfile.close()

            return field

    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None


def read_histogram(fname):
    """
    Reads a histogram contained in a csv file

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    hist , bin_edges : tupple
        The read data. None otherwise

    """
    try:
        with open(fname, "r", newline="") as csvfile:
            # first count the lines
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            nrows = sum(1 for row in reader)
            bin_edges = np.zeros(nrows + 1, dtype=float)
            hist = np.zeros(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            for i, row in enumerate(reader):
                bin_edges[i] = float(row["bin_edge_left"])
                bin_edges[i + 1] = float(row["bin_edge_right"])
                hist[i] = float(row["value"])

            return hist, bin_edges
    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None, None


def read_quantiles(fname):
    """
    Reads quantiles contained in a csv file

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    quantiles, values : tupple
        The read data. None otherwise

    """
    try:
        with open(fname, "r", newline="") as csvfile:
            # first count the lines
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            nrows = sum(1 for row in reader)
            quantiles = np.zeros(nrows, dtype=float)
            values = np.zeros(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            for i, row in enumerate(reader):
                quantiles[i] = float(row["quantile"])
                values[i] = float(row["value"])

            return quantiles, values
    except (EnvironmentError, ValueError) as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None, None


def read_excess_gates(fname):
    """
    Reads a csv files containing the position of gates exceeding
    a given percentile of frequency of occurrence

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi, rad1_rng,
    rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi, rad2_rng : tupple
        A tupple with the data read. None otherwise

    """
    try:
        with open(fname, "r", newline="") as csvfile:
            # first count the lines
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            nrows = sum(1 for row in reader)
            ray_ind = np.empty(nrows, dtype=int)
            rng_ind = np.empty(nrows, dtype=int)
            ele = np.empty(nrows, dtype=float)
            azi = np.empty(nrows, dtype=float)
            rng = np.empty(nrows, dtype=float)
            nsamples = np.empty(nrows, dtype=int)
            occurrence = np.empty(nrows, dtype=int)
            freq_occu = np.empty(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            for i, row in enumerate(reader):
                ray_ind[i] = int(row["ray_ind"])
                rng_ind[i] = int(row["rng_ind"])
                ele[i] = float(row["ele"])
                azi[i] = float(row["azi"])
                rng[i] = float(row["rng"])
                nsamples[i] = int(row["nsamples"])
                occurrence[i] = int(row["occurrence"])
                freq_occu[i] = float(row["freq_occu"])

            csvfile.close()

            return (ray_ind, rng_ind, ele, azi, rng, nsamples, occurrence, freq_occu)
    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None, None, None, None, None, None, None, None


def read_colocated_gates(fname):
    """
    Reads a csv files containing the position of colocated gates

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi, rad1_rng,
    rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi, rad2_rng : tupple
        A tupple with the data read. None otherwise

    """
    try:
        with open(fname, "r", newline="") as csvfile:
            # first count the lines
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            nrows = sum(1 for row in reader)
            rad1_ray_ind = np.empty(nrows, dtype=int)
            rad1_rng_ind = np.empty(nrows, dtype=int)
            rad1_ele = np.empty(nrows, dtype=float)
            rad1_azi = np.empty(nrows, dtype=float)
            rad1_rng = np.empty(nrows, dtype=float)
            rad2_ray_ind = np.empty(nrows, dtype=int)
            rad2_rng_ind = np.empty(nrows, dtype=int)
            rad2_ele = np.empty(nrows, dtype=float)
            rad2_azi = np.empty(nrows, dtype=float)
            rad2_rng = np.empty(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            for i, row in enumerate(reader):
                rad1_ray_ind[i] = int(row["rad1_ray_ind"])
                rad1_rng_ind[i] = int(row["rad1_rng_ind"])
                rad1_ele[i] = float(row["rad1_ele"])
                rad1_azi[i] = float(row["rad1_azi"])
                rad1_rng[i] = float(row["rad1_rng"])
                rad2_ray_ind[i] = int(row["rad2_ray_ind"])
                rad2_rng_ind[i] = int(row["rad2_rng_ind"])
                rad2_ele[i] = float(row["rad2_ele"])
                rad2_azi[i] = float(row["rad2_azi"])
                rad2_rng[i] = float(row["rad2_rng"])

            csvfile.close()

            return (
                rad1_ray_ind,
                rad1_rng_ind,
                rad1_ele,
                rad1_azi,
                rad1_rng,
                rad2_ray_ind,
                rad2_rng_ind,
                rad2_ele,
                rad2_azi,
                rad2_rng,
            )
    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None, None, None, None, None, None, None, None, None, None


def read_colocated_data(fname):
    """
    Reads a CSV file containing colocated data using pandas.

    Parameters
    ----------
    fname : str
        Path of time series file.

    Returns
    -------
    tuple
        A tuple with the data read:
        (rad1_time, rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi, rad1_rng,
        rad1_val, rad2_time, rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi,
        rad2_rng, rad2_val). Returns None for all if unable to read the file.
    """
    try:
        # Read the CSV file, skipping lines starting with "#"
        df = pd.read_csv(fname, comment="#", na_values="--")

        # Convert rad1_time and rad2_time to datetime
        df["rad1_time"] = pd.to_datetime(df["rad1_time"], format="%Y%m%d%H%M%S")
        df["rad2_time"] = pd.to_datetime(df["rad2_time"], format="%Y%m%d%H%M%S")

        # Create masked arrays
        def to_masked_array(column, dtype=None):
            return np.ma.masked_invalid(df[column].to_numpy(dtype=dtype))

        return (
            to_masked_array("rad1_time"),
            to_masked_array("rad1_ray_ind", dtype=int),
            to_masked_array("rad1_rng_ind", dtype=int),
            to_masked_array("rad1_ele", dtype=float),
            to_masked_array("rad1_azi", dtype=float),
            to_masked_array("rad1_rng", dtype=float),
            to_masked_array("rad1_val", dtype=float),
            to_masked_array("rad2_time"),
            to_masked_array("rad2_ray_ind", dtype=int),
            to_masked_array("rad2_rng_ind", dtype=int),
            to_masked_array("rad2_ele", dtype=float),
            to_masked_array("rad2_azi", dtype=float),
            to_masked_array("rad2_rng", dtype=float),
            to_masked_array("rad2_val", dtype=float),
        )

    except Exception as e:
        warn(str(e))
        warn("Unable to read file " + fname)
        return (None,) * 14


def read_colocated_data_time_avg(fname):
    """
    Reads a CSV file containing time-averaged colocated data using pandas.

    Parameters
    ----------
    fname : str
        Path of time series file.

    Returns
    -------
    tuple
        A tuple with the data read:
        (rad1_time, rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi, rad1_rng,
        rad1_dBZavg, rad1_PhiDPavg, rad1_Flagavg, rad2_time, rad2_ray_ind,
        rad2_rng_ind, rad2_ele, rad2_azi, rad2_rng, rad2_dBZavg, rad2_PhiDPavg,
        rad2_Flagavg). Returns None for all if unable to read the file.
    """
    try:
        # Read the CSV file, skipping lines starting with "#"
        df = pd.read_csv(fname, comment="#", na_values="--")

        # Convert rad1_time and rad2_time to datetime
        df["rad1_time"] = pd.to_datetime(df["rad1_time"], format="%Y%m%d%H%M%S")
        df["rad2_time"] = pd.to_datetime(df["rad2_time"], format="%Y%m%d%H%M%S")

        # Create masked arrays
        def to_masked_array(column, dtype=None):
            return np.ma.masked_invalid(df[column].to_numpy(dtype=dtype))

        return (
            to_masked_array("rad1_time"),
            to_masked_array("rad1_ray_ind", dtype=int),
            to_masked_array("rad1_rng_ind", dtype=int),
            to_masked_array("rad1_ele", dtype=float),
            to_masked_array("rad1_azi", dtype=float),
            to_masked_array("rad1_rng", dtype=float),
            to_masked_array("rad1_dBZavg", dtype=float),
            to_masked_array("rad1_PhiDPavg", dtype=float),
            to_masked_array("rad1_Flagavg", dtype=float),
            to_masked_array("rad2_time"),
            to_masked_array("rad2_ray_ind", dtype=int),
            to_masked_array("rad2_rng_ind", dtype=int),
            to_masked_array("rad2_ele", dtype=float),
            to_masked_array("rad2_azi", dtype=float),
            to_masked_array("rad2_rng", dtype=float),
            to_masked_array("rad2_dBZavg", dtype=float),
            to_masked_array("rad2_PhiDPavg", dtype=float),
            to_masked_array("rad2_Flagavg", dtype=float),
        )

    except Exception as e:
        warn(str(e))
        warn("Unable to read file " + fname)
        return (None,) * 18


def read_timeseries(fname):
    """
    Reads a time series contained in a CSV file using pandas.

    Parameters
    ----------
    fname : str
        Path of time series file

    Returns
    -------
    date, value : tuple
        A list of datetime objects and a numpy masked array of values.
    """
    try:
        # Read CSV while skipping commented lines
        df = pd.read_csv(
            fname,
            comment='#',
            parse_dates=['date'],
            usecols=['date', 'value']
        )

        # Drop duplicate timestamps
        df = df.drop_duplicates(subset='date')

        masked_values = pd.to_numeric(df['value'], errors='coerce')
        # Mask fill values
        fill_val = get_fillvalue()
        masked_values = np.ma.masked_values(masked_values.to_numpy(), fill_val)

        return df['date'].tolist(), masked_values

    except Exception as e:
        warn(str(e))
        warn('Unable to read file ' + fname)
        return None, None
    
def read_ts_cum(fname):
    """
    Reads a time series of precipitation accumulation contained in a csv file

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    date, np_radar, radar_value, np_sensor, sensor_value : tupple
        The data read

    """
    try:
        with open(fname, "r", newline="") as csvfile:
            # first count the lines
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            nrows = sum(1 for row in reader)
            np_radar = np.zeros(nrows, dtype=int)
            radar_value = np.ma.empty(nrows, dtype=float)
            np_sensor = np.zeros(nrows, dtype=int)
            sensor_value = np.ma.empty(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            date = list()
            for i, row in enumerate(reader):
                date.append(
                    datetime.datetime.strptime(row["date"], "%Y-%m-%d %H:%M:%S")
                )
                np_radar[i] = int(row["np_radar"])
                radar_value[i] = float(row["radar_value"])
                np_sensor[i] = int(row["np_sensor"])
                sensor_value[i] = float(row["sensor_value"])

            radar_value = np.ma.masked_values(radar_value, get_fillvalue())
            sensor_value = np.ma.masked_values(sensor_value, get_fillvalue())

            csvfile.close()

            return date, np_radar, radar_value, np_sensor, sensor_value
    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None, None, None, None, None


def read_ml_ts(fname):
    """
    Reads a melting layer time series contained in a csv file

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    dt_ml, ml_top_avg, ml_top_std, thick_avg, thick_std, nrays_valid,
    nrays_total : tupple
        The read data. None otherwise

    """
    try:
        with open(fname, "r", newline="") as csvfile:
            # first count the lines
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            nrows = sum(1 for row in reader)

            dt_ml = np.empty(nrows, dtype=datetime.datetime)
            ml_top_avg = np.ma.empty(nrows, dtype=float)
            ml_top_std = np.ma.empty(nrows, dtype=float)
            thick_avg = np.ma.empty(nrows, dtype=float)
            thick_std = np.ma.empty(nrows, dtype=float)
            nrays_valid = np.zeros(nrows, dtype=int)
            nrays_total = np.zeros(nrows, dtype=int)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))

            for i, row in enumerate(reader):
                dt_ml[i] = datetime.datetime.strptime(
                    row["date-time [UTC]"], "%Y-%m-%d %H:%M:%S"
                )
                ml_top_avg[i] = float(row["mean ml top height [m MSL]"])
                ml_top_std[i] = float(row["std ml top height [m MSL]"])
                thick_avg[i] = float(row["mean ml thickness [m]"])
                thick_std[i] = float(row["std ml thickness [m]"])
                nrays_valid[i] = int(row["N valid rays"])
                nrays_total[i] = int(row["rays total"])

            ml_top_avg = np.ma.masked_values(ml_top_avg, get_fillvalue())
            ml_top_std = np.ma.masked_values(ml_top_std, get_fillvalue())
            thick_avg = np.ma.masked_values(thick_avg, get_fillvalue())
            thick_std = np.ma.masked_values(thick_std, get_fillvalue())

            return (
                dt_ml,
                ml_top_avg,
                ml_top_std,
                thick_avg,
                thick_std,
                nrays_valid,
                nrays_total,
            )
    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None, None, None, None, None, None, None


def read_monitoring_ts(quantiles, fname, sort_by_date=False):
    """
    Reads a monitoring time series contained in a CSV file using pandas.

    Parameters
    ----------
    quantiles : 1D ndarray
        The computed quantiles.
    fname : str
        Path of time series file.
    sort_by_date : bool
        If True, the read data is sorted by date prior to exit.

    Returns
    -------
    date, np_t, quantile_data, geometric_mean_dB, linear_mean_dB : tuple
        - date: Array of timestamps.
        - np_t: Array of NP values.
        - quantile_data: Dictionary with quantile values, where keys are quantile names.
        - geometric_mean_dB: Mean reflectivity in dB.
        - linear_mean_dB: Mean reflectivity computed in linear units and converted back to dB.
    """
    try:
        with open(fname, "r+") as csvfile:
            lock_file(csvfile)

            # Read the CSV file, ignoring comment lines
            df = pd.read_csv(
                csvfile,
                comment="#",
                parse_dates=["date"],
                date_parser=lambda x: datetime.datetime.strptime(x, "%Y%m%d%H%M%S"),
            )

            if sort_by_date:
                df.sort_values(by="date", inplace=True)

            df = df.drop_duplicates(subset="date", keep="last")
            unlock_file(csvfile)

        # Extract date and NP values
        date = df["date"].dt.to_pydatetime()
        np_t = df["NP"].to_numpy(dtype=int)

        # Extract all quantile values dynamically
        quantile_data = {}
        for q in quantiles:
            column_name = f"quantile_{q}"
            if column_name in df.columns:
                quantile_data[q] = np.ma.masked_values(
                    df[column_name].to_numpy(dtype=float), get_fillvalue()
                )

        # Read mean values
        geometric_mean_dB = (
            np.ma.masked_values(df["geometric_mean_dB"].to_numpy(dtype=float), get_fillvalue())
            if "geometric_mean_dB" in df.columns else np.ma.masked
        )
        linear_mean_dB = (
            np.ma.masked_values(df["linear_mean_dB"].to_numpy(dtype=float), get_fillvalue())
            if "linear_mean_dB" in df.columns else np.ma.masked
        )

        return date, np_t, quantile_data, geometric_mean_dB, linear_mean_dB

    except Exception as e:
        warn(str(e))
        warn(f"Unable to read file {fname}")
        return None, None, None, None, None


def read_monitoring_ts_old(fname):
    """
    Reads an old format of the monitoring time series contained in a text file

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    date , np_t, central_quantile, low_quantile, high_quantile : tupple
        The read data in the current format. None otherwise

    """
    try:
        with open(fname, "r", newline="") as csvfile:
            # first count the lines
            reader = csv.DictReader(
                (row for row in csvfile if not row.startswith("#")),
                fieldnames=["Date [YYJJJ]", "Npoints", "Value"],
            )
            nrows = sum(1 for row in reader)
            np_t = np.zeros(nrows, dtype=int)
            central_quantile = np.ma.empty(nrows, dtype=float)
            low_quantile = np.ma.masked_all(nrows, dtype=float)
            high_quantile = np.ma.masked_all(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                (row for row in csvfile if not row.startswith("#")),
                fieldnames=["Date [YYJJJ]", "Npoints", "Value"],
            )

            date = list()
            for i, row in enumerate(reader):
                date.append(datetime.datetime.strptime(row["Date [YYJJJ]"], "%y%j"))
                np_t[i] = int(row["Npoints"])
                central_quantile[i] = float(row["Value"])

            central_quantile = np.ma.masked_invalid(central_quantile)

            csvfile.close()

            return date, np_t, central_quantile, low_quantile, high_quantile
    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None, None, None, None, None


def read_intercomp_scores_ts(fname, sort_by_date=False):
    """
    Reads a radar intercomparison scores csv file

    Parameters
    ----------
    fname : str
        path of time series file
    sort_by_date : bool
        if True, the read data is sorted by date prior to exit

    Returns
    -------
    date_vec, np_vec, meanbias_vec, medianbias_vec, quant25bias_vec,
    quant75bias_vec, modebias_vec, corr_vec, slope_vec, intercep_vec,
    intercep_slope1_vec : tupple
        The read data. None otherwise

    """
    try:
        with open(fname, "r", newline="") as csvfile:
            if FCNTL_AVAIL:
                while True:
                    try:
                        fcntl.flock(csvfile, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        break
                    except OSError as e:
                        if e.errno == errno.EAGAIN:
                            time.sleep(0.1)
                        elif e.errno == errno.EBADF:
                            warn(
                                "WARNING: No file locking is possible (NFS mount?), "
                                + "expect strange issues with multiprocessing..."
                            )
                            break
                        else:
                            raise
            else:
                while True:
                    try:
                        # Attempt to acquire an exclusive lock on the file
                        msvcrt.locking(
                            os.open(csvfile, os.O_RDWR | os.O_CREAT), msvcrt.LK_NBLCK, 0
                        )
                        break
                    except OSError as e:
                        if e.errno == 13:  # Permission denied
                            time.sleep(0.1)
                        elif e.errno == 13:  # No such file or directory
                            warn(
                                "WARNING: No file locking is possible (NFS mount?), "
                                + "expect strange issues with multiprocessing..."
                            )
                            break
                        else:
                            raise

            # first count the lines
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            nrows = sum(1 for row in reader)

            np_vec = np.zeros(nrows, dtype=int)
            meanbias_vec = np.ma.empty(nrows, dtype=float)
            medianbias_vec = np.ma.empty(nrows, dtype=float)
            quant25bias_vec = np.ma.empty(nrows, dtype=float)
            quant75bias_vec = np.ma.empty(nrows, dtype=float)
            modebias_vec = np.ma.empty(nrows, dtype=float)
            corr_vec = np.ma.empty(nrows, dtype=float)
            slope_vec = np.ma.empty(nrows, dtype=float)
            intercep_vec = np.ma.empty(nrows, dtype=float)
            intercep_slope1_vec = np.ma.empty(nrows, dtype=float)
            date_vec = np.empty(nrows, dtype=datetime.datetime)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(row for row in csvfile if not row.startswith("#"))
            for i, row in enumerate(reader):
                date_vec[i] = datetime.datetime.strptime(row["date"], "%Y%m%d%H%M%S")
                np_vec[i] = int(row["NP"])
                meanbias_vec[i] = float(row["mean_bias"])
                medianbias_vec[i] = float(row["median_bias"])
                quant25bias_vec[i] = float(row["quant25_bias"])
                quant75bias_vec[i] = float(row["quant75_bias"])
                modebias_vec[i] = float(row["mode_bias"])
                corr_vec[i] = float(row["corr"])
                slope_vec[i] = float(row["slope_of_linear_regression"])
                intercep_vec[i] = float(row["intercep_of_linear_regression"])
                intercep_slope1_vec[i] = float(
                    row["intercep_of_linear_regression_of_slope_1"]
                )

            fcntl.flock(csvfile, fcntl.LOCK_UN)
            csvfile.close()

            meanbias_vec = np.ma.masked_values(meanbias_vec, get_fillvalue())
            medianbias_vec = np.ma.masked_values(medianbias_vec, get_fillvalue())
            quant25bias_vec = np.ma.masked_values(quant25bias_vec, get_fillvalue())
            quant75bias_vec = np.ma.masked_values(quant75bias_vec, get_fillvalue())
            modebias_vec = np.ma.masked_values(modebias_vec, get_fillvalue())
            corr_vec = np.ma.masked_values(corr_vec, get_fillvalue())
            slope_vec = np.ma.masked_values(slope_vec, get_fillvalue())
            intercep_vec = np.ma.masked_values(intercep_vec, get_fillvalue())
            intercep_slope1_vec = np.ma.masked_values(
                intercep_slope1_vec, get_fillvalue()
            )

            if sort_by_date:
                ind = np.argsort(date_vec)
                date_vec = date_vec[ind]
                np_vec = np_vec[ind]
                meanbias_vec = meanbias_vec[ind]
                medianbias_vec = medianbias_vec[ind]
                quant25bias_vec = quant25bias_vec[ind]
                quant75bias_vec = quant75bias_vec[ind]
                modebias_vec = modebias_vec[ind]
                corr_vec = corr_vec[ind]
                slope_vec = slope_vec[ind]
                intercep_vec = intercep_vec[ind]
                intercep_slope1_vec = intercep_slope1_vec[ind]

            return (
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
            )
    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return (None, None, None, None, None, None, None, None, None, None, None)


def read_intercomp_scores_ts_old(fname):
    """
    Reads a radar intercomparison scores csv file in old format

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    date_vec, np_vec, meanbias_vec, medianbias_vec, quant25bias_vec,
    quant75bias_vec, modebias_vec, corr_vec, slope_vec, intercep_vec,
    intercep_slope1_vec : tupple
        The read data. None otherwise

    """
    try:
        with open(fname, "r", newline="") as csvfile:
            # first count the lines
            reader = csv.DictReader(
                (row for row in csvfile if not row.startswith("#")),
                fieldnames=[
                    "date",
                    "NP",
                    "mode_bias",
                    "median_bias",
                    "mean_bias",
                    "corr",
                    "slope_of_linear_regression",
                    "intercep_of_linear_regression",
                ],
                dialect="excel-tab",
            )
            nrows = sum(1 for row in reader)

            np_vec = np.zeros(nrows, dtype=int)
            meanbias_vec = np.ma.empty(nrows, dtype=float)
            medianbias_vec = np.ma.empty(nrows, dtype=float)
            quant25bias_vec = np.ma.masked_all(nrows, dtype=float)
            quant75bias_vec = np.ma.masked_all(nrows, dtype=float)
            modebias_vec = np.ma.empty(nrows, dtype=float)
            corr_vec = np.ma.empty(nrows, dtype=float)
            slope_vec = np.ma.empty(nrows, dtype=float)
            intercep_vec = np.ma.empty(nrows, dtype=float)
            intercep_slope1_vec = np.ma.masked_all(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                (row for row in csvfile if not row.startswith("#")),
                fieldnames=[
                    "date",
                    "NP",
                    "mode_bias",
                    "median_bias",
                    "mean_bias",
                    "corr",
                    "slope_of_linear_regression",
                    "intercep_of_linear_regression",
                ],
                dialect="excel-tab",
            )
            date_vec = list()
            for i, row in enumerate(reader):
                date_vec.append(datetime.datetime.strptime(row["date"], "%y%j"))
                np_vec[i] = int(row["NP"])
                meanbias_vec[i] = float(row["mean_bias"])
                medianbias_vec[i] = float(row["median_bias"])
                modebias_vec[i] = float(row["mode_bias"])
                corr_vec[i] = float(row["corr"])
                slope_vec[i] = float(row["slope_of_linear_regression"])
                intercep_vec[i] = float(row["intercep_of_linear_regression"])

            meanbias_vec = np.ma.masked_invalid(meanbias_vec)
            medianbias_vec = np.ma.masked_invalid(medianbias_vec)
            modebias_vec = np.ma.masked_invalid(modebias_vec)
            corr_vec = np.ma.masked_invalid(corr_vec)
            slope_vec = np.ma.masked_invalid(slope_vec)
            intercep_vec = np.ma.masked_invalid(intercep_vec)

            csvfile.close()

            return (
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
            )
    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return (None, None, None, None, None, None, None, None, None, None, None)


def read_intercomp_scores_ts_old_v0(fname, corr_min=0.6, np_min=9):
    """
    Reads a radar intercomparison scores csv file in the oldest format

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    date_vec, np_vec, meanbias_vec, medianbias_vec, quant25bias_vec,
    quant75bias_vec, modebias_vec, corr_vec, slope_vec, intercep_vec,
    intercep_slope1_vec : tupple
        The read data. None otherwise

    """
    try:
        with open(fname, "r", newline="") as csvfile:
            # first count the lines
            reader = csv.DictReader(
                (row for row in csvfile if not row.startswith("#")),
                fieldnames=["date", "ele", "NP", "mean_bias", "corr"],
                dialect="excel-tab",
            )
            nrows = sum(1 for row in reader)

            np_aux_vec = np.zeros(nrows, dtype=int)
            meanbias_aux_vec = np.ma.empty(nrows, dtype=float)
            corr_aux_vec = np.ma.empty(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                (row for row in csvfile if not row.startswith("#")),
                fieldnames=["date", "ele", "NP", "mean_bias", "corr"],
                dialect="excel-tab",
            )
            date_aux_vec = list()
            for i, row in enumerate(reader):
                date_aux_vec.append(datetime.datetime.strptime(row["date"], "%y%j"))
                np_aux_vec[i] = int(row["NP"])
                meanbias_aux_vec[i] = float(row["mean_bias"])
                corr_aux_vec[i] = float(row["corr"])

            meanbias_aux_vec = np.ma.masked_invalid(meanbias_aux_vec)
            corr_aux_vec = np.ma.masked_invalid(corr_aux_vec)

            csvfile.close()

            date_aux_vec = np.asarray(date_aux_vec)
            date_vec, unique_ind = np.unique(date_aux_vec, return_index=True)
            nelements = len(date_vec)

            np_vec = np.zeros(nelements, dtype=int)
            meanbias_vec = np.ma.empty(nelements, dtype=float)
            corr_vec = np.ma.empty(nelements, dtype=float)
            for i in range(nelements - 1):
                ind_aux = np.arange(unique_ind[i], unique_ind[i + 1])
                np_aux = np_aux_vec[ind_aux]
                meanbias_aux = meanbias_aux_vec[ind_aux]
                corr_aux = corr_aux_vec[ind_aux]

                ind_valid = np.where(
                    np.logical_and(corr_aux > corr_min, np_aux > np_min)
                )
                np_aux = np_aux[ind_valid]
                meanbias_aux = meanbias_aux[ind_valid]
                corr_aux = corr_aux[ind_valid]

                np_vec[i] = np.sum(np_aux, dtype=int)
                if np_vec[i] == 0:
                    continue
                meanbias_vec[i] = np.sum(np_aux * meanbias_aux) / np_vec[i]
                corr_vec[i] = np.sum(np_aux * corr_aux) / np_vec[i]

            # last date
            ind_aux = np.arange(unique_ind[-1], len(date_aux_vec))
            np_aux = np_aux_vec[ind_aux]
            meanbias_aux = meanbias_aux_vec[ind_aux]
            corr_aux = corr_aux_vec[ind_aux]

            ind_valid = np.where(np.logical_and(corr_aux > corr_min, np_aux > np_min))
            np_aux = np_aux[ind_valid]
            meanbias_aux = meanbias_aux[ind_valid]
            corr_aux = corr_aux[ind_valid]
            np_vec[-1] = np.sum(np_aux, dtype=int)
            if np_vec[-1] > 0:
                meanbias_vec[-1] = np.sum(np_aux * meanbias_aux) / np_vec[-1]
                corr_vec[-1] = np.sum(np_aux * corr_aux) / np_vec[-1]

            medianbias_vec = np.ma.masked_all(nelements, dtype=float)
            quant25bias_vec = np.ma.masked_all(nelements, dtype=float)
            quant75bias_vec = np.ma.masked_all(nelements, dtype=float)
            modebias_vec = np.ma.masked_all(nelements, dtype=float)
            slope_vec = np.ma.masked_all(nelements, dtype=float)
            intercep_vec = np.ma.masked_all(nelements, dtype=float)
            intercep_slope1_vec = np.ma.masked_all(nelements, dtype=float)

            return (
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
            )
    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return (None, None, None, None, None, None, None, None, None, None, None)


def read_selfconsistency(fname):
    """
    Reads a self-consistency table with Zdr, Kdp/Zh columns

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    zdr, kdpzh : arrays
        The read values

    """
    try:
        with open(fname, "r", newline="", encoding="utf-8", errors="ignore") as csvfile:
            # first count the lines
            reader = csv.DictReader(
                csvfile, fieldnames=["zdr", "kdpzh"], dialect="excel-tab"
            )
            nrows = sum(1 for row in reader)

            zdr_kdpzh_table = np.empty((2, nrows), dtype="float32")

            # now read the data
            csvfile.seek(0)

            reader = csv.DictReader(
                csvfile, fieldnames=["zdr", "kdpzh"], dialect="excel-tab"
            )
            for i, row in enumerate(reader):
                zdr_kdpzh_table[0, i] = float(row["zdr"])
                zdr_kdpzh_table[1, i] = float(row["kdpzh"])

            csvfile.close()

            return zdr_kdpzh_table
    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None


def read_antenna_pattern(fname, linear=False, twoway=False):
    """
    Read antenna pattern from file

    Parameters
    ----------
    fname : str
        path of the antenna pattern file
    linear : boolean
        if true the antenna pattern is given in linear units
    twoway : boolean
        if true the attenuation is two-way

    Returns
    -------
    pattern : dict
        dictionary with the fields angle and attenuation

    """

    try:
        pfile = open(fname, "r")
    except EnvironmentError as ee:
        warn(str(ee))
        raise Exception("ERROR: Could not find|open antenna file '" + fname + "'")

    # Get array length
    linenum = 0
    for line in pfile:
        if line.startswith("# CRC:"):
            break
        linenum += 1

    pattern = {"angle": np.empty(linenum), "attenuation": np.empty(linenum)}

    pfile.seek(0)
    cnt = 0
    for line in pfile:
        if line.startswith("# CRC:"):
            break
        line = line.strip()
        fields = line.split(",")
        pattern["angle"][cnt] = float(fields[0])
        pattern["attenuation"][cnt] = float(fields[1])
        cnt += 1

    pfile.close()

    if twoway:
        pattern["attenuation"] = 2.0 * pattern["attenuation"]

    if linear:
        pattern["attenuation"] = 10.0 ** (pattern["attenuation"] / 10.0)

    return pattern


def read_xls(xls_file):
    """Reads an excel file such as those used for CPC vlaidation

    Parameters
    ----------
    fname : str
        Full path of the excel file to be read

    Returns
    -------
    The excel file as a Pandas dataframe
    """

    data = pd.read_excel(xls_file, sheet_name=None)
    keys = list(data.keys())
    hourly_keys = []
    for k in keys:
        if "Data Hourly" in k:
            hourly_keys.append(k)
    out = pd.concat([data[k] for k in hourly_keys])
    return out
