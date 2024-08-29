"""
pyrad.graph.plots_aux
=====================

Auxiliary plotting functions

.. autosummary::
    :toctree: generated/

    generate_complex_range_Doppler_title
    generate_angle_Doppler_title
    generate_complex_Doppler_title
    generate_fixed_rng_span_title
    generate_fixed_rng_title
    generate_dda_map_title
    generate_dda_latitude_slice_title
    generate_dda_longitude_slice_title
    get_colobar_label
    get_field_name
    get_norm

"""

import numpy as np

import pyart

import matplotlib as mpl
import matplotlib.cm

mpl.use("Agg")

# Increase a bit font size
mpl.rcParams.update({"font.size": 16})
mpl.rcParams.update({"font.family": "sans-serif"})


def generate_complex_range_Doppler_title(radar, field, ray, datetime_format=None):
    """
    creates the fixed range plot title

    Parameters
    ----------
    radar : radar
        The radar object
    field : str
        name of the field
    stat : str
        The statistic computed
    datetime_forat : str or None
        The date time format to use

    Returns
    -------
    titl : str
        The plot title

    """
    begin_time = pyart.graph.common.generate_radar_time_begin(radar)
    if datetime_format:
        time_str = begin_time.strftime(datetime_format)
    else:
        time_str = begin_time.isoformat() + "Z"
    l1 = "%s azi%.1f-ele%.1f deg. %s " % (
        pyart.graph.common.generate_radar_name(radar),
        radar.azimuth["data"][ray],
        radar.elevation["data"][ray],
        time_str,
    )
    field_name = pyart.graph.common.generate_field_name(radar, field)
    return l1 + "\n" + field_name


def generate_angle_Doppler_title(
    radar, field, ang, ind_rng, along_azi=True, datetime_format=None
):
    """
    creates the angle-Doppler plot title

    Parameters
    ----------
    radar : radar
        The radar object
    field : str
        name of the field
    ang : float
        The fixed angle
    ind_rng : int
        the index of the fixed range
    along_azi : bool
        If true the plot is performed along azimuth, otherwise it is performed
        along elevation
    datetime_forat : str or None
        The date time format to use

    Returns
    -------
    titl : str
        The plot title

    """
    begin_time = pyart.graph.common.generate_radar_time_begin(radar)
    if datetime_format:
        time_str = begin_time.strftime(datetime_format)
    else:
        time_str = begin_time.isoformat() + "Z"
    if along_azi:
        ang_type = "ele"
    else:
        ang_type = "azi"
    l1 = "%s %s%.1f deg-rng%.1f m. %s " % (
        pyart.graph.common.generate_radar_name(radar),
        ang_type,
        ang,
        radar.range["data"][ind_rng],
        time_str,
    )
    field_name = pyart.graph.common.generate_field_name(radar, field)
    return l1 + "\n" + field_name


def generate_complex_Doppler_title(radar, field, ray, rng, datetime_format=None):
    """
    creates the fixed range plot title

    Parameters
    ----------
    radar : radar
        The radar object
    field : str
        name of the field
    stat : str
        The statistic computed
    datetime_forat : str or None
        The date time format to use

    Returns
    -------
    titl : str
        The plot title

    """
    begin_time = pyart.graph.common.generate_radar_time_begin(radar)
    if datetime_format:
        time_str = begin_time.strftime(datetime_format)
    else:
        time_str = begin_time.isoformat() + "Z"
    l1 = "%s azi%.1f-ele%.1f deg rng%.1f km. %s " % (
        pyart.graph.common.generate_radar_name(radar),
        radar.azimuth["data"][ray],
        radar.elevation["data"][ray],
        radar.range["data"][rng] / 1000.0,
        time_str,
    )
    field_name = pyart.graph.common.generate_field_name(radar, field)
    return l1 + "\n" + field_name


def generate_fixed_rng_span_title(radar, field, stat, datetime_format=None):
    """
    creates the fixed range plot title

    Parameters
    ----------
    radar : radar
        The radar object
    field : str
        name of the field
    stat : str
        The statistic computed
    datetime_forat : str or None
        The date time format to use

    Returns
    -------
    titl : str
        The plot title

    """
    begin_time = pyart.graph.common.generate_radar_time_begin(radar)
    if datetime_format:
        time_str = begin_time.strftime(datetime_format)
    else:
        time_str = begin_time.isoformat() + "Z"
    l1 = "%s %.1f-%.1f m %s. %s " % (
        pyart.graph.common.generate_radar_name(radar),
        np.min(radar.range["data"]),
        np.max(radar.range["data"]),
        stat,
        time_str,
    )
    field_name = pyart.graph.common.generate_field_name(radar, field)
    return l1 + "\n" + field_name


def generate_fixed_rng_title(radar, field, fixed_rng, datetime_format=None):
    """
    creates the fixed range plot title

    Parameters
    ----------
    radar : radar
        The radar object
    field : str
        name of the field
    fixed_rng : float
        The fixed range [m]
    datetime_forat : str or None
        The date time format to use

    Returns
    -------
    titl : str
        The plot title

    """
    begin_time = pyart.graph.common.generate_radar_time_begin(radar)
    if datetime_format:
        time_str = begin_time.strftime(datetime_format)
    else:
        time_str = begin_time.isoformat() + "Z"
    l1 = "%s %.1f m. %s " % (
        pyart.graph.common.generate_radar_name(radar),
        fixed_rng,
        time_str,
    )
    field_name = pyart.graph.common.generate_field_name(radar, field)
    return l1 + "\n" + field_name


def generate_dda_map_title(grid, field, level, datetime_format=None):
    """
    creates the dda map plot title

    Parameters
    ----------
    grid : grid
        The grid object
    field : str
        name of the background field
    level : int
        Verical level plotted.
    datetime_format : str or None
        The date time format to use

    Returns
    -------
    titl : str
        The plot title

    """
    begin_time = pyart.graph.common.generate_radar_time_begin(grid)
    if datetime_format:
        time_str = begin_time.strftime(datetime_format)
    else:
        time_str = begin_time.isoformat() + "Z"
    radar_names = ""

    for mdata in grid.metadata["additional_radars"]:
        if isinstance(mdata["radar_name"], bytes):
            mdata["radar_name"] = mdata["radar_name"].decode("utf-8")
        radar_names += "-" + mdata["radar_name"]
    height = grid.z["data"][level] / 1000.0
    l1 = f"DDA: {radar_names} {height:.1f} km {time_str}"
    field_name = "Hor. wind vectors (u,v) with " + field
    return l1 + "\n" + field_name


def generate_dda_latitude_slice_title(
    grid, field, level, datetime_format=None, wind_vectors="hor"
):
    """
    creates the dda latitude slice plot title

    Parameters
    ----------
    grid : grid
        The grid object
    field : str
        name of the background field
    level : int
        latitudinal level plotted.
    datetime_format : str or None
        The date time format to use
    wind_vectors : str
        'hor' if horizontal wind vectors are displayed (u and v) or 'ver'
        if vertical wind vectors are displayed (v and w)

    Returns
    -------
    titl : str
        The plot title

    """
    begin_time = pyart.graph.common.generate_radar_time_begin(grid)
    if datetime_format:
        time_str = begin_time.strftime(datetime_format)
    else:
        time_str = begin_time.isoformat() + "Z"
    radar_names = ""

    for mdata in grid.metadata["additional_radars"]:
        if isinstance(mdata["radar_name"], bytes):
            mdata["radar_name"] = mdata["radar_name"].decode("utf-8")
        radar_names += "-" + mdata["radar_name"]
    disp = grid.x["data"][level] / 1000.0
    if disp >= 0:
        direction = "east"
    else:
        direction = "west"
        disp = -disp

    l1 = f"DDA: {radar_names} {disp:.1f} km {direction} of origin {time_str}"
    if wind_vectors == "hor":
        field_name = "Hor. wind vectors (u,v)"
    elif wind_vectors == "ver":
        field_name = "Vert. wind vectors (v,w)"
    field_name += " with " + field
    return l1 + "\n" + field_name


def generate_dda_longitude_slice_title(
    grid, field, level, datetime_format=None, wind_vectors="hor"
):
    """
    creates the dda longitude slice plot title

    Parameters
    ----------
    grid : grid
        The grid object
    field : str
        name of the background field
    level : int
        longitudinal level plotted.
    datetime_format : str or None
        The date time format to use
    wind_vectors : str
        'hor' if horizontal wind vectors are displayed (u and v) or 'ver'
        if vertical wind vectors are displayed (v and w)

    Returns
    -------
    titl : str
        The plot title

    """
    begin_time = pyart.graph.common.generate_radar_time_begin(grid)
    if datetime_format:
        time_str = begin_time.strftime(datetime_format)
    else:
        time_str = begin_time.isoformat() + "Z"
    radar_names = ""

    for mdata in grid.metadata["additional_radars"]:
        if isinstance(mdata["radar_name"], bytes):
            mdata["radar_name"] = mdata["radar_name"].decode("utf-8")
        radar_names += "-" + mdata["radar_name"]
    disp = grid.x["data"][level] / 1000.0
    if disp >= 0:
        direction = "north"
    else:
        direction = "south"
        disp = -disp
    l1 = f"DDA: {radar_names} {disp:.1f} km {direction} of origin {time_str}"
    if wind_vectors == "hor":
        field_name = "Hor. wind vectors (u,v)"
    elif wind_vectors == "ver":
        field_name = "Vert. wind vectors (v,w)"
    field_name += " with " + field
    return l1 + "\n" + field_name


def get_colobar_label(field_dict, field_name):
    """
    creates the colorbar label using field metadata

    Parameters
    ----------
    field_dict : dict
        dictionary containing field metadata
    field_name : str
        name of the field

    Returns
    -------
    label : str
        colorbar label

    """
    if "standard_name" in field_dict:
        standard_name = field_dict["standard_name"]
    elif "long_name" in field_dict:
        standard_name = field_dict["long_name"]
    else:
        standard_name = field_name

    if "units" in field_dict:
        units = field_dict["units"]
    else:
        units = "?"

    return pyart.graph.common.generate_colorbar_label(standard_name, units)


def get_field_name(field_dict, field):
    """
    Return a nice field name for a particular field

    Parameters
    ----------
    field_dict : dict
        dictionary containing field metadata
    field : str
        name of the field

    Returns
    -------
    field_name : str
        the field name

    """
    if "standard_name" in field_dict:
        field_name = field_dict["standard_name"]
    elif "long_name" in field_dict:
        field_name = field_dict["long_name"]
    else:
        field_name = str(field)
    field_name = field_name.replace("_", " ")
    field_name = field_name[0].upper() + field_name[1:]

    return field_name


def get_norm(field_name, field_dict={}, isxarray=False):
    """
    Computes the normalization of the colormap, and gets the ticks and labels
    of the colorbar from the metadata of the field. Returns None if the
    required parameters are not present in the metadata. If field dict is
    not None the metadata is obtained directly from the field. Otherwise it is
    obtained from the Py-ART config file

    Parameters
    ----------
    field_name : str
        name of the field
    field_dict : dict or None
        dictionary containing the field and its metadata.
    isxarray : bool
        whether or not the norm will be used with xarray's plotting functions
        default is false, which means that matplotlib plotting functions will
        be used should be set to true when plotting Grid objects which are
        handled as xarray by pyart
    Returns
    -------
    norm : list
        the colormap index
    ticks : list
        the list of ticks in the colorbar
    labels : list
        the list of labels corresponding to each tick

    """
    norm = None
    ticks = None
    ticklabs = None

    ref_dict = pyart.config.get_metadata(field_name)
    cmap = mpl.colormaps.get_cmap(pyart.config.get_field_colormap(field_name))

    if field_dict is not None and "boundaries" in field_dict:
        if isxarray:
            ncolors = len(field_dict["boundaries"]) - 1
        else:
            ncolors = cmap.N
        norm = mpl.colors.BoundaryNorm(
            boundaries=field_dict["boundaries"], ncolors=ncolors
        )
    elif "boundaries" in ref_dict:
        if isxarray:
            ncolors = len(ref_dict["boundaries"]) - 1
        else:
            ncolors = cmap.N
        norm = mpl.colors.BoundaryNorm(
            boundaries=ref_dict["boundaries"], ncolors=ncolors
        )

    if field_dict is not None and "ticks" in field_dict:
        ticks = field_dict["ticks"]
        if "labels" in field_dict:
            ticklabs = field_dict["labels"]
    elif "ticks" in ref_dict:
        ticks = ref_dict["ticks"]
        if "labels" in ref_dict:
            ticklabs = ref_dict["labels"]
    return norm, ticks, ticklabs
