"""
pyrad.io.read_data_icon
========================

Functions for reading icon data

.. autosummary::
    :toctree: generated/

    icon2radar_data
    icon2radar_coord
    get_icon_fields
    read_icon_data
    read_icon_coord
    _ncvar_to_dict
    _prepare_for_interpolation
    _put_radar_in_swiss_coord


"""

from warnings import warn
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from scipy.spatial import cKDTree
import netCDF4

import pyart
from pyart.config import get_metadata, get_field_name

from ..io.io_aux import get_fieldname_icon

# from memory_profiler import profile

# import time


def icon2radar_data(
    radar,
    icon_coord,
    icon_data,
    time_index=0,
    slice_xy=True,
    slice_z=False,
    field_names=["temperature"],
    dtype=np.float32,
):
    """
    get the icon value corresponding to each radar gate using nearest
    neighbour interpolation

    Parameters
    ----------
    radar : Radar
        the radar object containing the information on the position of the
        radar gates
    icon_coord : dict
        dictionary containing the icon coordinates
    icon_data : dict
        dictionary containing the icon data
    time_index : int
        index of the forecasted data
    slice_xy : boolean
        if true the horizontal plane of the icon field is cut to the
        dimensions of the radar field
    slice_z : boolean
        if true the vertical plane of the icon field is cut to the dimensions
        of the radar field
    field_names : str
        names of icon fields to convert (default temperature)
    dtype : numpy data type object
        the data type of the output data

    Returns
    -------
    icon_fields : list of dict
        list of dictionary with the icon fields and metadata

    """
    # debugging
    # start_time = time.time()

    x_radar, y_radar, z_radar = _put_radar_in_swiss_coord(radar)

    (
        x_icon,
        y_icon,
        z_icon,
        ind_xmin,
        ind_ymin,
        ind_zmin,
        ind_xmax,
        ind_ymax,
        ind_zmax,
    ) = _prepare_for_interpolation(
        x_radar, y_radar, z_radar, icon_coord, slice_xy=slice_xy, slice_z=slice_z
    )

    icon_fields = []
    for field in field_names:
        if field not in icon_data:
            warn("icon field " + field + " data not available")
        else:
            values = icon_data[field]["data"][
                time_index,
                ind_zmin : ind_zmax + 1,
                ind_ymin : ind_ymax + 1,
                ind_xmin : ind_xmax + 1,
            ].flatten()
            # find interpolation function
            tree_options = {"compact_nodes": False, "balanced_tree": False}
            interp_func = NearestNDInterpolator(
                (z_icon, y_icon, x_icon), values, tree_options=tree_options
            )

            del values

            # interpolate
            data_interp = interp_func((z_radar, y_radar, x_radar))

            # put field
            field_dict = get_metadata(field)
            field_dict["data"] = data_interp.astype(dtype)
            icon_fields.append({field: field_dict})

            del data_interp

    if not icon_fields:
        warn("icon data not available")
        return None

    return icon_fields


def icon2radar_coord(radar, icon_coord, slice_xy=True, slice_z=False, field_name=None):
    """
    Given the radar coordinates find the nearest icon model pixel

    Parameters
    ----------
    radar : Radar
        the radar object containing the information on the position of the
        radar gates
    icon_coord : dict
        dictionary containing the icon coordinates
    slice_xy : boolean
        if true the horizontal plane of the icon field is cut to the
        dimensions of the radar field
    slice_z : boolean
        if true the vertical plane of the icon field is cut to the dimensions
        of the radar field
    field_name : str
        name of the field

    Returns
    -------
    icon_ind_field : dict
        dictionary containing a field of icon indices and metadata

    """
    # debugging
    # start_time = time.time()

    # parse the field parameters
    if field_name is None:
        field_name = get_field_name("icon_index")

    x_radar, y_radar, z_radar = _put_radar_in_swiss_coord(radar)

    (x_icon, y_icon, z_icon, ind_xmin, ind_ymin, ind_zmin, ind_xmax, ind_ymax, _) = (
        _prepare_for_interpolation(
            x_radar, y_radar, z_radar, icon_coord, slice_xy=slice_xy, slice_z=slice_z
        )
    )

    print("Generating tree")
    # default scipy compact_nodes and balanced_tree = True
    tree = cKDTree(
        np.transpose((z_icon, y_icon, x_icon)), compact_nodes=False, balanced_tree=False
    )
    print("Tree generated")
    _, ind_vec = tree.query(
        np.transpose((z_radar.flatten(), y_radar.flatten(), x_radar.flatten())), k=1
    )

    # put the index in the original icon coordinates
    nx_icon = len(icon_coord["x"]["data"])
    ny_icon = len(icon_coord["y"]["data"])

    nx = ind_xmax - ind_xmin + 1
    ny = ind_ymax - ind_ymin + 1

    ind_z = (ind_vec / (nx * ny)).astype(int) + ind_zmin
    ind_y = ((ind_vec - nx * ny * ind_z) / nx).astype(int) + ind_ymin
    ind_x = ((ind_vec - nx * ny * ind_z) % nx).astype(int) + ind_xmin
    ind_icon = (ind_x + nx_icon * ind_y + nx_icon * ny_icon * ind_z).astype(int)

    icon_ind_field = get_metadata(field_name)
    icon_ind_field["data"] = ind_icon.reshape(radar.nrays, radar.ngates)

    # debugging
    # print(" generating icon indices takes %s seconds " %
    #      (time.time() - start_time))

    return icon_ind_field


def get_icon_fields(icon_data, icon_ind, time_index=0, field_names=["temperature"]):
    """
    Get the icon data corresponding to each radar gate
    using a precomputed look up table of the nearest neighbour

    Parameters
    ----------
    icon_data : dict
        dictionary containing the icon data and metadata
    icon_ind : dict
        dictionary containing a field of icon indices and metadata
    time_index : int
        index of the forecasted data
    field_names : str
        names of icon parameters (default temperature)

    Returns
    -------
    icon_fields : list of dict
        dictionary with the icon fields and metadata

    """
    nrays, ngates = np.shape(icon_ind["data"])
    icon_fields = []
    for field in field_names:
        if field not in icon_data:
            warn("icon field " + field + " data not available")
        else:
            values = icon_data[field]["data"][time_index, :, :, :].flatten()

            # put field
            field_dict = get_metadata(field)
            field_dict["data"] = (
                values[icon_ind["data"].flatten()].reshape(nrays, ngates).astype(float)
            )
            icon_fields.append({field: field_dict})

    if not icon_fields:
        warn("icon data not available")
        return None

    return icon_fields


# @profile
def read_icon_data(fname, field_names=["temperature"], celsius=True):
    """
    Reads icon data from a netcdf file

    Parameters
    ----------
    fname : str
        name of the file to read
    field_names : str
        name of the variable to read
    celsius : Boolean
        if True and variable temperature converts data from Kelvin
        to Centigrade

    Returns
    -------
    icon_data : dictionary
        dictionary with the data and metadata

    """
    # read the data
    ncobj = netCDF4.Dataset(fname)
    ncvars = ncobj.variables

    # 4.1 Global attribute -> move to metadata dictionary
    metadata = dict([(k, getattr(ncobj, k)) for k in ncobj.ncattrs()])

    # read data for requested fields
    icon_data = dict()
    found = False
    for field in field_names:
        icon_name = get_fieldname_icon(field)
        if icon_name not in ncvars:
            warn(field + " data not present in icon file " + fname)
        else:
            var_data = _ncvar_to_dict(ncvars[icon_name], dtype="float16")

            # remove dimension ensemble member of icon-ch1-eps
            if var_data["data"].ndim == 5:
                var_data["data"] = np.squeeze(var_data["data"], axis=1)

            if field == "temperature" and celsius:
                var_data["data"] -= 273.15
                var_data["units"] = "degrees Celsius"
            if field == "vertical_wind_shear":
                var_data["data"] *= 1000.0
                var_data["units"] = "meters_per_second_per_km"
            icon_data.update({field: var_data})
            found = True
            del var_data
    if not found:
        warn("No field available in icon file " + fname)
        ncobj.close()
        return None

    # 4.2 put variables in dictionary
    x_1 = _ncvar_to_dict(ncvars["x_1"])
    y_1 = _ncvar_to_dict(ncvars["y_1"])
    lon_1 = _ncvar_to_dict(ncvars["lon_1"])
    lat_1 = _ncvar_to_dict(ncvars["lat_1"])
    z_1 = _ncvar_to_dict(ncvars["z_1"])
    z_bnds_1 = _ncvar_to_dict(ncvars["z_bnds_1"])
    time_data = _ncvar_to_dict(ncvars["time"])

    # close object
    ncobj.close()

    icon_data.update(
        {
            "metadata": metadata,
            "time": time_data,
            "x": x_1,
            "y": y_1,
            "z": z_1,
            "z_bnds": z_bnds_1,
            "lon": lon_1,
            "lat": lat_1,
        }
    )

    return icon_data


def read_icon_coord(fname, zmin=None):
    """
    Reads icon coordinates from a netcdf file

    Parameters
    ----------
    fname : str
        name of the file to read

    Returns
    -------
    icon_coord : dictionary
        dictionary with the data and metadata

    """
    # read the data
    try:
        ncobj = netCDF4.Dataset(fname)
        ncvars = ncobj.variables

        # 4.1 Global attribute -> move to metadata dictionary
        metadata = dict([(k, getattr(ncobj, k)) for k in ncobj.ncattrs()])

        # 4.2 put variables in dictionary
        x_1 = _ncvar_to_dict(ncvars["x_1"])
        y_1 = _ncvar_to_dict(ncvars["y_1"])
        lon_1 = _ncvar_to_dict(ncvars["lon_1"])
        lat_1 = _ncvar_to_dict(ncvars["lat_1"])
        z_1 = _ncvar_to_dict(ncvars["z_1"])
        z_bnds_1 = _ncvar_to_dict(ncvars["z_bnds_1"])
        hfl = _ncvar_to_dict(ncvars["HFL"])
        hsurf = _ncvar_to_dict(ncvars["HSURF"])
        fr_land = _ncvar_to_dict(ncvars["FR_LAND"])

        # close object
        ncobj.close()

        if zmin is not None:
            z_1["data"] = z_1["data"][z_1["data"] >= zmin]
            z_bnds_1["data"] = z_bnds_1["data"][z_bnds_1["data"] >= zmin]

        icon_coord = {
            "metadata": metadata,
            "x": x_1,
            "y": y_1,
            "z": z_1,
            "z_bnds": z_bnds_1,
            "lon": lon_1,
            "lat": lat_1,
            "hfl": hfl,
            "hsurf": hsurf,
            "fr_land": fr_land,
        }

        return icon_coord
    except EnvironmentError:
        warn("Unable to read file " + fname)
        return None


def _ncvar_to_dict(ncvar, dtype=np.float32):
    """Convert a NetCDF Dataset variable to a dictionary."""
    # copy all attributes
    d = dict((k, getattr(ncvar, k)) for k in ncvar.ncattrs())
    d.update({"data": ncvar[:]})
    if "_FillValue" in d:
        d["data"] = np.ma.asarray(d["data"], dtype=dtype)
        d["data"] = np.ma.masked_values(d["data"], float(d["_FillValue"]))
    else:
        d["data"] = np.asarray(d["data"], dtype=dtype)

    return d


def _prepare_for_interpolation(
    x_radar, y_radar, z_radar, icon_coord, slice_xy=True, slice_z=False
):
    """
    prepares the icon 3D volume for interpolation:
        1. if set slices the icon data to the area (or volume)
    covered by the radar
        2. creates the x, y, z grid for the interpolation

    Parameters
    ----------
    x_radar, y_radar, z_radar : arrays
        The Swiss coordinates of the radar
    icon_coord : dict
        dictionary containing the icon coordinates
    slice_xy : boolean
        if true the horizontal plane of the icon field is cut to the
        dimensions of the radar field
    slice_z : boolean
        if true the vertical plane of the icon field is cut to the dimensions
        of the radar field

    Returns
    -------
    x_icon, y_icon, z_icon : 1D arrays
        arrays containing the flatten swiss coordinates of the icon data in
        the area of interest
    ind_xmin, ind_ymin, ind_zmin, ind_xmax, ind_ymax, ind_zmax : ints
        the minimum and maximum indices of each dimension

    """
    nx_icon = len(icon_coord["x"]["data"])
    ny_icon = len(icon_coord["y"]["data"])
    nz_icon = len(icon_coord["z"]["data"])

    if slice_xy:
        # get the icon data within the radar range
        xmin = np.min(x_radar)
        xmax = np.max(x_radar)
        ymin = np.min(y_radar)
        ymax = np.max(y_radar)

        ind_xmin = np.where(icon_coord["x"]["data"] < xmin)[0]
        if ind_xmin.size == 0:
            ind_xmin = 0
        else:
            ind_xmin = ind_xmin[-1]

        ind_xmax = np.where(icon_coord["x"]["data"] > xmax)[0]
        if ind_xmax.size == 0:
            ind_xmax = nx_icon - 1
        else:
            ind_xmax = ind_xmax[0]

        ind_ymin = np.where(icon_coord["y"]["data"] < ymin)[0]
        if ind_ymin.size == 0:
            ind_ymin = 0
        else:
            ind_ymin = ind_ymin[-1]

        ind_ymax = np.where(icon_coord["y"]["data"] > ymax)[0]
        if ind_ymax.size == 0:
            ind_ymax = ny_icon - 1
        else:
            ind_ymax = ind_ymax[0]
    else:
        ind_xmin = 0
        ind_xmax = nx_icon - 1
        ind_ymin = 0
        ind_ymax = ny_icon - 1

    if slice_z:
        zmin = np.min(z_radar)
        zmax = np.max(z_radar)

        ind_z, _, _ = np.where(icon_coord["hfl"]["data"] < zmin)
        if ind_z.size == 0:
            ind_zmin = 0
        else:
            ind_zmin = np.min(ind_z)
        ind_z, _, _ = np.where(icon_coord["hfl"]["data"] > zmax)
        if ind_z.size == 0:
            ind_zmax = nz_icon - 1
        else:
            ind_zmax = np.max(ind_z)
    else:
        ind_zmin = 0
        ind_zmax = nz_icon - 1

    nx = ind_xmax - ind_xmin + 1
    ny = ind_ymax - ind_ymin + 1
    nz = ind_zmax - ind_zmin + 1

    x_icon = icon_coord["x"]["data"][ind_xmin : ind_xmax + 1]
    y_icon = icon_coord["y"]["data"][ind_ymin : ind_ymax + 1]
    z_icon = icon_coord["hfl"]["data"][
        ind_zmin : ind_zmax + 1, ind_ymin : ind_ymax + 1, ind_xmin : ind_xmax + 1
    ]

    x_icon = (np.broadcast_to(x_icon.reshape(1, 1, nx), (nz, ny, nx))).flatten()
    y_icon = (np.broadcast_to(y_icon.reshape(1, ny, 1), (nz, ny, nx))).flatten()
    z_icon = z_icon.flatten()

    return (
        x_icon,
        y_icon,
        z_icon,
        ind_xmin,
        ind_ymin,
        ind_zmin,
        ind_xmax,
        ind_ymax,
        ind_zmax,
    )


def _put_radar_in_swiss_coord(radar):
    """
    puts the Cartesian grid of the radar coordinates in Swiss coordinates

    Parameters
    ----------
    radar : Radar
        the radar object containing the information on the position of the
        radar gates

    Returns
    -------
    x_radar, y_radar, z_radar : 2D arrays
        arrays containing swiss coordinates of the radar [in m]

    """
    x0, y0, _ = pyart.core.wgs84_to_swissCH1903(
        radar.longitude["data"][0],
        radar.latitude["data"][0],
        radar.altitude["data"][0],
        no_altitude_transform=True,
    )

    x_radar = radar.gate_x["data"] + x0
    y_radar = radar.gate_y["data"] + y0
    z_radar = radar.gate_altitude["data"]

    return x_radar, y_radar, z_radar
