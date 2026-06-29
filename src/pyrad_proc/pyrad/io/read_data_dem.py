"""
pyrad.io.read_data_dem
========================

Functions for reading data derived from Digital Elevation Models (DEM)
"""

import pathlib
import numpy as np
import pyart

from scipy.interpolate import RegularGridInterpolator
from pyart.config import get_metadata
from .read_data_icon import _put_radar_in_swiss_coord
from ..util import warn

try:
    import rasterio
    from pyproj import CRS

    _RASTERIO_AVAILABLE = True
except ImportError:
    _RASTERIO_AVAILABLE = False
    warn("rasterio/pyproj unavailable, you will not be able to read DEM data...")


def dem2radar_data(radar, dem_data, slice_xy=True, field_name="visibility"):
    x_radar, y_radar, _ = _put_radar_in_swiss_coord(radar)

    x_dem, y_dem, ind_xmin, ind_ymin, ind_xmax, ind_ymax = _prepare_for_interpolation(
        x_radar, y_radar, dem_data, slice_xy=slice_xy
    )

    if field_name not in dem_data:
        warn("DEM field " + field_name + " data not available")
        return None

    values = dem_data[field_name]["data"][
        ind_xmin : ind_xmax + 1, ind_ymin : ind_ymax + 1
    ]

    values = np.ma.filled(values, np.nan)
    interp_func = RegularGridInterpolator((x_dem, y_dem), values, bounds_error=False)

    data_interp = interp_func((x_radar, y_radar))
    data_interp = np.ma.masked_invalid(data_interp)

    field_dict = get_metadata(field_name)
    field_dict["data"] = data_interp.astype(float)

    return field_dict


def read_dem(fname, field_name="terrain_altitude", fill_value=None, projparams=None):
    extension = pathlib.Path(fname).suffix.lower()

    if isinstance(projparams, int):
        proj4 = CRS.from_epsg(projparams).to_proj4()
        projparams = _proj4_str_to_dict(proj4)

    if extension in [".tif", ".tiff", ".gtif"]:
        dem = read_geotiff_data(fname, fill_value)
    elif extension in [".asc", ".dem", ".txt"]:
        dem = read_ascii_data(fname, fill_value)
    elif extension in [".rst"]:
        dem = read_idrisi_data(fname, fill_value)
    else:
        warn(
            f"WARNING: unable to read file {fname}, extension must be .tif .tiff "
            + ".gtif, .asc .dem .txt .rst"
        )
        return None

    if dem is None:
        return None

    if "projparams" in dem and dem["projparams"] is not None:
        projparams = dem["projparams"]

    field_dict = get_metadata(field_name)
    field_dict["data"] = dem["data"][::-1, :][None, :, :]
    field_dict["units"] = dem["metadata"].get("value units", "meters")

    x = get_metadata("x")
    y = get_metadata("y")
    z = get_metadata("z")

    orig_lat = get_metadata("origin_latitude")
    orig_lon = get_metadata("origin_longitude")
    orig_alt = get_metadata("origin_altitude")

    x["data"] = (
        np.arange(dem["metadata"]["columns"]) * dem["metadata"]["resolution"]
        + dem["metadata"]["resolution"] / 2.0
        + dem["metadata"]["min. X"]
    )

    y["data"] = (
        np.arange(dem["metadata"]["rows"]) * dem["metadata"]["resolution"]
        + dem["metadata"]["resolution"] / 2.0
        + dem["metadata"]["min. Y"]
    )

    z["data"] = np.array([0])
    orig_lat["data"] = [y["data"][0]]
    orig_lon["data"] = [x["data"][0]]
    orig_alt["data"] = [0]

    if projparams is None:
        warn(
            "WARNING: No proj could be read from file and no projparams were provided, "
            + "assuming the projection is CH1903.",
            use_debug=False,
        )
        projparams = _get_lv1903_proj4()

    time = get_metadata("grid_time")
    time["data"] = np.array([0.0])
    time["units"] = "seconds since 2000-01-01T00:00:00Z"

    dem_data = pyart.core.Grid(
        time,
        {field_name: field_dict},
        dem["metadata"],
        orig_lat,
        orig_lon,
        orig_alt,
        x,
        y,
        z,
        projection=projparams,
    )

    return dem_data


def _crs_to_projparams(crs):
    if crs is None:
        return None

    try:
        proj4 = CRS.from_user_input(crs).to_proj4()
        return _proj4_str_to_dict(proj4)
    except Exception:
        return None


def read_geotiff_data(fname, fill_value=None):
    if not _RASTERIO_AVAILABLE:
        warn("rasterio is required to use read_geotiff_data but is not installed")
        return None

    try:
        with rasterio.open(fname) as src:
            rasterarray = src.read(1)
            transform = src.transform

            nodata = src.nodata
            if fill_value is None:
                fill_value = nodata if nodata is not None else np.nan

            metadata = {
                "resolution": abs(transform.a),
                "min. X": transform.c,
                "max. Y": transform.f,
                "columns": src.width,
                "rows": src.height,
                "min. Y": transform.f + src.height * transform.e,
                "max. X": transform.c + src.width * transform.a,
                "value units": "meters",
                "flag value": fill_value,
            }

            if fill_value is None or not np.isfinite(fill_value):
                rasterarray = np.ma.masked_invalid(rasterarray)
            else:
                rasterarray = np.ma.masked_equal(rasterarray, fill_value)

            return {
                "projparams": _crs_to_projparams(src.crs),
                "metadata": metadata,
                "data": rasterarray,
            }

    except Exception as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None


def read_ascii_data(fname, fill_value=None):
    try:
        with open(fname, "r") as file:
            lines = [file.readline().strip().split() for _ in range(6)]

        metadata = {
            "columns": int(lines[0][1]),
            "rows": int(lines[1][1]),
            "min. X": float(lines[2][1]),
            "min. Y": float(lines[3][1]),
            "resolution": float(lines[4][1]),
            "flag value": float(lines[5][1]),
            "max. X": float(lines[2][1]) + float(lines[4][1]) * int(lines[0][1]),
            "max. Y": float(lines[3][1]) + float(lines[4][1]) * int(lines[1][1]),
            "value units": "m",
        }

        if fill_value is None:
            fill_value = metadata["flag value"]

        raster_array = np.loadtxt(fname, skiprows=6)
        raster_array = np.ma.masked_equal(raster_array, fill_value)

        return {"metadata": metadata, "data": raster_array}

    except EnvironmentError as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None


def read_idrisi_data(fname, fill_value=None):
    if not _RASTERIO_AVAILABLE:
        warn("rasterio is required to use read_idrisi_data but is not installed")
        return None

    try:
        if fill_value is None:
            fill_value = -99.0

        with rasterio.open(fname) as src:
            rasterarray = src.read(1)

            if fill_value is None or not np.isfinite(fill_value):
                rasterarray = np.ma.masked_invalid(rasterarray)
            else:
                rasterarray = np.ma.masked_equal(rasterarray, fill_value)

            projparams = _crs_to_projparams(src.crs)

        metadata = read_idrisi_metadata(fname)

        if metadata is None:
            return None

        return {
            "projparams": projparams,
            "metadata": metadata,
            "data": rasterarray,
        }

    except Exception as ee:
        warn(str(ee))
        warn("Unable to read file " + fname)
        return None


def read_idrisi_metadata(fname):
    fname_rdc = fname.replace(".rst", ".rdc")

    try:
        metadata = dict()
        with open(fname_rdc, "r", newline="") as txtfile:
            for line in txtfile:
                strs = line.split(":")
                if len(strs) < 2:
                    continue
                metadata.update(
                    {strs[0].strip(): pyart.aux_io.convert_data(strs[1].strip())}
                )

        return metadata

    except EnvironmentError:
        warn("Unable to read file " + fname_rdc)
        return None


def _prepare_for_interpolation(x_radar, y_radar, dem_coord, slice_xy=True):
    nx_dem = len(dem_coord["x"]["data"])
    ny_dem = len(dem_coord["y"]["data"])

    if slice_xy:
        xmin = np.min(x_radar)
        xmax = np.max(x_radar)
        ymin = np.min(y_radar)
        ymax = np.max(y_radar)

        ind_xmin = np.where(dem_coord["x"]["data"] < xmin)[0]
        ind_xmin = 0 if ind_xmin.size == 0 else ind_xmin[-1]

        ind_xmax = np.where(dem_coord["x"]["data"] > xmax)[0]
        ind_xmax = nx_dem - 1 if ind_xmax.size == 0 else ind_xmax[0]

        ind_ymin = np.where(dem_coord["y"]["data"] < ymin)[0]
        ind_ymin = 0 if ind_ymin.size == 0 else ind_ymin[-1]

        ind_ymax = np.where(dem_coord["y"]["data"] > ymax)[0]
        ind_ymax = ny_dem - 1 if ind_ymax.size == 0 else ind_ymax[0]
    else:
        ind_xmin = 0
        ind_xmax = nx_dem - 1
        ind_ymin = 0
        ind_ymax = ny_dem - 1

    x_dem = dem_coord["x"]["data"][ind_xmin : ind_xmax + 1]
    y_dem = dem_coord["y"]["data"][ind_ymin : ind_ymax + 1]

    return x_dem, y_dem, ind_xmin, ind_ymin, ind_xmax, ind_ymax


def _proj4_str_to_dict(proj4str):
    return dict(
        item.split("=")
        for item in proj4str.strip(" ").split("+")
        if len(item.split("=")) == 2
    )


def _get_lv1903_proj4():
    proj4 = CRS.from_epsg(21781).to_proj4()
    lv1903 = _proj4_str_to_dict(proj4)

    if "no_defs" not in lv1903:
        lv1903["no_defs"] = 0

    return lv1903
