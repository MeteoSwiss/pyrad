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
import os
import matplotlib as mpl
import matplotlib.cm
from warnings import warn

mpl.use("Agg")

# Increase a bit font size
mpl.rcParams.update({"font.size": 16})
mpl.rcParams.update({"font.family": "sans-serif"})

_CARTOPY_AVAILABLE = False
try:
    import cartopy
    from cartopy.io.img_tiles import GoogleTiles
    from PIL import Image

    # Define ESRI terrain tiles
    class ShadedReliefESRI(GoogleTiles):
        def __init__(self, cache):
            super().__init__(self, cache=cache)
            self.desired_tile_form = "RGB"

        # shaded relief
        def _image_url(self, tile):
            x, y, z = tile
            return (
                f"https://server.arcgisonline.com/ArcGIS/rest/services/Elevation/"
                f"World_Hillshade/MapServer/tile/{z}/{y}/{x}"
            )

    class OTM(GoogleTiles):
        def __init__(self, cache):
            super().__init__(self, cache=cache)
            self.desired_tile_form = "RGB"

        # OpenTopoMap
        def _image_url(self, tile):
            x, y, z = tile
            return f"https://a.tile.opentopomap.org/{z}/{x}/{y}.png"

    class OTM_BW(OTM):  # BW OTM
        """OpenTopoMap tiler converted to grayscale."""

        def get_image(self, tile):
            img, extent, origin = super().get_image(tile)

            # Convert numpy array -> PIL Image if needed
            if isinstance(img, np.ndarray):
                img = Image.fromarray(img)

            # Now safe to convert
            img = img.convert("L").convert("RGB")

            return img, extent, origin

    _CARTOPY_AVAILABLE = True
except ImportError:
    pass

try:
    import geopandas

    _GEOPANDAS_AVAILABLE = True
except ImportError:
    warn(
        "geopandas not available, you won't be able to display geojson files",
        use_debug=False,
    )
    _GEOPANDAS_AVAILABLE = False


def parse_cartomap_style(cartomap):
    """
    Parse a cartomap string with inline styling options.

    The input string is expected to follow a pipe-separated format:

        "filename|key=value|key=value|..."

    where the first element is the file path (e.g. shapefile or GeoJSON),
    and subsequent elements define optional styling parameters.

    Supported style keys include:
        - "color" or "c" : line color
        - "lw" or "linewidth" : line width (float)
        - "ls" or "linestyle" : line style (e.g. "-", "--", ":")
        - "alpha" : transparency (float between 0 and 1)

    Missing or invalid values are ignored and replaced with defaults.

    Parameters
    ----------
    cartomap : str
        Cartomap specification string. The first element is the filename,
        followed by optional styling parameters separated by "|".

        Example:
            "borders.shp|color=#ffffff|lw=1.5|ls=--|alpha=0.6"

    Returns
    -------
    filename : str
        Path to the cartomap file.
    style : dict
        Dictionary containing parsed styling parameters with keys:
        "color", "linewidth", "linestyle", and "alpha".

    Invalid values (e.g. non-numeric linewidth) are silently ignored.
    """
    parts = [p.strip() for p in cartomap.split("|")]

    filename = parts[0]

    # Defaults
    style = {
        "color": "black",
        "lw": 1.2,
        "ls": "-",
        "alpha": 1.0,
    }

    for part in parts[1:]:
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        key = key.strip().lower()
        value = value.strip()

        if key in ("lw", "linewidth"):
            try:
                style["linewidth"] = float(value)
            except ValueError:
                pass
        elif key in ("ls", "linestyle"):
            style["linestyle"] = value
        elif key in ("color", "c"):
            style["color"] = value
        elif key == "alpha":
            try:
                style["alpha"] = float(value)
            except ValueError:
                pass

    return filename, style


def embellish_plot(ax, plot_config):
    """
    Add map features and background layers to a Cartopy axis.

    This function enhances a given axis by adding cartographic elements such as:
    background tiles, country borders, coastlines, rivers, and user-provided
    vector layers (e.g. shapefiles or GeoJSON files).

    The behavior is controlled via the `plot_config` dictionary.

    Parameters
    ----------
    ax : cartopy.mpl.geoaxes.GeoAxes
        Cartopy axis to which map features will be added.
    plot_config : dict
        Dictionary containing plotting options. Supported keys include:

        - "maps" : list of str
            List of map layers to add. Elements can be:
                * predefined keywords:
                    "relief", "OTM", "OTM_BW",
                    "countries", "provinces", "urban_areas",
                    "roads", "railroads", "coastlines",
                    "lakes", "lakes_europe",
                    "rivers", "rivers_europe"
                * file paths (".shp", ".geojson", ".json") optionally
                  with inline styling:
                    "file.shp|color=red|lw=1.5|ls=--|alpha=0.5"

        - "resolution" : {"10m", "50m", "110m"}, optional
            Resolution of Natural Earth features. Default is "10m".

        - "bg_alpha" : float, optional
            Transparency of background tiles (e.g. relief or OTM).
            Default is 1.0.

        - "background_zoom" : int, optional
            Zoom level for background tile providers. Default is 8.


    Returns
    -------
    None
        The function modifies the input axis in place.
    """
    maps_list = plot_config.get("maps", [])
    resolution = plot_config.get("resolution", "10m")
    bg_alpha = plot_config.get("bg_alpha", 1.0)
    background_zoom = plot_config.get("background_zoom", 8)

    if "relief" in maps_list and "OTM" in maps_list:
        warn(
            "Plotting both 'relief' and 'OTM' is not supported, choosing only relief",
            use_debug=False,
        )
        maps_list.remove("OTM")

    if "relief" in maps_list or "OTM" in maps_list or "OTM_BW" in maps_list:
        if "PYRAD_CACHE" in os.environ:
            cache_dir = os.environ["PYRAD_CACHE"]
        else:
            cache_dir = os.path.join(os.path.expanduser("~"), "pyrad_cache")

        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir)

    if "relief" in maps_list:
        tiler = ShadedReliefESRI(cache=cache_dir)

    if "OTM_BW" in maps_list:
        tiler = OTM_BW(cache=cache_dir)

    if "OTM" in maps_list:
        tiler = OTM(cache=cache_dir)

    for cartomap in maps_list:
        if cartomap in ["relief", "OTM", "OTM_BW"]:
            ax.add_image(tiler, background_zoom, alpha=bg_alpha)

        elif cartomap == "countries":
            # add countries
            countries = cartopy.feature.NaturalEarthFeature(
                category="cultural",
                name="admin_0_countries",
                scale=resolution,
                facecolor="none",
            )
            ax.add_feature(countries, edgecolor="black")
        elif cartomap == "provinces":
            # Create a feature for States/Admin 1 regions at
            # 1:resolution from Natural Earth
            states_provinces = cartopy.feature.NaturalEarthFeature(
                category="cultural",
                name="admin_1_states_provinces_lines",
                scale=resolution,
                facecolor="none",
            )
            ax.add_feature(states_provinces, edgecolor="gray")
        elif cartomap == "urban_areas" and resolution in ("10m", "50m"):
            urban_areas = cartopy.feature.NaturalEarthFeature(
                category="cultural", name="urban_areas", scale=resolution
            )
            ax.add_feature(
                urban_areas, edgecolor="brown", facecolor="brown", alpha=0.25
            )
        elif cartomap == "roads" and resolution == "10m":
            roads = cartopy.feature.NaturalEarthFeature(
                category="cultural", name="roads", scale=resolution
            )
            ax.add_feature(roads, edgecolor="red", facecolor="none")
        elif cartomap == "railroads" and resolution == "10m":
            railroads = cartopy.feature.NaturalEarthFeature(
                category="cultural", name="railroads", scale=resolution
            )
            ax.add_feature(
                railroads, edgecolor="green", facecolor="none", linestyle=":"
            )
        elif cartomap == "coastlines":
            ax.coastlines(resolution=resolution)
        elif cartomap == "lakes":
            # add lakes
            lakes = cartopy.feature.NaturalEarthFeature(
                category="physical", name="lakes", scale=resolution
            )
            ax.add_feature(lakes, edgecolor="blue", facecolor="blue", alpha=0.25)
        elif resolution == "10m" and cartomap == "lakes_europe":
            lakes_europe = cartopy.feature.NaturalEarthFeature(
                category="physical", name="lakes_europe", scale=resolution
            )
            ax.add_feature(lakes_europe, edgecolor="blue", facecolor="blue", alpha=0.25)
        elif cartomap == "rivers":
            # add rivers
            rivers = cartopy.feature.NaturalEarthFeature(
                category="physical",
                name="rivers_lake_centerlines",
                scale=resolution,
            )
            ax.add_feature(rivers, edgecolor="blue", facecolor="none")
        elif resolution == "10m" and cartomap == "rivers_europe":
            rivers_europe = cartopy.feature.NaturalEarthFeature(
                category="physical", name="rivers_europe", scale=resolution
            )
            ax.add_feature(rivers_europe, edgecolor="blue", facecolor="none")
        elif ".shp" in cartomap or ".geojson" in cartomap or ".json" in cartomap:
            if _GEOPANDAS_AVAILABLE:
                cartomap_file, style = parse_cartomap_style(cartomap)
                try:
                    gdf = geopandas.read_file(cartomap_file)

                    if gdf.crs is None:
                        xmin, _, _, _ = gdf.total_bounds
                        if xmin > 1e6:
                            gdf = gdf.set_crs(epsg=2056)
                        else:
                            gdf = gdf.set_crs(epsg=21781)

                    gdf = gdf.to_crs(epsg=4326)

                    ax.add_geometries(
                        gdf.geometry,
                        crs=cartopy.crs.PlateCarree(),
                        facecolor="none",
                        edgecolor=style["color"],
                        linewidth=style["linewidth"],
                        linestyle=style["linestyle"],
                        alpha=style["alpha"],
                    )
                except Exception as e:
                    warn(f"Failed to load cartomap {cartomap}: {e}")
            else:
                warn(f"Geopandas not available, not able to display map {cartomap}")
        else:
            warn(
                "cartomap "
                + cartomap
                + " for resolution "
                + resolution
                + " not available"
            )


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
