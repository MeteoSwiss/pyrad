"""
pyrad.graph.plots_grid
======================

Functions to plot data in a Cartesian grid format

.. autosummary::
    :toctree: generated/

    plot_surface
    plot_surface_raw
    plot_surface_contour
    plot_latitude_slice
    plot_longitude_slice
    plot_cross_section
    plot_dda_map
    plot_dda_slice
"""

from warnings import warn
from copy import deepcopy

import numpy as np

import matplotlib.pyplot as plt

import pyart

from .plots_aux import get_norm

try:
    import cartopy
    from cartopy.io.img_tiles import Stamen

    _CARTOPY_AVAILABLE = True
except ImportError:
    _CARTOPY_AVAILABLE = False

try:
    import pydda

    _PYDDA_AVAILABLE = True
except ImportError:
    _PYDDA_AVAILABLE = False

import matplotlib as mpl

mpl.use("Agg")

# Increase a bit font size
mpl.rcParams.update({"font.size": 16})
mpl.rcParams.update({"font.family": "sans-serif"})


def plot_surface(
    grid,
    field_name,
    level,
    prdcfg,
    fname_list,
    titl=None,
    alpha=None,
    ax=None,
    fig=None,
    display=None,
    save_fig=True,
    use_basemap=False,
):
    """
    plots a surface from gridded data

    Parameters
    ----------
    grid : Grid object
        object containing the gridded data to plot
    field_name : str
        name of the radar field to plot
    level : int
        level index
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    titl : str
        Plot title
    alpha : float or None
        Set the alpha transparency of the grid plot. Useful for
        overplotting radar over other datasets.
    ax : Axis
        Axis to plot on. if fig is None a new axis will be created
    fig : Figure
        Figure to add the colorbar to. If none a new figure will be created
    display : GridMapDisplay object
        The display used
    save_fig : bool
        if true save the figure. If false it does not close the plot and
        returns the handle to the figure

    Returns
    -------
    fname_list : list of str or
    fig, ax, display : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    dpi = prdcfg["gridMapImageConfig"].get("dpi", 72)
    vmin = prdcfg.get("vmin", None)
    vmax = prdcfg.get("vmax", None)
    mask_outside = prdcfg.get("mask_outside", False)

    norm, ticks, ticklabs = get_norm(
        field_name, field_dict=grid.fields[field_name], isxarray=False
    )

    xsize = prdcfg["gridMapImageConfig"]["xsize"]
    ysize = prdcfg["gridMapImageConfig"]["ysize"]
    lonstep = prdcfg["gridMapImageConfig"].get("lonstep", 0.5)
    latstep = prdcfg["gridMapImageConfig"].get("latstep", 0.5)
    min_lon = prdcfg["gridMapImageConfig"].get("lonmin", 2.5)
    max_lon = prdcfg["gridMapImageConfig"].get("lonmax", 12.5)
    min_lat = prdcfg["gridMapImageConfig"].get("latmin", 43.5)
    max_lat = prdcfg["gridMapImageConfig"].get("latmax", 49.5)
    embellish = prdcfg["gridMapImageConfig"].get("embellish", True)
    colorbar_flag = prdcfg["gridMapImageConfig"].get("colorbar_flag", True)
    exact_limits = prdcfg["gridMapImageConfig"].get("exact_limits", True)

    if exact_limits:
        lon_lines = np.arange(min_lon, max_lon + lonstep, lonstep)
        lat_lines = np.arange(min_lat, max_lat + latstep, latstep)
    else:
        lon_lines = np.arange(np.floor(min_lon), np.ceil(max_lon) + 1, lonstep)
        lat_lines = np.arange(np.floor(min_lat), np.ceil(max_lat) + 1, latstep)

    if use_basemap or not _CARTOPY_AVAILABLE:
        resolution = prdcfg["gridMapImageConfig"].get("mapres", "l")
        if resolution not in ("c", "l", "i", "h", "f"):
            warn("Unknown map resolution: " + resolution)
            resolution = "l"

        if resolution == "c":
            area_thresh = 10000
        elif resolution == "l":
            area_thresh = 1000
        elif resolution == "i":
            area_thresh = 100
        elif resolution == "h":
            area_thresh = 10
        elif resolution == "f":
            area_thresh = 1
    else:
        resolution = prdcfg["gridMapImageConfig"].get("mapres", "110m")
        # Map from basemap to cartopy notation
        if resolution == "l":
            resolution = "110m"
        elif resolution == "i":
            resolution = "50m"
        elif resolution == "h":
            resolution = "10m"

        if resolution not in ("110m", "50m", "10m"):
            warn("Unknown map resolution: " + resolution)
            resolution = "110m"
        background_zoom = prdcfg["gridMapImageConfig"].get("background_zoom", 8)
        maps_list = prdcfg["gridMapImageConfig"].get("maps", [])

    if fig is None:
        fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)

        if use_basemap or not _CARTOPY_AVAILABLE:
            display = pyart.graph.GridMapDisplayBasemap(grid)
            display.plot_basemap(
                lat_lines=lat_lines,
                lon_lines=lon_lines,
                resolution=resolution,
                area_thresh=area_thresh,
                auto_range=False,
                min_lon=min_lon,
                max_lon=max_lon,
                min_lat=min_lat,
                max_lat=max_lat,
            )
            display.plot_grid(
                field_name,
                level=level,
                norm=norm,
                ticks=ticks,
                title=titl,
                ticklabs=ticklabs,
                mask_outside=mask_outside,
                vmin=vmin,
                vmax=vmax,
                alpha=alpha,
                fig=fig,
            )
        else:
            projection = cartopy.crs.PlateCarree()

            display = pyart.graph.GridMapDisplay(grid)
            display.plot_grid(
                field_name,
                level=level,
                norm=norm,
                ticks=ticks,
                ticklabs=ticklabs,
                lat_lines=lat_lines,
                projection=projection,
                lon_lines=lon_lines,
                vmin=vmin,
                embellish=False,
                add_grid_lines=True,
                vmax=vmax,
                mask_outside=mask_outside,
                alpha=alpha,
                title=titl,
                colorbar_flag=colorbar_flag,
            )

            # fig = plt.gcf()
            # ax.set_extent([min_lon, max_lon, min_lat, max_lat])
            # display.plot_crosshairs(lon=lon, lat=lat)
    else:
        if use_basemap or not _CARTOPY_AVAILABLE:
            display.plot_grid(
                field_name,
                level=level,
                norm=norm,
                ticks=ticks,
                lat_lines=lat_lines,
                lon_lines=lon_lines,
                title=titl,
                ticklabs=ticklabs,
                colorbar_flag=False,
                embellish=False,
                vmin=vmin,
                vmax=vmax,
                mask_outside=mask_outside,
                alpha=alpha,
                ax=ax,
                fig=fig,
            )
        else:
            display.plot_grid(
                field_name,
                level=level,
                norm=norm,
                ticks=ticks,
                projection=ax.projection,
                lat_lines=lat_lines,
                lon_lines=lon_lines,
                ticklabs=ticklabs,
                colorbar_flag=False,
                embellish=False,
                vmin=vmin,
                vmax=vmax,
                mask_outside=mask_outside,
                alpha=alpha,
                title=titl,
                ax=ax,
                fig=fig,
            )

    if embellish:
        ax = plt.gca()
        if "relief" in maps_list:
            tiler = Stamen("terrain-background")
            projection = tiler.crs
            fig.delaxes(ax)
            ax = fig.add_subplot(111, projection=projection)
            warn(
                "The projection of the image is set to that of the "
                + "background map, i.e. "
                + str(projection),
                UserWarning,
            )

        for cartomap in maps_list:
            if cartomap == "relief":
                ax.add_image(tiler, background_zoom)
            if cartomap == "countries":
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
                ax.add_feature(
                    lakes_europe, edgecolor="blue", facecolor="blue", alpha=0.25
                )
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
            else:
                warn(
                    "cartomap "
                    + cartomap
                    + " for resolution "
                    + resolution
                    + " not available"
                )

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list

    return (fig, ax, display)


def plot_surface_raw(
    grid,
    field_name,
    level,
    prdcfg,
    fname_list,
    titl=None,
    alpha=None,
    ax=None,
    fig=None,
    display=None,
    save_fig=True,
):
    """
    plots a surface from gridded data within reprojecting the data into a map

    Parameters
    ----------
    grid : Grid object
        object containing the gridded data to plot
    field_name : str
        name of the radar field to plot
    level : int
        level index
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    titl : str
        Plot title
    alpha : float or None
        Set the alpha transparency of the grid plot. Useful for
        overplotting radar over other datasets.
    ax : Axis
        Axis to plot on. if fig is None a new axis will be created
    fig : Figure
        Figure to add the colorbar to. If none a new figure will be created
    display : GridMapDisplay object
        The display used
    save_fig : bool
        if true save the figure. If false it does not close the plot and
        returns the handle to the figure

    Returns
    -------
    fname_list : list of str or
    fig, ax, display : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    dpi = prdcfg["gridMapImageConfig"].get("dpi", 72)
    vmin = prdcfg.get("vmin", None)
    vmax = prdcfg.get("vmax", None)

    norm, ticks, ticklabs = get_norm(
        field_name, field_dict=grid.fields[field_name], isxarray=False
    )

    xsize = prdcfg["gridMapImageConfig"]["xsize"]
    ysize = prdcfg["gridMapImageConfig"]["ysize"]
    colorbar_flag = prdcfg["gridMapImageConfig"].get("colorbar_flag", True)

    if fig is None:
        fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
        ax = fig.add_subplot(111)

        display = pyart.graph.GridMapDisplay(grid)
        display.plot_grid_raw(
            field_name,
            level=level,
            norm=norm,
            ticks=ticks,
            ticklabs=ticklabs,
            vmin=vmin,
            vmax=vmax,
            alpha=alpha,
            title=titl,
            ax=ax,
            fig=fig,
            colorbar_flag=colorbar_flag,
        )
        # display.plot_crosshairs(lon=lon, lat=lat)
    else:
        display.plot_grid_raw(
            field_name,
            level=level,
            norm=norm,
            ticks=ticks,
            ticklabs=ticklabs,
            colorbar_flag=False,
            vmin=vmin,
            vmax=vmax,
            alpha=alpha,
            title=titl,
            ax=ax,
            fig=fig,
        )

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list

    return (fig, ax, display)


def plot_surface_contour(
    grid,
    field_name,
    level,
    prdcfg,
    fname_list,
    contour_values=None,
    linewidths=1.5,
    colors="k",
    ax=None,
    fig=None,
    display=None,
    save_fig=True,
    use_basemap=False,
):
    """
    plots a contour plot from gridded data

    Parameters
    ----------
    grid : Grid object
        object containing the gridded data to plot
    field_name : str
        name of the radar field to plot
    level : int
        level index
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    contour_values : float array
        list of contours to plot
    linewidths : float
        width of the contour lines
    colors : color string or sequence of colors
        The contour colours
    ax : Axis
        Axis to plot on. if fig is None a new axis will be created
    fig : Figure
        Figure to add the colorbar to. If none a new figure will be created
    display : GridMapDisplay object
        The display used
    save_fig : bool
        if true save the figure if false it does not close the plot and
        returns the handle to the figure
    use_basemap : Bool
        If true uses basemap, otherwise uses cartopy.

    Returns
    -------
    fname_list : list of str or
    fig, ax : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    # get contour intervals
    if contour_values is None:
        field_dict = pyart.config.get_metadata(field_name)
        if "boundaries" in field_dict:
            vmin = field_dict["boundaries"][0]
            vmax = field_dict["boundaries"][-1]
            num = len(field_dict["boundaries"])
        else:
            vmin, vmax = pyart.config.get_field_limits(field_name)
            num = 10

        contour_values = np.linspace(vmin, vmax, num=num)

    dpi = prdcfg["gridMapImageConfig"].get("dpi", 72)

    xsize = prdcfg["gridMapImageConfig"]["xsize"]
    ysize = prdcfg["gridMapImageConfig"]["ysize"]
    lonstep = prdcfg["gridMapImageConfig"].get("lonstep", 0.5)
    latstep = prdcfg["gridMapImageConfig"].get("latstep", 0.5)
    min_lon = prdcfg["gridMapImageConfig"].get("lonmin", 2.5)
    max_lon = prdcfg["gridMapImageConfig"].get("lonmax", 12.5)
    min_lat = prdcfg["gridMapImageConfig"].get("latmin", 43.5)
    max_lat = prdcfg["gridMapImageConfig"].get("latmax", 49.5)
    exact_limits = prdcfg["gridMapImageConfig"].get("exact_limits", 0)

    if fig is None:
        if exact_limits:
            lon_lines = np.arange(min_lon, max_lon + lonstep, lonstep)
            lat_lines = np.arange(min_lat, max_lat + latstep, latstep)
        else:
            lon_lines = np.arange(np.floor(min_lon), np.ceil(max_lon) + 1, lonstep)
            lat_lines = np.arange(np.floor(min_lat), np.ceil(max_lat) + 1, latstep)

        fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
        ax = fig.add_subplot(111)

        if use_basemap or not _CARTOPY_AVAILABLE:
            resolution = prdcfg["gridMapImageConfig"].get("mapres", "l")
            if resolution not in ("c", "l", "i", "h", "f"):
                warn("Unknown map resolution: " + resolution)
                resolution = "l"

            if resolution == "c":
                area_thresh = 10000
            elif resolution == "l":
                area_thresh = 1000
            elif resolution == "i":
                area_thresh = 100
            elif resolution == "h":
                area_thresh = 10
            elif resolution == "f":
                area_thresh = 1

            display = pyart.graph.GridMapDisplayBasemap(grid)
            display.plot_basemap(
                lat_lines=lat_lines,
                lon_lines=lon_lines,
                resolution=resolution,
                auto_range=False,
                area_thresh=area_thresh,
                min_lon=min_lon,
                max_lon=max_lon,
                min_lat=min_lat,
                max_lat=max_lat,
                ax=ax,
            )

            lons, lats = grid.get_point_longitude_latitude(edges=False)
            data = grid.fields[field_name]["data"][level, :, :]

            basemap = display.get_basemap()
            basemap.contour(
                lons,
                lats,
                data,
                contour_values,
                colors=colors,
                linewidths=linewidths,
                latlon=True,
            )
            ax.set_title(display.generate_grid_title(field_name, level))
        else:
            resolution = prdcfg["gridMapImageConfig"].get("mapres", "110m")
            # Map from basemap to cartopy notation
            if resolution == "l":
                resolution = "110m"
            elif resolution == "i":
                resolution = "50m"
            elif resolution == "h":
                resolution = "10m"
            if resolution not in ("110m", "50m", "10m"):
                warn("Unknown map resolution: " + resolution)
                resolution = "110m"
            background_zoom = prdcfg["gridMapImageConfig"].get("background_zoom", 8)
            maps_list = prdcfg["gridMapImageConfig"].get("maps", [])

            display = pyart.graph.GridMapDisplay(grid)
            fig, ax = display.plot_grid_contour(
                field_name,
                level=level,
                ax=ax,
                fig=fig,
                lat_lines=lat_lines,
                lon_lines=lon_lines,
                contour_values=contour_values,
                linewidths=linewidths,
                colors=colors,
                resolution=resolution,
                background_zoom=background_zoom,
                maps_list=maps_list,
            )
    else:
        if use_basemap or not _CARTOPY_AVAILABLE:
            lons, lats = grid.get_point_longitude_latitude(edges=False)
            data = grid.fields[field_name]["data"][level, :, :]

            basemap = display.get_basemap()
            basemap.contour(
                lons,
                lats,
                data,
                contour_values,
                colors=colors,
                linewidths=linewidths,
                latlon=True,
            )
        else:
            lons, lats = grid.get_point_longitude_latitude(edges=False)
            data = grid.fields[field_name]["data"][level, :, :]

            ax.contour(
                lons,
                lats,
                data,
                contour_values,
                colors=colors,
                linewidths=linewidths,
                transform=cartopy.crs.PlateCarree(),
            )
            ax.set_extent([min_lon, max_lon, min_lat, max_lat])

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list

    return (fig, ax)


def plot_latitude_slice(grid, field_name, lon, lat, prdcfg, fname_list):
    """
    plots a latitude slice from gridded data

    Parameters
    ----------
    grid : Grid object
        object containing the gridded data to plot
    field_name : str
        name of the radar field to plot
    lon, lat : float
        coordinates of the slice to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    dpi = 72
    if "dpi" in prdcfg["xsecImageConfig"]:
        dpi = prdcfg["xsecImageConfig"]["dpi"]

    norm, ticks, ticklabs = get_norm(
        field_name, field_dict=grid.fields[field_name], isxarray=True
    )

    xsize = prdcfg["xsecImageConfig"].get("xsize", 10.0)
    ysize = prdcfg["xsecImageConfig"].get("ysize", 5.0)
    xmin = prdcfg["xsecImageConfig"].get("xmin", None)
    xmax = prdcfg["xsecImageConfig"].get("xmax", None)
    ymin = prdcfg["xsecImageConfig"].get("ymin", None)
    ymax = prdcfg["xsecImageConfig"].get("ymax", None)
    vmin = prdcfg.get("vmin", None)
    vmax = prdcfg.get("vmax", None)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    ax = fig.add_subplot(111, aspect="equal")
    display = pyart.graph.GridMapDisplay(grid)
    display.plot_latitude_slice(
        field_name,
        lon=lon,
        lat=lat,
        norm=norm,
        colorbar_orient="horizontal",
        ticks=ticks,
        ticklabs=ticklabs,
        ax=ax,
        fig=fig,
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)


def plot_longitude_slice(grid, field_name, lon, lat, prdcfg, fname_list):
    """
    plots a longitude slice from gridded data

    Parameters
    ----------
    grid : Grid object
        object containing the gridded data to plot
    field_name : str
        name of the radar field to plot
    lon, lat : float
        coordinates of the slice to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    dpi = 72
    if "dpi" in prdcfg["xsecImageConfig"]:
        dpi = prdcfg["xsecImageConfig"]["dpi"]

    norm, ticks, ticklabs = get_norm(
        field_name, field_dict=grid.fields[field_name], isxarray=True
    )

    xsize = prdcfg["xsecImageConfig"].get("xsize", 10.0)
    ysize = prdcfg["xsecImageConfig"].get("ysize", 5.0)
    xmin = prdcfg["xsecImageConfig"].get("xmin", None)
    xmax = prdcfg["xsecImageConfig"].get("xmax", None)
    ymin = prdcfg["xsecImageConfig"].get("ymin", None)
    ymax = prdcfg["xsecImageConfig"].get("ymax", None)
    vmin = prdcfg.get("vmin", None)
    vmax = prdcfg.get("vmax", None)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    ax = fig.add_subplot(111, aspect="equal")
    display = pyart.graph.GridMapDisplay(grid)
    display.plot_longitude_slice(
        field_name,
        lon=lon,
        lat=lat,
        norm=norm,
        colorbar_orient="horizontal",
        ticks=ticks,
        ticklabs=ticklabs,
        ax=ax,
        fig=fig,
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)


def plot_cross_section(grid, field_name, coord1, coord2, prdcfg, fname_list):
    """
    plots a croos section crossing two points in the grid

    Parameters
    ----------
    grid : Grid object
        object containing the gridded data to plot
    field_name : str
        name of the radar field to plot
    coord1 : tupple of floats
        lat, lon of the first point
    coord2 : tupple of floats
        lat, lon of the second point
    fname_list : list of str
        list of names of the files where to store the plot

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    dpi = 72
    if "dpi" in prdcfg["xsecImageConfig"]:
        dpi = prdcfg["xsecImageConfig"]["dpi"]

    norm, ticks, ticklabs = get_norm(
        field_name, field_dict=grid.fields[field_name], isxarray=True
    )

    xsize = prdcfg["xsecImageConfig"].get("xsize", 10.0)
    ysize = prdcfg["xsecImageConfig"].get("ysize", 5.0)
    vmin = prdcfg.get("vmin", None)
    vmax = prdcfg.get("vmax", None)
    # xmin = prdcfg['xsecImageConfig'].get('xmin', None)
    # xmax = prdcfg['xsecImageConfig'].get('xmax', None)
    # ymin = prdcfg['xsecImageConfig'].get('ymin', None)
    # ymax = prdcfg['xsecImageConfig'].get('ymax', None)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    ax = fig.add_subplot(111, aspect="equal")
    display = pyart.graph.GridMapDisplay(grid)
    display.plot_cross_section(
        field_name,
        start=coord1,
        end=coord2,
        norm=norm,
        colorbar_orient="vertical",
        ticks=ticks,
        ticklabs=ticklabs,
        fig=fig,
        ax=ax,
        axislabels_flag=True,
        vmin=vmin,
        vmax=vmax,
    )
    # ax.set_ylim(
    #    [prdcfg['xsecImageConfig']['ymin'], prdcfg['xsecImageConfig']['ymax']])

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)


def plot_dda_map(
    grid,
    bg_field_name,
    level,
    prdcfg,
    fname_list,
    titl=None,
    alpha=None,
    ax=None,
    fig=None,
    display=None,
    save_fig=True,
    display_type="quiver",
):
    """
    This procedure plots a horizontal cross section of winds from wind fields
    generated by PyDDA.

    Parameters
    ----------
    grid : Grid object
        object containing the gridded data to plot
    bg_field_name : str
        name of the background radar field to plot (behind  the wind vectors)
    level : int
        level index
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    titl : str
        Plot title
    alpha : float or None
        Set the alpha transparency of the grid plot. Useful for
        overplotting radar over other datasets.
    ax : Axis
        Axis to plot on. if fig is None a new axis will be created
    fig : Figure
        Figure to add the colorbar to. If none a new figure will be created
    display : GridMapDisplay object
        The display used
    save_fig : bool
        if true save the figure. If false it does not close the plot and
        returns the handle to the figure
    display_type : str
        Display type for the wind vectors, can be either 'quiver', 'barbs' or
        'streamline'

    Returns
    -------
    fname_list : list of str or
    fig, ax, display : tupple
        list of names of the saved plots or handle of the figure an axes

    """

    if not _PYDDA_AVAILABLE:
        warn("PyDDA package not available. Unable to display wind fields")
        return None

    dpi = prdcfg["gridMapImageConfig"].get("dpi", 72)
    vmin = prdcfg.get("vmin", None)
    vmax = prdcfg.get("vmax", None)
    u_vel_contours = prdcfg.get("u_vel_contours", None)
    v_vel_contours = prdcfg.get("v_vel_contours", None)
    w_vel_contours = prdcfg.get("w_vel_contours", None)
    vector_spacing_km = prdcfg.get("vector_spacing_km", 10.0)
    quiver_len = prdcfg.get("quiver_len", 10.0)
    quiver_width = prdcfg.get("quiver_width", 0.01)
    streamline_width = prdcfg.get("streamline_width", None)
    streamline_arrowsize = prdcfg.get("streamline_arrowsize", None)
    display_type = prdcfg.get("display_type", "quiver")

    norm, ticks, ticklabs = get_norm(
        bg_field_name, field_dict=grid.fields[bg_field_name]
    )

    xsize = prdcfg["gridMapImageConfig"]["xsize"]
    ysize = prdcfg["gridMapImageConfig"]["ysize"]
    lonstep = prdcfg["gridMapImageConfig"].get("lonstep", 0.5)
    latstep = prdcfg["gridMapImageConfig"].get("latstep", 0.5)
    min_lon = prdcfg["gridMapImageConfig"].get("lonmin", 2.5)
    max_lon = prdcfg["gridMapImageConfig"].get("lonmax", 12.5)
    min_lat = prdcfg["gridMapImageConfig"].get("latmin", 43.5)
    max_lat = prdcfg["gridMapImageConfig"].get("latmax", 49.5)
    embellish = prdcfg["gridMapImageConfig"].get("embellish", True)
    exact_limits = prdcfg["gridMapImageConfig"].get("exact_limits", 0)
    colorbar_flag = prdcfg["gridMapImageConfig"].get("colorbar_flag", True)
    colorbar_contour_flag = prdcfg["gridMapImageConfig"].get(
        "colorbar_contour_flag", False
    )
    cmap = pyart.config.get_field_colormap(bg_field_name)

    if exact_limits:
        lon_lines = np.arange(min_lon, max_lon + lonstep, lonstep)
        lat_lines = np.arange(min_lat, max_lat + latstep, latstep)
    else:
        lon_lines = np.arange(np.floor(min_lon), np.ceil(max_lon) + 1, lonstep)
        lat_lines = np.arange(np.floor(min_lat), np.ceil(max_lat) + 1, latstep)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)

    resolution = prdcfg["gridMapImageConfig"].get("mapres", "110m")
    # Map from basemap to cartopy notation
    if resolution == "l":
        resolution = "110m"
    elif resolution == "i":
        resolution = "50m"
    elif resolution == "h":
        resolution = "10m"

    if resolution not in ("110m", "50m", "10m"):
        warn("Unknown map resolution: " + resolution)
        resolution = "110m"

    maps_list = prdcfg["gridMapImageConfig"].get("maps", [])

    if display_type == "quiver":
        ax = pydda.vis.plot_horiz_xsection_quiver_map(
            [pydda.io.read_from_pyart_grid(deepcopy(grid))],
            background_field=bg_field_name,
            level=level,
            show_lobes=True,
            bg_grid_no=0,
            vmin=vmin,
            vmax=vmax,
            u_vel_contours=u_vel_contours,
            v_vel_contours=v_vel_contours,
            w_vel_contours=w_vel_contours,
            quiver_spacing_x_km=vector_spacing_km,
            quiver_spacing_y_km=vector_spacing_km,
            quiverkey_len=quiver_len,
            quiver_width=quiver_width,
            colorbar_flag=colorbar_flag,
            colorbar_contour_flag=colorbar_contour_flag,
            u_field="eastward_wind_component",
            v_field="northward_wind_component",
            w_field="vertical_wind_component",
            title_flag=False,
            cmap=cmap,
        )

    elif display_type == "barbs":
        ax = pydda.vis.plot_horiz_xsection_barbs_map(
            [pydda.io.read_from_pyart_grid(deepcopy(grid))],
            background_field=bg_field_name,
            level=level,
            show_lobes=True,
            bg_grid_no=0,
            vmin=vmin,
            vmax=vmax,
            u_vel_contours=u_vel_contours,
            v_vel_contours=v_vel_contours,
            w_vel_contours=w_vel_contours,
            barb_spacing_x_km=vector_spacing_km,
            barb_spacing_y_km=vector_spacing_km,
            colorbar_flag=colorbar_flag,
            colorbar_contour_flag=colorbar_contour_flag,
            u_field="eastward_wind_component",
            v_field="northward_wind_component",
            w_field="vertical_wind_component",
            title_flag=False,
            cmap=cmap,
        )

    elif display_type == "streamline":
        ax = pydda.vis.plot_horiz_xsection_streamlines_map(
            [pydda.io.read_from_pyart_grid(deepcopy(grid))],
            background_field=bg_field_name,
            level=level,
            show_lobes=True,
            bg_grid_no=0,
            vmin=vmin,
            vmax=vmax,
            u_vel_contours=u_vel_contours,
            v_vel_contours=v_vel_contours,
            w_vel_contours=w_vel_contours,
            linewidth=streamline_width,
            arrowsize=streamline_arrowsize,
            colorbar_flag=colorbar_flag,
            colorbar_contour_flag=colorbar_contour_flag,
            u_field="eastward_wind_component",
            v_field="northward_wind_component",
            w_field="vertical_wind_component",
            title_flag=False,
            cmap=cmap,
        )

    # Edit parameters a posteriori
    # since they cannot be given to pyDDA
    ax.set_title(titl)
    ax.collections[0].set_alpha(alpha)
    if norm:
        ax.collections[0].set_norm(norm)
    if ticks:
        ax.collections[0].colorbar.set_ticks(ticks)
    if ticklabs:
        ax.collections[0].colorbar.set_ticklabels(ticklabs)

    ax.set_xticks(lon_lines)
    ax.set_yticks(lat_lines)
    ax.set_xlim([lon_lines[0], lon_lines[-1]])
    ax.set_ylim([lat_lines[0], lat_lines[-1]])

    if embellish and _CARTOPY_AVAILABLE:
        for cartomap in maps_list:
            if cartomap == "countries":
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
                ax.add_feature(
                    lakes_europe, edgecolor="blue", facecolor="blue", alpha=0.25
                )
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
            else:
                warn(
                    "cartomap "
                    + cartomap
                    + " for resolution "
                    + resolution
                    + " not available"
                )

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list
    return (fig, ax, display)


def plot_dda_slice(
    grid,
    bg_field_name,
    slice_type,
    level,
    prdcfg,
    fname_list,
    titl=None,
    alpha=None,
    ax=None,
    fig=None,
    display=None,
    save_fig=True,
    display_type="quiver",
    wind_vectors="hor",
):
    """
    This procedure plots a cross section of wind vectors from wind fields
    generated by PyDDA in the X-Z plane at a specified latitude or longitude

        Parameters
    ----------
    grid : Grid object
        object containing the gridded data to plot
    bg_field_name : str
        name of the background radar field to plot (behind  the wind vectors)
    slice_type : str
        type of slice, can be either "latitude" or "longitude"
    level : int
        slicing level in latitudinal or longitudinal direction
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    titl : str
        Plot title
    alpha : float or None
        Set the alpha transparency of the grid plot. Useful for
        overplotting radar over other datasets.
    ax : Axis
        Axis to plot on. if fig is None a new axis will be created
    fig : Figure
        Figure to add the colorbar to. If none a new figure will be created
    display : GridMapDisplay object
        The display used
    save_fig : bool
        if true save the figure. If false it does not close the plot and
        returns the handle to the figure
    display_type : str
        Display type for the wind vectors, can be either 'quiver', 'barbs' or
        'streamline'
    wind_vectors : str
        'hor' if horizontal wind vectors are displayed (u and v) or 'ver'
        if vertical wind vectors are displayed (u and w) for latitude slice
        and (v and w) for longitude slice

    Returns
    -------
    fname_list : list of str or
    fig, ax, display : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    if not _PYDDA_AVAILABLE:
        warn("PyDDA package not available. Unable to display wind fields")
        return None
    dpi = prdcfg["xsecImageConfig"].get("dpi", 72)
    vmin = prdcfg.get("vmin", None)
    vmax = prdcfg.get("vmax", None)
    u_vel_contours = prdcfg.get("u_vel_contours", None)
    v_vel_contours = prdcfg.get("v_vel_contours", None)
    w_vel_contours = prdcfg.get("w_vel_contours", None)
    vector_spacing_km = prdcfg.get("vector_spacing_km", 10.0)
    quiver_len = prdcfg.get("quiver_len", 10.0)
    quiver_width = prdcfg.get("quiver_width", 0.01)
    streamline_width = prdcfg.get("streamline_width", None)
    streamline_arrowsize = prdcfg.get("streamline_arrowsize", None)
    display_type = prdcfg.get("display_type", "quiver")

    norm, ticks, ticklabs = get_norm(
        bg_field_name, field_dict=grid.fields[bg_field_name]
    )

    colorbar_flag = prdcfg["xsecImageConfig"].get("colorbar_flag", True)
    colorbar_contour_flag = prdcfg["xsecImageConfig"].get(
        "colorbar_contour_flag", False
    )
    cmap = pyart.config.get_field_colormap(bg_field_name)

    xsize = prdcfg["xsecImageConfig"].get("xsize", 10.0)
    ysize = prdcfg["xsecImageConfig"].get("ysize", 5.0)
    xmin = prdcfg["xsecImageConfig"].get("xmin", None)
    xmax = prdcfg["xsecImageConfig"].get("xmax", None)
    ymin = prdcfg["xsecImageConfig"].get("ymin", None)
    ymax = prdcfg["xsecImageConfig"].get("ymax", None)

    # This is a bit hackish
    # pyDDA does not support plotting vertical profiles of u and v
    # only (v and w) or (u and w), so we trick it by assigning the v variable
    # to w

    if wind_vectors == "hor":
        w_field = "eastward_wind_component"
    else:
        w_field = "vertical_wind_component"

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)

    if slice_type == "latitude":
        if display_type == "quiver":
            ax = pydda.vis.plot_yz_xsection_quiver(
                [pydda.io.read_from_pyart_grid(deepcopy(grid))],
                background_field=bg_field_name,
                level=level,
                bg_grid_no=0,
                vmin=vmin,
                vmax=vmax,
                u_vel_contours=u_vel_contours,
                v_vel_contours=v_vel_contours,
                w_vel_contours=w_vel_contours,
                quiver_spacing_y_km=vector_spacing_km,
                quiver_spacing_z_km=vector_spacing_km,
                quiver_width=quiver_width,
                quiverkey_len=quiver_len,
                colorbar_flag=colorbar_flag,
                colorbar_contour_flag=colorbar_contour_flag,
                u_field="eastward_wind_component",
                v_field="northward_wind_component",
                w_field=w_field,
                title_flag=False,
                cmap=cmap,
            )

        elif display_type == "barbs":
            ax = pydda.vis.plot_yz_xsection_barbs(
                [pydda.io.read_from_pyart_grid(deepcopy(grid))],
                background_field=bg_field_name,
                level=level,
                bg_grid_no=0,
                vmin=vmin,
                vmax=vmax,
                u_vel_contours=u_vel_contours,
                v_vel_contours=v_vel_contours,
                w_vel_contours=w_vel_contours,
                barb_spacing_y_km=vector_spacing_km,
                barb_spacing_z_km=vector_spacing_km,
                colorbar_flag=colorbar_flag,
                colorbar_contour_flag=colorbar_contour_flag,
                u_field="eastward_wind_component",
                v_field="northward_wind_component",
                w_field=w_field,
                title_flag=False,
                cmap=cmap,
            )

        elif display_type == "streamline":
            ax = pydda.vis.plot_yz_xsection_streamlines(
                [pydda.io.read_from_pyart_grid(deepcopy(grid))],
                background_field=bg_field_name,
                level=level,
                bg_grid_no=0,
                vmin=vmin,
                vmax=vmax,
                u_vel_contours=u_vel_contours,
                v_vel_contours=v_vel_contours,
                w_vel_contours=w_vel_contours,
                linewidth=streamline_width,
                arrowsize=streamline_arrowsize,
                colorbar_flag=colorbar_flag,
                colorbar_contour_flag=colorbar_contour_flag,
                u_field="eastward_wind_component",
                v_field="northward_wind_component",
                w_field=w_field,
                title_flag=False,
                cmap=cmap,
            )

    elif slice_type == "longitude":
        if display_type == "quiver":
            ax = pydda.vis.plot_xz_xsection_quiver(
                [pydda.io.read_from_pyart_grid(deepcopy(grid))],
                background_field=bg_field_name,
                level=level,
                bg_grid_no=0,
                vmin=vmin,
                vmax=vmax,
                u_vel_contours=u_vel_contours,
                v_vel_contours=v_vel_contours,
                w_vel_contours=w_vel_contours,
                quiver_width=quiver_width,
                quiver_spacing_x_km=vector_spacing_km,
                quiver_spacing_z_km=vector_spacing_km,
                quiverkey_len=quiver_len,
                colorbar_flag=colorbar_flag,
                colorbar_contour_flag=colorbar_contour_flag,
                u_field="eastward_wind_component",
                v_field="northward_wind_component",
                w_field=w_field,
                title_flag=False,
                cmap=cmap,
            )

        elif display_type == "barbs":
            ax = pydda.vis.plot_xz_xsection_barbs(
                [pydda.io.read_from_pyart_grid(deepcopy(grid))],
                background_field=bg_field_name,
                level=level,
                bg_grid_no=0,
                vmin=vmin,
                vmax=vmax,
                u_vel_contours=u_vel_contours,
                v_vel_contours=v_vel_contours,
                w_vel_contours=w_vel_contours,
                barb_spacing_x_km=vector_spacing_km,
                barb_spacing_z_km=vector_spacing_km,
                colorbar_flag=colorbar_flag,
                colorbar_contour_flag=colorbar_contour_flag,
                u_field="eastward_wind_component",
                v_field="northward_wind_component",
                w_field=w_field,
                title_flag=False,
                cmap=cmap,
            )

        elif display_type == "streamline":
            ax = pydda.vis.plot_xz_xsection_streamlines(
                [pydda.io.read_from_pyart_grid(deepcopy(grid))],
                background_field=bg_field_name,
                level=level,
                bg_grid_no=0,
                vmin=vmin,
                vmax=vmax,
                u_vel_contours=u_vel_contours,
                v_vel_contours=v_vel_contours,
                w_vel_contours=w_vel_contours,
                linewidth=streamline_width,
                arrowsize=streamline_arrowsize,
                colorbar_flag=colorbar_flag,
                colorbar_contour_flag=colorbar_contour_flag,
                u_field="eastward_wind_component",
                v_field="northward_wind_component",
                w_field=w_field,
                title_flag=False,
                cmap=cmap,
            )

    # Edit parameters a posteriori
    # since they cannot be given to pyDDA
    ax.set_title(titl)
    ax.collections[0].set_alpha(alpha)
    if norm:
        ax.collections[0].set_norm(norm)
    if ticks:
        ax.collections[0].colorbar.set_ticks(ticks)
    if ticklabs:
        ax.collections[0].colorbar.set_ticklabels(ticklabs)

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list

    return (fig, ax, display)
