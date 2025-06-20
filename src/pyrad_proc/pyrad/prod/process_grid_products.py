"""
pyrad.prod.process_grid_products
================================

Functions for obtaining Pyrad products from gridded datasets

.. autosummary::
    :toctree: generated/

    generate_grid_time_avg_products
    generate_sparse_grid_products
    generate_grid_products

"""

from ..util import warn
from copy import deepcopy

import numpy as np

import pyart

from ..io.io_aux import get_fieldname_pyart
from ..io.io_aux import get_save_dir, make_filename
from ..io.write_data import write_histogram, write_ts_stats

from ..graph.plots_grid import plot_surface, plot_surface_raw
from ..graph.plots_grid import plot_longitude_slice, plot_latitude_slice
from ..graph.plots_grid import plot_cross_section, plot_surface_contour
from ..graph.plots_grid import plot_dda_map, plot_dda_slice

from ..graph.plots_aux import get_colobar_label, get_field_name
from ..graph.plots_aux import generate_dda_map_title
from ..graph.plots_aux import generate_dda_latitude_slice_title
from ..graph.plots_aux import generate_dda_longitude_slice_title
from ..graph.plots import plot_histogram, plot_pos

from ..util.radar_utils import compute_histogram


def generate_grid_time_avg_products(dataset, prdcfg):
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

    return generate_grid_products(dataset, prdcfg)


def generate_sparse_grid_products(dataset, prdcfg):
    """
    generates products defined by sparse points. Accepted product types:
        'SURFACE_IMAGE': Generates a surface image
            User defined parameters:
                'field_limits': list of floats
                    The limits of the surface to plot [deg]
                    lon0, lon1, lat0, lat1

    Parameters
    ----------
    dataset : dictionary containing the points and their values

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    no return

    """
    dssavedir = prdcfg["dsname"]
    if "dssavename" in prdcfg:
        dssavedir = prdcfg["dssavename"]

    if prdcfg["type"] == "SURFACE_IMAGE":
        field_name = get_fieldname_pyart(prdcfg["voltype"])
        if field_name not in dataset["radar_out"]["fields"]:
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
            prdcfg["prdname"],
            timeinfo=prdcfg["timeinfo"],
        )

        fname_list = make_filename(
            "surface",
            prdcfg["dstype"],
            prdcfg["voltype"],
            prdcfg["imgformat"],
            timeinfo=prdcfg["timeinfo"],
            runinfo=prdcfg["runinfo"],
        )

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir + fname

        cb_label = get_colobar_label(
            dataset["radar_out"]["fields"][field_name], field_name
        )

        titl = (
            prdcfg["timeinfo"].strftime("%Y-%m-%dT%H:%M%SZ")
            + "\n"
            + get_field_name(dataset["radar_out"]["fields"][field_name], field_name)
        )

        if "field_limits" in prdcfg:
            field_limits = prdcfg["field_limits"]
        else:
            field_limits = dataset["radar_out"]["field_limits"]

        # get colobar limits
        vmin, vmax = pyart.config.get_field_limits(field_name)
        plot_pos(
            dataset["radar_out"]["lat"],
            dataset["radar_out"]["lon"],
            dataset["radar_out"]["fields"][field_name]["data"],
            fname_list,
            cb_label=cb_label,
            titl=titl,
            limits=field_limits,
            vmin=vmin,
            vmax=vmax,
        )

        print("----- save to " + " ".join(fname_list))

        return fname_list

    warn(" Unsupported product type: " + prdcfg["type"])
    return None


def generate_grid_products(dataset, prdcfg):
    """
    generates grid products. Accepted product types:
        'CROSS_SECTION': Plots a cross-section of gridded data
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                coord1, coord2: dict
                    The two lat-lon coordinates marking the limits. They have
                    the keywords 'lat' and 'lon' [degree]. The altitude limits
                    are defined by the parameters in 'xsecImageConfig' in the
                    'loc' configuration file
        'HISTOGRAM': Computes a histogram of the radar volum data
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                step: float or None
                    the data quantization step. If none it will be obtained
                    from the Py-ART configuration file
                vmin, vmax: float or None
                    The minimum and maximum values. If None they will be
                    obtained from the Py-ART configuration file
                mask_val: float or None
                    A value to mask.
                write_data: Bool
                    If true the histogram data is written in a csv file
        'LATITUDE_SLICE': Plots a cross-section of gridded data over a
            constant latitude.
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                lon, lat: floats
                    The starting point of the cross-section. The ending point
                    is defined by the parameters in 'xsecImageConfig' in the
                    'loc' configuration file
        'LONGITUDE_SLICE': Plots a cross-ection of gridded data over a
            constant longitude.
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                lon, lat: floats
                    The starting point of the cross-section. The ending point
                    is defined by the parameters in 'xsecImageConfig' in the
                    'loc' configuration file
        'SAVEALL': Saves a gridded data object including all or a list of
            user-defined fields in a netcdf file
            User defined parameters:
                datatypes: list of str or None
                    The list of data types to save. If it is None, all fields
                    in the radar object will be saved
        'SAVEALL_GRID' : Same as before but can be used in a mixed GRID/VOL
            dataset, as there is no ambiguity with SAVEALL for VOL datasets
        'SAVEVOL': Saves on field of a gridded data object in a netcdf file.
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                file_type: str
                    The type of file used to save the data. Can be 'nc' or
                    'h5'. Default 'nc'
                physical: Bool
                    If True the data will be saved in physical units (floats).
                    Otherwise it will be quantized and saved as binary.
                    Default True
                compression: str
                    For ODIM file formats, the type of compression. Can be any
                    of the allowed compression types for hdf5 files. Default
                    gzip
                compression_opts: any
                    The compression options allowed by the hdf5. Depends on
                    the type of compression. Default 6 (The gzip compression
                    level).
        'SAVEVOL_GRID' : Same as before but can be used in a mixed GRID/VOL
            dataset, as there is no ambiguity with SAVEVOL for VOL datasets
        'STATS': Computes statistics over the whole images and stores them in
            a file.
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                stat: str
                    The statistic used. Can be mean, median, min, max
        'SURFACE_RAW': Plots a surface image of gridded data without
            projecting it into a map
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                level: int
                    The altitude level to plot. The rest of the parameters are
                    defined by the parameters in 'ppiImageConfig' and
                    'ppiMapImageConfig' in the 'loc' configuration file
        'SURFACE_IMAGE': Plots a surface image of gridded data.
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                level: int
                    The altitude level to plot. The rest of the parameters are
                    defined by the parameters in 'ppiImageConfig' and
                    'ppiMapImageConfig' in the 'loc' configuration file
        'SURFACE_CONTOUR': Plots a surface image of contour gridded data.
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                level: int
                    The altitude level to plot. The rest of the parameters are
                    defined by the parameters in 'ppiImageConfig' and
                    'ppiMapImageConfig' in the 'loc' configuration file
                contour_values : float array or None
                    The contour values. If None the values are taken from the
                    'boundaries' keyword in the field description in the
                    Py-ART config file. If 'boundaries' is not set the
                    countours are 10 values linearly distributed from vmin to
                    vmax
                linewidths : float
                    width of the contour lines
                colors : color string or sequence of colors
                    The contour colours
        'SURFACE_CONTOUR_OVERPLOT':
            Plots a surface image of gridded data with a contour overplotted.
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                level: int
                    The altitude level to plot. The rest of the parameters are
                    defined by the parameters in 'ppiImageConfig' and
                    'ppiMapImageConfig' in the 'loc' configuration file
                contour_values : float array or None
                    The contour values. If None the values are taken from the
                    'boundaries' keyword in the field description in the
                    Py-ART config file. If 'boundaries' is not set the
                    countours are 10 values linearly distributed from vmin to
                    vmax
                linewidths : float
                    width of the contour lines
                colors : color string or sequence of colors
                    The contour colours
        'SURFACE_OVERPLOT':
            Plots on the same surface two images, one on top of the other.
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                level: int
                    The altitude level to plot. The rest of the parameters are
                    defined by the parameters in 'ppiImageConfig' and
                    'ppiMapImageConfig' in the 'loc' configuration file
                contour_values : float array or None
                    The contour values. If None the values are taken from the
                    'boundaries' keyword in the field description in the
                    Py-ART config file. If 'boundaries' is not set the
                    countours are 10 values linearly distributed from vmin to
                    vmax
        'DDA_MAP':
            Plots horizontal wind vectors obtained from a DDA analysis on a 2D map.
            The pyDDA package is required
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                level: int
                    The altitude level to plot. The rest of the parameters are
                    defined by the parameters in 'ppiImageConfig' and
                    'ppiMapImageConfig' in the 'loc' configuration file
                display_type: str
                    Display method for the wind vectors, can be either
                    'streamline', 'quiver' or 'barbs'
                bg_ref_rad: int
                    Which radar to use as reference to display the background
                    voltype.
                u_vel_contours: list of int
                    The contours to use for plotting contours of u. Set to None to not
                    display such contours.
                v_vel_contours: list of int
                    The contours to use for plotting contours of v. Set to None to not
                    display such contours.
                w_vel_contours: list of int
                    The contours to use for plotting contours of w. Set to None to not
                    display such contours.
                vector_spacing_km: float
                    Spacing in km between wind vectors in x and y axis
                    (only used for barbs and quiver plots)
                quiver_len: float
                    Length to use for the quiver key in m/s.
                    (only used for quiver plots)
                streamline_arrowsize: float
                    Factor scale arrow size for streamlines.
                    (only used for streamline plots)
                linewidth: float
                    Linewidths for streamlines.
                    (only used for streamline plots)
        'DDA_LONGITUDE_SLICE':
            Plots horizontal wind vectors obtained from a DDA analysis on a vertical transect
            along a given longitude
            The pyDDA package is required
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                level: int
                    The longitude level to plot. The rest of the parameters are
                    defined by the parameters in 'xsecImageConfig' in the 'loc'
                    configuration file
                display_type: str
                    Display method for the wind vectors, can be either
                    'streamline', 'quiver' or 'barbs'
                bg_ref_rad: int
                    Which radar to use as reference to display the background
                    voltype.
                u_vel_contours: list of int
                    The contours to use for plotting contours of u. Set to None to not
                    display such contours.
                v_vel_contours: list of int
                    The contours to use for plotting contours of v. Set to None to not
                    display such contours.
                w_vel_contours: list of int
                    The contours to use for plotting contours of w. Set to None to not
                    display such contours.
                vector_spacing_km: float
                    Spacing in km between wind vectors in x and y axis
                    (only used for barbs and quiver plots)
                quiver_len: float
                    Length to use for the quiver key in m/s.
                    (only used for quiver plots)
                streamline_arrowsize: float
                    Factor scale arrow size for streamlines.
                    (only used for streamline plots)
                linewidth: float
                    Linewidths for streamlines.
                    (only used for streamline plots)
        'DDA_LATITUDE_SLICE':
            Plots horizontal wind vectors obtained from a DDA analysis on a vertical transect
            along a given latitude
            The pyDDA package is required
            User defined parameters:
                voltype: str
                    name of the pyrad variable to use, it must be available in the dataset
                level: int
                    The latitude level to plot. The rest of the parameters are
                    defined by the parameters in 'xsecImageConfig' in the 'loc'
                    configuration file
                display_type: str
                    Display method for the wind vectors, can be either
                    'streamline', 'quiver' or 'barbs'
                bg_ref_rad: int
                    Which radar to use as reference to display the background
                    voltype.
                u_vel_contours: list of int
                    The contours to use for plotting contours of u. Set to None to not
                    display such contours.
                v_vel_contours: list of int
                    The contours to use for plotting contours of v. Set to None to not
                    display such contours.
                w_vel_contours: list of int
                    The contours to use for plotting contours of w. Set to None to not
                    display such contours.
                vector_spacing_km: float
                    Spacing in km between wind vectors in x and y axis
                    (only used for barbs and quiver plots)
                quiver_len: float
                    Length to use for the quiver key in m/s.
                    (only used for quiver plots)
                streamline_arrowsize: float
                    Factor scale arrow size for streamlines.
                    (only used for streamline plots)
                linewidth: float
                    Linewidths for streamlines.
                    (only used for streamline plots)
    Parameters
    ----------
    dataset : grid
        grid object

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    None or name of generated files

    """

    dssavedir = prdcfg["dsname"]
    if "dssavename" in prdcfg:
        dssavedir = prdcfg["dssavename"]

    if prdcfg["type"] == "STATS":
        field_name = get_fieldname_pyart(prdcfg["voltype"])
        if field_name not in dataset["radar_out"].fields:
            warn(
                " Field type "
                + field_name
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        # user defined values
        stat = prdcfg.get("stat", "mean")

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=None,
        )

        fname = make_filename(
            "stats",
            prdcfg["dstype"],
            prdcfg["voltype"],
            ["csv"],
            prdcfginfo=stat,
            timeinfo=None,
            runinfo=prdcfg["runinfo"],
        )[0]

        fname = savedir + fname

        if stat == "mean":
            value = np.ma.masked_all(1)
            value[0] = np.ma.mean(dataset["radar_out"].fields[field_name]["data"])
        elif stat == "median":
            value = np.ma.masked_all(1)
            value[0] = np.ma.median(dataset["radar_out"].fields[field_name]["data"])
        elif stat == "min":
            value = np.ma.masked_all(1)
            value[0] = np.ma.min(dataset["radar_out"].fields[field_name]["data"])
        elif stat == "max":
            value = np.ma.masked_all(1)
            value[0] = np.ma.max(dataset["radar_out"].fields[field_name]["data"])
        else:
            warn("Unsupported statistic " + stat)
            return None

        write_ts_stats(prdcfg["timeinfo"], value, fname, stat=stat)

        print("----- save to " + fname)

        return fname

    if prdcfg["type"] == "SURFACE_RAW":
        field_name = get_fieldname_pyart(prdcfg["voltype"])
        if field_name not in dataset["radar_out"].fields:
            warn(
                " Field type "
                + field_name
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        # user defined values
        level = prdcfg.get("level", 0)

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=prdcfg["timeinfo"],
        )

        fname_list = make_filename(
            "surface",
            prdcfg["dstype"],
            prdcfg["voltype"],
            prdcfg["imgformat"],
            prdcfginfo="l" + str(level),
            timeinfo=prdcfg["timeinfo"],
            runinfo=prdcfg["runinfo"],
        )

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir + fname

        plot_surface_raw(dataset["radar_out"], field_name, level, prdcfg, fname_list)

        print("----- save to " + " ".join(fname_list))

        return fname_list

    if prdcfg["type"] == "SURFACE_IMAGE":
        field_name = get_fieldname_pyart(prdcfg["voltype"])
        if field_name not in dataset["radar_out"].fields:
            warn(
                " Field type "
                + field_name
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        # user defined values
        level = prdcfg.get("level", 0)

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=prdcfg["timeinfo"],
        )

        fname_list = make_filename(
            "surface",
            prdcfg["dstype"],
            prdcfg["voltype"],
            prdcfg["imgformat"],
            prdcfginfo="l" + str(level),
            timeinfo=prdcfg["timeinfo"],
            runinfo=prdcfg["runinfo"],
        )

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir + fname

        plot_surface(dataset["radar_out"], field_name, level, prdcfg, fname_list)

        print("----- save to " + " ".join(fname_list))

        return fname_list

    if prdcfg["type"] == "SURFACE_CONTOUR":
        field_name = get_fieldname_pyart(prdcfg["voltype"])
        if field_name not in dataset["radar_out"].fields:
            warn(
                " Field type "
                + field_name
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        # user defined values
        contour_values = prdcfg.get("contour_values", None)
        linewidths = prdcfg.get("linewidths", 1.5)
        colors = prdcfg.get("colors", "k")
        level = prdcfg.get("level", 0)

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=prdcfg["timeinfo"],
        )

        fname_list = make_filename(
            "surface",
            prdcfg["dstype"],
            prdcfg["voltype"],
            prdcfg["imgformat"],
            prdcfginfo="l" + str(level),
            timeinfo=prdcfg["timeinfo"],
            runinfo=prdcfg["runinfo"],
        )

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir + fname

        plot_surface_contour(
            dataset["radar_out"],
            field_name,
            level,
            prdcfg,
            fname_list,
            contour_values=contour_values,
            linewidths=linewidths,
            colors=colors,
        )

        print("----- save to " + " ".join(fname_list))

        return fname_list

    if prdcfg["type"] == "SURFACE_CONTOUR_OVERPLOT":
        field_name = get_fieldname_pyart(prdcfg["voltype"])
        if field_name not in dataset["radar_out"].fields:
            warn(
                " Field type "
                + field_name
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        contour_name = get_fieldname_pyart(prdcfg["contourtype"])
        if contour_name not in dataset["radar_out"].fields:
            warn(
                "Contour type "
                + contour_name
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        contour_values = prdcfg.get("contour_values", None)
        linewidths = prdcfg.get("linewidths", 1.5)
        colors = prdcfg.get("colors", "k")
        level = prdcfg.get("level", 0)
        contour_level = prdcfg.get("contour_level", level)

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=prdcfg["timeinfo"],
        )

        fname_list = make_filename(
            "surface-contour",
            prdcfg["dstype"],
            prdcfg["voltype"] + "-" + prdcfg["contourtype"],
            prdcfg["imgformat"],
            prdcfginfo="l" + str(level),
            timeinfo=prdcfg["timeinfo"],
            runinfo=prdcfg["runinfo"],
        )

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir + fname

        titl = (
            pyart.graph.common.generate_grid_title(
                dataset["radar_out"], field_name, level
            )
            + " - "
            + pyart.graph.common.generate_field_name(dataset["radar_out"], contour_name)
        )

        fig, ax, display = plot_surface(
            dataset["radar_out"],
            field_name,
            level,
            prdcfg,
            fname_list,
            titl=titl,
            save_fig=False,
        )

        fname_list = plot_surface_contour(
            dataset["radar_out"],
            contour_name,
            contour_level,
            prdcfg,
            fname_list,
            contour_values=contour_values,
            linewidths=linewidths,
            colors=colors,
            ax=ax,
            fig=fig,
            display=display,
            save_fig=True,
        )

        print("----- save to " + " ".join(fname_list))

        return fname_list

    if prdcfg["type"] == "SURFACE_OVERPLOT":
        field_name_btm = get_fieldname_pyart(prdcfg["voltype_btm"])
        if field_name_btm not in dataset["radar_out"].fields:
            warn(
                " Field type "
                + field_name_btm
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        field_name_top = get_fieldname_pyart(prdcfg["voltype_top"])
        if field_name_top not in dataset["radar_out"].fields:
            warn(
                "Contour type "
                + field_name_top
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        level_btm = prdcfg.get("level_btm", 0)
        level_top = prdcfg.get("level_top", 0)
        alpha = prdcfg.get("alpha", None)

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=prdcfg["timeinfo"],
        )

        fname_list = make_filename(
            "surface_overplot",
            prdcfg["dstype"],
            prdcfg["voltype_btm"] + "-" + prdcfg["voltype_top"],
            prdcfg["imgformat"],
            prdcfginfo="l" + str(level_btm),
            timeinfo=prdcfg["timeinfo"],
            runinfo=prdcfg["runinfo"],
        )

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir + fname

        titl = (
            pyart.graph.common.generate_grid_title(
                dataset["radar_out"], field_name_btm, level_btm
            )
            + " - "
            + pyart.graph.common.generate_field_name(
                dataset["radar_out"], field_name_top
            )
        )

        fig, ax, display = plot_surface(
            dataset["radar_out"],
            field_name_btm,
            level_btm,
            prdcfg,
            fname_list,
            titl=titl,
            save_fig=False,
        )

        fname_list = plot_surface(
            dataset["radar_out"],
            field_name_top,
            level_top,
            prdcfg,
            fname_list,
            titl=titl,
            alpha=alpha,
            ax=ax,
            fig=fig,
            display=display,
            save_fig=True,
        )

        print("----- save to " + " ".join(fname_list))

        return fname_list

    if prdcfg["type"] == "LATITUDE_SLICE":
        field_name = get_fieldname_pyart(prdcfg["voltype"])
        if field_name not in dataset["radar_out"].fields:
            warn(
                " Field type "
                + field_name
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        # user defined values
        lon = prdcfg.get("lon", dataset["radar_out"].origin_longitude["data"][0])
        lat = prdcfg.get("lat", dataset["radar_out"].origin_latitude["data"][0])

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=prdcfg["timeinfo"],
        )

        fname_list = make_filename(
            "lat_slice",
            prdcfg["dstype"],
            prdcfg["voltype"],
            prdcfg["imgformat"],
            prdcfginfo="lat" + "{:.2f}".format(lat),
            timeinfo=prdcfg["timeinfo"],
            runinfo=prdcfg["runinfo"],
        )

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir + fname

        plot_latitude_slice(
            dataset["radar_out"], field_name, lon, lat, prdcfg, fname_list
        )

        print("----- save to " + " ".join(fname_list))

        return fname_list

    if prdcfg["type"] == "LONGITUDE_SLICE":
        field_name = get_fieldname_pyart(prdcfg["voltype"])
        if field_name not in dataset["radar_out"].fields:
            warn(
                " Field type "
                + field_name
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        # user defined values
        lon = prdcfg.get("lon", dataset["radar_out"].origin_longitude["data"][0])
        lat = prdcfg.get("lat", dataset["radar_out"].origin_latitude["data"][0])

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=prdcfg["timeinfo"],
        )

        fname_list = make_filename(
            "lon_slice",
            prdcfg["dstype"],
            prdcfg["voltype"],
            prdcfg["imgformat"],
            prdcfginfo="lon" + "{:.2f}".format(lon),
            timeinfo=prdcfg["timeinfo"],
            runinfo=prdcfg["runinfo"],
        )

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir + fname

        plot_longitude_slice(
            dataset["radar_out"], field_name, lon, lat, prdcfg, fname_list
        )

        print("----- save to " + " ".join(fname_list))

        return fname_list

    if prdcfg["type"] == "CROSS_SECTION":
        field_name = get_fieldname_pyart(prdcfg["voltype"])
        if field_name not in dataset["radar_out"].fields:
            warn(
                " Field type "
                + field_name
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        # user defined values
        lon1 = dataset["radar_out"].point_longitude["data"][0, 0, 0]
        lat1 = dataset["radar_out"].point_latitude["data"][0, 0, 0]

        lon2 = dataset["radar_out"].point_longitude["data"][0, -1, -1]
        lat2 = dataset["radar_out"].point_latitude["data"][0, -1, -1]
        if "coord1" in prdcfg:
            if "lon" in prdcfg["coord1"]:
                lon1 = prdcfg["coord1"]["lon"]
            if "lat" in prdcfg["coord1"]:
                lat1 = prdcfg["coord1"]["lat"]
        if "coord2" in prdcfg:
            if "lon" in prdcfg["coord2"]:
                lon2 = prdcfg["coord2"]["lon"]
            if "lat" in prdcfg["coord2"]:
                lat2 = prdcfg["coord2"]["lat"]

        coord1 = (lat1, lon1)
        coord2 = (lat2, lon2)

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=prdcfg["timeinfo"],
        )

        fname_list = make_filename(
            "lonlat",
            prdcfg["dstype"],
            prdcfg["voltype"],
            prdcfg["imgformat"],
            prdcfginfo="lon-lat1_"
            + "{:.2f}".format(lon1)
            + "-"
            + "{:.2f}".format(lat1)
            + "_lon-lat2_"
            + "{:.2f}".format(lon2)
            + "-"
            + "{:.2f}".format(lat2),
            timeinfo=prdcfg["timeinfo"],
            runinfo=prdcfg["runinfo"],
        )

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir + fname

        plot_cross_section(
            dataset["radar_out"], field_name, coord1, coord2, prdcfg, fname_list
        )

        print("----- save to " + " ".join(fname_list))

        return fname_list

    if prdcfg["type"] == "HISTOGRAM":
        field_name = get_fieldname_pyart(prdcfg["voltype"])
        if field_name not in dataset["radar_out"].fields:
            warn(
                " Field type "
                + field_name
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        step = prdcfg.get("step", None)
        vmin = prdcfg.get("vmin", None)
        vmax = prdcfg.get("vmax", None)
        mask_val = prdcfg.get("mask_val", None)
        write_data = prdcfg.get("write_data", 0)
        binwidth_equal = prdcfg.get("binwidth_equal", 0)

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=prdcfg["timeinfo"],
        )

        fname_list = make_filename(
            "histogram",
            prdcfg["dstype"],
            prdcfg["voltype"],
            prdcfg["imgformat"],
            timeinfo=prdcfg["timeinfo"],
            runinfo=prdcfg["runinfo"],
        )

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir + fname

        values = dataset["radar_out"].fields[field_name]["data"]
        if mask_val is not None:
            values = np.ma.masked_values(values, mask_val)
        bin_edges, values = compute_histogram(
            values, field_name, step=step, vmin=vmin, vmax=vmax
        )

        titl = (
            pyart.graph.common.generate_grid_time_begin(
                dataset["radar_out"]
            ).isoformat()
            + "Z"
            + "\n"
            + get_field_name(dataset["radar_out"].fields[field_name], field_name)
        )

        labelx = get_colobar_label(dataset["radar_out"].fields[field_name], field_name)

        plot_histogram(
            bin_edges,
            values,
            fname_list,
            labelx=labelx,
            labely="Number of Samples",
            titl=titl,
            binwidth_equal=binwidth_equal,
        )

        print("----- save to " + " ".join(fname_list))

        if write_data:
            fname = (
                savedir
                + make_filename(
                    "histogram",
                    prdcfg["dstype"],
                    prdcfg["voltype"],
                    ["csv"],
                    timeinfo=prdcfg["timeinfo"],
                    runinfo=prdcfg["runinfo"],
                )[0]
            )

            hist, _ = np.histogram(values, bins=bin_edges)
            write_histogram(
                bin_edges, hist, fname, datatype=prdcfg["voltype"], step=step
            )
            print("----- save to " + fname)

            return fname

        return fname_list

    if prdcfg["type"] == "SAVEVOL" or prdcfg["type"] == "SAVEVOL_GRID":
        field_name = get_fieldname_pyart(prdcfg["voltype"])
        if field_name not in dataset["radar_out"].fields:
            warn(
                " Field type "
                + field_name
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        file_type = prdcfg.get("file_type", "nc")
        prdcfg.get("physical", True)
        compression = prdcfg.get("compression", "gzip")
        compression_opts = prdcfg.get("compression_opts", 6)

        new_dataset = deepcopy(dataset["radar_out"])
        new_dataset.fields = dict()
        new_dataset.add_field(field_name, dataset["radar_out"].fields[field_name])

        # This key is required for pyDDA processing but poses
        # problem when saving to CFRadial
        if "additional_radars" in new_dataset.metadata:
            new_dataset.metadata.pop("additional_radars")

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=prdcfg["timeinfo"],
        )

        fname = make_filename(
            "savevol",
            prdcfg["dstype"],
            prdcfg["voltype"],
            [file_type],
            timeinfo=prdcfg["timeinfo"],
            runinfo=prdcfg["runinfo"],
        )[0]

        fname = savedir + fname

        if file_type == "nc":
            pyart.io.write_grid(
                fname, new_dataset, write_point_x_y_z=True, write_point_lon_lat_alt=True
            )
        elif file_type == "h5":
            pyart.aux_io.write_odim_grid_h5(
                fname,
                new_dataset,
                compression=compression,
                compression_opts=compression_opts,
            )
        else:
            warn("Data could not be saved. " + "Unknown saving file type " + file_type)
            return None

        print("saved file: " + fname)

        return fname

    if prdcfg["type"] == "SAVEALL" or prdcfg["type"] == "SAVEALL_GRID":
        datatypes = prdcfg.get("datatypes", None)

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=prdcfg["timeinfo"],
        )

        fname = make_filename(
            "savevol",
            prdcfg["dstype"],
            "all_fields",
            ["nc"],
            timeinfo=prdcfg["timeinfo"],
            runinfo=prdcfg["runinfo"],
        )[0]

        fname = savedir + fname

        field_names = None
        if datatypes is not None:
            field_names = []
            for datatype in datatypes:
                field_names.append(get_fieldname_pyart(datatype))

        if field_names is not None:
            new_dataset = deepcopy(dataset["radar_out"])
            new_dataset.fields = dict()
            for field_name in field_names:
                if field_name not in dataset["radar_out"].fields:
                    warn(field_name + " not in grid object")
                else:
                    new_dataset.add_field(
                        field_name, dataset["radar_out"].fields[field_name]
                    )
        else:
            new_dataset = dataset["radar_out"]

        pyart.io.write_grid(
            fname, new_dataset, write_point_x_y_z=True, write_point_lon_lat_alt=True
        )
        print("saved file: " + fname)

        return fname

    if prdcfg["type"] == "DDA_MAP":
        display_type = prdcfg.get("display_type", "quiver")
        if display_type not in ["quiver", "barbs", "streamline"]:
            warn("Invalid display_type, must be " + '"streamline", "quiver" or "barbs"')
            return None

        bg_field_name = get_fieldname_pyart(prdcfg["voltype"])
        bg_ref_rad = prdcfg.get("bg_ref_rad", 0)
        # get bg field name given which radar serves as reference
        if bg_ref_rad != 0:
            bg_field_name += "_radar{:d}".format(bg_ref_rad)
        if bg_field_name not in dataset["radar_out"].fields:
            warn(
                " Field type "
                + bg_field_name
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        # user defined values
        level = prdcfg.get("level", 0)
        alpha = prdcfg.get("alpha", 1)

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=prdcfg["timeinfo"],
        )

        fname_list = make_filename(
            "dda_map",
            prdcfg["dstype"],
            prdcfg["voltype"],
            prdcfg["imgformat"],
            prdcfginfo="l" + str(level),
            timeinfo=prdcfg["timeinfo"],
            runinfo=prdcfg["runinfo"],
        )

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir + fname

        titl = generate_dda_map_title(dataset["radar_out"], bg_field_name, level)

        plot_dda_map(
            dataset["radar_out"],
            bg_field_name,
            level,
            prdcfg,
            fname_list,
            alpha=alpha,
            display_type=display_type,
            titl=titl,
        )

        print("----- save to " + " ".join(fname_list))

    if prdcfg["type"] == "DDA_LONGITUDE_SLICE":
        display_type = prdcfg.get("display_type", "quiver")
        if display_type not in ["quiver", "barbs", "streamline"]:
            warn("Invalid display_type, must be " + '"streamline", "quiver" or "barbs"')
            return None

        bg_field_name = get_fieldname_pyart(prdcfg["voltype"])
        bg_ref_rad = prdcfg.get("bg_ref_rad", 0)
        # get bg field name given which radar serves as reference
        if bg_ref_rad != 0:
            bg_field_name += "_radar{:d}".format(bg_ref_rad)

        if bg_field_name not in dataset["radar_out"].fields:
            warn(
                " Field type "
                + field_name
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        # user defined values
        lon = prdcfg.get("lon", dataset["radar_out"].origin_longitude["data"][0])
        lat = prdcfg.get("lat", dataset["radar_out"].origin_latitude["data"][0])
        alpha = prdcfg.get("alpha", 1)
        wind_vectors = prdcfg.get("wind_vectors", "hor")

        # Find indexes where to slice with a pyart display object
        # this is a bit clumsy
        display = pyart.graph.GridMapDisplay(dataset["radar_out"])
        level, _ = display._find_nearest_grid_indices(lon, lat)

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=prdcfg["timeinfo"],
        )

        fname_list = make_filename(
            "dda_lon_slice",
            prdcfg["dstype"],
            prdcfg["voltype"],
            prdcfg["imgformat"],
            prdcfginfo="l" + str(level),
            timeinfo=prdcfg["timeinfo"],
            runinfo=prdcfg["runinfo"],
        )

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir + fname

        titl = generate_dda_longitude_slice_title(
            dataset["radar_out"], bg_field_name, level
        )

        plot_dda_slice(
            dataset["radar_out"],
            bg_field_name,
            "longitude",
            level,
            prdcfg,
            fname_list,
            alpha=alpha,
            display_type=display_type,
            titl=titl,
            wind_vectors=wind_vectors,
        )

        print("----- save to " + " ".join(fname_list))

    if prdcfg["type"] == "DDA_LATITUDE_SLICE":
        display_type = prdcfg.get("display_type", "quiver")
        if display_type not in ["quiver", "barbs", "streamline"]:
            warn("Invalid display_type, must be " + '"streamline", "quiver" or "barbs"')
            return None

        bg_field_name = get_fieldname_pyart(prdcfg["voltype"])
        bg_ref_rad = prdcfg.get("bg_ref_rad", 0)
        # get bg field name given which radar serves as reference
        if bg_ref_rad != 0:
            bg_field_name += "_radar{:d}".format(bg_ref_rad)

        if bg_field_name not in dataset["radar_out"].fields:
            warn(
                " Field type "
                + field_name
                + " not available in data set. Skipping product "
                + prdcfg["type"]
            )
            return None

        # user defined values
        lon = prdcfg.get("lon", dataset["radar_out"].origin_longitude["data"][0])
        lat = prdcfg.get("lat", dataset["radar_out"].origin_latitude["data"][0])
        alpha = prdcfg.get("alpha", 1)
        wind_vectors = prdcfg.get("wind_vectors", "hor")

        # Find indexes where to slice with a pyart display object
        # this is a bit clumsy
        display = pyart.graph.GridMapDisplay(dataset["radar_out"])
        _, level = display._find_nearest_grid_indices(lon, lat)

        savedir = get_save_dir(
            prdcfg["basepath"],
            prdcfg["procname"],
            dssavedir,
            prdcfg["prdname"],
            timeinfo=prdcfg["timeinfo"],
        )

        fname_list = make_filename(
            "dda_lat_slice",
            prdcfg["dstype"],
            prdcfg["voltype"],
            prdcfg["imgformat"],
            prdcfginfo="l" + str(level),
            timeinfo=prdcfg["timeinfo"],
            runinfo=prdcfg["runinfo"],
        )

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir + fname

        titl = generate_dda_latitude_slice_title(
            dataset["radar_out"], bg_field_name, level
        )

        plot_dda_slice(
            dataset["radar_out"],
            bg_field_name,
            "latitude",
            level,
            prdcfg,
            fname_list,
            alpha=alpha,
            display_type=display_type,
            titl=titl,
            wind_vectors=wind_vectors,
        )

        print("----- save to " + " ".join(fname_list))

    else:
        return None
