"""
pyrad.graph.plots
=================

Functions to plot Pyrad datasets

.. autosummary::
    :toctree: generated/

    plot_pos
    plot_pos_map
    plot_density
    plot_scatter
    plot_centroids
    plot_quantiles
    plot_histogram
    plot_histogram2
    plot_antenna_pattern
    plot_selfconsistency
    plot_selfconsistency_instrument
    plot_selfconsistency_instrument2
    plot_scatter_comp
    plot_sun_hits
    _plot_sunscan
    _plot_time_range

"""

from ..util.radar_utils import compute_quantiles_from_hist
from mpl_toolkits.axes_grid1 import make_axes_locatable
from .plots_aux import get_colobar_label, get_field_name, get_norm
import pyart
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker
from ..util import warn

import numpy as np

try:
    import cartopy
    from cartopy.io.img_tiles import Stamen

    _CARTOPY_AVAILABLE = True
except ImportError:
    _CARTOPY_AVAILABLE = False

import matplotlib as mpl

mpl.use("Agg")

# Increase a bit font size
mpl.rcParams.update({"font.size": 16})
mpl.rcParams.update({"font.family": "sans-serif"})


def plot_pos(
    lat,
    lon,
    alt,
    fname_list,
    ax=None,
    fig=None,
    save_fig=True,
    sort_altitude="No",
    dpi=72,
    alpha=1.0,
    cb_label="height [m MSL]",
    titl="Position",
    xlabel="Lon [Deg]",
    ylabel="Lat [Deg]",
    limits=None,
    vmin=None,
    vmax=None,
):
    """
    plots a trajectory on a Cartesian surface

    Parameters
    ----------
    lat, lon, alt : float array
        Points coordinates
    fname_list : list of str
        list of names of the files where to store the plot
    fig : Figure
        Figure to add the colorbar to. If none a new figure will be created
    ax : Axis
        Axis to plot on. if fig is None a new axis will be created
    save_fig : bool
        if true save the figure if false it does not close the plot and
        returns the handle to the figure
    sort_altitude : str
        String indicating whether to sort the altitude data. Can be 'No',
        'Lowest_on_top' or 'Highest_on_top'
    dpi : int
        Pixel density
    alpha : float
        Transparency
    cb_label : str
        Color bar label
    titl : str
        Plot title
    xlabel, ylabel : str
        The labels of the X and Y axis
    limits : tupple or None
        The limits of the field to plot
    vmin, vmax : float
        The limits of the color scale

    Returns
    -------
    fname_list : list of str or
    fig, ax : tupple
        list of names of the saved plots or handle of the figure an axes
    """
    if sort_altitude in ("Lowest_on_top", "Highest_on_top"):
        ind = np.argsort(alt)
        if sort_altitude == "Lowest_on_top":
            ind = ind[::-1]
        lat = lat[ind]
        lon = lon[ind]
        alt = alt[ind]

    if vmin is None:
        vmin = alt.min()
    if vmax is None:
        vmax = alt.max()
    marker = "x"
    col = alt
    cmap = "viridis"
    norm = plt.Normalize(vmin, vmax)

    if fig is None:
        fig = plt.figure(figsize=[10, 8], dpi=dpi)
        ax = fig.add_subplot(111, aspect="equal")
    else:
        ax.autoscale(False)

    cax = ax.scatter(lon, lat, c=col, marker=marker, alpha=alpha, cmap=cmap, norm=norm)

    if limits is not None:
        ax.set_xlim(limits[0], limits[1])
        ax.set_ylim(limits[2], limits[3])

    # plot colorbar
    cb = fig.colorbar(cax, orientation="horizontal")
    cb.set_label(cb_label)

    ax.set_title(titl)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Turn on the grid
    ax.grid()

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

        return fname_list

    return (fig, ax)


def plot_pos_map(
    lat,
    lon,
    alt,
    fname_list,
    ax=None,
    fig=None,
    save_fig=True,
    sort_altitude="No",
    dpi=72,
    alpha=1.0,
    cb_label="height [m MSL]",
    titl="Position",
    xlabel="Lon [Deg]",
    ylabel="Lat [Deg]",
    limits=None,
    vmin=None,
    vmax=None,
    lon_step=0.3,
    lat_step=0.1,
    background_zoom=8,
):
    """
    plots a trajectory on a map

    Parameters
    ----------
    lat, lon, alt : float array
        Points coordinates
    fname_list : list of str
        list of names of the files where to store the plot
    fig : Figure
        Figure to add the colorbar to. If none a new figure will be created
    ax : Axis
        Axis to plot on. if fig is None a new axis will be created
    save_fig : bool
        if true save the figure if false it does not close the plot and
        returns the handle to the figure
    sort_altitude : str
        String indicating whether to sort the altitude data. Can be 'No',
        'Lowest_on_top' or 'Highest_on_top'
    dpi : int
        Pixel density
    alpha : float
        Transparency
    cb_label : str
        Color bar label
    titl : str
        Plot title
    xlabel, ylabel : str
        The labels of the X and Y axis
    limits : tupple or None
        The limits of the field to plot
    vmin, vmax : float
        The limits of the color scale
    lon_step, lat_step : float
        The step interval of the latitude, longitude lines to plot
    background_zoom : int
        The zoom of the background image. A higher number will give more level
        of detail at the expense of speed.

    Returns
    -------
    fname_list : list of str or
    fig, ax : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    if not _CARTOPY_AVAILABLE:
        warn("Unable to plot trajectory on a map. Cartopy not available")
        return [""]

    if sort_altitude in ("Lowest_on_top", "Highest_on_top"):
        ind = np.argsort(alt)
        if sort_altitude == "Lowest_on_top":
            ind = ind[::-1]
        lat = lat[ind]
        lon = lon[ind]
        alt = alt[ind]

    if vmin is None:
        vmin = alt.min()
    if vmax is None:
        vmax = alt.max()
    marker = "x"
    col = alt
    cmap = "viridis"
    norm = plt.Normalize(vmin, vmax)

    # get background map instance
    stamen_terrain = Stamen("terrain-background")
    projection = cartopy.crs.PlateCarree()

    # set map limits and lat/lon lines
    if limits is None:
        limits = [np.min(lon), np.max(lon), np.min(lat), np.max(lat)]

    lon_lines = np.arange(limits[0], limits[1] + 1, lon_step)
    lat_lines = np.arange(limits[2], limits[3] + 1, lat_step)

    if fig is None:
        fig = plt.figure(figsize=[10, 8], dpi=dpi)

        # draw background
        ax = fig.add_subplot(111, projection=stamen_terrain.crs)
        ax.set_extent(limits, crs=projection)
        ax.add_image(stamen_terrain, background_zoom)

        # add countries
        countries = cartopy.feature.NaturalEarthFeature(
            category="cultural", name="admin_0_countries", scale="10m", facecolor="none"
        )
        ax.add_feature(countries, edgecolor="black")

        # draw grid lines and labels
        gl = ax.gridlines(xlocs=lon_lines, ylocs=lat_lines, draw_labels=True)
        gl.xlabels_top = False
        gl.ylabels_right = False

        ax.text(
            0.5,
            -0.2,
            xlabel,
            va="bottom",
            ha="center",
            rotation="horizontal",
            rotation_mode="anchor",
            transform=ax.transAxes,
        )

        ax.text(
            -0.1,
            0.55,
            ylabel,
            va="bottom",
            ha="center",
            rotation="vertical",
            rotation_mode="anchor",
            transform=ax.transAxes,
        )
    else:
        ax.autoscale(False)

    # plot data
    cax = ax.scatter(
        lon,
        lat,
        c=col,
        marker=marker,
        alpha=alpha,
        cmap=cmap,
        norm=norm,
        transform=projection,
    )

    # plot colorbar
    cb = fig.colorbar(cax, orientation="horizontal")
    cb.set_label(cb_label)

    ax.set_title(titl)

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

        return fname_list

    return (fig, ax)


def plot_density(
    hist_obj,
    hist_type,
    field_name,
    ind_sweep,
    prdcfg,
    fname_list,
    quantiles=(25.0, 50.0, 75.0),
    ref_value=0.0,
    vmin=None,
    vmax=None,
):
    """
    density plot (angle-values representation)

    Parameters
    ----------
    hist_obj : histogram object
        object containing the histogram data to plot
    hist_type : str
        type of histogram (instantaneous data or cumulative)
    field_name : str
        name of the radar field to plot
    ind_sweep : int
        sweep index to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    quantiles : array
        the quantile lines to plot
    ref_value : float
        the reference value
    vmin, vmax : float
        Minimum and maximum extent of the vertical axis

    Returns
    -------
    fname_list : list of str
        list of names of the created plots
    """

    hist_obj_aux = hist_obj.extract_sweeps([ind_sweep])

    # --------------------------------------------------
    # Get angular coordinate and sort field accordingly
    # --------------------------------------------------
    if hist_obj_aux.scan_type == "ppi":
        ang_raw = hist_obj_aux.azimuth["data"]
        ind_ang = np.argsort(ang_raw)
        ang = ang_raw[ind_ang]
        field = hist_obj_aux.fields[field_name]["data"][ind_ang, :]
        labelx = "azimuth angle (degrees)"
        force_full_angle_axis = True

    elif hist_obj_aux.scan_type == "rhi":
        ang_raw = hist_obj_aux.elevation["data"]
        ind_ang = np.argsort(ang_raw)
        ang = ang_raw[ind_ang]
        field = hist_obj_aux.fields[field_name]["data"][ind_ang, :]
        labelx = "elevation angle (degrees)"
        force_full_angle_axis = False

    else:
        field = hist_obj_aux.fields[field_name]["data"]
        ang = np.arange(hist_obj_aux.nrays)
        labelx = "ray number"
        force_full_angle_axis = False

    quantiles = np.asarray(quantiles)

    # --------------------------------------------------
    # Compute quantiles for each ray
    # --------------------------------------------------
    az_percentiles = np.ma.masked_all((len(ang), len(quantiles)))

    for ray in range(len(ang)):
        # If a ray has no samples, keep quantiles masked
        if np.ma.sum(field[ray, :]) <= 0:
            continue

        _, values_ray = compute_quantiles_from_hist(
            hist_obj.range["data"],
            field[ray, :],
            quantiles=quantiles,
        )
        az_percentiles[ray, :] = values_ray

    # --------------------------------------------------
    # Compute overall sweep quantiles
    # --------------------------------------------------
    _, values_sweep = compute_quantiles_from_hist(
        hist_obj.range["data"],
        np.ma.sum(field, axis=0),
        quantiles=quantiles,
    )

    # --------------------------------------------------
    # Mask empty histogram bins
    # --------------------------------------------------
    field = np.ma.masked_where(field == 0, field)

    # --------------------------------------------------
    # Title and labels
    # --------------------------------------------------
    if hist_type == "instant":
        titl = pyart.graph.common.generate_title(
            hist_obj_aux,
            field_name,
            ind_sweep,
        )
    else:
        titl = (
            "{:.1f}".format(hist_obj_aux.fixed_angle["data"][0])
            + " Deg. "
            + pyart.graph.common.generate_radar_time_begin(hist_obj).strftime(
                "%Y-%m-%d"
            )
            + "\n"
            + get_field_name(hist_obj.fields[field_name], field_name)
        )

    label = "Number of Points"
    labely = get_colobar_label(hist_obj_aux.fields[field_name], field_name)

    dpi = 72
    if "dpi" in prdcfg["ppiImageConfig"]:
        dpi = prdcfg["ppiImageConfig"]["dpi"]

    fig = plt.figure(
        figsize=[
            prdcfg["ppiImageConfig"]["xsize"],
            prdcfg["ppiImageConfig"]["ysize"],
        ],
        dpi=dpi,
    )
    ax = fig.add_subplot(111)

    # --------------------------------------------------
    # Build value-bin edges
    # --------------------------------------------------
    cmap = pyart.config.get_field_colormap(field_name)

    step = hist_obj.range["data"][1] - hist_obj.range["data"][0]
    val_edges = np.append(
        hist_obj.range["data"] - step / 2.0,
        hist_obj.range["data"][-1] + step / 2.0,
    )

    # --------------------------------------------------
    # Detect angular gaps
    # --------------------------------------------------
    if len(ang) > 1:
        dang = np.diff(ang)
        positive_dang = dang[dang > 0]

        if len(positive_dang) > 0:
            typical_step = np.ma.median(positive_dang)
        else:
            typical_step = 1.0

        gap_ind = np.where(dang > 3.0 * typical_step)[0]
    else:
        gap_ind = np.array([], dtype=int)

    segments = []
    start = 0

    for idx in gap_ind:
        segments.append((start, idx + 1))
        start = idx + 1

    segments.append((start, len(ang)))

    # --------------------------------------------------
    # Plot density using pcolormesh
    # --------------------------------------------------
    cax = None
    vmax_field = np.ma.max(field)

    for start, end in segments:
        ang_seg = ang[start:end]
        field_seg = field[start:end, :]

        if len(ang_seg) < 2:
            continue

        ang_step_seg = np.ma.median(np.diff(ang_seg))

        ang_edges = np.empty(len(ang_seg) + 1)
        ang_edges[1:-1] = 0.5 * (ang_seg[:-1] + ang_seg[1:])
        ang_edges[0] = ang_seg[0] - 0.5 * ang_step_seg
        ang_edges[-1] = ang_seg[-1] + 0.5 * ang_step_seg

        cax = ax.pcolormesh(
            ang_edges,
            val_edges,
            field_seg.T,
            cmap=cmap,
            vmin=0.0,
            vmax=vmax_field,
            shading="auto",
        )

    # --------------------------------------------------
    # Axis limits
    # --------------------------------------------------
    if force_full_angle_axis:
        ax.set_xlim(0.0, 360.0)
    else:
        if len(ang) > 1:
            ang_step = np.ma.median(np.diff(ang))
        else:
            ang_step = 1.0

        ax.set_xlim(
            np.ma.min(ang) - ang_step / 2.0,
            np.ma.max(ang) + ang_step / 2.0,
        )

    ax.set_ylim(bottom=vmin, top=vmax)
    ax.autoscale(False)

    # --------------------------------------------------
    # Plot reference line, split at angular gaps
    # --------------------------------------------------
    for start, end in segments:
        if end <= start:
            continue

        ax.plot(
            ang[start:end],
            np.zeros(end - start) + ref_value,
            "k--",
        )

    # --------------------------------------------------
    # Plot quantile lines, split at angular gaps
    # --------------------------------------------------
    quantile_colors = plt.cm.viridis_r(np.linspace(0, 1, len(quantiles)))

    for i, q in enumerate(quantiles):
        first_label = True

        for start, end in segments:
            if end <= start:
                continue

            label_sweep = f"{q:.1f}th sweep quantile" if first_label else None
            label_ray = f"{q:.1f}th ray quantile" if first_label else None

            ax.plot(
                ang[start:end],
                np.zeros(end - start) + values_sweep[i],
                color=quantile_colors[i],
                linestyle="--",
                label=label_sweep,
            )

            ax.plot(
                ang[start:end],
                az_percentiles[start:end, i],
                color=quantile_colors[i],
                linestyle="-",
                label=label_ray,
            )

            first_label = False

    # --------------------------------------------------
    # Labels, title, colorbar
    # --------------------------------------------------
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_title(titl)

    if cax is not None:
        cb = fig.colorbar(cax)
        cb.set_label(label)

    # --------------------------------------------------
    # Metadata box
    # --------------------------------------------------
    metadata = "npoints: " + str(np.ma.sum(field)) + "\n"

    for i, quant in enumerate(quantiles):
        val_quant_str = "--"
        if values_sweep[i] is not np.ma.masked:
            val_quant_str = f"{values_sweep[i]:.3f}"

        metadata += f"{quant} quant: {val_quant_str}\n"

    ax.text(
        0.05,
        0.05,
        metadata,
        horizontalalignment="left",
        verticalalignment="bottom",
        transform=ax.transAxes,
    )

    # --------------------------------------------------
    # Grid and save
    # --------------------------------------------------
    ax.grid()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")

    plt.close(fig)

    return fname_list


def plot_scatter(
    bin_edges1,
    bin_edges2,
    hist_2d,
    field_name1,
    field_name2,
    fname_list,
    prdcfg,
    metadata=None,
    lin_regr=None,
    lin_regr_slope1=None,
    rad1_name="RADAR001",
    rad2_name="RADAR002",
    titl="colocated radar gates",
    cmap=None,
    vmin=None,
    vmax=None,
    min_cnt=0,
    bins_transform="linear",
    marginals=True,
):
    """
    2D histogram (scatter density) plot with optional marginal distributions.

    Parameters
    ----------
    bin_edges1, bin_edges2 : array-like
        Bin edges of each field (x and y dimensions).
    hist_2d : ndarray (2D)
        The 2D histogram (counts per bin).
    field_name1, field_name2 : str
        Names of each field (used for axis labeling).
    fname_list : list of str
        List of filenames where the plot will be saved.
    prdcfg : dict
        Product configuration dictionary. Must contain the key
        "ppiImageConfig" with entries "xsize", "ysize", and optionally "dpi".
    metadata : str, optional
        String with metadata to display in the plot (top-left corner).
    lin_regr : tuple of float, optional
        Coefficients (slope, intercept) of a linear regression line to plot.
    lin_regr_slope1 : float, optional
        Intercept of a regression line with slope fixed to 1.
    rad1_name, rad2_name : str, optional
        Names of the radars providing the data.
    titl : str, optional
        Plot title.
    cmap : str or None, optional
        Name of the colormap. If None, the default Py-ART colormap for
        `field_name1` will be used.
    vmin, vmax : float, optional
        Minimum and maximum values for the axes limits.
    min_cnt : int, optional
        Minimum number of counts displayed in the histogram.
    bins_transform : {"linear", "log"}, optional
        Scaling of the histogram colorbar. If "log", a logarithmic scale is used.
        Default is "linear".
    marginals : bool, optional
        If True, marginal distributions (1D histograms) are plotted on the top
        and right side of the 2D histogram. Default is False.

    Returns
    -------
    fname_list : list of str
        List of filenames of the created plots.
    """
    # mask 0 data
    hist_2d = np.ma.masked_where(hist_2d == 0, hist_2d)

    label = "Number of Points"
    labelx = rad1_name + " " + field_name1
    labely = rad2_name + " " + field_name2

    dpi = 72
    if "dpi" in prdcfg["ppiImageConfig"]:
        dpi = prdcfg["ppiImageConfig"]["dpi"]

    fig = plt.figure(
        figsize=[
            prdcfg["ppiImageConfig"]["xsize"],
            prdcfg["ppiImageConfig"]["ysize"],
        ],
        dpi=dpi,
    )

    ax = fig.add_subplot(111)
    ax.set_aspect("equal", adjustable="box")
    ax.grid()
    if marginals:
        divider = make_axes_locatable(ax)

        ax_marg_x = divider.append_axes(
            "top",
            size="22%",
            pad=0.0,
            sharex=ax,
        )

        ax_marg_y = divider.append_axes(
            "right",
            size="22%",
            pad=0.0,
            sharey=ax,
        )

        cax_cb = divider.append_axes(
            "right",
            size="4%",
            pad=0.12,
        )

        ax_marg_x.tick_params(axis="x", labelbottom=False)
        ax_marg_y.tick_params(axis="y", labelleft=False)

    else:
        ax_marg_x = None
        ax_marg_y = None
        cax_cb = None

    if cmap is None:
        cmap = pyart.config.get_field_colormap(field_name1)

    X, Y = np.meshgrid(bin_edges1, bin_edges2)

    if bins_transform == "log":
        cax = ax.pcolormesh(
            X,
            Y,
            np.ma.transpose(hist_2d),
            cmap=cmap,
            norm=LogNorm(vmin=max(min_cnt, 1), vmax=np.max(hist_2d)),
        )
    else:
        cax = ax.pcolormesh(
            X,
            Y,
            np.ma.transpose(hist_2d),
            cmap=cmap,
            vmin=min_cnt,
            vmax=np.max(hist_2d),
        )

    step1 = bin_edges1[1] - bin_edges1[0]
    bin_centers1 = bin_edges1[:-1] + step1 / 2.0

    step2 = bin_edges2[1] - bin_edges2[0]
    bin_centers2 = bin_edges2[:-1] + step2 / 2.0

    # Marginal distributions
    if marginals:
        hist_for_marginals = np.ma.filled(hist_2d, 0)

        marginal_x = hist_for_marginals.sum(axis=1)
        marginal_y = hist_for_marginals.sum(axis=0)

        if marginal_x.max() > 0:
            marginal_x = marginal_x / marginal_x.max()
        if marginal_y.max() > 0:
            marginal_y = marginal_y / marginal_y.max()

        ax_marg_x.plot(
            bin_centers1,
            marginal_x,
            color="k",
            alpha=0.6,
            lw=1.2,
        )

        ax_marg_y.plot(
            marginal_y,
            bin_centers2,
            color="k",
            alpha=0.6,
            lw=1.2,
        )

        ax_marg_x.set_ylim(0, 1.05)
        ax_marg_y.set_xlim(0, 1.05)

        # Jointplot-like clean marginals
        for spine in ax_marg_x.spines.values():
            spine.set_visible(False)
        for spine in ax_marg_y.spines.values():
            spine.set_visible(False)

        ax_marg_x.tick_params(
            axis="both",
            which="both",
            bottom=False,
            top=False,
            left=False,
            right=False,
            labelbottom=False,
            labelleft=False,
        )

        ax_marg_y.tick_params(
            axis="both",
            which="both",
            bottom=False,
            top=False,
            left=False,
            right=False,
            labelbottom=False,
            labelleft=False,
        )

    # plot reference
    if bin_edges1.size == bin_edges2.size:
        ax.plot(bin_centers1, bin_centers2, "k--")

    # plot linear regression
    if lin_regr is not None:
        ax.plot(bin_centers1, lin_regr[0] * bin_centers1 + lin_regr[1], "r")

    if lin_regr_slope1 is not None:
        ax.plot(bin_centers1, bin_centers1 + lin_regr_slope1, "g")

    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    if marginals:
        ax_marg_x.set_title(titl, pad=10)
    else:
        ax.set_title(titl)

    if vmin is not None and vmax is not None:
        ax.set_xlim([vmin, vmax])
        ax.set_ylim([vmin, vmax])

    if marginals:
        cb = fig.colorbar(cax, cax=cax_cb)
    else:
        cb = fig.colorbar(cax, ax=ax)
    cb.set_label(label)

    if bins_transform == "log":
        cb.locator = ticker.LogLocator(base=10)
        formatter = ticker.FuncFormatter(lambda x, _: f"{x:g}")
        cb.formatter = formatter
        cb.ax.yaxis.set_minor_formatter(formatter)
        cb.minorticks_on()
        cb.update_ticks()

        cb.ax.tick_params(which="major", length=6, width=1.0)
        cb.ax.tick_params(which="minor", length=3, width=0.6, color="0.5")

    if metadata is not None:
        ax.text(
            0.05,
            0.95,
            metadata,
            horizontalalignment="left",
            verticalalignment="top",
            transform=ax.transAxes,
            color="#808080",
            alpha=0.8,
            fontsize=9,
        )

    if marginals:
        fig.canvas.draw()  # ensure positions are computed
        pos_main = ax.get_position()

        # Define a thinner width and attach it directly to the main axis
        new_width = pos_main.width * 0.22  # adjust (0.15–0.25 works well)

        ax_marg_y.set_position(
            [
                pos_main.x1,  # start exactly at right edge of main plot
                pos_main.y0,  # same vertical start
                new_width,  # thinner width
                pos_main.height,  # same height
            ]
        )
    else:
        fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")

    plt.close(fig)

    return fname_list


def plot_centroids(
    bin_edges1,
    bin_edges2,
    hist_2d,
    field_name1,
    field_name2,
    fname_list,
    prdcfg,
    titl="centroids",
    cmap=None,
    medoids_x=None,
    medoids_y=None,
    fmedoid_x=None,
    fmedoid_y=None,
):
    """
    2D histogram

    Parameters
    ----------
    bin_edges1, bin_edges2 : float array2
        the bins of each field
    hist_2d : ndarray 2D
        the 2D histogram
    field_name1, field_name2 : str
        the names of each field
    fname_list : list of str
        list of names of the files where to store the plot
    prdcfg : dict
        product configuration dictionary
    titl : str
        plot title
    cmap : str or None
        name of the colormap. If None it will be choosen the default for the
        field_name
    medoids_x, medoids_y : ndarray 1D or None
        intermediate medoids
    fmedoid_x, fmedoid_y : ndarray 1D or None
        final medoid

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    # mask 0 data
    hist_2d = np.ma.masked_where(hist_2d == 0, hist_2d)

    # display data
    label = "Number of Points"
    labelx = field_name1
    labely = field_name2

    dpi = 72
    if "dpi" in prdcfg["ppiImageConfig"]:
        dpi = prdcfg["ppiImageConfig"]["dpi"]

    fig = plt.figure(
        figsize=[prdcfg["ppiImageConfig"]["xsize"], prdcfg["ppiImageConfig"]["ysize"]],
        dpi=dpi,
    )
    ax = fig.add_subplot(111)

    if cmap is None:
        cmap = pyart.config.get_field_colormap(field_name1)

    X, Y = np.meshgrid(bin_edges1, bin_edges2)
    cax = ax.pcolormesh(
        X, Y, np.ma.transpose(hist_2d), cmap=cmap, vmin=0.0, vmax=np.max(hist_2d)
    )

    # plot intermediate medoids
    if medoids_x is not None:
        ax.scatter(medoids_x, medoids_y, c="w", marker="o")
    if fmedoid_x is not None:
        ax.plot(fmedoid_x, fmedoid_y, c="r", marker="D")

    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_title(titl)

    cb = fig.colorbar(cax)
    cb.set_label(label)

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def plot_quantiles(
    quant,
    value,
    fname_list,
    labelx="quantile",
    labely="value",
    titl="quantile",
    vmin=None,
    vmax=None,
    dpi=72,
):
    """
    plots quantiles

    Parameters
    ----------
    quant : array
        quantiles to be plotted
    value : array
        values of each quantile
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    titl : str
        The figure title
    vmin, vmax: float
        Lower/Upper limit of data values
    dpi : int
        dots per inch

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    fig, ax = plt.subplots(figsize=[10, 6], dpi=dpi)
    ax.plot(quant, value, "bx-")
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_ylim(bottom=vmin, top=vmax)
    ax.set_title(titl)

    # Turn on the grid
    ax.grid()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def plot_histogram(
    bin_edges,
    values,
    fname_list,
    labelx="bins",
    labely="Number of Samples",
    titl="histogram",
    binwidth_equal=False,
    dpi=72,
    labels=None,
    alpha=None,
    vert_line=None,
):
    """
    Compute and plot one or several histograms.

    Parameters
    ----------
    bin_edges : array-like or list of array-like
        Histogram bin edges. If a single array is provided, it is used for all
        value arrays. If a list is provided, each entry is used for the
        corresponding value array.
    values : array-like or list of array-like
        Data values. Can be a single array or a list of arrays to plot several
        distributions on the same axis.
    fname_list : list of str
        List of filenames where the plot will be saved.
    labelx : str, optional
        Label of the x-axis.
    labely : str, optional
        Label of the y-axis.
    titl : str, optional
        Figure title.
    binwidth_equal : bool, optional
        If True, bars are plotted with equal visual width regardless of the
        actual bin size.
    dpi : int, optional
        Dots per inch.
    labels : list of str, optional
        Labels for each distribution. If provided, a legend is added.
    alpha : float, optional
        Transparency of the histogram (only used for single distribution).
    vert_line : float or list of float, optional
        If provided, vertical line(s) are plotted at the specified x-value(s).
    Returns
    -------
    fname_list : list of str
        List of filenames of the created plots.
    """

    # Normalize inputs
    if not isinstance(values, (list, tuple)):
        values_list = [values]
    else:
        values_list = values

    n_hist = len(values_list)

    if isinstance(bin_edges, (list, tuple)) and len(bin_edges) == n_hist:
        bin_edges_list = bin_edges
    else:
        bin_edges_list = [bin_edges] * n_hist

    if labels is None:
        labels = [None] * n_hist

    # Alpha only relevant for single histogram
    if n_hist == 1:
        if alpha is None:
            alpha = 0.7
    else:
        alpha = None  # not used

    fig, ax = plt.subplots(figsize=[10, 6], dpi=dpi)

    for vals, bins, lab in zip(values_list, bin_edges_list, labels):
        hist, bins = np.histogram(vals, bins=bins)

        if n_hist == 1 and len(bins) < 100:
            # --- Standard filled histogram ---
            if binwidth_equal:
                x = np.arange(len(bins) - 1) + 0.5
                ax.bar(
                    x,
                    hist,
                    width=1,
                    alpha=alpha,
                    label=lab,
                    edgecolor="black",
                    linewidth=0.5,
                )
                ax.set_xticks(np.arange(len(bins)))
                ax.set_xticklabels(bins, rotation=45.0)
            else:
                ax.hist(
                    vals,
                    bins=bins,
                    alpha=alpha,
                    label=lab,
                    edgecolor="black",
                    linewidth=0.5,
                )

        else:
            # --- Line histogram (no fill) ---
            bin_centers = 0.5 * (bins[:-1] + bins[1:])
            ax.plot(
                bin_centers,
                hist,
                label=lab,
                linewidth=1.8,
            )

    ax.grid()
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_title(titl)

    if vert_line is not None:
        if isinstance(vert_line, (list, tuple)):
            for vline in vert_line:
                ax.axvline(vline, color="r", linestyle="--")
        else:
            ax.axvline(vert_line, color="r", linestyle="--")

    if any(label is not None for label in labels):
        ax.legend()

    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")

    plt.close(fig)

    return fname_list


def plot_histogram2(
    bin_centers,
    hist,
    fname_list,
    width=None,
    labelx="bins",
    labely="Number of Samples",
    titl="histogram",
    dpi=72,
    ax=None,
    fig=None,
    save_fig=True,
    color=None,
    alpha=None,
    invert_xaxis=False,
):
    """
    plots histogram

    Parameters
    ----------
    bin_centers : array
        histogram bin centers
    hist : array
        values for each bin
    fname_list : list of str
        list of names of the files where to store the plot
    width : scalar or array-like
        the width(s) of the bars. If None it is going to be estimated from the
        distances between centers
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    titl : str
        The figure title
    dpi : int
        dots per inch
    fig : Figure
        Figure to add the colorbar to. If none a new figure will be created
    ax : Axis
        Axis to plot on. if fig is None a new axis will be created
    save_fig : bool
        if true save the figure. If false it does not close the plot and
        returns the handle to the figure
    color : str
        color of the bars
    alpha : float
        parameter controling the transparency
    invert_xaxis : bool
        If true inverts the x axis

    Returns
    -------
    fname_list or fig, ax: list of str
        list of names of the created plots

    """
    if fig is None:
        fig = plt.figure(figsize=[10, 6], dpi=dpi)
        ax = fig.add_subplot(111)
    else:
        ax.autoscale(False)

    if width is None:
        width = bin_centers[1] - bin_centers[0]
    ax.bar(bin_centers, hist, width=width, color=color, alpha=alpha)
    if invert_xaxis:
        ax.invert_xaxis()

    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_title(titl)
    ax.grid()

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

        return fname_list

    return (fig, ax)


def plot_antenna_pattern(
    antpattern,
    fname_list,
    labelx="Angle [Deg]",
    linear=False,
    twoway=False,
    title="Antenna Pattern",
    ymin=None,
    ymax=None,
    dpi=72,
):
    """
    plots an antenna pattern

    Parameters
    ----------
    antpattern : dict
        dictionary with the angle and the attenuation
    value : float array
        values of the time series
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        The label of the X axis
    linear : boolean
        if true data is in linear units
    linear : boolean
        if true data represents the two way attenuation
    titl : str
        The figure title
    ymin, ymax: float
        Lower/Upper limit of y axis
    dpi : int
        dots per inch

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    waystr = "One-way "
    bw_threshold = -3.0
    if twoway:
        waystr = "Two-way "
        bw_threshold = -6

    linstr = "Att [dB]"
    if linear:
        linstr = "Att [-]"
        mainbeam = antpattern["angle"][
            antpattern["attenuation"] >= 10.0 ** (0.1 * bw_threshold)
        ]
    else:
        mainbeam = antpattern["angle"][antpattern["attenuation"] >= bw_threshold]
    beamwidth = np.abs(np.max(mainbeam) - np.min(mainbeam))

    labely = waystr + linstr

    fig, ax = plt.subplots(figsize=[10, 6], dpi=dpi)

    ax.plot(antpattern["angle"], antpattern["attenuation"])
    ax.set_title(title)
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_ylim(bottom=ymin, top=ymax)

    metadata = "3-dB beamwidth: " + "{:.2f}".format(float(beamwidth))
    ax.text(
        0.05,
        0.95,
        metadata,
        horizontalalignment="left",
        verticalalignment="top",
        transform=ax.transAxes,
    )

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def plot_selfconsistency(
    zdrkdp_table,
    fname_list,
    labelx="ZDR [dB]",
    labely="KDP/Zh [(deg*m3)/(km*mm6)]",
    title="Selfconsistency in rain",
    ymin=None,
    ymax=None,
    dpi=72,
    save_fig=True,
    ax=None,
    fig=None,
):
    """
    plots a ZDR-KDP/ZH selfconsistency in rain relation

    Parameters
    ----------
    antpattern : dict
        dictionary with the angle and the attenuation
    value : float array
        values of the time series
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        The label of the X axis
    linear : boolean
        if true data is in linear units
    linear : boolean
        if true data represents the two way attenuation
    titl : str
        The figure title
    ymin, ymax: float
        Lower/Upper limit of y axis
    dpi : int
        dots per inch

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    if fig is None:
        fig, ax = plt.subplots(figsize=[10, 6], dpi=dpi)
    else:
        ax.autoscale(False)

    ax.plot(zdrkdp_table[0], zdrkdp_table[1])
    ax.set_title(title)
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_ylim(bottom=ymin, top=ymax)

    # Make a tight layout
    fig.tight_layout()

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

        return fname_list

    return (fig, ax)


def plot_selfconsistency_instrument(
    zdr,
    kdp,
    zh,
    fname_list,
    bins_zdr_step=0.05,
    bins_zdr_min=0.0,
    bins_zdr_max=6.0,
    bins_kdpzh_step=0.1,
    bins_kdpzh_min=-2.0,
    bins_kdpzh_max=20.0,
    normalize=True,
    vmin=0.0,
    vmax=0.01,
    parametrization="None",
    zdr_kdpzh_dict=None,
    retrieve_relation=True,
    plot_theoretical=True,
    dpi=72,
):
    """
    plots the ZDR-KDP/ZH relationship obtained by an instrument. The
    theoretical curve and the retrieved curve

    Parameters
    ----------
    zdr, kdp, zh : 1D ndarray
        The valid values of ZDR [dB], KDP [deg/km] and Zh [mm6/m3] collected
        by the instrument
    fname_list : list of str
        list of names of the files where to store the plot
    bins_zdr_step : float
        The step of the ZDR axis of the histogram [dB]
    bins_zdr_min, bins_zdr_max : float
        The limits of the ZDR axis of the histogram (bins center) [dB]
    bins_kdpzh_step : float
        The step of the 1e5*KDP^a/ZH^b axis of the histogram
        [(deg*m3)/(km*mm6)]
    bins_kdpzh_min, bins_kdpzh_max : float
        The limits of the 1e5*KDP^a/ZH^b axis of the histogram (bins center)
        [(deg*m3)/(km*mm6)]
    normalize : Bool
        If True the occurrence density of ZH/KDP for each ZDR bin is going to
        be represented. Otherwise it will show the number of gates at each bin
    vmin, vmax : float
        min and max values of the colorbar
    parametrization : str
        The type of parametrization for the self-consistency curves. Can be
        'None', 'Gourley', 'Wolfensberger', 'Louf', 'Gorgucci' or 'Vaccarono'.
        'None' will use tables contained in zdr_kdpzh_dict. The parametrized
        curves are obtained from literature except for Wolfensberger that was
        derived from disdrometer data obtained by MeteoSwiss and EPFL. All
        parametrizations are valid for C-band only except that of Gourley.
    zdr_kdpzh_dict : dict
        dictionary containing a look up table relating ZDR with KDP/Zh for
        different elevations and the frequency band of the radar
    retrieve_relation : boolean
        if true a zdr-kdp/zh relationship is retrieved from the data
    plot_theoretical : bool
        if true the theoretical relationship is retrieved
    dpi : int
        dots per inch

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """

    # prepare bins for histogram
    bins_zdr_centers = np.arange(
        bins_zdr_min, bins_zdr_max + bins_zdr_step, bins_zdr_step
    )
    bins_zdr_edges = np.append(
        bins_zdr_centers - bins_zdr_step / 2.0,
        bins_zdr_centers[-1] + bins_zdr_step / 2.0,
    )
    bins_kdpzh_centers = np.arange(
        bins_kdpzh_min, bins_kdpzh_max + bins_kdpzh_step, bins_kdpzh_step
    )
    bins_kdpzh_edges = np.append(
        bins_kdpzh_centers - bins_kdpzh_step / 2.0,
        bins_kdpzh_centers[-1] + bins_kdpzh_step / 2.0,
    )

    labely = "10e5*KDP/ZH [(deg*m3)/(km*mm6)"

    kdpzh_th = None
    if plot_theoretical:
        if parametrization == "None":
            if zdr_kdpzh_dict is not None and "zdr_kdpzh" in zdr_kdpzh_dict:
                zdr_th = zdr_kdpzh_dict["zdr_kdpzh"][0][0]
                kdpzh_th = 1e5 * zdr_kdpzh_dict["zdr_kdpzh"][0][1]
            else:
                warn(
                    "Unable to plot theoretical self-consistency curve. "
                    "Relationship not provided"
                )
        elif parametrization == "Gourley":
            if zdr_kdpzh_dict is None or "freq_band" not in zdr_kdpzh_dict:
                warn(
                    "Unable to plot theoretical self-consistency curve. "
                    "Frequency band not provided"
                )
            else:
                zdr_th = np.arange(0.0, 3.5 + bins_zdr_step, bins_zdr_step)
                if zdr_kdpzh_dict["freq_band"] == "S":
                    kdpzh_th = (
                        3.696
                        - 1.963 * zdr_th
                        + 0.504 * zdr_th * zdr_th
                        - 0.051 * zdr_th * zdr_th * zdr_th
                    )
                elif zdr_kdpzh_dict["freq_band"] == "C":
                    kdpzh_th = (
                        6.746
                        - 2.970 * zdr_th
                        + 0.711 * zdr_th * zdr_th
                        - 0.079 * zdr_th * zdr_th * zdr_th
                    )
                elif zdr_kdpzh_dict["freq_band"] == "X":
                    kdpzh_th = (
                        11.74
                        - 4.020 * zdr_th
                        - 0.140 * zdr_th * zdr_th
                        + 0.130 * zdr_th * zdr_th * zdr_th
                    )
                else:
                    warn(
                        "Unable to plot theoretical self-consistency curve. "
                        "Unknown frequency band " + zdr_kdpzh_dict["freq_band"]
                    )
        elif parametrization == "Wolfensberger":
            if zdr_kdpzh_dict is None or "freq_band" not in zdr_kdpzh_dict:
                warn(
                    "Unable to plot theoretical self-consistency curve. "
                    "Frequency band not provided"
                )
            else:
                zdr_th = np.arange(0.0, 3.5 + bins_zdr_step, bins_zdr_step)
                if zdr_kdpzh_dict["freq_band"] == "C":
                    kdpzh_th = (
                        3.199 * np.ma.exp(-7.767e-1 * zdr_th) - 0.4436 * zdr_th + 3.464
                    )
                else:
                    warn(
                        "Unable to plot theoretical self-consistency curve. "
                        "Unknown frequency band " + zdr_kdpzh_dict["freq_band"]
                    )
        elif parametrization == "Louf":
            if zdr_kdpzh_dict is None or "freq_band" not in zdr_kdpzh_dict:
                warn(
                    "Unable to plot theoretical self-consistency curve. "
                    "Frequency band not provided"
                )
            else:
                zdr_th = np.arange(0.5, 3.5 + bins_zdr_step, bins_zdr_step)
                if zdr_kdpzh_dict["freq_band"] == "C":
                    kdpzh_th = (
                        6.607
                        - 4.577 * zdr_th
                        + 1.577 * zdr_th * zdr_th
                        - 0.23 * zdr_th * zdr_th * zdr_th
                    )
                else:
                    warn(
                        "Unable to plot theoretical self-consistency curve. "
                        "Unknown frequency band " + zdr_kdpzh_dict["freq_band"]
                    )
        elif parametrization == "Gorgucci":
            if zdr_kdpzh_dict is None or "freq_band" not in zdr_kdpzh_dict:
                warn(
                    "Unable to plot theoretical self-consistency curve. "
                    "Frequency band not provided"
                )
            else:
                zdr_th = np.arange(0.0, 3.5 + bins_zdr_step, bins_zdr_step)
                if zdr_kdpzh_dict["freq_band"] == "C":
                    zdr_lin = np.ma.power(10.0, 0.1 * zdr_th)
                    kdpzh_th = 18.2 * np.power(zdr_lin, -1.28)
                    labely = "10e5*KDP/ZH^0.95 [(deg*m3)/(km*mm6)"
                else:
                    warn(
                        "Unable to plot theoretical self-consistency curve. "
                        "Unknown frequency band " + zdr_kdpzh_dict["freq_band"]
                    )
        elif parametrization == "Vaccarono":
            if zdr_kdpzh_dict is None or "freq_band" not in zdr_kdpzh_dict:
                warn(
                    "Unable to plot theoretical self-consistency curve. "
                    "Frequency band not provided"
                )
            else:
                zdr_th = np.arange(0.0, 3.5 + bins_zdr_step, bins_zdr_step)
                if zdr_kdpzh_dict["freq_band"] == "C":
                    zdr_lin = np.ma.power(10.0, 0.1 * zdr_th)
                    kdpzh_th = 17.7 * np.power(zdr_lin, -2.09)
                    labely = "10e5*KDP^0.85/ZH^0.91 [(deg*m3)/(km*mm6)"
                else:
                    warn(
                        "Unable to plot theoretical self-consistency curve. "
                        "Unknown frequency band " + zdr_kdpzh_dict["freq_band"]
                    )

    if kdpzh_th is not None and parametrization == "Gorgucci":
        kdpzh = 1e5 * kdp / np.ma.power(zh, 0.95)
    elif kdpzh_th is not None and parametrization == "Vaccarono":
        kdpzh = 1e5 * np.ma.power(kdp, 0.85) / np.ma.power(zh, 0.95)
    else:
        kdpzh = 1e5 * kdp / zh

    if retrieve_relation:
        zdr_min = zdr.min()
        zdr_max = zdr.max()

        ind_min = np.where(bins_zdr_edges >= zdr_min)[0][0]
        ind_max = np.where(bins_zdr_edges < zdr_max)[0][-1]
        zdr_int = bins_zdr_edges[ind_min : ind_max + 1]
        zdr_par = bins_zdr_centers[ind_min:ind_max]
        kdpzh_par = np.ma.masked_all(zdr_par.size)
        for ind, zdr_left in enumerate(zdr_int[:-1]):
            zdr_right = zdr_int[ind + 1]
            ind_val = np.where(np.logical_and(zdr >= zdr_left, zdr < zdr_right))[0]
            if ind_val.size == 0:
                continue
            kdpzh_par[ind] = np.ma.median(kdpzh[ind_val])

    # prepare histogram
    zdr[zdr < bins_zdr_centers[0]] = bins_zdr_centers[0]
    zdr[zdr > bins_zdr_centers[-1]] = bins_zdr_centers[-1]

    kdpzh[kdpzh < bins_kdpzh_centers[0]] = bins_kdpzh_centers[0]
    kdpzh[kdpzh > bins_kdpzh_centers[-1]] = bins_kdpzh_centers[-1]

    hist2d, bins_zdr_edges, bins_kdpzh_edges = np.histogram2d(
        zdr, kdpzh, bins=[bins_zdr_edges, bins_kdpzh_edges]
    )

    hist2d = np.ma.masked_equal(hist2d, 0)

    if normalize:
        hist2d = hist2d / np.expand_dims(np.ma.sum(hist2d, axis=-1), axis=1)
        clabel = "Occurrence density\n(For each ZDR bin)"
    else:
        clabel = "Number of gates"
        vmin = None
        vmax = None

    if not retrieve_relation and not plot_theoretical:
        fname_list = _plot_time_range(
            bins_zdr_edges,
            bins_kdpzh_edges,
            hist2d,
            None,
            fname_list,
            titl="self-consistency",
            xlabel="ZDR [dB]",
            ylabel=labely,
            clabel=clabel,
            vmin=vmin,
            vmax=vmax,
            save_fig=True,
            dpi=dpi,
        )

        return fname_list

    fig, ax = _plot_time_range(
        bins_zdr_edges,
        bins_kdpzh_edges,
        hist2d,
        None,
        fname_list,
        titl="self-consistency",
        xlabel="ZDR [dB]",
        ylabel=labely,
        clabel=clabel,
        vmin=vmin,
        vmax=vmax,
        save_fig=False,
        dpi=dpi,
    )

    ax.autoscale(False)

    if retrieve_relation:
        ax.plot(zdr_par, kdpzh_par, "rx-", linewidth=3)

    if plot_theoretical and kdpzh_th is not None:
        ax.plot(zdr_th, kdpzh_th, "r", linewidth=3)

    for fname in fname_list:
        fig.savefig(fname, dpi=72)
    plt.close(fig)

    return fname_list


def plot_selfconsistency_instrument2(
    zdr,
    kdp,
    zh,
    fname_list,
    bins_zh_step=0.5,
    bins_zh_min=-30.0,
    bins_zh_max=80.0,
    normalize=True,
    vmin=0.0,
    vmax=0.01,
    parametrization="None",
    zdr_kdpzh_dict=None,
    dpi=72,
):
    """
    plots the ZDR-KDP/ZH relationship obtained by an instrument. The
    theoretical curve and the retrieved curve

    Parameters
    ----------
    zdr, kdp, zh : 1D ndarray
        The valid values of ZDR [dB], KDP [deg/km] and Zh [dBZ] collected by
        the instrument
    fname_list : list of str
        list of names of the files where to store the plot
    bins_zh_step : float
        The step of the ZH of the histogram [dB]
    bins_zh_min, bins_zh_max : float
        The limits of the ZH of the histograms (bins center) [dBZ]
    normalize : Bool
        If True the occurrence density of ZH/KDP for each ZDR bin is going to
        be represented. Otherwise it will show the number of gates at each bin
    vmin, vmax : float
        min and max values of the colorbar
    parametrization : str
        The type of parametrization for the self-consistency curves. Can be
        'None', 'Gourley', 'Wolfensberger', 'Louf', 'Gorgucci' or 'Vaccarono'.
        'None' will use tables contained in zdr_kdpzh_dict. The parametrized
        curves are obtained from literature except for Wolfensberger that was
        derived from disdrometer data obtained by MeteoSwiss and EPFL. All
        parametrizations are valid for C-band only except that of Gourley.
    zdr_kdpzh_dict : dict
        dictionary containing a look up table relating ZDR with KDP/Zh for
        different elevations and the frequency band of the radar
    dpi : int
        dots per inch

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """

    # prepare bins for histogram
    bins_zh_centers = np.arange(bins_zh_min, bins_zh_max + bins_zh_step, bins_zh_step)
    bins_zh_edges = np.append(
        bins_zh_centers - bins_zh_step / 2.0, bins_zh_centers[-1] + bins_zh_step / 2.0
    )

    if parametrization == "None":
        if zdr_kdpzh_dict is not None and "zdr_kdpzh" in zdr_kdpzh_dict:
            zdr_th_table = zdr_kdpzh_dict["zdr_kdpzh"][0][0]
            kdpzh_th_table = zdr_kdpzh_dict["zdr_kdpzh"][0][1]

            # sort values by ZDR
            ind_zdr_sorted = np.argsort(zdr)
            zh = zh[ind_zdr_sorted]
            kdp = kdp[ind_zdr_sorted]

            # get value of KDP/ZH as interpolation of look up table values
            kdpzh_th = np.interp(zdr[ind_zdr_sorted], zdr_th_table, kdpzh_th_table)

        else:
            warn("Unable to plot self-consistency curve. " "Relationship not provided")
            return None
    elif parametrization == "Gourley":
        if zdr_kdpzh_dict is None or "freq_band" not in zdr_kdpzh_dict:
            warn(
                "Unable to plot theoretical self-consistency curve. "
                "Frequency band not provided"
            )
            return None
        if zdr_kdpzh_dict["freq_band"] == "S":
            kdpzh_th = 1e-5 * (
                3.696 - 1.963 * zdr + 0.504 * zdr * zdr - 0.051 * zdr * zdr * zdr
            )
        elif zdr_kdpzh_dict["freq_band"] == "C":
            kdpzh_th = 1e-5 * (
                6.746 - 2.970 * zdr + 0.711 * zdr * zdr - 0.079 * zdr * zdr * zdr
            )
        elif zdr_kdpzh_dict["freq_band"] == "X":
            kdpzh_th = 1e-5 * (
                11.74 - 4.020 * zdr - 0.140 * zdr * zdr + 0.130 * zdr * zdr * zdr
            )
        else:
            warn(
                "Unable to plot theoretical self-consistency curve. "
                "Unknown frequency band " + zdr_kdpzh_dict["freq_band"]
            )
            return None
    elif parametrization == "Wolfensberger":
        if zdr_kdpzh_dict is None or "freq_band" not in zdr_kdpzh_dict:
            warn(
                "Unable to plot theoretical self-consistency curve. "
                "Frequency band not provided"
            )
            return None
        if zdr_kdpzh_dict["freq_band"] == "C":
            kdpzh_th = 1e-5 * (
                3.199 * np.ma.exp(-7.767e-1 * zdr) - 0.4436 * zdr + 3.464
            )
        else:
            warn(
                "Unable to plot theoretical self-consistency curve. "
                "Unknown frequency band " + zdr_kdpzh_dict["freq_band"]
            )
            return None
    elif parametrization == "Louf":
        if zdr_kdpzh_dict is None or "freq_band" not in zdr_kdpzh_dict:
            warn(
                "Unable to plot theoretical self-consistency curve. "
                "Frequency band not provided"
            )
            return None
        if zdr_kdpzh_dict["freq_band"] == "C":
            kdpzh_th = 1e-5 * (
                6.607 - 4.577 * zdr + 1.577 * zdr * zdr - 0.23 * zdr * zdr * zdr
            )
        else:
            warn(
                "Unable to plot theoretical self-consistency curve. "
                "Unknown frequency band " + zdr_kdpzh_dict["freq_band"]
            )
            return None
    elif parametrization == "Gorgucci":
        if zdr_kdpzh_dict is None or "freq_band" not in zdr_kdpzh_dict:
            warn(
                "Unable to plot theoretical self-consistency curve. "
                "Frequency band not provided"
            )
            return None
        if zdr_kdpzh_dict["freq_band"] == "C":
            zdr_lin = np.ma.power(10.0, 0.1 * zdr)
            kdpzh_th = 1e-5 * (18.2 * np.power(zdr_lin, -1.28))
        else:
            warn(
                "Unable to plot theoretical self-consistency curve. "
                "Unknown frequency band " + zdr_kdpzh_dict["freq_band"]
            )
            return None
    elif parametrization == "Vaccarono":
        if zdr_kdpzh_dict is None or "freq_band" not in zdr_kdpzh_dict:
            warn(
                "Unable to plot theoretical self-consistency curve. "
                "Frequency band not provided"
            )
            return None
        if zdr_kdpzh_dict["freq_band"] == "C":
            zdr_lin = np.ma.power(10.0, 0.1 * zdr)
            kdpzh_th = 1e-5 * (17.7 * np.power(zdr_lin, -2.09))
        else:
            warn(
                "Unable to plot theoretical self-consistency curve. "
                "Unknown frequency band " + zdr_kdpzh_dict["freq_band"]
            )
            return None

    if parametrization == "Gorgucci":
        zh_th = 10.0 * np.ma.log10(np.ma.power(kdp / kdpzh_th, 1 / 0.95))
    elif kdpzh_th is not None and parametrization == "Vaccarono":
        zh_th = 10.0 * np.ma.log10(
            np.ma.power(np.ma.power(kdp, 0.85) / kdpzh_th, 1 / 0.95)
        )
    else:
        zh_th = 10.0 * np.ma.log10(kdp / kdpzh_th)

    # prepare histogram
    zh[zh < bins_zh_centers[0]] = bins_zh_centers[0]
    zh[zh > bins_zh_centers[-1]] = bins_zh_centers[-1]

    zh_th[zh_th < bins_zh_centers[0]] = bins_zh_centers[0]
    zh_th[zh_th > bins_zh_centers[-1]] = bins_zh_centers[-1]

    hist2d, bins_zh_edges, _ = np.histogram2d(
        zh, zh_th, bins=[bins_zh_edges, bins_zh_edges]
    )
    hist2d = np.ma.masked_equal(hist2d, 0)

    if normalize:
        hist2d = hist2d / np.expand_dims(np.ma.sum(hist2d, axis=-1), axis=1)
        clabel = "Occurrence density\n(For each ZH bin)"
    else:
        clabel = "Number of gates"
        vmin = None
        vmax = None

    fig, ax = _plot_time_range(
        bins_zh_edges,
        bins_zh_edges,
        hist2d,
        None,
        fname_list,
        titl="self-consistency",
        xlabel="ZH measured [dBZ]",
        ylabel="ZH selfcons [dBZ]",
        clabel=clabel,
        vmin=vmin,
        vmax=vmax,
        save_fig=False,
        dpi=dpi,
    )

    ax.autoscale(False)
    ax.plot(bins_zh_centers, bins_zh_centers, "rx-", linewidth=3)

    for fname in fname_list:
        fig.savefig(fname, dpi=72)
    plt.close(fig)

    return fname_list


def plot_scatter_comp(
    values,
    fname_list,
    labelx="Sensor 1",
    labely=None,
    titl="Scatter",
    axis=None,
    metadata=None,
    dpi=72,
    ax=None,
    fig=None,
    save_fig=True,
    labels=None,
    write_stats=True,
    alpha=0.5,
    markersize=20,
):
    """
    Plot scatter comparisons of multiple arrays against the first array.

    Parameters
    ----------
    values : list of tuple(array-like, array-like)
        List of (x, y) pairs. Each tuple is plotted as one scatter series.
    fname_list : list of str
        List of names of the files where to store the plot.
    labelx : str
        Label of the X axis.
    labely : str or None
        Label of the Y axis. If None, defaults to "Estimated values".
    titl : str
        Figure title.
    axis : str
        Type of axis. If "equal", enforce equal axis limits and draw 1:1 line.
    metadata : str
        A string containing metadata.
    dpi : int
        Dots per inch.
    fig : Figure
        Figure to plot on. If None, a new figure will be created.
    ax : Axis
        Axis to plot on. If fig is None, a new axis will be created.
    save_fig : bool
        If True save the figure. If False, return figure and axis handles.
    labels : list of str or None
        Labels for arrays 2..N. Length must be len(values)-1 if provided.
    write_stats : bool
        If True, write stats for each series on the plot.
    alpha : float
        Transparency of scatter points.
    markersize : float
        Scatter marker size.

    Returns
    -------
    fname_list : list of str
        List of names of the created plots, if save_fig is True.
    or
    (fig, ax) : tuple
        Figure and axis handles, if save_fig is False.
    """

    if values is None or len(values) == 0:
        raise ValueError("values must be a non-empty list of (x, y) tuples")

    xy_pairs = []
    for i, pair in enumerate(values):
        if not isinstance(pair, (tuple, list)) or len(pair) != 2:
            raise ValueError(
                f"values[{i}] must be a tuple/list of two elements: (x, y)"
            )

        x = np.asanyarray(pair[0])
        y = np.asanyarray(pair[1])

        if len(x) != len(y):
            raise ValueError(
                f"x and y must have the same length for values[{i}]. "
                f"Got len(x)={len(x)}, len(y)={len(y)}."
            )

        xy_pairs.append((x, y))

    if labels is None:
        labels = [f"Series {i + 1}" for i in range(len(xy_pairs))]
    elif len(labels) != len(xy_pairs):
        raise ValueError("labels must have length len(values)")

    if labely is None:
        labely = "Estimated values"

    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=[7, 7], dpi=dpi)

    # Remove old stat texts only, keep metadata if present
    old_texts = list(ax.texts)
    for txt in old_texts:
        txt.remove()

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    markers = ["o", "s", "^", "D", "v", "P", "X", "<", ">"]

    all_valid = []

    for i, (x, y) in enumerate(xy_pairs):
        mask = np.isfinite(x) & np.isfinite(y)
        xv = x[mask]
        yv = y[mask]

        all_valid.append((xv, yv))

        color = colors[i % len(colors)]
        marker = markers[i % len(markers)]

        ax.scatter(
            xv,
            yv,
            s=markersize,
            alpha=alpha,
            color=color,
            marker=marker,
            edgecolors="none",
            label=labels[i],
        )

    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_title(titl)

    # Determine common axis range from all valid data
    all_data = []
    for xv, yv in all_valid:
        if xv.size > 0:
            all_data.append(xv)
        if yv.size > 0:
            all_data.append(yv)

    if len(all_data) > 0:
        data_min = min(np.nanmin(a) for a in all_data)
        data_max = max(np.nanmax(a) for a in all_data)
    else:
        data_min, data_max = 0.0, 1.0

    if axis == "equal":
        ax.axis([data_min, data_max, data_min, data_max])
        ax.plot([data_min, data_max], [data_min, data_max], "k--", linewidth=1)
        ax.set(adjustable="box", aspect="equal")

    if metadata is not None:
        ax.text(
            0.05,
            0.95,
            metadata,
            horizontalalignment="left",
            verticalalignment="top",
            transform=ax.transAxes,
        )

    if write_stats:
        stats_lines = []

        for i, (xv, yv) in enumerate(all_valid):
            if xv.size == 0:
                rmse = np.nan
                bias = np.nan
                corr = np.nan
            else:
                rmse = np.sqrt(np.nanmean((xv - yv) ** 2))
                bias = np.nanmean(yv - xv)

                if xv.size > 1 and np.nanstd(xv) > 0 and np.nanstd(yv) > 0:
                    corr = np.corrcoef(xv, yv)[0, 1]
                else:
                    corr = np.nan

            stats_lines.append(
                f"{labels[i]}: RMSE={rmse:.2f}, Bias={bias:.2f}, Corr={corr:.2f}"
            )

        ax.text(
            0.01,
            0.98,
            "\n".join(stats_lines),
            horizontalalignment="left",
            verticalalignment="top",
            transform=ax.transAxes,
            fontsize=9,
            bbox=dict(
                boxstyle="round",
                facecolor="white",
                alpha=0.7,
                edgecolor="none",
            ),
        )

    ax.grid()
    ax.legend(loc="best")
    fig.tight_layout()

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        return fname_list


def plot_sun_hits(field, field_name, fname_list, prdcfg):
    """
    plots the sun hits

    Parameters
    ----------
    radar : Radar object
        object containing the radar data to plot
    field_name : str
        name of the radar field to plot
    altitude : float
        the altitude [m MSL] to be plotted
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    azmin = prdcfg["sunhitsImageConfig"]["azmin"]
    azmax = prdcfg["sunhitsImageConfig"]["azmax"]
    elmin = prdcfg["sunhitsImageConfig"]["elmin"]
    elmax = prdcfg["sunhitsImageConfig"]["elmax"]

    field_dict = pyart.config.get_metadata(field_name)

    # display data
    dpi = 72
    if "dpi" in prdcfg["sunhitsImageConfig"]:
        dpi = prdcfg["sunhitsImageConfig"]["dpi"]

    fig = plt.figure(
        figsize=[
            prdcfg["sunhitsImageConfig"]["xsize"],
            prdcfg["sunhitsImageConfig"]["ysize"],
        ],
        dpi=dpi,
    )
    ax = fig.add_subplot(111)
    cmap = pyart.config.get_field_colormap(field_name)
    vmin, vmax = pyart.config.get_field_limits(field_name)
    titl = (
        prdcfg["timeinfo"].strftime("%Y-%m-%d")
        + "\n"
        + get_field_name(field_dict, field_name)
    )

    cax = ax.imshow(
        field,
        extent=(azmin, azmax, elmin, elmax),
        origin="lower",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation="none",
    )
    ax.set_xlabel("rad_az-sun_az (deg)")
    ax.set_ylabel("rad_el-sun_el (deg)")
    ax.set_title(titl)

    # plot the colorbar and set the label.
    label = get_colobar_label(field_dict, field_name)
    cb = fig.colorbar(cax)
    cb.set_label(label)

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def _plot_sunscan(
    rad_az,
    rad_el,
    rad_data,
    sun_hits,
    field_name,
    fname_list,
    titl="AZ-EL sunscan plot",
    xlabel="Azimuth (deg)",
    ylabel="Elevation (deg)",
    clabel=None,
    vmin=None,
    vmax=None,
    figsize=(10, 8),
    save_fig=True,
    dpi=72,
):
    """
    plots a AZ-EL plot of a sunscan

    Parameters
    ----------
    rad_time : 1D array
        array containing the x dimension (typically time)
    rad_range : 1D array
        array containing the y dimension (typically range)
    rad_data : 2D array
        array containing the data to plot
    sun_hits : dict
        dictionary containing the sun hits data
    field_name : str or None
        field name. Used to define plot characteristics
    fname_list : list of str
        list of names of the files where to store the plot
    titl : str
        Plot title
    xlabel, ylabel : str
        x- and y-axis labels
    clabel : str or None
        colorbar label
    vmin, vmax : float
        min and max values of the color bar
    figsize : list
        figure size [xsize, ysize]
    save_fig : bool
        If true the figure is saved in files fname_list. Otherwise the fig and
        ax is returned
    dpi : int
        dpi

    Returns
    -------
    fname_list : list of str
        list of names of the created plots or
    fig, ax : object
        handles to fig and ax objects

    """
    # display data
    norm = None
    cmap = None
    ticks = None
    ticklabs = None
    if field_name is not None:
        field_dict = pyart.config.get_metadata(field_name)
        if clabel is None:
            clabel = get_colobar_label(field_dict, field_name)

        cmap = pyart.config.get_field_colormap(field_name)

        norm, ticks, ticklabs = get_norm(field_name)
        if vmin is None or vmax is None:
            vmin = vmax = None
            if norm is None:  # if norm is set do not override with vmin/vmax
                vmin, vmax = pyart.config.get_field_limits(field_name)
        else:
            norm = None
    else:
        if clabel is None:
            clabel = "value"
        if vmin is None:
            vmin = np.ma.min(rad_data)
        if vmax is None:
            vmax = np.ma.max(rad_data)

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111)

    T, R = np.meshgrid(rad_az, rad_el)
    cax = ax.pcolormesh(
        T, R, np.ma.transpose(rad_data), cmap=cmap, vmin=vmin, vmax=vmax, norm=norm
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(titl)

    # Add scatter points of detected sun hits location
    ax.scatter(sun_hits["rad_az"], sun_hits["rad_el"], s=7.0, c="red")

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label(clabel)

    # Make a tight layout
    fig.tight_layout()

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

        return fname_list
    return fig, ax


def _plot_time_range(
    rad_time,
    rad_range,
    rad_data,
    field_name,
    fname_list,
    titl="Time-Range plot",
    xlabel="time (s from start time)",
    ylabel="range (Km)",
    clabel=None,
    vmin=None,
    vmax=None,
    figsize=(10, 8),
    save_fig=True,
    dpi=72,
    plot_max=False,
):
    """
    Plot a time-range plot.

    Parameters
    ----------
    rad_time : 1D array
        Array containing the x dimension, typically time.
    rad_range : 1D array
        Array containing the y dimension, typically range.
    rad_data : 2D array
        Array containing the data to plot. Expected shape is
        (ntime, nrange).
    field_name : str or None
        Field name. Used to define plot characteristics.
    fname_list : list of str
        Names of the files where the plot will be stored.
    titl : str, optional
        Plot title.
    xlabel, ylabel : str, optional
        X- and y-axis labels.
    clabel : str or None, optional
        Colorbar label.
    vmin, vmax : float or None, optional
        Minimum and maximum values of the colorbar.
    figsize : tuple, optional
        Figure size as (xsize, ysize).
    save_fig : bool, optional
        If True, save the figure to each file in ``fname_list``.
        Otherwise, return the figure and axes.
    dpi : int, optional
        Figure resolution in dots per inch.
    plot_max : bool, optional
        If True, mark the maximum finite, unmasked value with a red dot
        and annotate its time, range, and value.

    Returns
    -------
    fname_list : list of str
        Names of the created plots when ``save_fig`` is True.
    fig, ax : matplotlib objects
        Figure and axes handles when ``save_fig`` is False.
    """
    # Display characteristics
    norm = None
    cmap = None
    ticks = None
    ticklabs = None

    if field_name is not None:
        field_dict = pyart.config.get_metadata(field_name)

        if clabel is None:
            clabel = get_colobar_label(field_dict, field_name)

        cmap = pyart.config.get_field_colormap(field_name)
        norm, ticks, ticklabs = get_norm(field_name)

        if vmin is None or vmax is None:
            vmin = vmax = None

            # If norm is set, do not override it with vmin/vmax.
            if norm is None:
                vmin, vmax = pyart.config.get_field_limits(field_name)
        else:
            norm = None

    else:
        if clabel is None:
            clabel = "value"

        if vmin is None:
            vmin = np.ma.min(rad_data)

        if vmax is None:
            vmax = np.ma.max(rad_data)

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111)

    T, R = np.meshgrid(rad_time, rad_range)

    cax = ax.pcolormesh(
        T,
        R,
        np.ma.transpose(rad_data),
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        norm=norm,
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(titl)

    cb = fig.colorbar(cax)

    if ticks is not None:
        cb.set_ticks(ticks)

    if ticklabs:
        cb.set_ticklabels(ticklabs)

    cb.set_label(clabel)

    if plot_max:
        data = np.ma.masked_invalid(np.ma.asarray(rad_data))

        if data.count() > 0:
            # rad_data is expected to have shape (ntime, nrange).
            itime, irange = np.unravel_index(
                np.ma.argmax(data),
                data.shape,
            )

            max_time = rad_time[itime]
            max_range = rad_range[irange]
            max_value = data[itime, irange].item()

            ax.plot(
                max_time,
                max_range,
                marker="o",
                markersize=7,
                markerfacecolor="red",
                markeredgecolor="white",
                markeredgewidth=1.0,
                linestyle="none",
                zorder=10,
            )

            annotation = (
                f"max = {max_value:.2f}\n"
                f"x = {max_time:.2f}\n"
                f"y = {max_range:.2f}"
            )

            ax.annotate(
                annotation,
                xy=(max_time, max_range),
                xytext=(10, 10),
                textcoords="offset points",
                ha="left",
                va="bottom",
                fontsize=9,
                color="red",
                bbox={
                    "boxstyle": "round,pad=0.3",
                    "facecolor": "white",
                    "edgecolor": "red",
                    "alpha": 0.85,
                },
                arrowprops={
                    "arrowstyle": "->",
                    "color": "red",
                },
                zorder=11,
            )

    fig.tight_layout()

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi, bbox_inches="tight")

        plt.close(fig)
        return fname_list

    return fig, ax


def plot_bias_over_range(
    range_centers,
    meanbias,
    percentile_dict,
    fname_list,
    nsamples=None,
    reference_label="rad2 - rad1",
    field_name="field",
    titl=None,
    labelx="Range to radar 1 (m)",
    labely=None,
    mean_label="mean bias",
    range_marker=None,
    ref_value=0.0,
    plot_npoints=True,
    bias_vmin=None,
    bias_vmax=None,
    figsize=(12, 7),
    dpi=72,
    percentile_style="band",
):
    """
    Plot intercomparison bias over range.

    Percentiles can be displayed either as discrete lines or as shaded
    percentile envelopes. The default style, "band", is recommended for
    readability because it highlights the median and the percentile spread
    without overloading the figure.

    Parameters
    ----------
    range_centers : array-like
        Range bin centers in meters.

    meanbias : array-like
        Mean bias per range bin.

    percentile_dict : dict
        Dictionary mapping percentile value to percentile bias array.
        Example:
            {
                5: quant05_vec,
                25: quant25_vec,
                50: median_vec,
                75: quant75_vec,
                95: quant95_vec,
            }

    fname_list : list of str
        Output figure filenames.

    nsamples : array-like, optional
        Number of samples per range bin.

    reference_label : str
        Bias label, e.g. "310CHX - 144CHX".

    field_name : str
        Name of plotted field.

    titl : str, optional
        Plot title.

    labelx : str
        x-axis label.

    labely : str, optional
        y-axis label. If None, generated from reference_label.

    mean_label : str
        Legend label for the mean bias curve.

    range_marker : float, optional
        Vertical range marker in meters.

    ref_value : float
        Horizontal reference line value.

    plot_npoints : bool
        If True and nsamples is not None, plot sample count on secondary axis.

    bias_vmin, bias_vmax : float, optional
        y-axis limits.

    figsize : tuple
        Figure size.

    dpi : int
        Figure resolution.

    percentile_style : str
        How to display percentiles. Accepted values are:
            "band":
                Display percentile ranges as shaded envelopes.
                If available, 5-95 and 25-75 bands are plotted, and the
                50th percentile is shown as a discrete median line.
            "lines":
                Display all percentiles as thin, color-sorted lines.

        Default: "band".

    Returns
    -------
    fname_list : list of str
        List of saved figure filenames.
    """

    range_centers = np.asarray(range_centers, dtype=float)
    meanbias = np.asarray(meanbias, dtype=float)

    percentile_dict = {
        float(perc): np.asarray(values, dtype=float)
        for perc, values in percentile_dict.items()
    }

    fig, ax1 = plt.subplots(figsize=figsize)

    # ------------------------------------------------------------------
    # Percentile display
    # ------------------------------------------------------------------
    if percentile_style == "band":
        # Wide envelope: 5-95 %
        if 5.0 in percentile_dict and 95.0 in percentile_dict:
            ax1.fill_between(
                range_centers,
                percentile_dict[5.0],
                percentile_dict[95.0],
                alpha=0.12,
                label="5-95th percentile",
            )

        # Inner envelope: 25-75 %
        if 25.0 in percentile_dict and 75.0 in percentile_dict:
            ax1.fill_between(
                range_centers,
                percentile_dict[25.0],
                percentile_dict[75.0],
                alpha=0.25,
                label="25-75th percentile",
            )

        # Median
        if 50.0 in percentile_dict:
            ax1.plot(
                range_centers,
                percentile_dict[50.0],
                "-",
                linewidth=1.4,
                alpha=0.9,
                label="median",
            )

    elif percentile_style == "lines":
        # Color-sorted but visually discrete percentile lines
        percentiles = sorted(percentile_dict.keys())
        cmap = plt.get_cmap("viridis")

        for i, perc in enumerate(percentiles):
            color = cmap(i / max(len(percentiles) - 1, 1))

            linewidth = 1.4 if perc == 50.0 else 1.0
            alpha = 0.9 if perc == 50.0 else 0.45
            linestyle = "-" if perc == 50.0 else "--"

            ax1.plot(
                range_centers,
                percentile_dict[perc],
                linestyle=linestyle,
                linewidth=linewidth,
                alpha=alpha,
                color=color,
                label=f"{perc:g}th percentile",
            )

    else:
        raise ValueError("percentile_style must be either 'band' or 'lines'")

    # ------------------------------------------------------------------
    # Mean bias on top
    # ------------------------------------------------------------------
    ax1.plot(
        range_centers,
        meanbias,
        "-o",
        linewidth=2.0,
        markersize=4,
        label=mean_label,
        zorder=5,
    )

    ax1.axhline(ref_value, linestyle="--", linewidth=1)

    if range_marker is not None:
        ax1.axvline(range_marker, linestyle="--", linewidth=1)

    ax1.set_xlabel(labelx)

    if labely is None:
        labely = f"Bias ({reference_label}) [dB]"

    ax1.set_ylabel(labely)

    if bias_vmin is not None or bias_vmax is not None:
        ax1.set_ylim(bias_vmin, bias_vmax)

    if titl is None:
        titl = f"{field_name} bias over range"

    ax1.set_title(titl)
    ax1.grid(True, axis="y", alpha=0.3)

    lines, labels = ax1.get_legend_handles_labels()

    if plot_npoints and nsamples is not None:
        nsamples = np.asarray(nsamples)

        ax2 = ax1.twinx()
        ax2.plot(
            range_centers,
            nsamples,
            "k--",
            linewidth=1,
            alpha=0.5,
            label="N samples",
        )
        ax2.set_ylabel("Number of samples")

        lines2, labels2 = ax2.get_legend_handles_labels()
        lines.extend(lines2)
        labels.extend(labels2)

    ax1.legend(lines, labels)

    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(
            fname,
            dpi=dpi,
            bbox_inches="tight",
        )

    plt.close(fig)

    return fname_list
