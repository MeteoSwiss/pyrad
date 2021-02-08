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

from warnings import warn
from copy import deepcopy

import numpy as np

try:
    import cartopy
    from cartopy.io.img_tiles import Stamen
    _CARTOPY_AVAILABLE = True
except ImportError:
    _CARTOPY_AVAILABLE = False

import matplotlib as mpl
mpl.use('Agg')

# Increase a bit font size
mpl.rcParams.update({'font.size': 16})
mpl.rcParams.update({'font.family':  "sans-serif"})

import matplotlib.pyplot as plt

import pyart

from .plots_aux import get_colobar_label, get_field_name, get_norm

from ..util.radar_utils import compute_quantiles_from_hist


def plot_pos(lat, lon, alt, fname_list, ax=None, fig=None, save_fig=True,
             sort_altitude='No', dpi=72, alpha=1., cb_label='height [m MSL]',
             titl='Position', xlabel='Lon [Deg]', ylabel='Lat [Deg]',
             limits=None, vmin=None, vmax=None):
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
    if sort_altitude in ('Lowest_on_top', 'Highest_on_top'):
        ind = np.argsort(alt)
        if sort_altitude == 'Lowest_on_top':
            ind = ind[::-1]
        lat = lat[ind]
        lon = lon[ind]
        alt = alt[ind]

    if vmin is None:
        vmin = alt.min()
    if vmax is None:
        vmax = alt.max()
    marker = 'x'
    col = alt
    cmap = 'viridis'
    norm = plt.Normalize(vmin, vmax)

    if fig is None:
        fig = plt.figure(figsize=[10, 8], dpi=dpi)
        ax = fig.add_subplot(111, aspect='equal')
    else:
        ax.autoscale(False)

    cax = ax.scatter(
        lon, lat, c=col, marker=marker, alpha=alpha, cmap=cmap, norm=norm)

    if limits is not None:
        ax.set_xlim(limits[0], limits[1])
        ax.set_ylim(limits[2], limits[3])

    # plot colorbar
    cb = fig.colorbar(cax, orientation='horizontal')
    cb.set_label(cb_label)

    ax.set_title(titl)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Turn on the grid
    ax.grid()

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list

    return (fig, ax)


def plot_pos_map(lat, lon, alt, fname_list, ax=None, fig=None, save_fig=True,
                 sort_altitude='No', dpi=72, alpha=1.,
                 cb_label='height [m MSL]', titl='Position',
                 xlabel='Lon [Deg]', ylabel='Lat [Deg]', limits=None,
                 vmin=None, vmax=None, lon_step=0.3, lat_step=0.1,
                 background_zoom=8):
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
        warn('Unable to plot trajectory on a map. Cartopy not available')
        return ['']

    if sort_altitude in ('Lowest_on_top', 'Highest_on_top'):
        ind = np.argsort(alt)
        if sort_altitude == 'Lowest_on_top':
            ind = ind[::-1]
        lat = lat[ind]
        lon = lon[ind]
        alt = alt[ind]

    if vmin is None:
        vmin = alt.min()
    if vmax is None:
        vmax = alt.max()
    marker = 'x'
    col = alt
    cmap = 'viridis'
    norm = plt.Normalize(vmin, vmax)

    # get background map instance
    stamen_terrain = Stamen('terrain-background')
    projection = cartopy.crs.PlateCarree()

    # set map limits and lat/lon lines
    if limits is None:
        limits = [np.min(lon), np.max(lon), np.min(lat), np.max(lat)]

    lon_lines = np.arange(limits[0], limits[1]+1, lon_step)
    lat_lines = np.arange(limits[2], limits[3]+1, lat_step)

    if fig is None:
        fig = plt.figure(figsize=[10, 8], dpi=dpi)

        # draw background
        ax = fig.add_subplot(111, projection=stamen_terrain.crs)
        ax.set_extent(limits, crs=projection)
        ax.add_image(stamen_terrain, background_zoom)

        # add countries
        countries = cartopy.feature.NaturalEarthFeature(
            category='cultural',
            name='admin_0_countries',
            scale='10m',
            facecolor='none')
        ax.add_feature(countries, edgecolor='black')

        # draw grid lines and labels
        gl = ax.gridlines(xlocs=lon_lines, ylocs=lat_lines, draw_labels=True)
        gl.xlabels_top = False
        gl.ylabels_right = False

        ax.text(0.5, -0.2, xlabel, va='bottom', ha='center',
                rotation='horizontal', rotation_mode='anchor',
                transform=ax.transAxes)

        ax.text(-0.1, 0.55, ylabel, va='bottom', ha='center',
                rotation='vertical', rotation_mode='anchor',
                transform=ax.transAxes)
    else:
        ax.autoscale(False)

    # plot data
    cax = ax.scatter(
        lon, lat, c=col, marker=marker, alpha=alpha, cmap=cmap, norm=norm,
        transform=projection)

    # plot colorbar
    cb = fig.colorbar(cax, orientation='horizontal')
    cb.set_label(cb_label)

    ax.set_title(titl)

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list

    return (fig, ax)


def plot_density(hist_obj, hist_type, field_name, ind_sweep, prdcfg,
                 fname_list, quantiles=[25., 50., 75.], ref_value=0.,
                 vmin=None, vmax=None):
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
        Minim and maximum extend of the vertical axis

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    hist_obj_aux = hist_obj.extract_sweeps([ind_sweep])
    if hist_obj_aux.scan_type == 'ppi':
        ang = np.sort(hist_obj_aux.azimuth['data'])
        ind_ang = np.argsort(hist_obj_aux.azimuth['data'])
        ang_min = np.min(hist_obj_aux.azimuth['data'])
        ang_max = np.max(hist_obj_aux.azimuth['data'])
        ang_step = ang[1]-ang[0]
        field = hist_obj_aux.fields[field_name]['data'][ind_ang, :]
        labelx = 'azimuth angle (degrees)'
    elif hist_obj_aux.scan_type == 'rhi':
        ang = np.sort(hist_obj_aux.elevation['data'])
        ind_ang = np.argsort(hist_obj_aux.elevation['data'])
        ang_min = np.min(hist_obj_aux.elevation['data'])
        ang_max = np.max(hist_obj_aux.elevation['data'])
        ang_step = ang[1]-ang[0]
        field = hist_obj_aux.fields[field_name]['data'][ind_ang, :]
        labelx = 'elevation angle (degrees)'
    else:
        field = hist_obj_aux.fields[field_name]['data']
        ang = np.array(range(hist_obj_aux.nrays))
        ang_min = 0
        ang_max = hist_obj_aux.nrays-1
        ang_step = ang[1]-ang[0]
        labelx = 'ray number'

    # compute percentiles of the histogram
    az_percentile_ref = np.ma.masked_all(len(ang))
    az_percentile_low = deepcopy(az_percentile_ref)
    az_percentile_high = deepcopy(az_percentile_ref)
    for ray in range(len(ang)):
        quantiles, values_ray = compute_quantiles_from_hist(
            hist_obj.range['data'], field[ray, :], quantiles=quantiles)

        az_percentile_low[ray] = values_ray[0]
        az_percentile_ref[ray] = values_ray[1]
        az_percentile_high[ray] = values_ray[2]

    quantiles, values_sweep = compute_quantiles_from_hist(
        hist_obj.range['data'], np.ma.sum(field, axis=0),
        quantiles=quantiles)

    # mask 0 data
    field = np.ma.masked_where(field == 0, field)

    # display data
    if hist_type == 'instant':
        titl = pyart.graph.common.generate_title(
            hist_obj_aux, field_name, ind_sweep)
    else:
        titl = (
            '{:.1f}'.format(hist_obj_aux.fixed_angle['data'][0])+' Deg. ' +
            pyart.graph.common.generate_radar_time_begin(
                hist_obj).strftime('%Y-%m-%d') + '\n' +
            get_field_name(hist_obj.fields[field_name], field_name))
    label = 'Number of Points'
    labely = get_colobar_label(hist_obj_aux.fields[field_name], field_name)

    dpi = 72
    if 'dpi' in prdcfg['ppiImageConfig']:
        dpi = prdcfg['ppiImageConfig']['dpi']

    fig = plt.figure(figsize=[prdcfg['ppiImageConfig']['xsize'],
                              prdcfg['ppiImageConfig']['ysize']],
                     dpi=dpi)
    ax = fig.add_subplot(111)

    # define limits of field and color map
    cmap = pyart.config.get_field_colormap(field_name)
    step = hist_obj.range['data'][1]-hist_obj.range['data'][0]
    xmin = ang_min-ang_step/2.
    xmax = ang_max+ang_step/2.
    ymin = hist_obj.range['data'][0]-step/2.
    ymax = hist_obj.range['data'][-1]+step/2.

    cax = ax.imshow(
        np.ma.transpose(field), origin='lower', cmap=cmap, vmin=0.,
        vmax=np.max(field), extent=(xmin, xmax, ymin, ymax),
        aspect='auto', interpolation='none')
    ax.set_ylim(bottom=vmin, top=vmax)
    ax.autoscale(False)

    # plot reference
    ax.plot(ang, np.zeros(len(ang))+ref_value, 'k--')

    # plot quantiles
    ax.plot(ang, np.zeros(len(ang))+values_sweep[1], 'r')
    ax.plot(ang, np.zeros(len(ang))+values_sweep[0], 'r--')
    ax.plot(ang, np.zeros(len(ang))+values_sweep[2], 'r--')

    ax.plot(ang, az_percentile_ref, 'k')
    ax.plot(ang, az_percentile_low, 'k--')
    ax.plot(ang, az_percentile_high, 'k--')

    # ax.autoscale(enable=True, axis='both', tight=True)

    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_title(titl)

    cb = fig.colorbar(cax)
    cb.set_label(label)

    val_quant0_str = '--'
    if values_sweep[0] is not np.ma.masked:
        val_quant0_str = '{:.3f}'.format(values_sweep[0])
    val_quant1_str = '--'
    if values_sweep[1] is not np.ma.masked:
        val_quant1_str = '{:.3f}'.format(values_sweep[1])
    val_quant2_str = '--'
    if values_sweep[2] is not np.ma.masked:
        val_quant2_str = '{:.3f}'.format(values_sweep[2])

    metadata = (
        'npoints: '+str(np.ma.sum(field))+'\n' +
        str(quantiles[1])+' quant: '+val_quant1_str+'\n' +
        str(quantiles[0])+' quant: '+val_quant0_str+'\n' +
        str(quantiles[2])+' quant: '+val_quant2_str+'\n')

    ax.text(0.05, 0.05, metadata, horizontalalignment='left',
            verticalalignment='bottom', transform=ax.transAxes)

    # Turn on the grid
    ax.grid()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_scatter(bin_edges1, bin_edges2, hist_2d, field_name1, field_name2, fname_list,
                 prdcfg, metadata=None, lin_regr=None, lin_regr_slope1=None,
                 rad1_name='RADAR001', rad2_name='RADAR002'):
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
    metadata : str
        a string with metadata to write in the plot
    lin_regr : tupple with 2 values
        the coefficients for a linear regression
    lin_regr_slope1 : float
        the intercep point of a linear regression of slope 1
    rad1_name, rad2_name : str
        name of the radars which data is used

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    # mask 0 data
    hist_2d = np.ma.masked_where(hist_2d == 0, hist_2d)

    # display data
    titl = 'colocated radar gates'
    label = 'Number of Points'
    labelx = rad1_name+' '+field_name1
    labely = rad2_name+' '+field_name2

    dpi = 72
    if 'dpi' in prdcfg['ppiImageConfig']:
        dpi = prdcfg['ppiImageConfig']['dpi']

    fig = plt.figure(figsize=[prdcfg['ppiImageConfig']['xsize'],
                              prdcfg['ppiImageConfig']['ysize']],
                     dpi=dpi)
    ax = fig.add_subplot(111)

    cmap = pyart.config.get_field_colormap(field_name1)

    cax = ax.imshow(
        np.ma.transpose(hist_2d), origin='lower', cmap=cmap, vmin=0.,
        vmax=np.max(hist_2d),
        extent=(bin_edges1[0], bin_edges1[-1], bin_edges2[0], bin_edges2[-1]),
        aspect='auto', interpolation='none')
    ax.autoscale(False)

    # plot reference
    step1 = bin_edges1[1]-bin_edges1[0]
    bin_centers1 = bin_edges1[0:-1]+step1/2.

    step2 = bin_edges2[1]-bin_edges2[0]
    bin_centers2 = bin_edges2[0:-1]+step2/2.

    ax.plot(bin_centers1, bin_centers2, 'k--')

    # plot linear regression
    if lin_regr is not None:
        ax.plot(bin_centers1, lin_regr[0]*bin_centers1+lin_regr[1], 'r')
    if lin_regr_slope1 is not None:
        ax.plot(bin_centers1, bin_centers1+lin_regr_slope1, 'g')

    # ax.autoscale(enable=True, axis='both', tight=True)

    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_title(titl)

    cb = fig.colorbar(cax)
    cb.set_label(label)

    if metadata is not None:
        ax.text(0.05, 0.95, metadata, horizontalalignment='left',
                verticalalignment='top', transform=ax.transAxes)

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_quantiles(quant, value, fname_list, labelx='quantile', labely='value',
                   titl='quantile', vmin=None, vmax=None, dpi=72):
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
    ax.plot(quant, value, 'bx-')
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_ylim(bottom=vmin, top=vmax)
    ax.set_title(titl)

    # Turn on the grid
    ax.grid()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_histogram(bin_edges, values, fname_list, labelx='bins',
                   labely='Number of Samples', titl='histogram', dpi=72):
    """
    computes and plots histogram

    Parameters
    ----------
    bin_edges : array
        histogram bin edges
    values : array
        data values
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    titl : str
        The figure title
    dpi : int
        dots per inch

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    fig, ax = plt.subplots(figsize=[10, 6], dpi=dpi)
    ax.hist(values, bins=bin_edges)
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_title(titl)

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_histogram2(bin_centers, hist, fname_list, width=None, labelx='bins',
                    labely='Number of Samples', titl='histogram', dpi=72,
                    ax=None, fig=None, save_fig=True, color=None, alpha=None,
                    invert_xaxis=False):
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
        width = bin_centers[1]-bin_centers[0]
    ax.bar(bin_centers, hist, width=width, color=color, alpha=alpha)
    if invert_xaxis:
        ax.invert_xaxis()

    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_title(titl)

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list

    return (fig, ax)


def plot_antenna_pattern(antpattern, fname_list, labelx='Angle [Deg]',
                         linear=False, twoway=False, title='Antenna Pattern',
                         ymin=None, ymax=None, dpi=72):
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
    waystr = 'One-way '
    if twoway:
        waystr = 'Two-way '

    linstr = 'Att [dB]'
    if linear:
        linstr = 'Att [-]'
        mainbeam = (
            antpattern['angle'][antpattern['attenuation'] >= 10.**(-0.3)])
    else:
        mainbeam = antpattern['angle'][antpattern['attenuation'] >= -3.]
    beamwidth = np.abs(np.max(mainbeam)-np.min(mainbeam))

    labely = waystr+linstr

    fig, ax = plt.subplots(figsize=[10, 6], dpi=dpi)

    ax.plot(antpattern['angle'], antpattern['attenuation'])
    ax.set_title(title)
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_ylim(bottom=ymin, top=ymax)

    metadata = '3-dB beamwidth: '+'{:.2f}'.format(float(beamwidth))
    ax.text(0.05, 0.95, metadata, horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes)

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_selfconsistency(zdrkdp_table, fname_list, labelx='ZDR [dB]',
                         labely='KDP/Zh [(deg*m3)/(km*mm6)]',
                         title='Selfconsistency in rain', ymin=None,
                         ymax=None, dpi=72, save_fig=True, ax=None, fig=None):
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
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list

    return (fig, ax)


def plot_selfconsistency_instrument(zdr, kdp, zh, fname_list,
                                    bins_zdr_step=0.05, bins_zdr_min=0.,
                                    bins_zdr_max=6., bins_kdpzh_step=0.1,
                                    bins_kdpzh_min=-2., bins_kdpzh_max=20.,
                                    normalize=True, vmin=0., vmax=0.01,
                                    parametrization='None',
                                    zdr_kdpzh_dict=None,
                                    retrieve_relation=True,
                                    plot_theoretical=True, dpi=72):
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
        bins_zdr_min, bins_zdr_max+bins_zdr_step, bins_zdr_step)
    bins_zdr_edges = np.append(
        bins_zdr_centers-bins_zdr_step/2.,
        bins_zdr_centers[-1]+bins_zdr_step/2.)
    bins_kdpzh_centers = np.arange(
        bins_kdpzh_min, bins_kdpzh_max+bins_kdpzh_step, bins_kdpzh_step)
    bins_kdpzh_edges = np.append(
        bins_kdpzh_centers-bins_kdpzh_step/2.,
        bins_kdpzh_centers[-1]+bins_kdpzh_step/2.)

    labely = '10e5*KDP/ZH [(deg*m3)/(km*mm6)'

    kdpzh_th = None
    if plot_theoretical:
        if parametrization == 'None':
            if zdr_kdpzh_dict is not None and 'zdr_kdpzh' in zdr_kdpzh_dict:
                zdr_th = zdr_kdpzh_dict['zdr_kdpzh'][0][0]
                kdpzh_th = 1e5*zdr_kdpzh_dict['zdr_kdpzh'][0][1]
            else:
                warn('Unable to plot theoretical self-consistency curve. '
                     'Relationship not provided')
        elif parametrization == 'Gourley':
            if zdr_kdpzh_dict is None or 'freq_band' not in zdr_kdpzh_dict:
                warn('Unable to plot theoretical self-consistency curve. '
                     'Frequency band not provided')
            else:
                zdr_th = np.arange(0., 3.5+bins_zdr_step, bins_zdr_step)
                if zdr_kdpzh_dict['freq_band'] == 'S':
                    kdpzh_th = (
                        3.696-1.963*zdr_th+0.504*zdr_th*zdr_th -
                        0.051*zdr_th*zdr_th*zdr_th)
                elif zdr_kdpzh_dict['freq_band'] == 'C':
                    kdpzh_th = (
                        6.746-2.970*zdr_th+0.711*zdr_th*zdr_th -
                        0.079*zdr_th*zdr_th*zdr_th)
                elif zdr_kdpzh_dict['freq_band'] == 'X':
                    kdpzh_th = (
                        11.74-4.020*zdr_th-0.140*zdr_th*zdr_th +
                        0.130*zdr_th*zdr_th*zdr_th)
                else:
                    warn(
                        'Unable to plot theoretical self-consistency curve. '
                        'Unknown frequency band '+zdr_kdpzh_dict['freq_band'])
        elif parametrization == 'Wolfensberger':
            if zdr_kdpzh_dict is None or 'freq_band' not in zdr_kdpzh_dict:
                warn('Unable to plot theoretical self-consistency curve. '
                     'Frequency band not provided')
            else:
                zdr_th = np.arange(0., 3.5+bins_zdr_step, bins_zdr_step)
                if zdr_kdpzh_dict['freq_band'] == 'C':
                    kdpzh_th = (
                        3.199*np.ma.exp(-7.767e-1*zdr_th)-0.4436*zdr_th+3.464)
                else:
                    warn(
                        'Unable to plot theoretical self-consistency curve. '
                        'Unknown frequency band '+zdr_kdpzh_dict['freq_band'])
        elif parametrization == 'Louf':
            if zdr_kdpzh_dict is None or 'freq_band' not in zdr_kdpzh_dict:
                warn('Unable to plot theoretical self-consistency curve. '
                     'Frequency band not provided')
            else:
                zdr_th = np.arange(0.5, 3.5+bins_zdr_step, bins_zdr_step)
                if zdr_kdpzh_dict['freq_band'] == 'C':
                    kdpzh_th = (
                        6.607-4.577*zdr_th+1.577*zdr_th*zdr_th -
                        0.23*zdr_th*zdr_th*zdr_th)
                else:
                    warn(
                        'Unable to plot theoretical self-consistency curve. '
                        'Unknown frequency band '+zdr_kdpzh_dict['freq_band'])
        elif parametrization == 'Gorgucci':
            if zdr_kdpzh_dict is None or 'freq_band' not in zdr_kdpzh_dict:
                warn('Unable to plot theoretical self-consistency curve. '
                     'Frequency band not provided')
            else:
                zdr_th = np.arange(0., 3.5+bins_zdr_step, bins_zdr_step)
                if zdr_kdpzh_dict['freq_band'] == 'C':
                    zdr_lin = np.ma.power(10., 0.1*zdr_th)
                    kdpzh_th = 18.2*np.power(zdr_lin, -1.28)
                    labely = '10e5*KDP/ZH^0.95 [(deg*m3)/(km*mm6)'
                else:
                    warn(
                        'Unable to plot theoretical self-consistency curve. '
                        'Unknown frequency band '+zdr_kdpzh_dict['freq_band'])
        elif parametrization == 'Vaccarono':
            if zdr_kdpzh_dict is None or 'freq_band' not in zdr_kdpzh_dict:
                warn('Unable to plot theoretical self-consistency curve. '
                     'Frequency band not provided')
            else:
                zdr_th = np.arange(0., 3.5+bins_zdr_step, bins_zdr_step)
                if zdr_kdpzh_dict['freq_band'] == 'C':
                    zdr_lin = np.ma.power(10., 0.1*zdr_th)
                    kdpzh_th = 17.7*np.power(zdr_lin, -2.09)
                    labely = '10e5*KDP^0.85/ZH^0.91 [(deg*m3)/(km*mm6)'
                else:
                    warn(
                        'Unable to plot theoretical self-consistency curve. '
                        'Unknown frequency band '+zdr_kdpzh_dict['freq_band'])

    if kdpzh_th is not None and parametrization == 'Gorgucci':
        kdpzh = 1e5*kdp/np.ma.power(zh, 0.95)
    elif kdpzh_th is not None and parametrization == 'Vaccarono':
        kdpzh = 1e5*np.ma.power(kdp, 0.85)/np.ma.power(zh, 0.95)
    else:
        kdpzh = 1e5*kdp/zh

    if retrieve_relation:
        zdr_min = zdr.min()
        zdr_max = zdr.max()

        ind_min = np.where(bins_zdr_edges >= zdr_min)[0][0]
        ind_max = np.where(bins_zdr_edges < zdr_max)[0][-1]
        zdr_int = bins_zdr_edges[ind_min:ind_max+1]
        zdr_par = bins_zdr_centers[ind_min:ind_max]
        kdpzh_par = np.ma.masked_all(zdr_par.size)
        for ind, zdr_left in enumerate(zdr_int[:-1]):
            zdr_right = zdr_int[ind+1]
            ind_val = np.where(
                np.logical_and(zdr >= zdr_left, zdr < zdr_right))[0]
            if ind_val.size == 0:
                continue
            kdpzh_par[ind] = np.ma.median(kdpzh[ind_val])

    # prepare histogram
    zdr[zdr < bins_zdr_centers[0]] = bins_zdr_centers[0]
    zdr[zdr > bins_zdr_centers[-1]] = bins_zdr_centers[-1]

    kdpzh[kdpzh < bins_kdpzh_centers[0]] = bins_kdpzh_centers[0]
    kdpzh[kdpzh > bins_kdpzh_centers[-1]] = bins_kdpzh_centers[-1]

    hist2d, bins_zdr_edges, bins_kdpzh_edges = np.histogram2d(
        zdr, kdpzh, bins=[bins_zdr_edges, bins_kdpzh_edges])

    hist2d = np.ma.masked_equal(hist2d, 0)


    if normalize:
        hist2d = hist2d/np.expand_dims(np.ma.sum(hist2d, axis=-1), axis=1)
        clabel = 'Occurrence density\n(For each ZDR bin)'
    else:
        clabel = 'Number of gates'
        vmin = None
        vmax = None

    if not retrieve_relation and not plot_theoretical:
        fname_list = _plot_time_range(
            bins_zdr_edges, bins_kdpzh_edges, hist2d, None, fname_list,
            titl='self-consistency', xlabel='ZDR [dB]', ylabel=labely,
            clabel=clabel, vmin=vmin, vmax=vmax, save_fig=True, dpi=dpi)

        return fname_list

    fig, ax = _plot_time_range(
        bins_zdr_edges, bins_kdpzh_edges, hist2d, None, fname_list,
        titl='self-consistency', xlabel='ZDR [dB]', ylabel=labely,
        clabel=clabel, vmin=vmin, vmax=vmax, save_fig=False, dpi=dpi)

    ax.autoscale(False)

    if retrieve_relation:
        ax.plot(zdr_par, kdpzh_par, 'rx-', linewidth=3)

    if plot_theoretical and kdpzh_th is not None:
        ax.plot(zdr_th, kdpzh_th, 'r', linewidth=3)

    for fname in fname_list:
        fig.savefig(fname, dpi=72)
    plt.close(fig)

    return fname_list


def plot_selfconsistency_instrument2(zdr, kdp, zh, fname_list,
                                     bins_zh_step=0.5, bins_zh_min=-30.,
                                     bins_zh_max=80.,
                                     normalize=True, vmin=0., vmax=0.01,
                                     parametrization='None',
                                     zdr_kdpzh_dict=None, dpi=72):
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
    bins_zh_centers = np.arange(
        bins_zh_min, bins_zh_max+bins_zh_step, bins_zh_step)
    bins_zh_edges = np.append(
        bins_zh_centers-bins_zh_step/2.,
        bins_zh_centers[-1]+bins_zh_step/2.)

    if parametrization == 'None':
        if zdr_kdpzh_dict is not None and 'zdr_kdpzh' in zdr_kdpzh_dict:
            zdr_th_table = zdr_kdpzh_dict['zdr_kdpzh'][0][0]
            kdpzh_th_table = zdr_kdpzh_dict['zdr_kdpzh'][0][1]

            # sort values by ZDR
            ind_zdr_sorted = np.argsort(zdr)
            zh = zh[ind_zdr_sorted]
            kdp = kdp[ind_zdr_sorted]

            # get value of KDP/ZH as interpolation of look up table values
            kdpzh_th = np.interp(
                zdr[ind_zdr_sorted], zdr_th_table, kdpzh_th_table)

        else:
            warn('Unable to plot self-consistency curve. '
                 'Relationship not provided')
            return None
    elif parametrization == 'Gourley':
        if zdr_kdpzh_dict is None or 'freq_band' not in zdr_kdpzh_dict:
            warn('Unable to plot theoretical self-consistency curve. '
                 'Frequency band not provided')
            return None
        else:
            if zdr_kdpzh_dict['freq_band'] == 'S':
                kdpzh_th = 1e-5*(
                    3.696-1.963*zdr+0.504*zdr*zdr-0.051*zdr*zdr*zdr)
            elif zdr_kdpzh_dict['freq_band'] == 'C':
                kdpzh_th = 1e-5*(
                    6.746-2.970*zdr+0.711*zdr*zdr-0.079*zdr*zdr*zdr)
            elif zdr_kdpzh_dict['freq_band'] == 'X':
                kdpzh_th = 1e-5*(
                    11.74-4.020*zdr-0.140*zdr*zdr+0.130*zdr*zdr*zdr)
            else:
                warn(
                    'Unable to plot theoretical self-consistency curve. '
                    'Unknown frequency band '+zdr_kdpzh_dict['freq_band'])
                return None
    elif parametrization == 'Wolfensberger':
        if zdr_kdpzh_dict is None or 'freq_band' not in zdr_kdpzh_dict:
            warn('Unable to plot theoretical self-consistency curve. '
                 'Frequency band not provided')
            return None
        else:
            if zdr_kdpzh_dict['freq_band'] == 'C':
                kdpzh_th = 1e-5*(
                    3.199*np.ma.exp(-7.767e-1*zdr)-0.4436*zdr+3.464)
            else:
                warn(
                    'Unable to plot theoretical self-consistency curve. '
                    'Unknown frequency band '+zdr_kdpzh_dict['freq_band'])
                return None
    elif parametrization == 'Louf':
        if zdr_kdpzh_dict is None or 'freq_band' not in zdr_kdpzh_dict:
            warn('Unable to plot theoretical self-consistency curve. '
                 'Frequency band not provided')
            return None
        else:
            if zdr_kdpzh_dict['freq_band'] == 'C':
                kdpzh_th = 1e-5*(
                    6.607-4.577*zdr+1.577*zdr*zdr-0.23*zdr*zdr*zdr)
            else:
                warn(
                    'Unable to plot theoretical self-consistency curve. '
                    'Unknown frequency band '+zdr_kdpzh_dict['freq_band'])
                return None
    elif parametrization == 'Gorgucci':
        if zdr_kdpzh_dict is None or 'freq_band' not in zdr_kdpzh_dict:
            warn('Unable to plot theoretical self-consistency curve. '
                 'Frequency band not provided')
            return None
        else:
            if zdr_kdpzh_dict['freq_band'] == 'C':
                zdr_lin = np.ma.power(10., 0.1*zdr)
                kdpzh_th = 1e-5*(18.2*np.power(zdr_lin, -1.28))
            else:
                warn(
                    'Unable to plot theoretical self-consistency curve. '
                    'Unknown frequency band '+zdr_kdpzh_dict['freq_band'])
                return None
    elif parametrization == 'Vaccarono':
        if zdr_kdpzh_dict is None or 'freq_band' not in zdr_kdpzh_dict:
            warn('Unable to plot theoretical self-consistency curve. '
                 'Frequency band not provided')
            return None
        else:
            if zdr_kdpzh_dict['freq_band'] == 'C':
                zdr_lin = np.ma.power(10., 0.1*zdr)
                kdpzh_th = 1e-5*(17.7*np.power(zdr_lin, -2.09))
            else:
                warn(
                    'Unable to plot theoretical self-consistency curve. '
                    'Unknown frequency band '+zdr_kdpzh_dict['freq_band'])
                return None

    if parametrization == 'Gorgucci':
        zh_th = 10.*np.ma.log10(np.ma.power(kdp/kdpzh_th, 1/0.95))
    elif kdpzh_th is not None and parametrization == 'Vaccarono':
        zh_th = 10.*np.ma.log10(
            np.ma.power(np.ma.power(kdp, 0.85)/kdpzh_th, 1/0.95))
    else:
        zh_th = 10.*np.ma.log10(kdp/kdpzh_th)

    # prepare histogram
    zh[zh < bins_zh_centers[0]] = bins_zh_centers[0]
    zh[zh > bins_zh_centers[-1]] = bins_zh_centers[-1]

    zh_th[zh_th < bins_zh_centers[0]] = bins_zh_centers[0]
    zh_th[zh_th > bins_zh_centers[-1]] = bins_zh_centers[-1]

    hist2d, bins_zh_edges, _ = np.histogram2d(
        zh, zh_th, bins=[bins_zh_edges, bins_zh_edges])
    hist2d = np.ma.masked_equal(hist2d, 0)

    if normalize:
        hist2d = hist2d/np.expand_dims(np.ma.sum(hist2d, axis=-1), axis=1)
        clabel = 'Occurrence density\n(For each ZH bin)'
    else:
        clabel = 'Number of gates'
        vmin = None
        vmax = None

    fig, ax = _plot_time_range(
        bins_zh_edges, bins_zh_edges, hist2d, None, fname_list,
        titl='self-consistency', xlabel='ZH measured [dBZ]',
        ylabel='ZH selfcons [dBZ]',
        clabel=clabel, vmin=vmin, vmax=vmax, save_fig=False, dpi=dpi)

    ax.autoscale(False)
    ax.plot(bins_zh_centers, bins_zh_centers, 'rx-', linewidth=3)

    for fname in fname_list:
        fig.savefig(fname, dpi=72)
    plt.close(fig)

    return fname_list


def plot_scatter_comp(value1, value2, fname_list, labelx='Sensor 1',
                      labely='Sensor 2', titl='Scatter', axis=None,
                      metadata=None, dpi=72, ax=None, fig=None,
                      save_fig=True, point_format='bx'):
    """
    plots the scatter between two time series

    Parameters
    ----------
    value1 : float array
        values of the first time series
    value2 : float array
        values of the second time series
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    titl : str
        The figure title
    axis : str
        type of axis
    metadata : string
        a string containing metadata
    dpi : int
        dots per inch
    fig : Figure
        Figure to add the colorbar to. If none a new figure will be created
    ax : Axis
        Axis to plot on. if fig is None a new axis will be created
    save_fig : bool
        if true save the figure if false it does not close the plot and
        returns the handle to the figure
    point_format : str
        format of the scatter point

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    max_value = np.ma.max([np.max(value1), np.max(value2)])

    if fig is None:
        fig, ax = plt.subplots(figsize=[7, 7], dpi=dpi)

        ax.plot(value1, value2, point_format)
        ax.set_xlabel(labelx)
        ax.set_ylabel(labely)
        ax.set_title(titl)

        if axis == 'equal':
            ax.axis([0, max_value, 0, max_value])
            ax.plot([0, max_value], [0, max_value], 'k--')
            ax.set(adjustable='box-forced', aspect='equal')

        if metadata is not None:
            ax.text(0.05, 0.95, metadata, horizontalalignment='left',
                    verticalalignment='top', transform=ax.transAxes)

        # Make a tight layout
        fig.tight_layout()
    else:
        ax.autoscale(False)
        ax.plot(value1, value2, point_format)

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list
    return (fig, ax)


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
    azmin = prdcfg['sunhitsImageConfig']['azmin']
    azmax = prdcfg['sunhitsImageConfig']['azmax']
    elmin = prdcfg['sunhitsImageConfig']['elmin']
    elmax = prdcfg['sunhitsImageConfig']['elmax']

    field_dict = pyart.config.get_metadata(field_name)

    # display data
    dpi = 72
    if 'dpi' in prdcfg['sunhitsImageConfig']:
        dpi = prdcfg['sunhitsImageConfig']['dpi']

    fig = plt.figure(figsize=[prdcfg['sunhitsImageConfig']['xsize'],
                              prdcfg['sunhitsImageConfig']['ysize']],
                     dpi=dpi)
    ax = fig.add_subplot(111)
    cmap = pyart.config.get_field_colormap(field_name)
    vmin, vmax = pyart.config.get_field_limits(field_name)
    titl = (prdcfg['timeinfo'].strftime('%Y-%m-%d') + '\n' +
            get_field_name(field_dict, field_name))

    cax = ax.imshow(
        field, extent=(azmin, azmax, elmin, elmax),
        origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, interpolation='none')
    ax.set_xlabel('rad_az-sun_az (deg)')
    ax.set_ylabel('rad_el-sun_el (deg)')
    ax.set_title(titl)

    # plot the colorbar and set the label.
    label = get_colobar_label(field_dict, field_name)
    cb = fig.colorbar(cax)
    cb.set_label(label)

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def _plot_sunscan(rad_az, rad_el, rad_data, sun_hits, field_name, fname_list,
                  titl='AZ-EL sunscan plot', xlabel='Azimuth (deg)',
                  ylabel='Elevation (deg)', clabel=None, vmin=None, vmax=None,
                  figsize=[10, 8], save_fig=True, dpi=72):
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
            clabel = 'value'
        if vmin is None:
            vmin = np.ma.min(rad_data)
        if vmax is None:
            vmax = np.ma.max(rad_data)

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111)

    T, R = np.meshgrid(rad_az, rad_el)
    cax = ax.pcolormesh(
        T, R, np.ma.transpose(rad_data), cmap=cmap, vmin=vmin, vmax=vmax,
        norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(titl)

    # Add scatter points of detected sun hits location
    ax.scatter(sun_hits['rad_az'], sun_hits['rad_el'], s=7., c='red')

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
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list

def _plot_time_range(rad_time, rad_range, rad_data, field_name, fname_list,
                     titl='Time-Range plot',
                     xlabel='time (s from start time)', ylabel='range (Km)',
                     clabel=None, vmin=None, vmax=None, figsize=[10, 8],
                     save_fig=True, dpi=72):
    """
    plots a time-range plot

    Parameters
    ----------
    rad_time : 1D array
        array containing the x dimension (typically time)
    rad_range : 1D array
        array containing the y dimension (typically range)
    rad_data : 2D array
        array containing the data to plot
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
            clabel = 'value'
        if vmin is None:
            vmin = np.ma.min(rad_data)
        if vmax is None:
            vmax = np.ma.max(rad_data)

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111)

    T, R = np.meshgrid(rad_time, rad_range)
    cax = ax.pcolormesh(
        T, R, np.ma.transpose(rad_data), cmap=cmap, vmin=vmin, vmax=vmax,
        norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(titl)

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
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list

    return fig, ax
