"""
pyrad.graph.plots_vol
=====================

Functions to plot radar volume data

.. autosummary::
    :toctree: generated/

    plot_ray
    plot_ppi
    plot_ppi_map
    plot_rhi
    plot_bscope
    plot_time_range
    plot_fixed_rng
    plot_fixed_rng_span
    plot_fixed_rng_sun
    plot_cappi
    plot_traj
    plot_rhi_contour
    plot_ppi_contour
    plot_roi_contour
    plot_rhi_profile
    plot_along_coord
    plot_field_coverage

"""
from warnings import warn
from copy import deepcopy

import numpy as np


from netCDF4 import num2date

try:
    import cartopy
    from cartopy.io.img_tiles import Stamen
    _CARTOPY_AVAILABLE = True
except ImportError:
    _CARTOPY_AVAILABLE = False

try:
    import shapely
    _SHAPELY_AVAILABLE = True
except ImportError:
    warn('shapely not available')
    _SHAPELY_AVAILABLE = False

import matplotlib as mpl
mpl.use('Agg')

# Increase a bit font size
mpl.rcParams.update({'font.size': 16})
mpl.rcParams.update({'font.family':  "sans-serif"})

import matplotlib.pyplot as plt

import pyart

from .plots_aux import get_colobar_label, get_norm, generate_fixed_rng_title
from .plots_aux import generate_fixed_rng_span_title
from .plots_aux import generate_complex_range_Doppler_title
from .plots import plot_quantiles, plot_histogram, _plot_time_range
from .plots import _plot_sunscan

from ..util.radar_utils import compute_quantiles_sweep, find_ang_index
from ..util.radar_utils import compute_histogram_sweep


def plot_ray(radar, field_name, ind_ray, prdcfg, fname_list, titl=None,
             vmin=None, vmax=None, save_fig=True):
    """
    plots a ray

    Parameters
    ----------
    radar : Radar object
        object containing the radar data to plot
    field_name : str
        name of the radar field to plot
    ind_ray : int
        ray index to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    plot_type : str
        type of plot (PPI, QUANTILES or HISTOGRAM)
    titl : str
        Plot title
    vmin, vmax : float
        min and max values of the y axis
    save_fig : bool
        if true save the figure. If false it does not close the plot and
        returns the handle to the figure

    Returns
    -------
    fname_list : list of str or
    fig, ax : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    rng_km = radar.range['data']/1000.
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)

    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']
    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)

    if titl is None:
        titl = generate_complex_range_Doppler_title(
            radar, field_name, ind_ray)
    labely = get_colobar_label(radar.fields[field_name], field_name)

    ax = fig.add_subplot(111)

    ax.plot(rng_km, radar.fields[field_name]['data'][ind_ray, :], marker='x')

    ax.set_title(titl)
    ax.set_xlabel('Range (km)')
    ax.set_ylabel(labely)
    ax.set_ylim(bottom=vmin, top=vmax)
    ax.set_xlim([rng_km[0], rng_km[-1]])

    # Turn on the grid
    ax.grid()

    # Make a tight layout
    fig.tight_layout()

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list

    return (fig, ax)


def plot_ppi(radar, field_name, ind_el, prdcfg, fname_list, plot_type='PPI',
             titl=None, vmin=None, vmax=None, step=None, quantiles=None,
             save_fig=True):
    """
    plots a PPI

    Parameters
    ----------
    radar : Radar object
        object containing the radar data to plot
    field_name : str
        name of the radar field to plot
    ind_el : int
        sweep index to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    plot_type : str
        type of plot (PPI, QUANTILES or HISTOGRAM)
    titl : str
        Plot title
    vmin, vmax : float
        The minimum and maximum value. If None the scale is going to be
        obtained from the Py-ART config file.
    step : float
        step for histogram plotting
    quantiles : float array
        quantiles to plot
    save_fig : bool
        if true save the figure. If false it does not close the plot and
        returns the handle to the figure

    Returns
    -------
    fname_list : list of str or
    fig, ax : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    if plot_type == 'PPI':
        dpi = prdcfg['ppiImageConfig'].get('dpi', 72)

        norm = None
        ticks = None
        ticklabs = None
        if vmin is None or vmax is None:
            norm, ticks, ticklabs = get_norm(
                field_name, field_dict=radar.fields[field_name])
            vmin = None
            vmax = None

        xsize = prdcfg['ppiImageConfig']['xsize']
        ysize = prdcfg['ppiImageConfig']['ysize']
        fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)

        ax = fig.add_subplot(111, aspect='equal')

        display = pyart.graph.RadarDisplay(radar)
        display.plot_ppi(
            field_name, title=titl, sweep=ind_el, norm=norm, ticks=ticks,
            vmin=vmin, vmax=vmax, ticklabs=ticklabs, fig=fig, ax=ax)

        display.set_limits(
            ylim=[prdcfg['ppiImageConfig']['ymin'],
                  prdcfg['ppiImageConfig']['ymax']],
            xlim=[prdcfg['ppiImageConfig']['xmin'],
                  prdcfg['ppiImageConfig']['xmax']], ax=ax)
        if 'rngRing' in prdcfg['ppiImageConfig']:
            if prdcfg['ppiImageConfig']['rngRing'] > 0:
                display.plot_range_rings(np.arange(
                    0., radar.range['data'][-1]/1000.,
                    prdcfg['ppiImageConfig']['rngRing']), ax=ax)
        display.plot_cross_hair(5., ax=ax)

        # Turn on the grid
        ax.grid()

        # Make a tight layout
        fig.tight_layout()

        if save_fig:
            for fname in fname_list:
                fig.savefig(fname, dpi=dpi)
            plt.close(fig)

            return fname_list

        return (fig, ax)

    if plot_type == 'QUANTILES':
        quantiles, values = compute_quantiles_sweep(
            radar.fields[field_name]['data'],
            radar.sweep_start_ray_index['data'][ind_el],
            radar.sweep_end_ray_index['data'][ind_el], quantiles=quantiles)

        titl = pyart.graph.common.generate_title(radar, field_name, ind_el)
        labely = get_colobar_label(radar.fields[field_name], field_name)

        plot_quantiles(quantiles, values, fname_list, labelx='quantile',
                       labely=labely, titl=titl)

    elif plot_type == 'HISTOGRAM':
        bins, values = compute_histogram_sweep(
            radar.fields[field_name]['data'],
            radar.sweep_start_ray_index['data'][ind_el],
            radar.sweep_end_ray_index['data'][ind_el], field_name, step=step)

        titl = pyart.graph.common.generate_title(radar, field_name, ind_el)
        labelx = get_colobar_label(radar.fields[field_name], field_name)

        plot_histogram(bins, values, fname_list, labelx=labelx,
                       labely='Number of Samples', titl=titl)
    else:
        warn('Unknown plot type '+plot_type)

    return fname_list


def plot_ppi_map(radar, field_name, ind_el, prdcfg, fname_list,
                 save_fig=True):
    """
    plots a PPI on a geographic map

    Parameters
    ----------
    radar : Radar object
        object containing the radar data to plot
    field_name : str
        name of the radar field to plot
    ind_el : int
        sweep index to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    save_fig : bool
        if true save the figure. If false it does not close the plot and
        returns the handle to the figure

    Returns
    -------
    fname_list : list of str or
    fig, ax, display : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)

    norm, ticks, ticklabs = get_norm(
        field_name, field_dict=radar.fields[field_name])

    xsize = prdcfg['ppiMapImageConfig']['xsize']
    ysize = prdcfg['ppiMapImageConfig']['ysize']
    lonstep = prdcfg['ppiMapImageConfig'].get('lonstep', 0.5)
    latstep = prdcfg['ppiMapImageConfig'].get('latstep', 0.5)
    min_lon = prdcfg['ppiMapImageConfig'].get('lonmin', 2.5)
    max_lon = prdcfg['ppiMapImageConfig'].get('lonmax', 12.5)
    min_lat = prdcfg['ppiMapImageConfig'].get('latmin', 43.5)
    max_lat = prdcfg['ppiMapImageConfig'].get('latmax', 49.5)
    resolution = prdcfg['ppiMapImageConfig'].get('mapres', '110m')
    if resolution not in ('110m', '50m', '10m'):
        warn('Unknown map resolution: '+resolution)
        resolution = '110m'
    background_zoom = prdcfg['ppiMapImageConfig'].get('background_zoom', 8)

    lon_lines = np.arange(np.floor(min_lon), np.ceil(max_lon)+1, lonstep)
    lat_lines = np.arange(np.floor(min_lat), np.ceil(max_lat)+1, latstep)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    ax = fig.add_subplot(111)

    display_map = pyart.graph.RadarMapDisplay(radar)
    display_map.plot_ppi_map(
        field_name, sweep=ind_el, norm=norm, ticks=ticks, ticklabs=ticklabs,
        min_lon=min_lon, max_lon=max_lon, min_lat=min_lat, max_lat=max_lat,
        resolution=resolution, background_zoom=background_zoom,
        lat_lines=lat_lines, lon_lines=lon_lines,
        maps_list=prdcfg['ppiMapImageConfig']['maps'], ax=ax, fig=fig,
        colorbar_flag=True, alpha=1)

    ax = display_map.ax

    if 'rngRing' in prdcfg['ppiMapImageConfig']:
        if prdcfg['ppiMapImageConfig']['rngRing'] > 0:
            rng_rings = np.arange(
                0., radar.range['data'][-1]/1000.,
                prdcfg['ppiMapImageConfig']['rngRing'])
            for rng_ring in rng_rings:
                display_map.plot_range_ring(rng_ring, ax=ax)

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list

    return (fig, ax, display_map)


def plot_rhi(radar, field_name, ind_az, prdcfg, fname_list, plot_type='RHI',
             titl=None, vmin=None, vmax=None, step=None, quantiles=None,
             save_fig=True):
    """
    plots an RHI

    Parameters
    ----------
    radar : Radar object
        object containing the radar data to plot
    field_name : str
        name of the radar field to plot
    ind_az : int
        sweep index to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    plot_type : str
        type of plot (PPI, QUANTILES or HISTOGRAM)
    titl : str
        Plot title
    vmin, vmax : float
        The minimum and maximum value. If None the scale is going to be
        obtained from the Py-ART config file.
    step : float
        step for histogram plotting
    quantiles : float array
        quantiles to plot
    save_fig : bool
        if true save the figure. If false it does not close the plot and
        returns the handle to the figure

    Returns
    -------
    fname_list : list of str or
    fig, ax : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    if plot_type == 'RHI':
        dpi = prdcfg['ppiImageConfig'].get('dpi', 72)

        norm = None
        ticks = None
        ticklabs = None
        if vmin is None or vmax is None:
            norm, ticks, ticklabs = get_norm(
                field_name, field_dict=radar.fields[field_name])
            vmin = None
            vmax = None

        xsize = prdcfg['rhiImageConfig']['xsize']
        ysize = prdcfg['rhiImageConfig']['ysize']
        fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
        ax = fig.add_subplot(111, aspect='equal')
        display = pyart.graph.RadarDisplay(radar)
        display.plot_rhi(
            field_name, title=titl, sweep=ind_az, norm=norm, ticks=ticks,
            ticklabs=ticklabs, vmin=vmin, vmax=vmax,
            colorbar_orient='horizontal', reverse_xaxis=False, fig=fig, ax=ax)
        display.set_limits(
            ylim=[prdcfg['rhiImageConfig']['ymin'],
                  prdcfg['rhiImageConfig']['ymax']],
            xlim=[prdcfg['rhiImageConfig']['xmin'],
                  prdcfg['rhiImageConfig']['xmax']],
            ax=ax)
        display.plot_cross_hair(5., ax=ax)

        # Turn on the grid
        ax.grid()

        # Make a tight layout
        fig.tight_layout()

        if save_fig:
            for fname in fname_list:
                fig.savefig(fname, dpi=dpi)
            plt.close(fig)

            return fname_list

        return (fig, ax)

    if plot_type == 'QUANTILES':
        quantiles, values = compute_quantiles_sweep(
            radar.fields[field_name]['data'],
            radar.sweep_start_ray_index['data'][ind_az],
            radar.sweep_end_ray_index['data'][ind_az], quantiles=quantiles)

        if titl is None:
            titl = pyart.graph.common.generate_title(
                radar, field_name, ind_az)
        labely = get_colobar_label(radar.fields[field_name], field_name)

        plot_quantiles(quantiles, values, fname_list, labelx='quantile',
                       labely=labely, titl=titl)

    elif plot_type == 'HISTOGRAM':
        bins, values = compute_histogram_sweep(
            radar.fields[field_name]['data'],
            radar.sweep_start_ray_index['data'][ind_az],
            radar.sweep_end_ray_index['data'][ind_az], field_name, step=step)

        if titl is None:
            titl = pyart.graph.common.generate_title(
                radar, field_name, ind_az)
        labelx = get_colobar_label(radar.fields[field_name], field_name)

        plot_histogram(bins, values, fname_list, labelx=labelx,
                       labely='Number of Samples', titl=titl)
    else:
        warn('Unknown plot type '+plot_type)

    return fname_list


def plot_bscope(radar, field_name, ind_sweep, prdcfg, fname_list,
                vmin=None, vmax=None, ray_dim='ang', xaxis_rng=True):
    """
    plots a B-Scope (angle-range representation)

    Parameters
    ----------
    radar : Radar object
        object containing the radar data to plot
    field_name : str
        name of the radar field to plot
    ind_sweep : int
        sweep index to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    vmin, vmax : float
        Min and max values of the colorbar
    ray_dim : str
        the ray dimension. Can be 'ang' or 'time'
    xaxis : bool
        if true the range will be in the x-axis. Otherwise it will be in the
        y-axis.

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    norm = None
    ticks = None
    ticklabs = None
    if vmin is None or vmax is None:
        norm, ticks, ticklabs = get_norm(
            field_name, field_dict=radar.fields[field_name])

        if norm is None:  # if norm is set do not override with vmin/vmax
            vmin, vmax = pyart.config.get_field_limits(field_name)

    radar_aux = radar.extract_sweeps([ind_sweep])
    if ray_dim == 'ang':
        if radar_aux.scan_type == 'ppi':
            ray = np.sort(radar_aux.azimuth['data'])
            ind_ray = np.argsort(radar_aux.azimuth['data'])
            field = radar_aux.fields[field_name]['data'][ind_ray, :]
            ray_label = 'azimuth angle (degrees)'
        elif radar_aux.scan_type == 'rhi':
            ray = np.sort(radar_aux.elevation['data'])
            ind_ray = np.argsort(radar_aux.elevation['data'])
            field = radar_aux.fields[field_name]['data'][ind_ray, :]
            ray_label = 'elevation angle (degrees)'
        else:
            field = radar_aux.fields[field_name]['data']
            ray = np.array(range(radar_aux.nrays))
            ray_label = 'ray number'
    else:
        ray = np.sort(radar_aux.time['data'])
        start_time = ray[0]
        ray -= start_time
        ind_ray = np.argsort(radar_aux.time['data'])
        field = radar_aux.fields[field_name]['data'][ind_ray, :]
        sweep_start_time = num2date(
            start_time, radar_aux.time['units'], radar_aux.time['calendar'])
        ray_label = (
            'time [s from ' +
            sweep_start_time.strftime('%Y-%m-%d %H:%M:%S')+' UTC]')

    # display data
    titl = pyart.graph.common.generate_title(radar_aux, field_name, 0)
    label = get_colobar_label(radar_aux.fields[field_name], field_name)

    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)

    fig = plt.figure(figsize=[prdcfg['ppiImageConfig']['xsize'],
                              prdcfg['ppiImageConfig']['ysize']],
                     dpi=dpi)
    ax = fig.add_subplot(111)
    if radar_aux.ngates == 1:
        ax.plot(ray, field, 'bx', figure=fig)
        ax.set_xlabel(ray_label)
        ax.set_ylabel(label)
        ax.set_title(titl)
    else:
        cmap = pyart.config.get_field_colormap(field_name)

        rng_aux = radar_aux.range['data']/1000.
        rng_res = rng_aux[1]-rng_aux[0]
        rng_aux = np.append(rng_aux-rng_res/2., rng_aux[-1]+rng_res/2.)
        rng_label = 'Range (km)'

        ray_res = np.ma.median(ray[1:]-ray[:-1])
        ray_aux = np.append(ray-ray_res/2, ray[-1]+ray_res/2)

        if xaxis_rng:
            cax = ax.pcolormesh(
                rng_aux, ray_aux, field, cmap=cmap, vmin=vmin, vmax=vmax,
                norm=norm)
            ax.set_xlabel(rng_label)
            ax.set_ylabel(ray_label)
        else:
            cax = ax.pcolormesh(
                ray_aux, rng_aux, np.ma.transpose(field), cmap=cmap,
                vmin=vmin, vmax=vmax, norm=norm)
            ax.set_xlabel(ray_label)
            ax.set_ylabel(rng_label)
        ax.set_title(titl)

        cb = fig.colorbar(cax)
        if ticks is not None:
            cb.set_ticks(ticks)
        if ticklabs:
            cb.set_ticklabels(ticklabs)
        cb.set_label(label)

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_time_range(radar, field_name, ind_sweep, prdcfg, fname_list,
                    vmin=None, vmax=None, ylabel='range (Km)'):
    """
    plots a time-range plot

    Parameters
    ----------
    radar : Radar object
        object containing the radar data to plot
    field_name : str
        name of the radar field to plot
    ind_sweep : int
        sweep index to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    vmin, vmax : float
        Min and max values of the colorbar
    ylabel : str
        The y-axis label

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    radar_aux = radar.extract_sweeps([ind_sweep])
    field = radar_aux.fields[field_name]['data']

    # display data
    titl = pyart.graph.common.generate_title(radar_aux, field_name, ind_sweep)
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig'].get('xsize', 10)
    ysize = prdcfg['ppiImageConfig'].get('ysize', 8)

    rng_aux = deepcopy(radar_aux.range['data'])
    if ylabel == 'range (km)':
        rng_aux /= 1000.
    rng_res = rng_aux[1]-rng_aux[0]
    rng_aux = np.append(rng_aux-rng_res/2., rng_aux[-1]+rng_res/2.)

    time_res = np.mean(radar_aux.time['data'][1:]-radar_aux.time['data'][0:-1])
    time_aux = np.append(
        radar_aux.time['data'], radar_aux.time['data'][-1]+time_res)
    return _plot_time_range(
        time_aux, rng_aux, field, field_name, fname_list, titl=titl,
        ylabel=ylabel, vmin=vmin, vmax=vmax, figsize=[xsize, ysize], dpi=dpi)


def plot_fixed_rng(radar, field_name, prdcfg, fname_list, azi_res=None,
                   ele_res=None, ang_tol=1., vmin=None, vmax=None):
    """
    plots a fixed range plot

    Parameters
    ----------
    radar : radar object
        The radar object containing the fixed range data
    field_name : str
        The name of the field to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    azi_res, ele_res : float
        The nominal azimuth and elevation angle resolution [deg]
    ang_tol : float
        The tolerance between the nominal and the actual radar angle
    vmin, vmax : float
        Min and Max values of the color scale. If None it is going to be taken
        from the Py-ART config files

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    # Get radar azimuth angles within limits taking as reference
    # the first elevation angle
    fixed_rng = radar.range['data'][0]

    if radar.scan_type == 'ppi':
        ele_vec = np.sort(radar.fixed_angle['data'])
        azi_vec = np.sort(
            radar.azimuth['data'][radar.sweep_start_ray_index['data'][0]:
                                  radar.sweep_end_ray_index['data'][0]+1])
    else:
        ele_vec = np.sort(
            radar.elevation['data'][radar.sweep_start_ray_index['data'][0]:
                                    radar.sweep_end_ray_index['data'][0]+1])
        azi_vec = np.sort(radar.fixed_angle['data'])

    # put data in a regular 2D grid
    field_2D = np.ma.masked_all((azi_vec.size, ele_vec.size))
    sweep_start_inds = radar.sweep_start_ray_index['data']
    sweep_end_inds = radar.sweep_end_ray_index['data']

    if radar.scan_type == 'ppi':
        for j, ele in enumerate(ele_vec):
            field_1D = radar.fields[field_name]['data'][
                sweep_start_inds[j]:sweep_end_inds[j]+1]
            azi_1D = radar.azimuth['data'][
                sweep_start_inds[j]:sweep_end_inds[j]+1]

            for i, azi in enumerate(azi_vec):
                ind = find_ang_index(azi_1D, azi, ang_tol=ang_tol)

                if ind is None:
                    continue
                try:
                    field_2D[i, j] = field_1D[ind]
                except ValueError:
                    field_2D[i, j] = field_1D[ind][0]
    else:
        for i, azi in enumerate(azi_vec):
            field_1D = radar.fields[field_name]['data'][
                sweep_start_inds[i]:sweep_end_inds[i]+1]
            ele_1D = radar.elevation['data'][
                sweep_start_inds[i]:sweep_end_inds[i]+1]

            for j, ele in enumerate(ele_vec):
                ind = find_ang_index(ele_1D, ele, ang_tol=ang_tol)
                if ind is None:
                    continue
                field_2D[i, j] = field_1D[ind]

    # get limits of angle bins
    if radar.scan_type == 'ppi':
        if azi_res is None:
            azi_res = np.median(azi_vec[1:]-azi_vec[0:-1])
            if radar.ray_angle_res is not None:
                azi_res = np.min(
                    [radar.ray_angle_res['data'][0], azi_res])

        azi_vec = np.append(azi_vec-azi_res/2., azi_vec[-1]+azi_res/2.)

        if ele_res is None:
            ele_res = np.median(ele_vec[1:]-ele_vec[0:-1])
            if radar.instrument_parameters is not None:
                if 'radar_beam_width_h' in radar.instrument_parameters:
                    bwidth = radar.instrument_parameters[
                        'radar_beam_width_h']['data'][0]
                    ele_res = np.min([bwidth, ele_res])
                elif 'radar_beam_width_v' in radar.instrument_parameters:
                    bwidth = radar.instrument_parameters[
                        'radar_beam_width_v']['data'][0]
                    ele_res = np.min([bwidth, ele_res])

        ele_vec = np.append(ele_vec-ele_res/2., ele_vec[-1]+ele_res/2.)
    else:
        if ele_res is None:
            ele_res = np.median(ele_vec[1:]-ele_vec[0:-1])
            if radar.ray_angle_res is not None:
                ele_res = np.min(
                    [radar.ray_angle_res['data'][0], ele_res])

        ele_vec = np.append(ele_vec-ele_res/2., ele_vec[-1]+ele_res/2.)

        if azi_res is None:
            azi_res = np.median(azi_vec[1:]-azi_vec[0:-1])
            if radar.instrument_parameters is not None:
                if 'radar_beam_width_h' in radar.instrument_parameters:
                    bwidth = radar.instrument_parameters[
                        'radar_beam_width_h']['data'][0]
                    azi_res = np.min([bwidth, azi_res])
                elif 'radar_beam_width_v' in radar.instrument_parameters:
                    bwidth = radar.instrument_parameters[
                        'radar_beam_width_v']['data'][0]
                    azi_res = np.min([bwidth, azi_res])

        azi_vec = np.append(azi_vec-azi_res/2., azi_vec[-1]+azi_res/2.)

    titl = generate_fixed_rng_title(radar, field_name, fixed_rng)

    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig'].get('xsize', 10)
    ysize = prdcfg['ppiImageConfig'].get('ysize', 8)

    return _plot_time_range(
        azi_vec, ele_vec, field_2D, field_name, fname_list, titl=titl,
        xlabel='azimuth (deg)', ylabel='elevation (deg)',
        figsize=[xsize, ysize], vmin=vmin, vmax=vmax, dpi=dpi)


def plot_fixed_rng_span(radar, field_name, prdcfg, fname_list, azi_res=None,
                        ele_res=None, ang_tol=1., stat='max'):
    """
    plots a fixed range plot

    Parameters
    ----------
    radar : radar object
        The radar object containing the fixed range data
    field_name : str
        The name of the field to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    azi_res, ele_res : float
        The nominal azimuth and elevation angle resolution [deg]
    ang_tol : float
        The tolerance between the nominal and the actual radar angle

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    # Get radar azimuth angles within limits taking as reference
    # the first elevation angle
    if radar.scan_type == 'ppi':
        ele_vec = np.sort(radar.fixed_angle['data'])
        azi_vec = np.sort(
            radar.azimuth['data'][radar.sweep_start_ray_index['data'][0]:
                                  radar.sweep_end_ray_index['data'][0]+1])
    else:
        ele_vec = np.sort(
            radar.elevation['data'][radar.sweep_start_ray_index['data'][0]:
                                    radar.sweep_end_ray_index['data'][0]+1])
        azi_vec = np.sort(radar.fixed_angle['data'])

    # put data in a regular 2D grid
    field_2D = np.ma.masked_all((azi_vec.size, ele_vec.size))
    rng_2D = np.ma.masked_all((azi_vec.size, ele_vec.size))
    sweep_start_inds = radar.sweep_start_ray_index['data']
    sweep_end_inds = radar.sweep_end_ray_index['data']

    if radar.scan_type == 'ppi':
        for j, ele in enumerate(ele_vec):
            field = radar.fields[field_name]['data'][
                sweep_start_inds[j]:sweep_end_inds[j]+1, :]

            if stat == 'max':
                field_1D = np.ma.max(field, axis=-1)
                ind = np.ma.argmax(field, axis=-1)
                rng_1D = radar.range['data'][ind]
            elif stat == 'min':
                field_1D = np.ma.min(field, axis=-1)
                ind = np.ma.argmin(field, axis=-1)
                rng_1D = radar.range['data'][ind]
            elif stat == 'mean':
                field_1D = np.ma.mean(field, axis=-1)
                mid_rng = radar.range['data'][int(radar.ngates/2)]
                rng_1D = np.ma.zeros(np.shape(field_1D))+mid_rng
            elif stat == 'median':
                field_1D = np.ma.median(field, axis=-1)
                mid_rng = radar.range['data'][int(radar.ngates/2)]
                rng_1D = np.ma.zeros(np.shape(field_1D))+mid_rng

            azi_1D = radar.azimuth['data'][
                sweep_start_inds[j]:sweep_end_inds[j]+1]

            for i, azi in enumerate(azi_vec):
                ind = find_ang_index(azi_1D, azi, ang_tol=ang_tol)
                if ind is None:
                    continue
                field_2D[i, j] = field_1D[ind]
                rng_2D[i, j] = rng_1D[ind]
    else:
        for i, azi in enumerate(azi_vec):
            field = radar.fields[field_name]['data'][
                sweep_start_inds[i]:sweep_end_inds[i]+1, :]

            if stat == 'max':
                field_1D = np.ma.max(field, axis=-1)
                ind = np.ma.argmax(field, axis=-1)
                rng_1D = radar.range['data'][ind]
            elif stat == 'min':
                field_1D = np.ma.min(field, axis=-1)
                ind = np.ma.argmin(field, axis=-1)
                rng_1D = radar.range['data'][ind]
            elif stat == 'mean':
                field_1D = np.ma.mean(field, axis=-1)
                mid_rng = radar.range['data'][int(radar.ngates/2)]
                rng_1D = np.ma.zeros(np.shape(field_1D))+mid_rng
            elif stat == 'median':
                field_1D = np.ma.median(field, axis=-1)
                mid_rng = radar.range['data'][int(radar.ngates/2)]
                rng_1D = np.ma.zeros(np.shape(field_1D))+mid_rng

            ele_1D = radar.elevation['data'][
                sweep_start_inds[i]:sweep_end_inds[i]+1]

            for j, ele in enumerate(ele_vec):
                ind = find_ang_index(ele_1D, ele, ang_tol=ang_tol)
                if ind is None:
                    continue
                field_2D[i, j] = field_1D[ind]
                rng_2D[i, j] = rng_1D[ind]

    # get limits of angle bins
    if radar.scan_type == 'ppi':
        if azi_res is None:
            azi_res = np.median(azi_vec[1:]-azi_vec[0:-1])
            if radar.ray_angle_res is not None:
                azi_res = np.min(
                    [radar.ray_angle_res['data'][0], azi_res])

        azi_vec = np.append(azi_vec-azi_res/2., azi_vec[-1]+azi_res/2.)

        if ele_res is None:
            ele_res = np.median(ele_vec[1:]-ele_vec[0:-1])
            if radar.instrument_parameters is not None:
                if 'radar_beam_width_h' in radar.instrument_parameters:
                    bwidth = radar.instrument_parameters[
                        'radar_beam_width_h']['data'][0]
                    ele_res = np.min([bwidth, ele_res])
                elif 'radar_beam_width_v' in radar.instrument_parameters:
                    bwidth = radar.instrument_parameters[
                        'radar_beam_width_v']['data'][0]
                    ele_res = np.min([bwidth, ele_res])

        ele_vec = np.append(ele_vec-ele_res/2., ele_vec[-1]+ele_res/2.)
    else:
        if ele_res is None:
            ele_res = np.median(ele_vec[1:]-ele_vec[0:-1])
            if radar.ray_angle_res is not None:
                ele_res = np.min(
                    [radar.ray_angle_res['data'][0], ele_res])

        ele_vec = np.append(ele_vec-ele_res/2., ele_vec[-1]+ele_res/2.)

        if azi_res is None:
            azi_res = np.median(azi_vec[1:]-azi_vec[0:-1])
            if radar.instrument_parameters is not None:
                if 'radar_beam_width_h' in radar.instrument_parameters:
                    bwidth = radar.instrument_parameters[
                        'radar_beam_width_h']['data'][0]
                    azi_res = np.min([bwidth, azi_res])
                elif 'radar_beam_width_v' in radar.instrument_parameters:
                    bwidth = radar.instrument_parameters[
                        'radar_beam_width_v']['data'][0]
                    azi_res = np.min([bwidth, azi_res])

        azi_vec = np.append(azi_vec-azi_res/2., azi_vec[-1]+azi_res/2.)

    titl = generate_fixed_rng_span_title(radar, field_name, stat)

    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig'].get('xsize', 10)
    ysize = prdcfg['ppiImageConfig'].get('ysize', 8)

    fname_rng_list = []
    for fname in fname_list:
        fname_rng_list.append(
            fname.rsplit('.', 1)[0]+'_RNG.'+fname.rsplit('.', 1)[1])

    _plot_time_range(
        azi_vec, ele_vec, rng_2D, 'radar_range', fname_rng_list, titl=titl,
        xlabel='azimuth (deg)', ylabel='elevation (deg)',
        figsize=[xsize, ysize], dpi=dpi)

    return _plot_time_range(
        azi_vec, ele_vec, field_2D, field_name, fname_list, titl=titl,
        xlabel='azimuth (deg)', ylabel='elevation (deg)',
        figsize=[xsize, ysize], dpi=dpi)


def plot_fixed_rng_sun(radar, field_name, sun_hits, prdcfg, fname_list,
                       azi_res=None, ele_res=None, ang_tol=1., vmin=None,
                       vmax=None):
    """
    plots a fixed range plot

    Parameters
    ----------
    radar : radar object
        The radar object containing the fixed range data
    field_name : str
        The name of the field to plot
    sun_hits: dict
        dictionary containing the sun hits data
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    azi_res, ele_res : float
        The nominal azimuth and elevation angle resolution [deg]
    ang_tol : float
        The tolerance between the nominal and the actual radar angle
    vmin, vmax : float
        Min and Max values of the color scale. If None it is going to be taken
        from the Py-ART config files

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    # Get radar azimuth angles within limits taking as reference
    # the first elevation angle
    fixed_rng = radar.range['data'][0]

    if radar.scan_type == 'ppi':
        ele_vec = np.sort(radar.fixed_angle['data'])
        azi_vec = np.sort(
            radar.azimuth['data'][radar.sweep_start_ray_index['data'][0]:
                                  radar.sweep_end_ray_index['data'][0]+1])
    else:
        ele_vec = np.sort(
            radar.elevation['data'][radar.sweep_start_ray_index['data'][0]:
                                    radar.sweep_end_ray_index['data'][0]+1])
        azi_vec = np.sort(radar.fixed_angle['data'])

    # put data in a regular 2D grid
    field_2D = np.ma.masked_all((azi_vec.size, ele_vec.size))
    sweep_start_inds = radar.sweep_start_ray_index['data']
    sweep_end_inds = radar.sweep_end_ray_index['data']

    if radar.scan_type == 'ppi':
        for j, ele in enumerate(ele_vec):
            field_1D = radar.fields[field_name]['data'][
                sweep_start_inds[j]:sweep_end_inds[j]+1]
            azi_1D = radar.azimuth['data'][
                sweep_start_inds[j]:sweep_end_inds[j]+1]

            for i, azi in enumerate(azi_vec):
                ind = find_ang_index(azi_1D, azi, ang_tol=ang_tol)
                # print('IND: ',ind)
                if ind is None:
                    continue
                # print('FIELD_1D: ',field_1D[ind])
                field_2D[i, j] = field_1D[ind][0]
    else:
        for i, azi in enumerate(azi_vec):
            field_1D = radar.fields[field_name]['data'][
                sweep_start_inds[i]:sweep_end_inds[i]+1]
            ele_1D = radar.elevation['data'][
                sweep_start_inds[i]:sweep_end_inds[i]+1]

            for j, ele in enumerate(ele_vec):
                ind = find_ang_index(ele_1D, ele, ang_tol=ang_tol)
                if ind is None:
                    continue
                field_2D[i, j] = field_1D[ind]

    # get limits of angle bins
    if radar.scan_type == 'ppi':
        if azi_res is None:
            azi_res = np.median(azi_vec[1:]-azi_vec[0:-1])
            if radar.ray_angle_res is not None:
                azi_res = np.min(
                    [radar.ray_angle_res['data'][0], azi_res])

        azi_vec = np.append(azi_vec-azi_res/2., azi_vec[-1]+azi_res/2.)

        if ele_res is None:
            ele_res = np.median(ele_vec[1:]-ele_vec[0:-1])
            if radar.instrument_parameters is not None:
                if 'radar_beam_width_h' in radar.instrument_parameters:
                    bwidth = radar.instrument_parameters[
                        'radar_beam_width_h']['data'][0]
                    ele_res = np.min([bwidth, ele_res])
                elif 'radar_beam_width_v' in radar.instrument_parameters:
                    bwidth = radar.instrument_parameters[
                        'radar_beam_width_v']['data'][0]
                    ele_res = np.min([bwidth, ele_res])

        ele_vec = np.append(ele_vec-ele_res/2., ele_vec[-1]+ele_res/2.)
    else:
        if ele_res is None:
            ele_res = np.median(ele_vec[1:]-ele_vec[0:-1])
            if radar.ray_angle_res is not None:
                ele_res = np.min(
                    [radar.ray_angle_res['data'][0], ele_res])

        ele_vec = np.append(ele_vec-ele_res/2., ele_vec[-1]+ele_res/2.)

        if azi_res is None:
            azi_res = np.median(azi_vec[1:]-azi_vec[0:-1])
            if radar.instrument_parameters is not None:
                if 'radar_beam_width_h' in radar.instrument_parameters:
                    bwidth = radar.instrument_parameters[
                        'radar_beam_width_h']['data'][0]
                    azi_res = np.min([bwidth, azi_res])
                elif 'radar_beam_width_v' in radar.instrument_parameters:
                    bwidth = radar.instrument_parameters[
                        'radar_beam_width_v']['data'][0]
                    azi_res = np.min([bwidth, azi_res])

        azi_vec = np.append(azi_vec-azi_res/2., azi_vec[-1]+azi_res/2.)

    titl = generate_fixed_rng_title(radar, field_name, fixed_rng)

    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig'].get('xsize', 10)
    ysize = prdcfg['ppiImageConfig'].get('ysize', 8)

    return _plot_sunscan(
        azi_vec, ele_vec, field_2D, sun_hits, field_name, fname_list,
        titl=titl, xlabel='azimuth (deg)', ylabel='elevation (deg)',
        figsize=[xsize, ysize], vmin=vmin, vmax=vmax, dpi=dpi)


def plot_cappi(radar, field_name, altitude, prdcfg, fname_list,
               beamwidth=1., beam_spacing=1., save_fig=True):
    """
    plots a Constant Altitude Plan Position Indicator CAPPI

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
    beamwidth : float
        The radar beamwidth
    beam_spacing : float
        the ray angle resolution
    save_fig : bool
        if true save the figure. If false it does not close the plot and
        returns the handle to the figure

    Returns
    -------
    fname_list : list of str or
    fig, ax : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    norm, ticks, ticklabs = get_norm(
        field_name, field_dict=radar.fields[field_name])

    xmin = prdcfg['ppiImageConfig']['xmin']
    xmax = prdcfg['ppiImageConfig']['xmax']
    ymin = prdcfg['ppiImageConfig']['ymin']
    ymax = prdcfg['ppiImageConfig']['ymax']

    wfunc = prdcfg.get('wfunc', 'NEAREST')
    cappi_res = prdcfg.get('res', 500.)

    # number of grid points in cappi
    ny = int((ymax-ymin)*1000./cappi_res)+1
    nx = int((xmax-xmin)*1000./cappi_res)+1

    # parameters to determine the gates to use for each grid point
    if (radar.instrument_parameters is not None and
            'radar_beam_width_h' in radar.instrument_parameters):
        beamwidth = radar.instrument_parameters[
            'radar_beam_width_h']['data'][0]

    if radar.ray_angle_res is not None:
        beam_spacing = radar.ray_angle_res['data'][0]

    lat = float(radar.latitude['data'])
    lon = float(radar.longitude['data'])
    alt = 0.

    # cartesian mapping
    grid = pyart.map.grid_from_radars(
        (radar,), gridding_algo='map_to_grid', weighting_function=wfunc,
        roi_func='dist_beam', h_factor=1.0, nb=beamwidth, bsp=beam_spacing,
        min_radius=cappi_res/2.,
        grid_shape=(1, ny, nx),
        grid_limits=((altitude, altitude), (ymin*1000., ymax*1000.),
                     (xmin*1000., xmax*1000.)),
        grid_origin=(lat, lon), grid_origin_alt=alt,
        fields=[field_name])

    # display data
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)

    fig = plt.figure(figsize=[prdcfg['ppiImageConfig']['xsize'],
                              prdcfg['ppiImageConfig']['ysize']],
                     dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    cmap = pyart.config.get_field_colormap(field_name)

    vmin = vmax = None
    if norm is None:  # if norm is set do not override with vmin/vmax
        vmin, vmax = pyart.config.get_field_limits(field_name)

    titl = pyart.graph.common.generate_grid_title(grid, field_name, 0)

    cax = ax.imshow(
        grid.fields[field_name]['data'][0], extent=(xmin, xmax, ymin, ymax),
        origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, norm=norm,
        interpolation='none')
    ax.set_xlabel('East West distance from radar(km)')
    ax.set_ylabel('North South distance from radar(km)')
    ax.set_title(titl)

    # plot the colorbar and set the label.
    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    label = get_colobar_label(grid.fields[field_name], field_name)
    cb.set_label(label)

    # Make a tight layout
    fig.tight_layout()

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list

    return (fig, ax)


def plot_traj(rng_traj, azi_traj, ele_traj, time_traj, prdcfg, fname_list,
              rad_alt=None, rad_tstart=None, ax=None, fig=None,
              save_fig=True):
    """
    plots a trajectory on a Cartesian surface

    Parameters
    ----------
    rng_traj, azi_traj, ele_traj : float array
        antenna coordinates of the trajectory [m and deg]
    time_traj : datetime array
        trajectory time
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    rad_alt : float or None
        radar altitude [m MSL]
    rad_tstart : datetime object or None
        start time of the radar scan
    surface_alt : float
        surface altitude [m MSL]
    color_ref : str
        What the color code represents. Can be 'None', 'rel_altitude',
        'altitude' or 'time'
    fig : Figure
        Figure to add the colorbar to. If none a new figure will be created
    ax : Axis
        Axis to plot on. if fig is None a new axis will be created
    save_fig : bool
        if true save the figure if false it does not close the plot and
        returns the handle to the figure

    Returns
    -------
    fname_list : list of str or
    fig, ax : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    color_ref = prdcfg.get('color_ref', 'None')
    if 'altitude' not in prdcfg and color_ref == 'rel_altitude':
        warn('Unable to plot trajectory relative to surface altitude. ' +
             'Unknown surface altitude.')
        color_ref = 'None'
    if rad_tstart is None and color_ref == 'time':
        warn('Unable to plot trajectory relative to radar scan start time. ' +
             'Unknown radar scan start time.')
        color_ref = 'None'
    if rad_alt is None and color_ref in ('rel_altitude', 'altitude'):
        warn('Unable to plot trajectory altitude. ' +
             'Unknown radar altitude.')
        color_ref = 'None'

    x, y, z = pyart.core.antenna_to_cartesian(
        rng_traj/1000., azi_traj, ele_traj)

    if color_ref == 'rel_altitude':
        h = z+rad_alt
        h_rel = h-prdcfg['altitude']

        marker = 'x'
        col = h_rel
        cmap = 'coolwarm'
        norm = plt.Normalize(-2000., 2000.)
        cb_label = 'Altitude relative to CAPPI [m]'
        plot_cb = True
    elif color_ref == 'altitude':
        h = z+rad_alt

        marker = 'x'
        col = h
        cmap = 'Greys'
        norm = plt.Normalize(h.min(), h.max())
        cb_label = 'Altitude [m MSL]'
        plot_cb = True
    elif color_ref == 'time':
        td_vec = time_traj-rad_tstart
        tt_s = []
        for td in td_vec:
            tt_s.append(td.total_seconds())
        tt_s = np.asarray(tt_s)

        marker = 'x'
        col = tt_s
        cmap = 'Greys'
        norm = plt.Normalize(tt_s.min(), tt_s.max())
        cb_label = 'Time from start of radar scan [s]'
        plot_cb = True
    else:
        col = 'k'
        marker = 'x'
        cmap = None
        plot_cb = False

    # display data
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    if fig is None:
        fig = plt.figure(figsize=[prdcfg['ppiImageConfig']['xsize'],
                                  prdcfg['ppiImageConfig']['ysize']],
                         dpi=dpi)
        ax = fig.add_subplot(111, aspect='equal')
    else:
        ax.autoscale(False)

    cax = ax.scatter(
        x/1000., y/1000., c=col, marker=marker, alpha=0.5, cmap=cmap,
        norm=norm)

    # plot colorbar
    if plot_cb:
        cb = fig.colorbar(cax, orientation='horizontal')
        cb.set_label(cb_label)

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list

    return (fig, ax)


def plot_rhi_contour(radar, field_name, ind_az, prdcfg, fname_list,
                     contour_values=None, linewidths=1.5, ax=None, fig=None,
                     save_fig=True):
    """
    plots contour data on an RHI

    Parameters
    ----------
    radar : Radar object
        object containing the radar data to plot
    field_name : str
        name of the radar field to plot
    ind_az : int
        sweep index to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    contour_values : float array
        list of contours to plot
    linewidths : float
        width of the contour lines
    fig : Figure
        Figure to add the colorbar to. If none a new figure will be created
    ax : Axis
        Axis to plot on. if fig is None a new axis will be created
    save_fig : bool
        if true save the figure if false it does not close the plot and
        returns the handle to the figure

    Returns
    -------
    fname_list : list of str or
    fig, ax : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)

    # get contour intervals
    if contour_values is None:
        field_dict = pyart.config.get_metadata(field_name)
        if 'boundaries' in field_dict:
            vmin = field_dict['boundaries'][0]
            vmax = field_dict['boundaries'][-1]
            num = len(field_dict['boundaries'])
        else:
            vmin, vmax = pyart.config.get_field_limits(field_name)
            num = 10

        contour_values = np.linspace(vmin, vmax, num=num)

    # get data and position
    display = pyart.graph.RadarDisplay(radar)
    data = display._get_data(field_name, ind_az, None, True, None)

    x_edges, y_edges, z_edges = display._get_x_y_z(ind_az, True, True)
    delta_x = x_edges[1:, 1:]-x_edges[:-1, :-1]
    delta_y = y_edges[1:, 1:]-y_edges[:-1, :-1]
    delta_z = z_edges[1:, 1:]-z_edges[:-1, :-1]

    x = x_edges[:-1, :-1]+delta_x/2.
    y = y_edges[:-1, :-1]+delta_y/2.
    z = z_edges[:-1, :-1]+delta_z/2.

    R = np.sqrt(x ** 2 + y ** 2) * np.sign(x)

    # display data
    if fig is None:
        xsize = prdcfg['rhiImageConfig']['xsize']
        ysize = prdcfg['rhiImageConfig']['ysize']
        fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
        ax = fig.add_subplot(111, aspect='equal')

        ax.contour(R, z, data, contour_values, colors='k',
                   linewidths=linewidths)

        display._set_title(field_name, ind_az, None, ax)
        display._label_axes_rhi((None, None), ax)

        display.set_limits(
            ylim=[prdcfg['rhiImageConfig']['ymin'],
                  prdcfg['rhiImageConfig']['ymax']],
            xlim=[prdcfg['rhiImageConfig']['xmin'],
                  prdcfg['rhiImageConfig']['xmax']], ax=ax)
        display.plot_cross_hair(5., ax=ax)

        # Turn on the grid
        ax.grid()

        # Make a tight layout
        fig.tight_layout()
    else:
        ax.autoscale(False)
        ax.contour(R, z, data, contour_values, colors='k',
                   linewidths=linewidths)

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list

    return (fig, ax)


def plot_ppi_contour(radar, field_name, ind_el, prdcfg, fname_list,
                     contour_values=None, linewidths=1.5, ax=None, fig=None,
                     save_fig=True):
    """
    plots contour data on a PPI

    Parameters
    ----------
    radar : Radar object
        object containing the radar data to plot
    field_name : str
        name of the radar field to plot
    ind_el : int
        sweep index to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    contour_values : float array
        list of contours to plot
    linewidths : float
        width of the contour lines
    fig : Figure
        Figure to add the contour to. If none a new figure will be created
    ax : Axis
        Axis to plot on. if fig is None a new axis will be created
    save_fig : bool
        if true save the figure if false it does not close the plot and
        returns the handle to the figure

    Returns
    -------
    fname_list : list of str or
    fig, ax : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)

    # get contour intervals
    if contour_values is None:
        field_dict = pyart.config.get_metadata(field_name)
        if 'boundaries' in field_dict:
            vmin = field_dict['boundaries'][0]
            vmax = field_dict['boundaries'][-1]
            num = len(field_dict['boundaries'])
        else:
            vmin, vmax = pyart.config.get_field_limits(field_name)
            num = 10

        contour_values = np.linspace(vmin, vmax, num=num)

    # get data and position
    display = pyart.graph.RadarDisplay(radar)
    data = display._get_data(field_name, ind_el, None, True, None)

    x_edges, y_edges = display._get_x_y(ind_el, True, True)
    delta_x = x_edges[1:, 1:]-x_edges[:-1, :-1]
    delta_y = y_edges[1:, 1:]-y_edges[:-1, :-1]

    x = x_edges[:-1, :-1]+delta_x/2.
    y = y_edges[:-1, :-1]+delta_y/2.

    # display data
    if fig is None:
        xsize = prdcfg['ppiImageConfig']['xsize']
        ysize = prdcfg['ppiImageConfig']['ysize']
        fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
        ax = fig.add_subplot(111, aspect='equal')

        ax.contour(x, y, data, contour_values, colors='k',
                   linewidths=linewidths)

        display._set_title(field_name, ind_el, None, ax)
        display._label_axes_ppi((None, None), ax)

        display.set_limits(
            ylim=[prdcfg['ppiImageConfig']['ymin'],
                  prdcfg['ppiImageConfig']['ymax']],
            xlim=[prdcfg['ppiImageConfig']['xmin'],
                  prdcfg['ppiImageConfig']['xmax']], ax=ax)
        display.plot_cross_hair(5., ax=ax)

        # Turn on the grid
        ax.grid()

        # Make a tight layout
        fig.tight_layout()
    else:
        ax.autoscale(False)
        ax.contour(x, y, data, contour_values, colors='k',
                   linewidths=linewidths)

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list

    return (fig, ax)


def plot_roi_contour(roi_dict, prdcfg, fname_list, plot_center=True,
                     xlabel='Lon [Deg]', ylabel='Lat [Deg]',
                     titl='TRT cell position', ax=None,
                     fig=None, save_fig=True):
    """
    plots the contour of a region of interest on a map

    Parameters
    ----------
    roi_dict : dict
        dictionary containing lon_roi, lat_roi, the points defining the
        contour
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    plot_center : bool
        If True a marked with the center of the roi is plotted
    fig : Figure
        Figure to add the colorbar to. If none a new figure will be created
    ax : Axis
        Axis to plot on. if fig is None a new axis will be created
    save_fig : bool
        if true save the figure if false it does not close the plot and
        returns the handle to the figure

    Returns
    -------
    fname_list : list of str or
    fig, ax : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    if not _SHAPELY_AVAILABLE or not _CARTOPY_AVAILABLE:
        warn('Unable to plot ROI contour: Missing shapely and/or'
             ' cartopy modules')
        return None

    dpi = prdcfg['ppiMapImageConfig'].get('dpi', 72)

    # create polygon to plot
    polygon = shapely.geometry.Polygon(list(zip(
        roi_dict['lon'], roi_dict['lat'])))

    # display data
    if fig is None:
        xsize = prdcfg['ppiMapImageConfig']['xsize']
        ysize = prdcfg['ppiMapImageConfig']['ysize']
        lonstep = prdcfg['ppiMapImageConfig'].get('lonstep', 0.5)
        latstep = prdcfg['ppiMapImageConfig'].get('latstep', 0.5)
        min_lon = prdcfg['ppiMapImageConfig'].get('lonmin', 2.5)
        max_lon = prdcfg['ppiMapImageConfig'].get('lonmax', 12.5)
        min_lat = prdcfg['ppiMapImageConfig'].get('latmin', 43.5)
        max_lat = prdcfg['ppiMapImageConfig'].get('latmax', 49.5)
        resolution = prdcfg['ppiMapImageConfig'].get('mapres', '110m')
        if resolution not in ('110m', '50m', '10m'):
            warn('Unknown map resolution: '+resolution)
            resolution = '110m'
        background_zoom = prdcfg['ppiMapImageConfig'].get(
            'background_zoom', 8)

        lon_lines = np.arange(np.floor(min_lon), np.ceil(max_lon)+1, lonstep)
        lat_lines = np.arange(np.floor(min_lat), np.ceil(max_lat)+1, latstep)
        limits = [min_lon, max_lon, min_lat, max_lat]

        # get background map instance
        stamen_terrain = Stamen('terrain-background')
        projection = cartopy.crs.PlateCarree()

        fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)

        # draw background
        ax = fig.add_subplot(111, projection=stamen_terrain.crs)
        ax.set_extent(limits, crs=projection)
        ax.add_image(stamen_terrain, background_zoom)

        # add countries
        countries = cartopy.feature.NaturalEarthFeature(
            category='cultural',
            name='admin_0_countries',
            scale=resolution,
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

        ax.set_title(titl)
    else:
        ax.autoscale(False)

    ax.add_geometries(
        [polygon], cartopy.crs.PlateCarree(), facecolor='none',
        edgecolor='k')

    if 'lon_center' in roi_dict and 'lat_center' in roi_dict and plot_center:
        ax.scatter(
            [roi_dict['lon_center']], [roi_dict['lat_center']], c='k',
            marker='x', transform=cartopy.crs.PlateCarree())

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        return fname_list

    return (fig, ax)


def plot_rhi_profile(data_list, hvec, fname_list, labelx='Value',
                     labely='Height (m MSL)', labels=('Mean'),
                     title='RHI profile', colors=None, linestyles=None,
                     vmin=None, vmax=None, hmin=None, hmax=None, dpi=72):
    """
    plots an RHI profile

    Parameters
    ----------
    data_list : list of float array
        values of the profile
    hvec : float array
        height points of the profile
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    labels : array of str
        The label of the legend
    title : str
        The figure title
    colors : array of str
        Specifies the colors of each line
    linestyles : array of str
        Specifies the line style of each line
    vmin, vmax: float
        Lower/Upper limit of data values
    hmin, hmax: float
        Lower/Upper limit of altitude
    dpi : int
        dots per inch

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    fig, ax = plt.subplots(figsize=[10, 6], dpi=dpi)

    lab = None
    col = None
    lstyle = None
    for i, data in enumerate(data_list):
        if labels is not None:
            lab = labels[i]
        if colors is not None:
            col = colors[i]
        if linestyles is not None:
            lstyle = linestyles[i]
        ax.plot(
            data, hvec, label=lab, color=col, linestyle=lstyle, marker='x')

    ax.set_title(title)
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_xlim(left=vmin, right=vmax)
    ax.set_ylim(bottom=hmin, top=hmax)
    ax.legend(loc='best')

    # Turn on the grid
    ax.grid()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_along_coord(xval_list, yval_list, fname_list, labelx='coord',
                     labely='Value', labels=None,
                     title='Plot along coordinate', colors=None,
                     linestyles=None, ymin=None, ymax=None, dpi=72,
                     data_on_y=True, plot_legend=True):
    """
    plots data along a certain radar coordinate

    Parameters
    ----------
    xval_list : list of float arrays
        the x values, range, azimuth or elevation
    yval_list : list of float arrays
        the y values. Parameter to plot
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    labels : array of str
        The label of the legend
    title : str
        The figure title
    colors : array of str
        Specifies the colors of each line
    linestyles : array of str
        Specifies the line style of each line
    ymin, ymax: float
        Lower/Upper limit of y axis
    dpi : int
        dots per inch
    data_on_y : bool
        If True the data is in the y axis and the coordinate in the x axis.
        False swaps it
    plot_legend : bool
        If True a legend will be plotted

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    fig, ax = plt.subplots(figsize=[10, 6], dpi=dpi)

    lab = None
    col = None
    lstyle = None

    for i, xval in enumerate(xval_list):
        yval = yval_list[i]
        if labels is not None:
            lab = labels[i]
        if colors is not None:
            col = colors[i]
        if linestyles is not None:
            lstyle = linestyles[i]
        if data_on_y:
            ax.plot(xval, yval, label=lab, color=col, linestyle=lstyle)
        else:
            ax.plot(yval, xval, label=lab, color=col, linestyle=lstyle)

    ax.set_title(title)
    if data_on_y:
        ax.set_xlabel(labelx)
        ax.set_ylabel(labely)
        ax.set_ylim(bottom=ymin, top=ymax)
    else:
        ax.set_xlabel(labely)
        ax.set_ylabel(labelx)
        ax.set_xlim(left=ymin, right=ymax)
    if plot_legend:
        ax.legend(loc='best')

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_field_coverage(xval_list, yval_list, fname_list,
                        labelx='Azimuth (deg)', labely='Range extension [m]',
                        labels=None, title='Field coverage', ymin=None,
                        ymax=None, xmeanval=None, ymeanval=None,
                        labelmeanval=None, dpi=72):
    """
    plots a time series

    Parameters
    ----------
    xval_list : list of float arrays
        the x values, azimuth
    yval_list : list of float arrays
        the y values. Range extension
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    labels : array of str
        The label of the legend
    title : str
        The figure title
    ymin, ymax : float
        Lower/Upper limit of y axis
    xmeanval, ymeanval : float array
        the x and y values of a mean along elevation
    labelmeanval : str
        the label of the mean
    dpi : int
        dots per inch

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    fig, ax = plt.subplots(figsize=[10, 6], dpi=dpi)

    lab = None

    for i, xval in enumerate(xval_list):
        yval = yval_list[i]
        if labels is not None:
            lab = labels[i]
        ax.plot(xval, yval, label=lab, linestyle='None', marker='o',
                fillstyle='full')

    if xmeanval is not None and ymeanval is not None:
        ax.plot(xmeanval, ymeanval, label=labelmeanval, linestyle='-',
                color='r', marker='x')

    ax.set_title(title)
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_ylim(bottom=ymin, top=ymax)
    if labels is not None:
        ax.legend(loc='best')

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list
