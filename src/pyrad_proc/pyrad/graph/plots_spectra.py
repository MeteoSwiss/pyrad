"""
pyrad.graph.plots_spectra
=========================

Functions to plot spectral data

.. autosummary::
    :toctree: generated/

    plot_range_Doppler
    plot_angle_Doppler
    plot_time_Doppler
    plot_Doppler
    plot_complex_range_Doppler
    plot_complex_angle_Doppler
    plot_complex_time_Doppler
    plot_amp_phase_range_Doppler
    plot_amp_phase_angle_Doppler
    plot_amp_phase_time_Doppler
    plot_complex_Doppler
    plot_amp_phase_Doppler
    _create_irregular_grid
    _adapt_data_to_irregular_grid

"""

import numpy as np

import matplotlib as mpl
mpl.use('Agg')

# Increase a bit font size
mpl.rcParams.update({'font.size': 16})
mpl.rcParams.update({'font.family':  "sans-serif"})

import matplotlib.pyplot as plt

import pyart

from .plots_aux import get_colobar_label, get_norm
from .plots_aux import generate_complex_range_Doppler_title
from .plots_aux import generate_complex_Doppler_title
from .plots_aux import generate_angle_Doppler_title


def plot_range_Doppler(spectra, field_name, ray, prdcfg, fname_list,
                       xaxis_info='Doppler_velocity', titl=None,
                       clabel=None, vmin=None, vmax=None):
    """
    Makes a range-Doppler plot

    Parameters
    ----------
    spectra : radar spectra object
        object containing the spectra or the IQ data to plot
    field_name : str
        name of the field to plot
    ray : int
        ray index
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    xaxis_info : str
        Type of x-axis. Can be 'Doppler_velocity', 'Doppler_frequency' or
        'pulse_number'
    titl : str or None
        The plot title
    clabel : str or None
        The color bar label
    vmin, vmax : float or None
        The value limits

    Returns
    -------
    fname_list : list of str
        list of names of the saved plots

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']

    if xaxis_info == 'Doppler_velocity':
        xaxis = spectra.Doppler_velocity['data'][ray, :].compressed()
        xlabel = 'Doppler velocity (m/s)'
    elif xaxis_info == 'Doppler_frequency':
        xaxis = spectra.Doppler_frequency['data'][ray, :].compressed()
        xlabel = 'Doppler frequency (Hz)'
    else:
        xaxis = np.arange(spectra.npulses['data'][ray])
        xlabel = 'Pulse number'

    xres = np.abs(xaxis[1]-xaxis[0])

    yaxis = spectra.range['data']/1000.
    yres = np.abs(yaxis[1]-yaxis[0])
    ylabel = 'range (km)'

    xaxis_lim = np.append(xaxis-xres/2, xaxis[-1]+xres/2)
    yaxis_lim = np.append(yaxis-yres/2, yaxis[-1]+yres/2)

    field_2D = spectra.fields[field_name]['data'][ray, :, 0:xaxis.size]

    X, Y = np.meshgrid(xaxis_lim, yaxis_lim)

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

        norm, ticks, ticklabs = get_norm(
            field_name, field_dict=spectra.fields[field_name])
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
            vmin = np.ma.min(field_2D)
        if vmax is None:
            vmax = np.ma.max(field_2D)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_complex_range_Doppler_title(
            spectra, field_name, ray, datetime_format=None)

    ax = fig.add_subplot(111)
    cax = ax.pcolormesh(
        X, Y, field_2D, cmap=cmap, vmin=vmin, vmax=vmax,
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

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_angle_Doppler(spectra, field_name, ang, ind_rays, ind_rng, prdcfg,
                       fname_list, xaxis_info='Doppler_velocity',
                       yaxis_pos='centre', along_azi=True, titl=None,
                       clabel=None, vmin=None, vmax=None):
    """
    Makes an angle-Doppler plot

    Parameters
    ----------
    spectra : radar spectra object
        object containing the spectra or the IQ data to plot
    field_name : str
        name of the field to plot
    ang : float
        The fixed angle
    ind_rays : 1D int array
        The indices of the rays to plot
    ind_rng : int
        The index of the range to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    xaxis_info : str
        Type of x-axis. Can be 'Doppler_velocity', 'Doppler_frequency' or
        'pulse_number'
    yaxis_pos : str
        the position that the y point represents in the y-axis bin. Can be
        'start', end' or 'centre'
    along_azi : bool
        If true the plot is performed along azimuth. If false it is performed
        along elevation
    titl : str or None
        The plot title
    clabel : str or None
        The color bar label
    vmin, vmax : float or None
        The value limits

    Returns
    -------
    fname_list : list of str
        list of names of the saved plots

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']

    if xaxis_info == 'Doppler_velocity':
        xaxis = spectra.Doppler_velocity['data'][ind_rays, :]
        xlabel = 'Doppler velocity (m/s)'
    elif xaxis_info == 'Doppler_frequency':
        xaxis = spectra.Doppler_frequency['data'][ind_rays, :]
        xlabel = 'Doppler frequency (Hz)'
    else:
        xaxis = np.ma.masked_all((ind_rays.size, spectra.npulses_max))
        for ray, ind_ray in enumerate(ind_rays):
            npuls = spectra.npulses['data'][ind_ray]
            xaxis[ray, 0:npuls] = np.arange(npuls)
        xlabel = 'Pulse number'
    ylabel = 'Angle (deg)'

    if along_azi:
        yaxis = spectra.azimuth['data'][ind_rays]
    else:
        yaxis = spectra.elevation['data'][ind_rays]

    xaxis_lim, yaxis_lim, xaxis_size = _create_irregular_grid(
        xaxis, yaxis, yaxis_pos=yaxis_pos)

    data = spectra.fields[field_name]['data'][ind_rays, :, :]
    data = data[:, ind_rng, :]
    zpoints = _adapt_data_to_irregular_grid(
        data, xaxis_size, xaxis_lim.shape[0], xaxis_lim.shape[1])

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

        norm, ticks, ticklabs = get_norm(
            field_name, field_dict=spectra.fields[field_name])
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
            vmin = np.ma.min(zpoints)
        if vmax is None:
            vmax = np.ma.max(zpoints)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_angle_Doppler_title(
            spectra, field_name, ang, ind_rng, along_azi=along_azi,
            datetime_format=None)

    ax = fig.add_subplot(111)
    cax = ax.pcolor(
        xaxis_lim, yaxis_lim, zpoints, cmap=cmap, vmin=vmin, vmax=vmax,
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

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_time_Doppler(spectra, field_name, prdcfg, fname_list,
                      xaxis_info='Doppler_velocity', yaxis_pos='start',
                      titl=None, clabel=None, vmin=None, vmax=None, xmin=None,
                      xmax=None, ymin=None, ymax=None):
    """
    Makes a time-Doppler plot

    Parameters
    ----------
    spectra : radar spectra object
        object containing the spectra or the IQ data to plot
    field_name : str
        name of the field to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    xaxis_info : str
        Type of x-axis. Can be 'Doppler_velocity', 'Doppler_frequency' or
        'pulse_number'
    yaxis_pos : str
        the position that the y point represents in the y-axis bin. Can be
        'start', end' or 'centre'
    titl : str or None
        The plot title
    clabel : str or None
        The color bar label
    vmin, vmax : float or None
        The value limits
    xmin, xmax, ymin, ymax : float or None
        The axis limits

    Returns
    -------
    fname_list : list of str
        list of names of the saved plots

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']

    if xaxis_info == 'Doppler_velocity':
        xaxis = spectra.Doppler_velocity['data']
        xlabel = 'Doppler velocity (m/s)'
    elif xaxis_info == 'Doppler_frequency':
        xaxis = spectra.Doppler_frequency['data']
        xlabel = 'Doppler frequency (Hz)'
    else:
        xaxis = np.ma.masked_all((spectra.nrays, spectra.npulses_max))
        for ray, npuls in enumerate(spectra.npulses['data']):
            xaxis[ray, 0:npuls] = np.arange(npuls)
        xlabel = 'Pulse number'
    ylabel = 'time (s from start time)'

    xaxis_lim, yaxis_lim, xaxis_size = _create_irregular_grid(
        xaxis, spectra.time['data'], yaxis_pos=yaxis_pos)

    zpoints = _adapt_data_to_irregular_grid(
        np.ma.squeeze(spectra.fields[field_name]['data'], axis=1), xaxis_size,
        xaxis_lim.shape[0], xaxis_lim.shape[1])

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

        norm, ticks, ticklabs = get_norm(
            field_name, field_dict=spectra.fields[field_name])
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
            vmin = np.ma.min(zpoints)
        if vmax is None:
            vmax = np.ma.max(zpoints)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_complex_Doppler_title(
            spectra, field_name, 0, 0, datetime_format=None)

    ax = fig.add_subplot(111)
    cax = ax.pcolor(
        xaxis_lim, yaxis_lim, zpoints, cmap=cmap, vmin=vmin, vmax=vmax,
        norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(titl)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label(clabel)

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_Doppler(spectra, field_name, ray, rng, prdcfg, fname_list,
                 xaxis_info='Doppler_velocity', ylabel=None,
                 titl=None, vmin=None, vmax=None):
    """
    Makes a Doppler plot

    Parameters
    ----------
    spectra : radar spectra object
        object containing the spectra or the IQ data to plot
    field_name : str
        name of the field to plot
    ray, rng : int
        ray and rng index
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    xaxis_info : str
        Type of x-axis. Can be 'Doppler_velocity', 'Doppler_frequency' or
        'pulse_number'
    ylabel : str or None
        The label of the y-axis
    titl : str or None
        The plot title
    vmin, vmax : float or None
        The value limits

    Returns
    -------
    fname_list : list of str
        list of names of the saved plots

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']

    if xaxis_info == 'Doppler_velocity':
        xaxis = spectra.Doppler_velocity['data'][ray, :].compressed()
        xlabel = 'Doppler velocity (m/s)'
    elif xaxis_info == 'Doppler_frequency':
        xaxis = spectra.Doppler_frequency['data'][ray, :].compressed()
        xlabel = 'Doppler frequency (Hz)'
    else:
        xaxis = np.arange(spectra.npulses['data'][ray])
        xlabel = 'Pulse number'

    field = spectra.fields[field_name]['data'][ray, rng, 0:xaxis.size]

    # display data
    if field_name is not None:
        if ylabel is None:
            field_dict = pyart.config.get_metadata(field_name)
            ylabel = get_colobar_label(field_dict, field_name)
        if vmin is None or vmax is None:
            vmin, vmax = pyart.config.get_field_limits(field_name)
    else:
        if vmin is None:
            vmin = np.ma.min(field)
        if vmax is None:
            vmax = np.ma.max(field)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_complex_Doppler_title(
            spectra, field_name, ray, rng, datetime_format=None)

    ax = fig.add_subplot(111)
    ax.plot(xaxis, field, marker='x')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=vmin, top=vmax)
    ax.set_xlim([xaxis[0], xaxis[-1]])
    ax.set_title(titl)

    # Turn on the grid
    ax.grid()

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_complex_range_Doppler(spectra, field_name, ray, prdcfg, fname_list,
                               xaxis_info='Doppler_velocity', titl=None,
                               clabel=None, vmin=None, vmax=None):
    """
    Makes a complex range-Doppler plot. Plotting separately the real and the
    imaginary part

    Parameters
    ----------
    spectra : radar spectra object
        object containing the spectra or the IQ data to plot
    field_name : str
        name of the field to plot
    ray : int
        ray index
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    xaxis_info : str
        Type of x-axis. Can be 'Doppler_velocity', 'Doppler_frequency' or
        'pulse_number'
    titl : str or None
        The plot title
    clabel : str or None
        The label of color bar
    vmin, vmax : float or None
        The value limits

    Returns
    -------
    fname_list : list of str
        list of names of the saved plots

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']

    if xaxis_info == 'Doppler_velocity':
        xaxis = spectra.Doppler_velocity['data'][ray, :].compressed()
        xlabel = 'Doppler velocity (m/s)'
    elif xaxis_info == 'Doppler_frequency':
        xaxis = spectra.Doppler_frequency['data'][ray, :].compressed()
        xlabel = 'Doppler frequency (Hz)'
    else:
        xaxis = np.arange(spectra.npulses['data'][ray])
        xlabel = 'Pulse number'

    xres = np.abs(xaxis[1]-xaxis[0])

    yaxis = spectra.range['data']/1000.
    yres = np.abs(yaxis[1]-yaxis[0])
    ylabel = 'range (km)'

    xaxis_lim = np.append(xaxis-xres/2, xaxis[-1]+xres/2)
    yaxis_lim = np.append(yaxis-yres/2, yaxis[-1]+yres/2)

    field_2D = spectra.fields[field_name]['data'][ray, :, 0:xaxis.size]

    re_field_2D = np.real(field_2D)
    im_field_2D = np.imag(field_2D)

    X, Y = np.meshgrid(xaxis_lim, yaxis_lim)

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

        norm, ticks, ticklabs = get_norm(
            field_name, field_dict=spectra.fields[field_name])
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
            vmin = np.ma.min([np.ma.min(re_field_2D), np.ma.min(im_field_2D)])
        if vmax is None:
            vmax = np.ma.max([np.ma.max(re_field_2D), np.ma.max(im_field_2D)])

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_complex_range_Doppler_title(
            spectra, field_name, ray, datetime_format=None)

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    cax = ax.pcolormesh(
        X, Y, re_field_2D, cmap=cmap, vmin=vmin, vmax=vmax,
        norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('Real part')

    ax = fig.add_subplot(122)
    cax = ax.pcolormesh(
        X, Y, im_field_2D, cmap=cmap, vmin=vmin, vmax=vmax,
        norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('Imaginary part')

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label(clabel)

    # Make a tight layout
    fig.tight_layout()

    plt.subplots_adjust(top=0.8)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_complex_angle_Doppler(spectra, field_name, ang, ind_rays, ind_rng,
                               prdcfg, fname_list,
                               xaxis_info='Doppler_velocity',
                               yaxis_pos='centre', along_azi=True, titl=None,
                               clabel=None, vmin=None, vmax=None):
    """
    Makes an angle-Doppler plot of complex spectra

    Parameters
    ----------
    spectra : radar spectra object
        object containing the spectra or the IQ data to plot
    field_name : str
        name of the field to plot
    ang : float
        The fixed angle
    ind_rays : 1D int array
        The indices of the rays to plot
    ind_rng : int
        The index of the range to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    xaxis_info : str
        Type of x-axis. Can be 'Doppler_velocity', 'Doppler_frequency' or
        'pulse_number'
    yaxis_pos : str
        the position that the y point represents in the y-axis bin. Can be
        'start', end' or 'centre'
    along_azi : bool
        If true the plot is performed along azimuth. If false it is performed
        along elevation
    titl : str or None
        The plot title
    clabel : str or None
        The color bar label
    vmin, vmax : float or None
        The value limits

    Returns
    -------
    fname_list : list of str
        list of names of the saved plots

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']

    if xaxis_info == 'Doppler_velocity':
        xaxis = spectra.Doppler_velocity['data'][ind_rays, :]
        xlabel = 'Doppler velocity (m/s)'
    elif xaxis_info == 'Doppler_frequency':
        xaxis = spectra.Doppler_frequency['data'][ind_rays, :]
        xlabel = 'Doppler frequency (Hz)'
    else:
        xaxis = np.ma.masked_all((ind_rays.size, spectra.npulses_max))
        for ray, ind_ray in enumerate(ind_rays):
            npuls = spectra.npulses['data'][ind_ray]
            xaxis[ray, 0:npuls] = np.arange(npuls)
        xlabel = 'Pulse number'
    ylabel = 'Angle (deg)'

    if along_azi:
        yaxis = spectra.azimuth['data'][ind_rays]
    else:
        yaxis = spectra.elevation['data'][ind_rays]

    xaxis_lim, yaxis_lim, xaxis_size = _create_irregular_grid(
        xaxis, yaxis, yaxis_pos=yaxis_pos)

    data = spectra.fields[field_name]['data'][ind_rays, :, :]
    data = data[:, ind_rng, :]
    field_2D = _adapt_data_to_irregular_grid(
        data, xaxis_size, xaxis_lim.shape[0], xaxis_lim.shape[1])

    re_field_2D = np.real(field_2D)
    im_field_2D = np.imag(field_2D)

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

        norm, ticks, ticklabs = get_norm(
            field_name, field_dict=spectra.fields[field_name])
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
            vmin = np.ma.min([np.ma.min(re_field_2D), np.ma.min(im_field_2D)])
        if vmax is None:
            vmax = np.ma.max([np.ma.max(re_field_2D), np.ma.max(im_field_2D)])

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_angle_Doppler_title(
            spectra, field_name, ang, ind_rng, along_azi=along_azi,
            datetime_format=None)

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    cax = ax.pcolor(
        xaxis_lim, yaxis_lim, re_field_2D, cmap=cmap, vmin=vmin, vmax=vmax,
        norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('Real part')

    ax = fig.add_subplot(122)
    cax = ax.pcolor(
        xaxis_lim, yaxis_lim, im_field_2D, cmap=cmap, vmin=vmin, vmax=vmax,
        norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('Imaginary part')

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label(clabel)

    # Make a tight layout
    fig.tight_layout()

    plt.subplots_adjust(top=0.8)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_complex_time_Doppler(spectra, field_name, prdcfg, fname_list,
                              xaxis_info='Doppler_velocity',
                              yaxis_pos='start', titl=None, clabel=None,
                              vmin=None, vmax=None):
    """
    Makes a complex time-Doppler plot. Plotting separately the real and the
    imaginary part

    Parameters
    ----------
    spectra : radar spectra object
        object containing the spectra or the IQ data to plot
    field_name : str
        name of the field to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    xaxis_info : str
        Type of x-axis. Can be 'Doppler_velocity', 'Doppler_frequency' or
        'pulse_number'
    yaxis_pos : str
        the position that the y point represents in the y-axis bin. Can be
        'start', end' or 'centre'
    titl : str or None
        The plot title
    clabel : str or None
        The label of color bar
    vmin, vmax : float or None
        The value limits

    Returns
    -------
    fname_list : list of str
        list of names of the saved plots

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']

    if xaxis_info == 'Doppler_velocity':
        xaxis = spectra.Doppler_velocity['data']
        xlabel = 'Doppler velocity (m/s)'
    elif xaxis_info == 'Doppler_frequency':
        xaxis = spectra.Doppler_frequency['data']
        xlabel = 'Doppler frequency (Hz)'
    else:
        xaxis = np.ma.masked_all((spectra.nrays, spectra.npulses_max))
        for ray, npuls in enumerate(spectra.npulses['data']):
            xaxis[ray, 0:npuls] = np.arange(npuls)
        xlabel = 'Pulse number'
    ylabel = 'time (s from start time)'

    xaxis_lim, yaxis_lim, xaxis_size = _create_irregular_grid(
        xaxis, spectra.time['data'], yaxis_pos=yaxis_pos)

    field_2D = _adapt_data_to_irregular_grid(
        np.ma.squeeze(spectra.fields[field_name]['data'], axis=1), xaxis_size,
        xaxis_lim.shape[0], xaxis_lim.shape[1])

    re_field_2D = np.real(field_2D)
    im_field_2D = np.imag(field_2D)

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

        norm, ticks, ticklabs = get_norm(
            field_name, field_dict=spectra.fields[field_name])
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
            vmin = np.ma.min([np.ma.min(re_field_2D), np.ma.min(im_field_2D)])
        if vmax is None:
            vmax = np.ma.max([np.ma.max(re_field_2D), np.ma.max(im_field_2D)])

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_complex_Doppler_title(
            spectra, field_name, 0, 0, datetime_format=None)

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    cax = ax.pcolor(
        xaxis_lim, yaxis_lim, re_field_2D, cmap=cmap, vmin=vmin, vmax=vmax,
        norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('Real part')

    ax = fig.add_subplot(122)
    cax = ax.pcolor(
        xaxis_lim, yaxis_lim, im_field_2D, cmap=cmap, vmin=vmin, vmax=vmax,
        norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('Imaginary part')

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label(clabel)

    # Make a tight layout
    fig.tight_layout()

    plt.subplots_adjust(top=0.8)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_amp_phase_range_Doppler(spectra, field_name, ray, prdcfg, fname_list,
                                 xaxis_info='Doppler_velocity', titl=None,
                                 ampli_vmin=None, ampli_vmax=None,
                                 phase_vmin=None, phase_vmax=None):
    """
    Makes a complex range-Doppler plot plotting separately the module and the
    phase of the signal

    Parameters
    ----------
    spectra : radar spectra object
        object containing the spectra or the IQ data to plot
    field_name : str
        name of the field to plot
    ray : int
        ray index
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    xaxis_info : str
        Type of x-axis. Can be 'Doppler_velocity', 'Doppler_frequency' or
        'pulse_number'
    titl : str or None
        The plot title
    ampli_vmin, ampli_vmax, phase_vmin, phase_vmax : float or None
        The value limits

    Returns
    -------
    fname_list : list of str
        list of names of the saved plots

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']

    if xaxis_info == 'Doppler_velocity':
        xaxis = spectra.Doppler_velocity['data'][ray, :].compressed()
        xlabel = 'Doppler velocity (m/s)'
    elif xaxis_info == 'Doppler_frequency':
        xaxis = spectra.Doppler_frequency['data'][ray, :].compressed()
        xlabel = 'Doppler frequency (Hz)'
    else:
        xaxis = np.arange(spectra.npulses['data'][ray])
        xlabel = 'Pulse number'
    xres = np.abs(xaxis[1]-xaxis[0])

    yaxis = spectra.range['data']/1000.
    yres = np.abs(yaxis[1]-yaxis[0])
    ylabel = 'range (km)'

    xaxis_lim = np.append(xaxis-xres/2, xaxis[-1]+xres/2)
    yaxis_lim = np.append(yaxis-yres/2, yaxis[-1]+yres/2)

    field_2D = spectra.fields[field_name]['data'][ray, :, 0:xaxis.size]

    ampli_field_2D = np.ma.abs(field_2D)
    phase_field_2D = np.ma.angle(field_2D, deg=True)

    X, Y = np.meshgrid(xaxis_lim, yaxis_lim)

    # display data
    norm = None
    cmap = None
    ticks = None
    ticklabs = None
    if field_name is not None:
        cmap = pyart.config.get_field_colormap(field_name)

        norm, ticks, ticklabs = get_norm(
            field_name, field_dict=spectra.fields[field_name])
        if ampli_vmin is None or ampli_vmax is None:
            ampli_vmin = ampli_vmax = None
            if norm is None:  # if norm is set do not override with vmin/vmax
                ampli_vmin, ampli_vmax = pyart.config.get_field_limits(
                    field_name)
        else:
            norm = None

        if phase_vmin is None or phase_vmax is None:
            phase_vmin = phase_vmax = None
            if norm is None:  # if norm is set do not override with vmin/vmax
                phase_vmin, phase_vmax = (-180., 180.)
        else:
            norm = None
    else:
        if ampli_vmin is None:
            ampli_vmin = np.ma.min(ampli_field_2D)
        if ampli_vmax is None:
            ampli_vmax = np.ma.max(ampli_field_2D)
        if phase_vmin is None:
            phase_vmin = np.ma.min(phase_field_2D)
        if phase_vmax is None:
            phase_vmax = np.ma.max(phase_field_2D)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_complex_range_Doppler_title(
            spectra, field_name, ray, datetime_format=None)

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    cax = ax.pcolormesh(
        X, Y, ampli_field_2D, cmap=cmap, vmin=ampli_vmin, vmax=ampli_vmax,
        norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('Amplitude')

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label('Amplitude (-)')

    ax = fig.add_subplot(122)
    cax = ax.pcolormesh(
        X, Y, phase_field_2D, cmap=cmap, vmin=phase_vmin, vmax=phase_vmax,
        norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('Phase')

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label('Phase (deg)')

    # Make a tight layout
    fig.tight_layout()

    plt.subplots_adjust(top=0.8)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_amp_phase_angle_Doppler(spectra, field_name, ang, ind_rays, ind_rng,
                                 prdcfg, fname_list,
                                 xaxis_info='Doppler_velocity',
                                 yaxis_pos='centre', along_azi=True,
                                 titl=None, ampli_vmin=None, ampli_vmax=None,
                                 phase_vmin=None, phase_vmax=None,):
    """
    Makes an angle-Doppler plot of complex spectra

    Parameters
    ----------
    spectra : radar spectra object
        object containing the spectra or the IQ data to plot
    field_name : str
        name of the field to plot
    ang : float
        The fixed angle
    ind_rays : 1D int array
        The indices of the rays to plot
    ind_rng : int
        The index of the range to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    xaxis_info : str
        Type of x-axis. Can be 'Doppler_velocity', 'Doppler_frequency' or
        'pulse_number'
    yaxis_pos : str
        the position that the y point represents in the y-axis bin. Can be
        'start', end' or 'centre'
    along_azi : bool
        If true the plot is performed along azimuth. If false it is performed
        along elevation
    titl : str or None
        The plot title
    ampli_vmin, ampli_vmax, phase_vmin, phase_vmax : float or None
        The value limits

    Returns
    -------
    fname_list : list of str
        list of names of the saved plots

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']

    if xaxis_info == 'Doppler_velocity':
        xaxis = spectra.Doppler_velocity['data'][ind_rays, :]
        xlabel = 'Doppler velocity (m/s)'
    elif xaxis_info == 'Doppler_frequency':
        xaxis = spectra.Doppler_frequency['data'][ind_rays, :]
        xlabel = 'Doppler frequency (Hz)'
    else:
        xaxis = np.ma.masked_all((ind_rays.size, spectra.npulses_max))
        for ray, ind_ray in enumerate(ind_rays):
            npuls = spectra.npulses['data'][ind_ray]
            xaxis[ray, 0:npuls] = np.arange(npuls)
        xlabel = 'Pulse number'
    ylabel = 'Angle (deg)'

    if along_azi:
        yaxis = spectra.azimuth['data'][ind_rays]
    else:
        yaxis = spectra.elevation['data'][ind_rays]

    xaxis_lim, yaxis_lim, xaxis_size = _create_irregular_grid(
        xaxis, yaxis, yaxis_pos=yaxis_pos)

    data = spectra.fields[field_name]['data'][ind_rays, :, :]
    data = data[:, ind_rng, :]
    field_2D = _adapt_data_to_irregular_grid(
        data, xaxis_size, xaxis_lim.shape[0], xaxis_lim.shape[1])

    ampli_field_2D = np.ma.abs(field_2D)
    phase_field_2D = np.ma.angle(field_2D, deg=True)

    # display data
    norm = None
    cmap = None
    ticks = None
    ticklabs = None
    if field_name is not None:
        cmap = pyart.config.get_field_colormap(field_name)

        norm, ticks, ticklabs = get_norm(
            field_name, field_dict=spectra.fields[field_name])
        if ampli_vmin is None or ampli_vmax is None:
            ampli_vmin = ampli_vmax = None
            if norm is None:  # if norm is set do not override with vmin/vmax
                ampli_vmin, ampli_vmax = pyart.config.get_field_limits(
                    field_name)
        else:
            norm = None

        if phase_vmin is None or phase_vmax is None:
            phase_vmin = phase_vmax = None
            if norm is None:  # if norm is set do not override with vmin/vmax
                phase_vmin, phase_vmax = (-180., 180.)
        else:
            norm = None
    else:
        if ampli_vmin is None:
            ampli_vmin = np.ma.min(ampli_field_2D)
        if ampli_vmax is None:
            ampli_vmax = np.ma.max(ampli_field_2D)
        if phase_vmin is None:
            phase_vmin = np.ma.min(phase_field_2D)
        if phase_vmax is None:
            phase_vmax = np.ma.max(phase_field_2D)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_angle_Doppler_title(
            spectra, field_name, ang, ind_rng, along_azi=along_azi,
            datetime_format=None)

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    cax = ax.pcolor(
        xaxis_lim, yaxis_lim, ampli_field_2D, cmap=cmap, vmin=ampli_vmin,
        vmax=ampli_vmax, norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('Amplitude')

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label('Amplitude (-)')

    ax = fig.add_subplot(122)
    cax = ax.pcolor(
        xaxis_lim, yaxis_lim, phase_field_2D, cmap=cmap, vmin=phase_vmin,
        vmax=phase_vmax, norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('Phase')

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label('Phase (deg)')

    # Make a tight layout
    fig.tight_layout()

    plt.subplots_adjust(top=0.8)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_amp_phase_time_Doppler(spectra, field_name, prdcfg, fname_list,
                                xaxis_info='Doppler_velocity',
                                yaxis_pos='start', titl=None, ampli_vmin=None,
                                ampli_vmax=None, phase_vmin=None,
                                phase_vmax=None):
    """
    Makes a complex time-Doppler plot plotting separately the module and the
    phase of the signal

    Parameters
    ----------
    spectra : radar spectra object
        object containing the spectra or the IQ data to plot
    field_name : str
        name of the field to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    xaxis_info : str
        Type of x-axis. Can be 'Doppler_velocity', 'Doppler_frequency' or
        'pulse_number'
    yaxis_pos : str
        the position that the y point represents in the y-axis bin. Can be
        'start', end' or 'centre'
    titl : str or None
        The plot title
    ampli_vmin, ampli_vmax, phase_vmin, phase_vmax : float or None
        The value limits

    Returns
    -------
    fname_list : list of str
        list of names of the saved plots

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']

    if xaxis_info == 'Doppler_velocity':
        xaxis = spectra.Doppler_velocity['data']
        xlabel = 'Doppler velocity (m/s)'
    elif xaxis_info == 'Doppler_frequency':
        xaxis = spectra.Doppler_frequency['data']
        xlabel = 'Doppler frequency (Hz)'
    else:
        xaxis = np.ma.masked_all((spectra.nrays, spectra.npulses_max))
        for ray, npuls in enumerate(spectra.npulses['data']):
            xaxis[ray, 0:npuls] = np.arange(npuls)
        xlabel = 'Pulse number'
    ylabel = 'time (s from start time)'

    xaxis_lim, yaxis_lim, xaxis_size = _create_irregular_grid(
        xaxis, spectra.time['data'], yaxis_pos=yaxis_pos)

    field_2D = _adapt_data_to_irregular_grid(
        np.ma.squeeze(spectra.fields[field_name]['data'], axis=1), xaxis_size,
        xaxis_lim.shape[0], xaxis_lim.shape[1])

    ampli_field_2D = np.ma.abs(field_2D)
    phase_field_2D = np.ma.angle(field_2D, deg=True)

    # display data
    norm = None
    cmap = None
    ticks = None
    ticklabs = None
    if field_name is not None:
        cmap = pyart.config.get_field_colormap(field_name)

        norm, ticks, ticklabs = get_norm(
            field_name, field_dict=spectra.fields[field_name])
        if ampli_vmin is None or ampli_vmax is None:
            ampli_vmin = ampli_vmax = None
            if norm is None:  # if norm is set do not override with vmin/vmax
                ampli_vmin, ampli_vmax = pyart.config.get_field_limits(
                    field_name)
        else:
            norm = None

        if phase_vmin is None or phase_vmax is None:
            phase_vmin = phase_vmax = None
            if norm is None:  # if norm is set do not override with vmin/vmax
                phase_vmin, phase_vmax = (-180., 180.)
        else:
            norm = None
    else:
        if ampli_vmin is None:
            ampli_vmin = np.ma.min(ampli_field_2D)
        if ampli_vmax is None:
            ampli_vmax = np.ma.max(ampli_field_2D)
        if phase_vmin is None:
            phase_vmin = np.ma.min(phase_field_2D)
        if phase_vmax is None:
            phase_vmax = np.ma.max(phase_field_2D)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_complex_Doppler_title(
            spectra, field_name, 0, 0, datetime_format=None)

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    cax = ax.pcolor(
        xaxis_lim, yaxis_lim, ampli_field_2D, cmap=cmap, vmin=ampli_vmin,
        vmax=ampli_vmax, norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('Amplitude')

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label('Amplitude (-)')

    ax = fig.add_subplot(122)
    cax = ax.pcolor(
        xaxis_lim, yaxis_lim, phase_field_2D, cmap=cmap, vmin=phase_vmin,
        vmax=phase_vmax, norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('Phase')

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label('Phase (deg)')

    # Make a tight layout
    fig.tight_layout()

    plt.subplots_adjust(top=0.8)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_complex_Doppler(spectra, field_name, ray, rng, prdcfg, fname_list,
                         xaxis_info='Doppler_velocity', ylabel=None,
                         titl=None, vmin=None, vmax=None):
    """
    Makes a complex Doppler plot plotting separately the real and the
    imaginary parts

    Parameters
    ----------
    spectra : radar spectra object
        object containing the spectra or the IQ data to plot
    field_name : str
        name of the field to plot
    ray, rng : int
        ray and range index
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    xaxis_info : str
        Type of x-axis. Can be 'Doppler_velocity', 'Doppler_frequency' or
        'pulse_number'
    ylabel : str or None
        The label of the y-axis
    titl : str or None
        The plot title
    vmin, vmax : float or None
        The value limits

    Returns
    -------
    fname_list : list of str
        list of names of the saved plots

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']

    if xaxis_info == 'Doppler_velocity':
        xaxis = spectra.Doppler_velocity['data'][ray, :].compressed()
        xlabel = 'Doppler velocity (m/s)'
    elif xaxis_info == 'Doppler_frequency':
        xaxis = spectra.Doppler_frequency['data'][ray, :].compressed()
        xlabel = 'Doppler frequency (Hz)'
    else:
        xaxis = np.arange(spectra.npulses['data'][ray])
        xlabel = 'Pulse number'

    field = spectra.fields[field_name]['data'][ray, rng, 0:xaxis.size]

    re_field = np.real(field)
    im_field = np.imag(field)

    # display data
    if field_name is not None:
        if ylabel is None:
            field_dict = pyart.config.get_metadata(field_name)
            ylabel = get_colobar_label(field_dict, field_name)
        if vmin is None or vmax is None:
            vmin, vmax = pyart.config.get_field_limits(field_name)
    else:
        if vmin is None:
            vmin = np.ma.min([np.ma.min(re_field), np.ma.min(im_field)])
        if vmax is None:
            vmax = np.ma.max([np.ma.max(re_field), np.ma.max(im_field)])

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_complex_Doppler_title(
            spectra, field_name, ray, rng, datetime_format=None)

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    ax.plot(xaxis, re_field, marker='x')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=vmin, top=vmax)
    ax.set_xlim([xaxis[0], xaxis[-1]])
    ax.set_title('Real part')

    # Turn on the grid
    ax.grid()

    ax = fig.add_subplot(122)
    ax.plot(xaxis, im_field, marker='x')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=vmin, top=vmax)
    ax.set_xlim([xaxis[0], xaxis[-1]])
    ax.set_title('Imaginary part')

    # Turn on the grid
    ax.grid()

    # Make a tight layout
    fig.tight_layout()

    plt.subplots_adjust(top=0.8)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_amp_phase_Doppler(spectra, field_name, ray, rng, prdcfg, fname_list,
                           xaxis_info='Doppler_velocity', titl=None,
                           ampli_vmin=None, ampli_vmax=None, phase_vmin=None,
                           phase_vmax=None):
    """
    Makes a complex Doppler plot plotting separately the module and the phase
    of the signal

    Parameters
    ----------
    spectra : radar spectra object
        object containing the spectra or the IQ data to plot
    field_name : str
        name of the field to plot
    ray, rng : int
        ray and range index
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    xaxis_info : str
        Type of x-axis. Can be 'Doppler_velocity', 'Doppler_frequency' or
        'pulse_number'
    titl : str or None
        The plot title
    ampli_vmin, ampli_vmax, phase_vmin, phase_vmax : float or None
        The value limits

    Returns
    -------
    fname_list : list of str
        list of names of the saved plots

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']

    if xaxis_info == 'Doppler_velocity':
        xaxis = spectra.Doppler_velocity['data'][ray, :].compressed()
        xlabel = 'Doppler velocity (m/s)'
    elif xaxis_info == 'Doppler_frequency':
        xaxis = spectra.Doppler_frequency['data'][ray, :].compressed()
        xlabel = 'Doppler frequency (Hz)'
    else:
        xaxis = np.arange(spectra.npulses['data'][ray])
        xlabel = 'Pulse number'

    field = spectra.fields[field_name]['data'][ray, rng, 0:xaxis.size]

    ampli_field = np.ma.abs(field)
    phase_field = np.ma.angle(field, deg=True)

    # display data
    if field_name is not None:
        if ampli_vmin is None or ampli_vmax is None:
            ampli_vmin, ampli_vmax = pyart.config.get_field_limits(
                field_name)
        if phase_vmin is None or phase_vmax is None:
            phase_vmin, phase_vmax = (-180., 180.)
    else:
        if ampli_vmin is None:
            ampli_vmin = np.ma.min(ampli_field)
        if ampli_vmax is None:
            ampli_vmax = np.ma.max(ampli_field)
        if phase_vmin is None:
            phase_vmin = np.ma.min(phase_field)
        if phase_vmax is None:
            phase_vmax = np.ma.max(phase_field)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_complex_Doppler_title(
            spectra, field_name, ray, rng, datetime_format=None)

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    ax.plot(xaxis, ampli_field, marker='x')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Amplitude (-)')
    ax.set_ylim(bottom=ampli_vmin, top=ampli_vmax)
    ax.set_xlim([xaxis[0], xaxis[-1]])
    ax.set_title('Amplitude')

    # Turn on the grid
    ax.grid()

    ax = fig.add_subplot(122)
    ax.plot(xaxis, phase_field, marker='x')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Phase (deg)')
    ax.set_ylim(bottom=phase_vmin, top=phase_vmax)
    ax.set_xlim([xaxis[0], xaxis[-1]])
    ax.set_title('Phase')

    # Turn on the grid
    ax.grid()

    # Make a tight layout
    fig.tight_layout()

    plt.subplots_adjust(top=0.8)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def _create_irregular_grid(xaxis, yaxis, yaxis_pos='start'):
    """
    Create an irregular grid to be able to plot data with variable x-axis

    Parameters
    ----------
    xaxis : 2D float array
        array containing the x points
    yaxis : 1D float array
        array containing the y points
    yaxis_pos : str
        the position that the y point represents in the y-axis bin. Can be
        'start', end' or 'centre'

    Returns
    -------
    xaxis_lim, yaxis_lim : 2D float array
        Matrix containing the edges of the X and Y axis

    """
    nbinsy = yaxis.size

    # One dimensional y-axis bin limits
    yres = np.median(yaxis[1:]-yaxis[:-1])
    if yaxis_pos == 'start':
        yaxis_lim = np.append(yaxis, yaxis[-1]+yres)
    elif yaxis_pos == 'end':
        yaxis_lim = np.append(yaxis-yres, yaxis[-1])
    else:
        yaxis_lim = np.append(yaxis-yres/2., yaxis[-1]+yres/2.)

    # List of x-axis bin limits
    xres = np.abs(xaxis[:, 1]-xaxis[:, 0])
    xaxis_size = list()
    xaxis_list = list()
    for i in range(nbinsy):
        xaxis_aux = xaxis[i, :].compressed()
        xaxis_lim = np.append(
            xaxis_aux-xres[i]/2., xaxis_aux[-1]+xres[i]/2.)
        xaxis_size.append(xaxis_lim.size)
        xaxis_list.append(xaxis_lim)

    nxbin_lim = max(xaxis_size)

    # Create two dimensional grid of y-axis
    yaxis_lim = np.repeat(yaxis_lim, 2)[1: -1]
    nybin_lim = yaxis_lim.size
    yaxis_lim = np.broadcast_to(
        np.reshape(yaxis_lim, (1, nybin_lim)), (nxbin_lim, nybin_lim))

    # Create two dimensional grid of x-axis
    xaxis_lim = np.ma.masked_all((nxbin_lim, nbinsy))
    for i in range(nbinsy):
        xaxis_lim[0:xaxis_size[i], i] = xaxis_list[i]
    xaxis_lim = np.repeat(xaxis_lim, 2, axis=1)

    return xaxis_lim, yaxis_lim, xaxis_size


def _adapt_data_to_irregular_grid(data, xaxis_size, nxbin_lim, nybin_lim):
    """
    Adapts data to irregular grid to allow plotting

    Parameters
    ----------
    data : 2D float array
        The data to plot
    xaxis_size : 1D-array
        The size of each x-axis
    nxbin_lim, nybin_lim : float
        number of gate limits at each axis

    Returns
    -------
    zpoints : 2D float array
        Matrix containing the data points.

    """
    nbinsy = data.shape[0]

    # adapt Z data
    zpoints = np.ma.masked_all((nxbin_lim-1, nybin_lim-1), dtype=data.dtype)
    for i in range(nbinsy-1):
        zpoints[0:xaxis_size[i]-1, 2*i] = data[i, 0:xaxis_size[i]-1]
        zpoints[0:xaxis_size[i]-1, 2*i+1] = data[i, 0:xaxis_size[i]-1]
    zpoints[0:xaxis_size[-1]-1, -1] = data[-1, 0:xaxis_size[-1]-1]

    return zpoints
