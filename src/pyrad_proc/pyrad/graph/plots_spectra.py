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

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pyart

from .plots_aux import (
    generate_angle_Doppler_title,
    generate_complex_Doppler_title,
    generate_complex_range_Doppler_title,
    get_colobar_label,
    get_norm,
)

mpl.use("Agg")

# Increase a bit font size
mpl.rcParams.update({"font.size": 16})
mpl.rcParams.update({"font.family": "sans-serif"})


def plot_range_Doppler(
    spectra,
    field_name,
    ray,
    prdcfg,
    fname_list,
    xaxis_info="Doppler_velocity",
    titl=None,
    clabel=None,
    vmin=None,
    vmax=None,
    plot_mean_velocity=False,
    mean_velocity_field_name=None,
):
    """
    Make a range-Doppler plot.

    Optionally, the power-weighted mean Doppler velocity can be overlaid as a
    function of range. For each range gate, it is calculated as

        mean_velocity = sum(velocity * linear_power) / sum(linear_power)

    over all valid Doppler bins.

    If the weighting field is expressed in logarithmic units such as dB or
    dBZ, it is converted to linear units before calculating the weighted
    average.

    Parameters
    ----------
    spectra : radar spectra object
        Object containing the spectra or IQ data to plot.
    field_name : str
        Name of the field to plot.
    ray : int
        Ray index.
    prdcfg : dict
        Dictionary containing the product configuration.
    fname_list : list of str
        Names of the files where the plot will be stored.
    xaxis_info : str, optional
        Type of x-axis. Can be ``'Doppler_velocity'``,
        ``'Doppler_frequency'`` or ``'pulse_number'``.
    titl : str or None, optional
        Plot title.
    clabel : str or None, optional
        Colorbar label.
    vmin, vmax : float or None, optional
        Color scale limits.
    plot_mean_velocity : bool, optional
        If True, overlay the power-weighted mean Doppler velocity at each
        range gate. This is only supported when
        ``xaxis_info='Doppler_velocity'``.
    mean_velocity_field_name : str or None, optional
        Field used as the weighting power spectrum. If None, ``field_name`` is
        used. This allows, for example, plotting one spectral variable while
        calculating mean velocity from spectral reflectivity or spectral
        power.

    Returns
    -------
    fname_list : list of str
        Names of the saved plots.

    """
    dpi = prdcfg["spectraImageConfig"].get("dpi", 72)
    xsize = prdcfg["spectraImageConfig"]["xsize"]
    ysize = prdcfg["spectraImageConfig"]["ysize"]
    ymin = prdcfg["spectraImageConfig"]["ymin"]
    ymax = prdcfg["spectraImageConfig"]["ymax"]
    velmin = prdcfg["spectraImageConfig"]["velmin"]
    velmax = prdcfg["spectraImageConfig"]["velmax"]

    if xaxis_info == "Doppler_velocity":
        xaxis = np.ma.asarray(spectra.Doppler_velocity["data"][ray, :]).compressed()
        xlabel = "Doppler velocity (m/s)"
    elif xaxis_info == "Doppler_frequency":
        xaxis = np.ma.asarray(spectra.Doppler_frequency["data"][ray, :]).compressed()
        xlabel = "Doppler frequency (Hz)"
    else:
        xaxis = np.arange(spectra.npulses["data"][ray])
        xlabel = "Pulse number"

    if xaxis.size < 2:
        raise ValueError(
            "At least two Doppler bins or pulses are required to create "
            "the range-Doppler plot."
        )

    xres = np.abs(xaxis[1] - xaxis[0])

    yaxis = spectra.range["data"] / 1000.0

    if yaxis.size < 2:
        raise ValueError(
            "At least two range gates are required to create the " "range-Doppler plot."
        )

    yres = np.abs(yaxis[1] - yaxis[0])
    ylabel = "range (km)"

    xaxis_lim = np.append(
        xaxis - xres / 2.0,
        xaxis[-1] + xres / 2.0,
    )
    yaxis_lim = np.append(
        yaxis - yres / 2.0,
        yaxis[-1] + yres / 2.0,
    )

    field_2D = np.ma.asarray(spectra.fields[field_name]["data"][ray, :, 0 : xaxis.size])

    X, Y = np.meshgrid(xaxis_lim, yaxis_lim)

    # Display configuration
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
            field_name,
            field_dict=spectra.fields[field_name],
        )

        if vmin is None or vmax is None:
            vmin = None
            vmax = None

            # Do not combine explicit limits with a BoundaryNorm.
            if norm is None:
                vmin, vmax = pyart.config.get_field_limits(field_name)
        else:
            norm = None
    else:
        if clabel is None:
            clabel = "value"

        if vmin is None:
            vmin = np.ma.min(field_2D)

        if vmax is None:
            vmax = np.ma.max(field_2D)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)

    if titl is None:
        titl = generate_complex_range_Doppler_title(
            spectra,
            field_name,
            ray,
            datetime_format=None,
        )

    ax = fig.add_subplot(111)

    cax = ax.pcolormesh(
        X,
        Y,
        field_2D,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        norm=norm,
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(titl)

    if velmin is not None and velmax is not None:
        ax.set_xlim([velmin, velmax])

    if ymin is not None and ymax is not None:
        ax.set_ylim([ymin, ymax])

    cb = fig.colorbar(cax)

    if ticks is not None:
        cb.set_ticks(ticks)

    if ticklabs:
        cb.set_ticklabels(ticklabs)

    cb.set_label(clabel)

    # Overlay the power-weighted mean Doppler velocity.
    if plot_mean_velocity:
        if xaxis_info != "Doppler_velocity":
            raise ValueError(
                "The mean Doppler velocity can only be plotted when "
                "xaxis_info='Doppler_velocity'."
            )

        if mean_velocity_field_name is None:
            mean_velocity_field_name = field_name

        if mean_velocity_field_name not in spectra.fields:
            raise ValueError(
                f"Weighting field '{mean_velocity_field_name}' is not "
                "available in the spectra object."
            )

        weighting_field = spectra.fields[mean_velocity_field_name]
        spectral_weights = np.ma.asarray(
            weighting_field["data"][ray, :, 0 : xaxis.size],
            dtype=float,
        )

        # Ensure that invalid numerical values are masked.
        spectral_weights = np.ma.masked_invalid(spectral_weights)

        units = str(weighting_field.get("units", "")).lower()
        field_name_lower = mean_velocity_field_name.lower()

        is_logarithmic = (
            "db" in units
            or "decibel" in units
            or "dbz" in units
            or "db" in field_name_lower
        )

        if is_logarithmic:
            # Convert dB quantities to linear power before integration.
            linear_weights = np.ma.power(10.0, spectral_weights / 10.0)
        else:
            # Linear power or reflectivity values can be used directly.
            linear_weights = spectral_weights.copy()

        # Negative linear powers are not physically meaningful as weights.
        linear_weights = np.ma.masked_less_equal(linear_weights, 0.0)

        velocity_2D = np.broadcast_to(
            xaxis[np.newaxis, :],
            linear_weights.shape,
        )

        weighted_velocity_sum = np.ma.sum(
            linear_weights * velocity_2D,
            axis=1,
        )
        power_sum = np.ma.sum(linear_weights, axis=1)

        mean_velocity = np.ma.divide(
            weighted_velocity_sum,
            power_sum,
        )

        mean_velocity = np.ma.masked_where(
            np.ma.filled(power_sum, 0.0) <= 0.0,
            mean_velocity,
        )
        mean_velocity = np.ma.masked_invalid(mean_velocity)

        line_color = prdcfg["spectraImageConfig"].get(
            "mean_velocity_color",
            "white",
        )
        line_width = prdcfg["spectraImageConfig"].get(
            "mean_velocity_linewidth",
            1.0,
        )
        line_style = prdcfg["spectraImageConfig"].get(
            "mean_velocity_linestyle",
            "--",
        )

        ax.plot(
            mean_velocity,
            yaxis,
            color=line_color,
            linewidth=line_width,
            linestyle=line_style,
            zorder=10,
        )

    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")

    plt.close(fig)

    return fname_list


def plot_angle_Doppler(
    spectra,
    field_name,
    ang,
    ind_rays,
    ind_rng,
    prdcfg,
    fname_list,
    xaxis_info="Doppler_velocity",
    yaxis_pos="centre",
    along_azi=True,
    titl=None,
    clabel=None,
    vmin=None,
    vmax=None,
):
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
    dpi = prdcfg["spectraImageConfig"].get("dpi", 72)
    xsize = prdcfg["spectraImageConfig"]["xsize"]
    ysize = prdcfg["spectraImageConfig"]["ysize"]
    velmin = prdcfg["spectraImageConfig"]["velmin"]
    velmax = prdcfg["spectraImageConfig"]["velmax"]

    if xaxis_info == "Doppler_velocity":
        xaxis = spectra.Doppler_velocity["data"][ind_rays, :]
        xlabel = "Doppler velocity (m/s)"
    elif xaxis_info == "Doppler_frequency":
        xaxis = spectra.Doppler_frequency["data"][ind_rays, :]
        xlabel = "Doppler frequency (Hz)"
    else:
        xaxis = np.ma.masked_all((ind_rays.size, spectra.npulses_max))
        for ray, ind_ray in enumerate(ind_rays):
            npuls = spectra.npulses["data"][ind_ray]
            xaxis[ray, 0:npuls] = np.arange(npuls)
        xlabel = "Pulse number"
    ylabel = "Angle (deg)"

    if along_azi:
        yaxis = spectra.azimuth["data"][ind_rays]
    else:
        yaxis = spectra.elevation["data"][ind_rays]

    xaxis_lim, yaxis_lim, xaxis_size = _create_irregular_grid(
        xaxis, yaxis, yaxis_pos=yaxis_pos
    )

    data = spectra.fields[field_name]["data"][ind_rays, :, :]
    data = data[:, ind_rng, :]
    zpoints = _adapt_data_to_irregular_grid(
        data, xaxis_size, xaxis_lim.shape[0], xaxis_lim.shape[1]
    )

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
            field_name, field_dict=spectra.fields[field_name]
        )
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
            vmin = np.ma.min(zpoints)
        if vmax is None:
            vmax = np.ma.max(zpoints)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_angle_Doppler_title(
            spectra, field_name, ang, ind_rng, along_azi=along_azi, datetime_format=None
        )

    ax = fig.add_subplot(111)
    cax = ax.pcolor(
        xaxis_lim, yaxis_lim, zpoints, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(titl)
    if not (velmin is None or velmax is None):
        ax.set_xlim([velmin, velmax])

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label(clabel)

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def plot_time_Doppler(
    spectra,
    field_name,
    prdcfg,
    fname_list,
    xaxis_info="Doppler_velocity",
    yaxis_pos="start",
    titl=None,
    clabel=None,
    vmin=None,
    vmax=None,
    xmin=None,
    xmax=None,
    ymin=None,
    ymax=None,
):
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
    dpi = prdcfg["spectraImageConfig"].get("dpi", 72)
    xsize = prdcfg["spectraImageConfig"]["xsize"]
    ysize = prdcfg["spectraImageConfig"]["ysize"]

    if xaxis_info == "Doppler_velocity":
        xaxis = spectra.Doppler_velocity["data"]
        xlabel = "Doppler velocity (m/s)"
    elif xaxis_info == "Doppler_frequency":
        xaxis = spectra.Doppler_frequency["data"]
        xlabel = "Doppler frequency (Hz)"
    else:
        xaxis = np.ma.masked_all((spectra.nrays, spectra.npulses_max))
        for ray, npuls in enumerate(spectra.npulses["data"]):
            xaxis[ray, 0:npuls] = np.arange(npuls)
        xlabel = "Pulse number"
    ylabel = "time (s from start time)"

    xaxis_lim, yaxis_lim, xaxis_size = _create_irregular_grid(
        xaxis, spectra.time["data"], yaxis_pos=yaxis_pos
    )

    zpoints = _adapt_data_to_irregular_grid(
        np.ma.squeeze(spectra.fields[field_name]["data"], axis=1),
        xaxis_size,
        xaxis_lim.shape[0],
        xaxis_lim.shape[1],
    )

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
            field_name, field_dict=spectra.fields[field_name]
        )
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
            vmin = np.ma.min(zpoints)
        if vmax is None:
            vmax = np.ma.max(zpoints)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_complex_Doppler_title(
            spectra, field_name, 0, 0, datetime_format=None
        )

    ax = fig.add_subplot(111)
    cax = ax.pcolor(
        xaxis_lim, yaxis_lim, zpoints, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm
    )
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
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def plot_Doppler(
    spectra,
    field_name,
    ray,
    rng,
    prdcfg,
    fname_list,
    xaxis_info="Doppler_velocity",
    ylabel=None,
    titl=None,
    vmin=None,
    vmax=None,
):
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
    dpi = prdcfg["spectraImageConfig"].get("dpi", 72)
    xsize = prdcfg["spectraImageConfig"]["xsize"]
    ysize = prdcfg["spectraImageConfig"]["ysize"]
    velmin = prdcfg["spectraImageConfig"]["velmin"]
    velmax = prdcfg["spectraImageConfig"]["velmax"]

    if xaxis_info == "Doppler_velocity":
        xaxis = spectra.Doppler_velocity["data"][ray, :].compressed()
        xlabel = "Doppler velocity (m/s)"
    elif xaxis_info == "Doppler_frequency":
        xaxis = spectra.Doppler_frequency["data"][ray, :].compressed()
        xlabel = "Doppler frequency (Hz)"
    else:
        xaxis = np.arange(spectra.npulses["data"][ray])
        xlabel = "Pulse number"

    field = spectra.fields[field_name]["data"][ray, rng, 0 : xaxis.size]

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
            spectra, field_name, ray, rng, datetime_format=None
        )

    ax = fig.add_subplot(111)
    ax.plot(xaxis, field, marker="x")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=vmin, top=vmax)
    ax.set_xlim([xaxis[0], xaxis[-1]])
    ax.set_title(titl)

    if not (velmin is None or velmax is None):
        ax.set_xlim([velmin, velmax])

    # Turn on the grid
    ax.grid()

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def plot_complex_range_Doppler(
    spectra,
    field_name,
    ray,
    prdcfg,
    fname_list,
    xaxis_info="Doppler_velocity",
    titl=None,
    clabel=None,
    vmin=None,
    vmax=None,
):
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
    dpi = prdcfg["spectraImageConfig"].get("dpi", 72)
    xsize = prdcfg["spectraImageConfig"]["xsize"]
    ysize = prdcfg["spectraImageConfig"]["ysize"]
    ymin = prdcfg["spectraImageConfig"]["ymin"]
    ymax = prdcfg["spectraImageConfig"]["ymax"]
    velmin = prdcfg["spectraImageConfig"]["velmin"]
    velmax = prdcfg["spectraImageConfig"]["velmax"]

    if xaxis_info == "Doppler_velocity":
        xaxis = spectra.Doppler_velocity["data"][ray, :].compressed()
        xlabel = "Doppler velocity (m/s)"
    elif xaxis_info == "Doppler_frequency":
        xaxis = spectra.Doppler_frequency["data"][ray, :].compressed()
        xlabel = "Doppler frequency (Hz)"
    else:
        xaxis = np.arange(spectra.npulses["data"][ray])
        xlabel = "Pulse number"

    xres = np.abs(xaxis[1] - xaxis[0])

    yaxis = spectra.range["data"] / 1000.0
    yres = np.abs(yaxis[1] - yaxis[0])
    ylabel = "range (km)"

    xaxis_lim = np.append(xaxis - xres / 2, xaxis[-1] + xres / 2)
    yaxis_lim = np.append(yaxis - yres / 2, yaxis[-1] + yres / 2)

    field_2D = spectra.fields[field_name]["data"][ray, :, 0 : xaxis.size]

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
            field_name, field_dict=spectra.fields[field_name]
        )
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
            vmin = np.ma.min([np.ma.min(re_field_2D), np.ma.min(im_field_2D)])
        if vmax is None:
            vmax = np.ma.max([np.ma.max(re_field_2D), np.ma.max(im_field_2D)])

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_complex_range_Doppler_title(
            spectra, field_name, ray, datetime_format=None
        )

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    cax = ax.pcolormesh(X, Y, re_field_2D, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Real part")
    if not (velmin is None or velmax is None):
        ax.set_xlim([velmin, velmax])
    if not (ymin is None or ymax is None):
        ax.set_ylim([ymin, ymax])

    ax = fig.add_subplot(122)
    cax = ax.pcolormesh(X, Y, im_field_2D, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Imaginary part")
    if not (velmin is None or velmax is None):
        ax.set_xlim([velmin, velmax])
    if not (ymin is None or ymax is None):
        ax.set_ylim([ymin, ymax])

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
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def plot_complex_angle_Doppler(
    spectra,
    field_name,
    ang,
    ind_rays,
    ind_rng,
    prdcfg,
    fname_list,
    xaxis_info="Doppler_velocity",
    yaxis_pos="centre",
    along_azi=True,
    titl=None,
    clabel=None,
    vmin=None,
    vmax=None,
):
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
    dpi = prdcfg["spectraImageConfig"].get("dpi", 72)
    xsize = prdcfg["spectraImageConfig"]["xsize"]
    ysize = prdcfg["spectraImageConfig"]["ysize"]
    velmin = prdcfg["spectraImageConfig"]["velmin"]
    velmax = prdcfg["spectraImageConfig"]["velmax"]

    if xaxis_info == "Doppler_velocity":
        xaxis = spectra.Doppler_velocity["data"][ind_rays, :]
        xlabel = "Doppler velocity (m/s)"
    elif xaxis_info == "Doppler_frequency":
        xaxis = spectra.Doppler_frequency["data"][ind_rays, :]
        xlabel = "Doppler frequency (Hz)"
    else:
        xaxis = np.ma.masked_all((ind_rays.size, spectra.npulses_max))
        for ray, ind_ray in enumerate(ind_rays):
            npuls = spectra.npulses["data"][ind_ray]
            xaxis[ray, 0:npuls] = np.arange(npuls)
        xlabel = "Pulse number"
    ylabel = "Angle (deg)"

    if along_azi:
        yaxis = spectra.azimuth["data"][ind_rays]
    else:
        yaxis = spectra.elevation["data"][ind_rays]

    xaxis_lim, yaxis_lim, xaxis_size = _create_irregular_grid(
        xaxis, yaxis, yaxis_pos=yaxis_pos
    )

    data = spectra.fields[field_name]["data"][ind_rays, :, :]
    data = data[:, ind_rng, :]
    field_2D = _adapt_data_to_irregular_grid(
        data, xaxis_size, xaxis_lim.shape[0], xaxis_lim.shape[1]
    )

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
            field_name, field_dict=spectra.fields[field_name]
        )
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
            vmin = np.ma.min([np.ma.min(re_field_2D), np.ma.min(im_field_2D)])
        if vmax is None:
            vmax = np.ma.max([np.ma.max(re_field_2D), np.ma.max(im_field_2D)])

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_angle_Doppler_title(
            spectra, field_name, ang, ind_rng, along_azi=along_azi, datetime_format=None
        )

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    cax = ax.pcolor(
        xaxis_lim, yaxis_lim, re_field_2D, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Real part")
    if not (velmin is None or velmax is None):
        ax.set_xlim([velmin, velmax])

    ax = fig.add_subplot(122)
    cax = ax.pcolor(
        xaxis_lim, yaxis_lim, im_field_2D, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Imaginary part")
    if not (velmin is None or velmax is None):
        ax.set_xlim([velmin, velmax])

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
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def plot_complex_time_Doppler(
    spectra,
    field_name,
    prdcfg,
    fname_list,
    xaxis_info="Doppler_velocity",
    yaxis_pos="start",
    titl=None,
    clabel=None,
    vmin=None,
    vmax=None,
):
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
    dpi = prdcfg["spectraImageConfig"].get("dpi", 72)
    xsize = prdcfg["spectraImageConfig"]["xsize"]
    ysize = prdcfg["spectraImageConfig"]["ysize"]
    velmin = prdcfg["spectraImageConfig"]["velmin"]
    velmax = prdcfg["spectraImageConfig"]["velmax"]

    if xaxis_info == "Doppler_velocity":
        xaxis = spectra.Doppler_velocity["data"]
        xlabel = "Doppler velocity (m/s)"
    elif xaxis_info == "Doppler_frequency":
        xaxis = spectra.Doppler_frequency["data"]
        xlabel = "Doppler frequency (Hz)"
    else:
        xaxis = np.ma.masked_all((spectra.nrays, spectra.npulses_max))
        for ray, npuls in enumerate(spectra.npulses["data"]):
            xaxis[ray, 0:npuls] = np.arange(npuls)
        xlabel = "Pulse number"
    ylabel = "time (s from start time)"

    xaxis_lim, yaxis_lim, xaxis_size = _create_irregular_grid(
        xaxis, spectra.time["data"], yaxis_pos=yaxis_pos
    )

    field_2D = _adapt_data_to_irregular_grid(
        np.ma.squeeze(spectra.fields[field_name]["data"], axis=1),
        xaxis_size,
        xaxis_lim.shape[0],
        xaxis_lim.shape[1],
    )

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
            field_name, field_dict=spectra.fields[field_name]
        )
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
            vmin = np.ma.min([np.ma.min(re_field_2D), np.ma.min(im_field_2D)])
        if vmax is None:
            vmax = np.ma.max([np.ma.max(re_field_2D), np.ma.max(im_field_2D)])

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_complex_Doppler_title(
            spectra, field_name, 0, 0, datetime_format=None
        )

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    cax = ax.pcolor(
        xaxis_lim, yaxis_lim, re_field_2D, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Real part")
    if not (velmin is None or velmax is None):
        ax.set_xlim([velmin, velmax])

    ax = fig.add_subplot(122)
    cax = ax.pcolor(
        xaxis_lim, yaxis_lim, im_field_2D, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Imaginary part")
    if not (velmin is None or velmax is None):
        ax.set_xlim([velmin, velmax])

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
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def plot_amp_phase_range_Doppler(
    spectra,
    field_name,
    ray,
    prdcfg,
    fname_list,
    xaxis_info="Doppler_velocity",
    titl=None,
    ampli_vmin=None,
    ampli_vmax=None,
    phase_vmin=None,
    phase_vmax=None,
):
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
    dpi = prdcfg["spectraImageConfig"].get("dpi", 72)
    xsize = prdcfg["spectraImageConfig"]["xsize"]
    ysize = prdcfg["spectraImageConfig"]["ysize"]
    ymin = prdcfg["spectraImageConfig"]["ymin"]
    ymax = prdcfg["spectraImageConfig"]["ymax"]
    velmin = prdcfg["spectraImageConfig"]["velmin"]
    velmax = prdcfg["spectraImageConfig"]["velmax"]

    if xaxis_info == "Doppler_velocity":
        xaxis = spectra.Doppler_velocity["data"][ray, :].compressed()
        xlabel = "Doppler velocity (m/s)"
    elif xaxis_info == "Doppler_frequency":
        xaxis = spectra.Doppler_frequency["data"][ray, :].compressed()
        xlabel = "Doppler frequency (Hz)"
    else:
        xaxis = np.arange(spectra.npulses["data"][ray])
        xlabel = "Pulse number"
    xres = np.abs(xaxis[1] - xaxis[0])

    yaxis = spectra.range["data"] / 1000.0
    yres = np.abs(yaxis[1] - yaxis[0])
    ylabel = "range (km)"

    xaxis_lim = np.append(xaxis - xres / 2, xaxis[-1] + xres / 2)
    yaxis_lim = np.append(yaxis - yres / 2, yaxis[-1] + yres / 2)

    field_2D = spectra.fields[field_name]["data"][ray, :, 0 : xaxis.size]

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
            field_name, field_dict=spectra.fields[field_name]
        )
        if ampli_vmin is None or ampli_vmax is None:
            ampli_vmin = ampli_vmax = None
            if norm is None:  # if norm is set do not override with vmin/vmax
                ampli_vmin, ampli_vmax = pyart.config.get_field_limits(field_name)
        else:
            norm = None

        if phase_vmin is None or phase_vmax is None:
            phase_vmin = phase_vmax = None
            if norm is None:  # if norm is set do not override with vmin/vmax
                phase_vmin, phase_vmax = (-180.0, 180.0)
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
            spectra, field_name, ray, datetime_format=None
        )

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    cax = ax.pcolormesh(
        X, Y, ampli_field_2D, cmap=cmap, vmin=ampli_vmin, vmax=ampli_vmax, norm=norm
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Amplitude")
    if not (velmin is None or velmax is None):
        ax.set_xlim([velmin, velmax])
    if not (ymin is None or ymax is None):
        ax.set_ylim([ymin, ymax])

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label("Amplitude (-)")

    ax = fig.add_subplot(122)
    cax = ax.pcolormesh(
        X, Y, phase_field_2D, cmap=cmap, vmin=phase_vmin, vmax=phase_vmax, norm=norm
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Phase")
    if not (velmin is None or velmax is None):
        ax.set_xlim([velmin, velmax])
    if not (ymin is None or ymax is None):
        ax.set_ylim([ymin, ymax])

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label("Phase (deg)")

    # Make a tight layout
    fig.tight_layout()

    plt.subplots_adjust(top=0.8)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def plot_amp_phase_angle_Doppler(
    spectra,
    field_name,
    ang,
    ind_rays,
    ind_rng,
    prdcfg,
    fname_list,
    xaxis_info="Doppler_velocity",
    yaxis_pos="centre",
    along_azi=True,
    titl=None,
    ampli_vmin=None,
    ampli_vmax=None,
    phase_vmin=None,
    phase_vmax=None,
):
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
    dpi = prdcfg["spectraImageConfig"].get("dpi", 72)
    xsize = prdcfg["spectraImageConfig"]["xsize"]
    ysize = prdcfg["spectraImageConfig"]["ysize"]
    velmin = prdcfg["spectraImageConfig"]["velmin"]
    velmax = prdcfg["spectraImageConfig"]["velmax"]

    if xaxis_info == "Doppler_velocity":
        xaxis = spectra.Doppler_velocity["data"][ind_rays, :]
        xlabel = "Doppler velocity (m/s)"
    elif xaxis_info == "Doppler_frequency":
        xaxis = spectra.Doppler_frequency["data"][ind_rays, :]
        xlabel = "Doppler frequency (Hz)"
    else:
        xaxis = np.ma.masked_all((ind_rays.size, spectra.npulses_max))
        for ray, ind_ray in enumerate(ind_rays):
            npuls = spectra.npulses["data"][ind_ray]
            xaxis[ray, 0:npuls] = np.arange(npuls)
        xlabel = "Pulse number"
    ylabel = "Angle (deg)"

    if along_azi:
        yaxis = spectra.azimuth["data"][ind_rays]
    else:
        yaxis = spectra.elevation["data"][ind_rays]

    xaxis_lim, yaxis_lim, xaxis_size = _create_irregular_grid(
        xaxis, yaxis, yaxis_pos=yaxis_pos
    )

    data = spectra.fields[field_name]["data"][ind_rays, :, :]
    data = data[:, ind_rng, :]
    field_2D = _adapt_data_to_irregular_grid(
        data, xaxis_size, xaxis_lim.shape[0], xaxis_lim.shape[1]
    )

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
            field_name, field_dict=spectra.fields[field_name]
        )
        if ampli_vmin is None or ampli_vmax is None:
            ampli_vmin = ampli_vmax = None
            if norm is None:  # if norm is set do not override with vmin/vmax
                ampli_vmin, ampli_vmax = pyart.config.get_field_limits(field_name)
        else:
            norm = None

        if phase_vmin is None or phase_vmax is None:
            phase_vmin = phase_vmax = None
            if norm is None:  # if norm is set do not override with vmin/vmax
                phase_vmin, phase_vmax = (-180.0, 180.0)
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
            spectra, field_name, ang, ind_rng, along_azi=along_azi, datetime_format=None
        )

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    cax = ax.pcolor(
        xaxis_lim,
        yaxis_lim,
        ampli_field_2D,
        cmap=cmap,
        vmin=ampli_vmin,
        vmax=ampli_vmax,
        norm=norm,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Amplitude")
    if not (velmin is None or velmax is None):
        ax.set_xlim([velmin, velmax])

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label("Amplitude (-)")

    ax = fig.add_subplot(122)
    cax = ax.pcolor(
        xaxis_lim,
        yaxis_lim,
        phase_field_2D,
        cmap=cmap,
        vmin=phase_vmin,
        vmax=phase_vmax,
        norm=norm,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Phase")
    if not (velmin is None or velmax is None):
        ax.set_xlim([velmin, velmax])

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label("Phase (deg)")

    # Make a tight layout
    fig.tight_layout()

    plt.subplots_adjust(top=0.8)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def plot_amp_phase_time_Doppler(
    spectra,
    field_name,
    prdcfg,
    fname_list,
    xaxis_info="Doppler_velocity",
    yaxis_pos="start",
    titl=None,
    ampli_vmin=None,
    ampli_vmax=None,
    phase_vmin=None,
    phase_vmax=None,
):
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
    dpi = prdcfg["spectraImageConfig"].get("dpi", 72)
    xsize = prdcfg["spectraImageConfig"]["xsize"]
    ysize = prdcfg["spectraImageConfig"]["ysize"]
    velmin = prdcfg["spectraImageConfig"]["velmin"]
    velmax = prdcfg["spectraImageConfig"]["velmax"]

    if xaxis_info == "Doppler_velocity":
        xaxis = spectra.Doppler_velocity["data"]
        xlabel = "Doppler velocity (m/s)"
    elif xaxis_info == "Doppler_frequency":
        xaxis = spectra.Doppler_frequency["data"]
        xlabel = "Doppler frequency (Hz)"
    else:
        xaxis = np.ma.masked_all((spectra.nrays, spectra.npulses_max))
        for ray, npuls in enumerate(spectra.npulses["data"]):
            xaxis[ray, 0:npuls] = np.arange(npuls)
        xlabel = "Pulse number"
    ylabel = "time (s from start time)"

    xaxis_lim, yaxis_lim, xaxis_size = _create_irregular_grid(
        xaxis, spectra.time["data"], yaxis_pos=yaxis_pos
    )

    field_2D = _adapt_data_to_irregular_grid(
        np.ma.squeeze(spectra.fields[field_name]["data"], axis=1),
        xaxis_size,
        xaxis_lim.shape[0],
        xaxis_lim.shape[1],
    )

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
            field_name, field_dict=spectra.fields[field_name]
        )
        if ampli_vmin is None or ampli_vmax is None:
            ampli_vmin = ampli_vmax = None
            if norm is None:  # if norm is set do not override with vmin/vmax
                ampli_vmin, ampli_vmax = pyart.config.get_field_limits(field_name)
        else:
            norm = None

        if phase_vmin is None or phase_vmax is None:
            phase_vmin = phase_vmax = None
            if norm is None:  # if norm is set do not override with vmin/vmax
                phase_vmin, phase_vmax = (-180.0, 180.0)
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
            spectra, field_name, 0, 0, datetime_format=None
        )

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    cax = ax.pcolor(
        xaxis_lim,
        yaxis_lim,
        ampli_field_2D,
        cmap=cmap,
        vmin=ampli_vmin,
        vmax=ampli_vmax,
        norm=norm,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Amplitude")
    if not (velmin is None or velmax is None):
        ax.set_xlim([velmin, velmax])

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label("Amplitude (-)")

    ax = fig.add_subplot(122)
    cax = ax.pcolor(
        xaxis_lim,
        yaxis_lim,
        phase_field_2D,
        cmap=cmap,
        vmin=phase_vmin,
        vmax=phase_vmax,
        norm=norm,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Phase")
    if not (velmin is None or velmax is None):
        ax.set_xlim([velmin, velmax])

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label("Phase (deg)")

    # Make a tight layout
    fig.tight_layout()

    plt.subplots_adjust(top=0.8)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def plot_complex_Doppler(
    spectra,
    field_name,
    ray,
    rng,
    prdcfg,
    fname_list,
    xaxis_info="Doppler_velocity",
    ylabel=None,
    titl=None,
    vmin=None,
    vmax=None,
):
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
    dpi = prdcfg["spectraImageConfig"].get("dpi", 72)
    xsize = prdcfg["spectraImageConfig"]["xsize"]
    ysize = prdcfg["spectraImageConfig"]["ysize"]
    velmin = prdcfg["spectraImageConfig"]["velmin"]
    velmax = prdcfg["spectraImageConfig"]["velmax"]

    if xaxis_info == "Doppler_velocity":
        xaxis = spectra.Doppler_velocity["data"][ray, :].compressed()
        xlabel = "Doppler velocity (m/s)"
    elif xaxis_info == "Doppler_frequency":
        xaxis = spectra.Doppler_frequency["data"][ray, :].compressed()
        xlabel = "Doppler frequency (Hz)"
    else:
        xaxis = np.arange(spectra.npulses["data"][ray])
        xlabel = "Pulse number"

    field = spectra.fields[field_name]["data"][ray, rng, 0 : xaxis.size]

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
            spectra, field_name, ray, rng, datetime_format=None
        )

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    ax.plot(xaxis, re_field, marker="x")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=vmin, top=vmax)
    ax.set_xlim([xaxis[0], xaxis[-1]])
    ax.set_title("Real part")
    if not (velmin is None or velmax is None):
        ax.set_xlim([velmin, velmax])

    # Turn on the grid
    ax.grid()

    ax = fig.add_subplot(122)
    ax.plot(xaxis, im_field, marker="x")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=vmin, top=vmax)
    ax.set_xlim([xaxis[0], xaxis[-1]])
    ax.set_title("Imaginary part")
    if not (velmin is None or velmax is None):
        ax.set_xlim([velmin, velmax])

    # Turn on the grid
    ax.grid()

    # Make a tight layout
    fig.tight_layout()

    plt.subplots_adjust(top=0.8)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def plot_amp_phase_Doppler(
    spectra,
    field_name,
    ray,
    rng,
    prdcfg,
    fname_list,
    xaxis_info="Doppler_velocity",
    titl=None,
    ampli_vmin=None,
    ampli_vmax=None,
    phase_vmin=None,
    phase_vmax=None,
):
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
    dpi = prdcfg["spectraImageConfig"].get("dpi", 72)
    xsize = prdcfg["spectraImageConfig"]["xsize"]
    ysize = prdcfg["spectraImageConfig"]["ysize"]
    velmin = prdcfg["spectraImageConfig"]["velmin"]
    velmax = prdcfg["spectraImageConfig"]["velmax"]

    if xaxis_info == "Doppler_velocity":
        xaxis = spectra.Doppler_velocity["data"][ray, :].compressed()
        xlabel = "Doppler velocity (m/s)"
    elif xaxis_info == "Doppler_frequency":
        xaxis = spectra.Doppler_frequency["data"][ray, :].compressed()
        xlabel = "Doppler frequency (Hz)"
    else:
        xaxis = np.arange(spectra.npulses["data"][ray])
        xlabel = "Pulse number"

    field = spectra.fields[field_name]["data"][ray, rng, 0 : xaxis.size]

    ampli_field = np.ma.abs(field)
    phase_field = np.ma.angle(field, deg=True)

    # display data
    if field_name is not None:
        if ampli_vmin is None or ampli_vmax is None:
            ampli_vmin, ampli_vmax = pyart.config.get_field_limits(field_name)
        if phase_vmin is None or phase_vmax is None:
            phase_vmin, phase_vmax = (-180.0, 180.0)
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
            spectra, field_name, ray, rng, datetime_format=None
        )

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    ax.plot(xaxis, ampli_field, marker="x")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Amplitude (-)")
    ax.set_ylim(bottom=ampli_vmin, top=ampli_vmax)
    ax.set_xlim([xaxis[0], xaxis[-1]])
    ax.set_title("Amplitude")
    if not (velmin is None or velmax is None):
        ax.set_xlim([velmin, velmax])

    # Turn on the grid
    ax.grid()

    ax = fig.add_subplot(122)
    ax.plot(xaxis, phase_field, marker="x")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Phase (deg)")
    ax.set_ylim(bottom=phase_vmin, top=phase_vmax)
    ax.set_xlim([xaxis[0], xaxis[-1]])
    ax.set_title("Phase")

    # Turn on the grid
    ax.grid()

    # Make a tight layout
    fig.tight_layout()

    plt.subplots_adjust(top=0.8)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return fname_list


def _create_irregular_grid(xaxis, yaxis, yaxis_pos="start"):
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
    yres = np.median(yaxis[1:] - yaxis[:-1])
    if yaxis_pos == "start":
        yaxis_lim = np.append(yaxis, yaxis[-1] + yres)
    elif yaxis_pos == "end":
        yaxis_lim = np.append(yaxis - yres, yaxis[-1])
    else:
        yaxis_lim = np.append(yaxis - yres / 2.0, yaxis[-1] + yres / 2.0)

    # List of x-axis bin limits
    xres = np.abs(xaxis[:, 1] - xaxis[:, 0])
    xaxis_size = list()
    xaxis_list = list()
    for i in range(nbinsy):
        xaxis_aux = xaxis[i, :].compressed()
        xaxis_lim = np.append(xaxis_aux - xres[i] / 2.0, xaxis_aux[-1] + xres[i] / 2.0)
        xaxis_size.append(xaxis_lim.size)
        xaxis_list.append(xaxis_lim)

    nxbin_lim = max(xaxis_size)

    # Create two dimensional grid of y-axis
    yaxis_lim = np.repeat(yaxis_lim, 2)[1:-1]
    nybin_lim = yaxis_lim.size
    yaxis_lim = np.broadcast_to(
        np.reshape(yaxis_lim, (1, nybin_lim)), (nxbin_lim, nybin_lim)
    )

    # Create two dimensional grid of x-axis
    xaxis_lim = np.ma.masked_all((nxbin_lim, nbinsy))
    for i in range(nbinsy):
        xaxis_lim[0 : xaxis_size[i], i] = xaxis_list[i]
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
    zpoints = np.ma.masked_all((nxbin_lim - 1, nybin_lim - 1), dtype=data.dtype)
    for i in range(nbinsy - 1):
        zpoints[0 : xaxis_size[i] - 1, 2 * i] = data[i, 0 : xaxis_size[i] - 1]
        zpoints[0 : xaxis_size[i] - 1, 2 * i + 1] = data[i, 0 : xaxis_size[i] - 1]
    zpoints[0 : xaxis_size[-1] - 1, -1] = data[-1, 0 : xaxis_size[-1] - 1]

    return zpoints
