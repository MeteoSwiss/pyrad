"""
pyrad.graph.plot_timeseries
===========================

Functions to plot Pyrad datasets

.. autosummary::
    :toctree: generated/

    plot_timeseries
    plot_timeseries_comp
    plot_monitoring_ts
    plot_intercomp_scores_ts
    plot_ml_ts
    plot_sun_retrieval_ts

"""

from warnings import warn

import numpy as np
import datetime

import matplotlib as mpl
mpl.use('Agg')

# Increase a bit font size
mpl.rcParams.update({'font.size': 16})
mpl.rcParams.update({'font.family':  "sans-serif"})

import matplotlib.dates as mdates
import matplotlib.pyplot as plt

import pyart


def plot_timeseries(tvec, data_list, fname_list, labelx='Time [UTC]',
                    labely='Value', labels=['Sensor'], title='Time Series',
                    period=0, timeformat=None, colors=None, linestyles=None,
                    markers=None, ymin=None, ymax=None, dpi=72):
    """
    plots a time series

    Parameters
    ----------
    tvec : datetime object
        time of the time series
    data_list : list of float array
        values of the time series
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
    period : float
        measurement period in seconds used to compute accumulation. If 0 no
        accumulation is computed
    timeformat : str
        Specifies the tvec and time format on the x axis
    colors : array of str
        Specifies the colors of each line
    linestyles : array of str
        Specifies the line style of each line
    markers: array of str
        Specify the markers to be used for each line
    ymin, ymax: float
        Lower/Upper limit of y axis
    dpi : int
        dots per inch

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    History
    --------
    201?.??.?? -fvj- creation
    2017.08.21 -jgr- modified margins and grid + minor graphical updates
    2018.03.05 -jgr- added x-limit of x axis to avoid unwanted error messages

    """
    if period > 0:
        for i, data in enumerate(data_list):
            data *= (period/3600.)
            data_list[i] = np.ma.cumsum(data)

    fig, ax = plt.subplots(figsize=[10, 6], dpi=dpi)

    lab = None
    col = None
    lstyle = '--'
    marker = 'o'

    for i, data in enumerate(data_list):
        if labels is not None:
            lab = labels[i]
        if colors is not None:
            col = colors[i]
        if linestyles is not None:
            lstyle = linestyles[i]
        if markers is not None:
            marker = markers[i]
        ax.plot(tvec, data, label=lab, color=col, linestyle=lstyle,
                marker=marker)

    ax.set_title(title)
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_ylim(bottom=ymin, top=ymax)
    ax.set_xlim([tvec[0], tvec[-1]])

    # Turn on the grid
    ax.grid()

    if timeformat is not None:
        ax.xaxis.set_major_formatter(mdates.DateFormatter(timeformat))

    # rotates and right aligns the x labels, and moves the bottom of the
    # axes up to make room for them
    fig.autofmt_xdate()

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_timeseries_comp(date1, value1, date2, value2, fname_list,
                         labelx='Time [UTC]', labely='Value',
                         label1='Sensor 1', label2='Sensor 2',
                         titl='Time Series Comparison', period1=0, period2=0,
                         ymin=None, ymax=None, dpi=72):
    """
    plots 2 time series in the same graph

    Parameters
    ----------
    date1 : datetime object
        time of the first time series
    value1 : float array
        values of the first time series
    date2 : datetime object
        time of the second time series
    value2 : float array
        values of the second time series
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    label1, label2 : str
        legend label for each time series
    titl : str
        The figure title
     period1, period2 : float
        measurement period in seconds used to compute accumulation. If 0 no
        accumulation is computed
    dpi : int
        dots per inch
    ymin, ymax : float
        The limits of the Y-axis. None will keep the default limit.

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    History
    --------
    201?.??.?? -fvj- created
    2017.08.21 -jgr- changed some graphical aspects

    """
    if (period1 > 0) and (period2 > 0):
        # TODO: document this and check (sometimes artefacts)
        value1 *= (period1/3600.)
        value1 = np.ma.cumsum(value1)

        value2 *= (period2/3600.)
        value2 = np.ma.cumsum(value2)

    fig, ax = plt.subplots(figsize=[10, 6.5], dpi=dpi)
    ax.plot(date1, value1, 'b', label=label1, linestyle='--', marker='o')
    ax.plot(date2, value2, 'r', label=label2, linestyle='--', marker='s')
    ax.legend(loc='best')
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_title(titl)

    ax.grid()

    ax.set_ylim(bottom=ymin, top=ymax)
    ax.set_xlim([date2[0], date2[-1]])

    # rotates and right aligns the x labels, and moves the bottom of the
    # axes up to make room for them
    fig.autofmt_xdate()

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_monitoring_ts(date, np_t, cquant, lquant, hquant, field_name,
                       fname_list, ref_value=None, vmin=None, vmax=None,
                       np_min=0, labelx='Time [UTC]', labely='Value',
                       titl='Time Series', dpi=72, plot_until_year_end = False):
    """
    plots a time series of monitoring data

    Parameters
    ----------
    date : datetime object
        time of the time series
    np_t : int array
        number of points
    cquant, lquant, hquant : float array
        values of the central, low and high quantiles
    field_name : str
        name of the field
    fname_list : list of str
        list of names of the files where to store the plot
    ref_value : float
        the reference value
    vmin, vmax : float
        The limits of the y axis
    np_min : int
        minimum number of points to consider the sample plotable
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    titl : str
        The figure title
    dpi : int
        dots per inch
    plot_until_year_end : bool
        if true will always set the maximum xtick to the end of the year
    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    vmin_pyart, vmax_pyart = pyart.config.get_field_limits(field_name)
    if vmin is None:
        vmin = vmin_pyart
    if vmax is None:
        vmax = vmax_pyart

    # plot only valid data (but keep first and last date)
    date2 = np.array(date)
    isvalid = np.logical_not(np.ma.getmaskarray(cquant))
    if np_min > 0:
        has_np = np_t > np_min
        isvalid = np.logical_and(isvalid, has_np)

    cquant_plt = cquant[isvalid]
    lquant_plt = lquant[isvalid]
    hquant_plt = hquant[isvalid]
    date_plt = date2[isvalid]
    if not isvalid[0]:
        cquant_plt = np.ma.append(np.ma.masked, cquant_plt)
        lquant_plt = np.ma.append(np.ma.masked, lquant_plt)
        hquant_plt = np.ma.append(np.ma.masked, hquant_plt)
        date_plt = np.ma.append(date2[0], date_plt)
    if not isvalid[-1]:
        cquant_plt = np.ma.append(cquant_plt, np.ma.masked)
        lquant_plt = np.ma.append(lquant_plt, np.ma.masked)
        hquant_plt = np.ma.append(hquant_plt, np.ma.masked)
        date_plt = np.ma.append(date_plt, date2[-1])

    fig = plt.figure(figsize=[15, 13], dpi=dpi)

    ax = fig.add_subplot(2, 1, 1)
    ax.plot(date_plt, cquant_plt, 'x-')
    ax.plot(date_plt, lquant_plt, 'rx-')
    ax.plot(date_plt, hquant_plt, 'rx-')
    if ref_value is not None:
        ax.plot(date_plt, np.zeros(len(date_plt))+ref_value, 'k--')
    ax.set_ylabel(labely)
    ax.set_title(titl)
    ax.set_ylim([vmin, vmax])
    if plot_until_year_end:
        t0 = date_plt[0]
        tend = datetime.datetime(year = t0.year, 
                                 month = 12,
                                 day = 31, 
                                 hour = 23, 
                                 minute = 59)
        ax.set_xlim([t0, tend])
    else:
        # tight x axis
        ax.autoscale(enable=True, axis='x', tight=True)
    
    ax.grid(True)

    ax = fig.add_subplot(2, 1, 2)
    ax.plot(date, np_t, 'x-')

    if np_min is not None:
        ax.plot(date, np.zeros(len(date))+np_min, 'k--')

    if plot_until_year_end:
         ax.set_xlim([t0, tend])
    else:
        # tight x axis
        ax.autoscale(enable=True, axis='x', tight=True)

    ax.set_ylabel('Number of Samples')
    ax.set_xlabel(labelx)

    # rotates and right aligns the x labels, and moves the bottom of the
    # axes up to make room for them
    fig.autofmt_xdate()

    
    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_intercomp_scores_ts(date_vec, np_vec, meanbias_vec, medianbias_vec,
                             quant25bias_vec, quant75bias_vec, modebias_vec,
                             corr_vec, slope_vec, intercep_vec,
                             intercep_slope1_vec, fname_list, ref_value=0.,
                             np_min=0, corr_min=0.,
                             labelx='Time UTC',
                             titl='RADAR001-RADAR002 intercomparison',
                             dpi=72):
    """
    plots a time series of radar intercomparison scores

    Parameters
    ----------
    date_vec : datetime object
        time of the time series
    np_vec : int array
        number of points
    meanbias_vec, medianbias_vec, modebias_vec : float array
        mean, median and mode bias
    quant25bias_vec, quant75bias_vec: 25th and 75th percentile of the bias
    corr_vec : float array
        correlation
    slope_vec, intercep_vec : float array
        slope and intercep of a linear regression
    intercep_slope1_vec : float
        the intercep point of a inear regression of slope 1
    ref_value : float
        the reference value
    np_min : int
        The minimum number of points to consider the result valid
    corr_min : float
        The minimum correlation to consider the results valid
    labelx : str
        The label of the X axis
    titl : str
        The figure title

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    # plot only valid data (but keep first and last date)
    date2 = np.array(date_vec)
    isvalid = np.logical_not(np.ma.getmaskarray(meanbias_vec))
    isvalid_corr = np.logical_not(np.ma.getmaskarray(corr_vec))
    if np_min > 0:
        has_np = np_vec > np_min
        isvalid = np.logical_and(isvalid, has_np)
    if corr_min > 0:
        has_corr_min = corr_vec > corr_min
        isvalid = np.logical_and(isvalid, has_corr_min)

    meanbias_plt = meanbias_vec[isvalid]
    medianbias_plt = medianbias_vec[isvalid]
    quant25bias_plt = quant25bias_vec[isvalid]
    quant75bias_plt = quant75bias_vec[isvalid]
    modebias_plt = modebias_vec[isvalid]
    intercep_plt = intercep_slope1_vec[isvalid]
    corr_plt = corr_vec[isvalid_corr]
    date_corr = date2[isvalid_corr]
    date_plt = date2[isvalid]
    if not isvalid[0]:
        meanbias_plt = np.ma.append(np.ma.masked, meanbias_plt)
        medianbias_plt = np.ma.append(np.ma.masked, medianbias_plt)
        quant25bias_plt = np.ma.append(np.ma.masked, quant25bias_plt)
        quant75bias_plt = np.ma.append(np.ma.masked, quant75bias_plt)
        modebias_plt = np.ma.append(np.ma.masked, modebias_plt)
        intercep_plt = np.ma.append(np.ma.masked, intercep_plt)
        date_plt = np.ma.append(date2[0], date_plt)
    if not isvalid[-1]:
        meanbias_plt = np.ma.append(meanbias_plt, np.ma.masked)
        medianbias_plt = np.ma.append(medianbias_plt, np.ma.masked)
        quant25bias_plt = np.ma.append(quant25bias_plt, np.ma.masked)
        quant75bias_plt = np.ma.append(quant75bias_plt, np.ma.masked)
        modebias_plt = np.ma.append(modebias_plt, np.ma.masked)
        intercep_plt = np.ma.append(intercep_plt, np.ma.masked)
        date_plt = np.ma.append(date_plt, date2[-1])

    if not isvalid_corr[0]:
        corr_plt = np.ma.append(np.ma.masked, corr_plt)
        date_corr = np.ma.append(date2[0], date_corr)
    if not isvalid_corr[-1]:
        corr_plt = np.ma.append(corr_plt, np.ma.masked)
        date_corr = np.ma.append(date_corr, date2[-1])

    fig = plt.figure(figsize=[10, 20], dpi=dpi)

    ax = fig.add_subplot(4, 1, 1)
    ax.plot(date_plt, medianbias_plt, 'bx-', label='median')
    ax.plot(date_plt, meanbias_plt, 'rx-', label='mean')
    ax.plot(date_plt, modebias_plt, 'gx-', label='mode')
    ax.plot(date_plt, intercep_plt, 'yx-', label='intercep of slope 1 LR')
    if ref_value is not None:
        ax.plot(date_plt, np.zeros(len(date_plt))+ref_value, 'k--')
    # plt.legend(loc='best')
    ax.set_ylabel('bias [dB]')
    ax.set_title(titl)
    ax.set_ylim([-5., 5.])

    # tight x axis
    ax.autoscale(enable=True, axis='x', tight=True)
    ax.grid(True)

    ax = fig.add_subplot(4, 1, 2)
    ax.plot(date_plt, medianbias_plt, 'bx-', label='median')
    ax.plot(date_plt, quant25bias_plt, 'rx-', label='25-percentile')
    ax.plot(date_plt, quant75bias_plt, 'rx-', label='75-percentile')
    if ref_value is not None:
        ax.plot(date_plt, np.zeros(len(date_plt))+ref_value, 'k--')
    # plt.legend(loc='best')
    ax.set_ylabel('bias [dB]')
    ax.set_ylim([-5., 5.])

    # tight x axis
    ax.autoscale(enable=True, axis='x', tight=True)
    ax.grid(True)

    ax = fig.add_subplot(4, 1, 3)
    ax.plot(date_corr, corr_plt, 'bx-')

    if corr_min > 0:
        ax.plot(date_corr, np.zeros(len(date_corr))+corr_min, 'k--')

    ax.set_ylabel('correlation')
    ax.set_ylim([0., 1.])

    # tight x axis
    ax.autoscale(enable=True, axis='x', tight=True)
    ax.grid(True)

    ax = fig.add_subplot(4, 1, 4)
    ax.plot(date2, np_vec, 'bx-')

    if np_min > 0:
        ax.plot(date2, np.zeros(len(date2))+np_min, 'k--')

    ax.set_ylabel('Number of Samples')
    ax.set_xlabel(labelx)

    # tight x axis
    ax.autoscale(enable=True, axis='x', tight=True)

    # rotates and right aligns the x labels, and moves the bottom of the
    # axes up to make room for them
    fig.autofmt_xdate()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_ml_ts(dt_ml_arr, ml_top_avg_arr, ml_top_std_arr, thick_avg_arr,
               thick_std_arr, nrays_valid_arr, nrays_total_arr, fname_list,
               labelx='Time UTC', titl='Melting layer time series', dpi=72):
    """
    plots a time series of melting layer data

    Parameters
    ----------
    dt_ml_arr : datetime object
        time of the time series
    np_vec : int array
        number of points
    meanbias_vec, medianbias_vec, modebias_vec : float array
        mean, median and mode bias
    quant25bias_vec, quant75bias_vec: 25th and 75th percentile of the bias
    corr_vec : float array
        correlation
    slope_vec, intercep_vec : float array
        slope and intercep of a linear regression
    intercep_slope1_vec : float
        the intercep point of a inear regression of slope 1
    ref_value : float
        the reference value
    np_min : int
        The minimum number of points to consider the result valid
    corr_min : float
        The minimum correlation to consider the results valid
    labelx : str
        The label of the X axis
    titl : str
        The figure title

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    fig = plt.figure(figsize=[10, 15], dpi=dpi)

    ax = fig.add_subplot(3, 1, 1)
    ax.plot(dt_ml_arr, ml_top_avg_arr, 'bx-', label='avg')
    ax.plot(dt_ml_arr, ml_top_avg_arr+ml_top_std_arr, 'rx-', label='avg+std')
    ax.plot(dt_ml_arr, ml_top_avg_arr-ml_top_std_arr, 'rx-', label='avg-std')
    # plt.legend(loc='best')
    ax.set_ylabel('Top height [m MSL]')
    ax.set_title(titl)
    ax.set_ylim([0., 6000.])
    ax.set_xlim([dt_ml_arr[0], dt_ml_arr[-1]])

    # tight x axis
    ax.autoscale(enable=True, axis='x', tight=True)
    ax.grid(True)

    ax = fig.add_subplot(3, 1, 2)
    ax.plot(dt_ml_arr, thick_avg_arr, 'bx-', label='avg')
    ax.plot(dt_ml_arr, thick_avg_arr+thick_std_arr, 'rx-', label='avg+std')
    ax.plot(dt_ml_arr, thick_avg_arr-thick_std_arr, 'rx-', label='avg-std')
    # plt.legend(loc='best')
    ax.set_ylabel('Thickness [m]')
    ax.set_ylim([0., 3000.])
    ax.set_xlim([dt_ml_arr[0], dt_ml_arr[-1]])

    # tight x axis
    ax.autoscale(enable=True, axis='x', tight=True)
    ax.grid(True)

    ax = fig.add_subplot(3, 1, 3)
    ax.plot(dt_ml_arr, nrays_valid_arr, 'bx-', label='N valid rays')
    ax.plot(dt_ml_arr, nrays_total_arr, 'rx-', label='rays total')
    # plt.legend(loc='best')
    ax.set_ylabel('Rays')
    ax.set_xlabel(labelx)
    ax.set_ylim([0, np.max(nrays_total_arr)+5])
    ax.set_xlim([dt_ml_arr[0], dt_ml_arr[-1]])

    # tight x axis
    ax.autoscale(enable=True, axis='x', tight=True)
    ax.grid(True)

    # rotates and right aligns the x labels, and moves the bottom of the
    # axes up to make room for them
    fig.autofmt_xdate()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_sun_retrieval_ts(sun_retrieval, data_type, fname_list, labelx='Date',
                          titl='Sun retrieval Time Series', dpi=72):
    """
    plots sun retrieval time series series

    Parameters
    ----------
    sun_retrieval : tuple
        tuple containing the retrieved parameters
    data_type : str
        parameter to be plotted
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        the x label
    titl : str
        the title of the plot
    dpi : int
        dots per inch

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    value_std = None
    ref = None
    date = sun_retrieval[1]
    if data_type == 'nhits_h':
        value = sun_retrieval[2]
        labely = 'Number of sun hits H channel'
        vmin = 0
        vmax = np.max(sun_retrieval[2])+1
    elif data_type == 'el_width_h':
        value = sun_retrieval[3]
        labely = 'Elevation beamwidth H channel (Deg)'
        vmin = 0.
        vmax = 4.
    elif data_type == 'az_width_h':
        value = sun_retrieval[4]
        labely = 'Azimuth beamwidth H channel (Deg)'
        vmin = 0.
        vmax = 4.
    elif data_type == 'el_bias_h':
        value = sun_retrieval[5]
        ref = np.zeros(len(value))
        labely = 'Elevation pointing bias H channel (Deg)'
        vmin = -2.
        vmax = 2.
    elif data_type == 'az_bias_h':
        value = sun_retrieval[6]
        ref = np.zeros(len(value))
        labely = 'Azimuth pointing bias H channel (Deg)'
        vmin = -2.
        vmax = 2.
    elif data_type == 'dBm_sun_est':
        value = sun_retrieval[7]
        value_std = sun_retrieval[8]
        labely = 'Sun Power H channel (dBm)'
        vmin = -110.
        vmax = -90.
    elif data_type == 'rx_bias_h':
        value = (10.*np.ma.log10(sun_retrieval[9]) -
                 10.*np.ma.log10(sun_retrieval[21]))
        value_std = sun_retrieval[8]
        ref = np.zeros(len(value))
        labely = 'Receiver bias H channel (dB)'
        vmin = -5.
        vmax = 5.
    elif data_type == 'sf_h':
        value = 10.*np.ma.log10(sun_retrieval[9])
        # value_std = sun_retrieval[8]
        ref = 10.*np.ma.log10(sun_retrieval[21])
        labely = 'Observed solar flux H channel (dB(sfu))'
        vmin = 15.
        vmax = 30.
    elif data_type == 'nhits_v':
        value = sun_retrieval[10]
        labely = 'Number of sun hits V channel'
        vmin = 0
        vmax = np.max(sun_retrieval[10])+1
    elif data_type == 'el_width_v':
        value = sun_retrieval[11]
        labely = 'Elevation beamwidth V channel (Deg)'
        vmin = 0.
        vmax = 4.
    elif data_type == 'az_width_v':
        value = sun_retrieval[12]
        labely = 'Azimuth beamwidth V channel (Deg)'
        vmin = 0.
        vmax = 4.
    elif data_type == 'el_bias_v':
        value = sun_retrieval[13]
        ref = np.zeros(len(value))
        labely = 'Elevation pointing bias V channel (Deg)'
        vmin = -2.
        vmax = 2.
    elif data_type == 'az_bias_v':
        value = sun_retrieval[14]
        ref = np.zeros(len(value))
        labely = 'Azimuth pointing bias V channel (Deg)'
        vmin = -2.
        vmax = 2.
    elif data_type == 'dBmv_sun_est':
        value = sun_retrieval[15]
        value_std = sun_retrieval[16]
        labely = 'Sun Power V channel (dBm)'
        vmin = -110.
        vmax = -90.
    elif data_type == 'rx_bias_v':
        value = (10.*np.ma.log10(sun_retrieval[17]) -
                 10.*np.ma.log10(sun_retrieval[21]))
        value_std = sun_retrieval[16]
        ref = np.zeros(len(value))
        labely = 'Receiver bias V channel (dB)'
        vmin = -5.
        vmax = 5.
    elif data_type == 'sf_v':
        value = 10.*np.ma.log10(sun_retrieval[17])
        # value_std = sun_retrieval[16]
        ref = 10.*np.ma.log10(sun_retrieval[21])
        labely = 'Observed solar flux V channel (dB(sfu))'
        vmin = 15.
        vmax = 30.
    elif data_type == 'nhits_zdr':
        value = sun_retrieval[18]
        labely = 'Number of sun hits ZDR'
        vmin = 0
        vmax = np.max(sun_retrieval[18])+1
    elif data_type == 'ZDR_sun_est':
        value = sun_retrieval[19]
        value_std = sun_retrieval[20]
        ref = np.zeros(len(value))
        labely = 'Sun ZDR (dB)'
        vmin = -2.
        vmax = 2.

    mask = np.ma.getmaskarray(value)
    if mask.all():
        warn('Unable to create figure '+' '.join(fname_list) +
             '. No valid data')
        return None

    # plot only valid data (but keep first and last date)
    isvalid = np.logical_not(mask)
    date2 = np.array(date)

    value_plt = value[isvalid]
    date_plt = date2[isvalid]
    if not isvalid[0]:
        value_plt = np.ma.append(np.ma.masked, value_plt)
        date_plt = np.ma.append(date2[0], date_plt)
    if not isvalid[-1]:
        value_plt = np.ma.append(value_plt, np.ma.masked)
        date_plt = np.ma.append(date_plt, date2[-1])

    fig, ax = plt.subplots(figsize=[10, 6], dpi=dpi)
    ax.plot(date_plt, value_plt, 'x-')
    if value_std is not None:
        value_std_plt = value_std[isvalid]
        if not isvalid[0]:
            value_std_plt = np.ma.append(np.ma.masked, value_std_plt)
        if not isvalid[-1]:
            value_std_plt = np.ma.append(value_std_plt, np.ma.masked)

        ax.plot(date_plt, value_plt+value_std_plt, 'rx-')
        ax.plot(date_plt, value_plt-value_std_plt, 'rx-')
    if ref is not None:
        ref_plt = ref[isvalid]
        if not isvalid[0]:
            ref_plt = np.ma.append(ref[0], ref_plt)
        if not isvalid[-1]:
            ref_plt = np.ma.append(ref_plt, ref[-1])
        ax.plot(date_plt, ref_plt, 'k--')
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_title(titl)
    ax.set_ylim([vmin, vmax])
    ax.set_xlim([date_plt[0], date_plt[-1]])
    # tight x axis
    ax.autoscale(enable=True, axis='x', tight=True)
    ax.grid(True)

    # rotates and right aligns the x labels, and moves the bottom of the
    # axes up to make room for them
    fig.autofmt_xdate()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list
