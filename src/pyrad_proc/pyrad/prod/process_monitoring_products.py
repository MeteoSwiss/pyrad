"""
pyrad.prod.process_monitoring_products
======================================

Functions for obtaining Pyrad products from monitoring datasets

.. autosummary::
    :toctree: generated/

    generate_monitoring_products

"""

from copy import deepcopy
from warnings import warn
import os

import numpy as np

import pyart

from ..io.io_aux import get_fieldname_pyart
from ..io.io_aux import get_save_dir, make_filename
from ..io.io_aux import generate_field_name_str

from ..io.read_data_other import read_monitoring_ts

from ..io.write_data import write_monitoring_ts, write_alarm_msg, send_msg
from ..io.write_data import write_histogram

from ..graph.plots import plot_histogram2, plot_density
from ..graph.plots_timeseries import plot_monitoring_ts
from ..graph.plots_aux import get_colobar_label, get_field_name

from ..util.radar_utils import compute_quantiles_from_hist


def generate_monitoring_products(dataset, prdcfg):
    """
    generates a monitoring product. With the parameter 'hist_type' the user
    may define if the product is computed for each radar volume ('instant') or
    at the end of the processing period ('cumulative'). Default is
    'cumulative'.
    Accepted product types:
        'ANGULAR_DENSITY': For a specified elevation angle, plots a 2D
            histogram with the azimuth angle in the X-axis and the data values
            in the Y-axis. The reference values and the user defined quantiles
            are also plot on the same figure
            User defined parameters:
                anglenr: int
                    The elevation angle number to plot
                quantiles: list of floats
                    The quantiles to plot. Default 25., 50., 75.
                ref_value: float
                    The reference value
                vmin, vmax : floats or None
                    The minimum and maximum values of the data points. If not
                    specified they are obtained from the Py-ART config file
        'CUMUL_VOL_TS': Plots time series of the average of instantaneous
            quantiles stored in a csv file.
            User defined parameters:
                quantiles: list of 3 floats
                    the quantiles to compute. Default 25., 50., 75.
                ref_value: float
                    The reference value. Default 0
                sort_by_date: Bool
                    If true when reading the csv file containing the
                    statistics the data is sorted by date. Default False
                rewrite: Bool
                    If true the csv file containing the statistics is
                    rewritten
                add_data_in_fname: Bool
                    If true and the data used is cumulative the year is
                    written in the csv file name and the plot file name
                npoints_min: int
                    Minimum number of points to use the data point in the
                    plotting and to send an alarm. Default 0
                vmin, vmax: float or None
                    Limits of the Y-axis (data value). If None the limits
                    are obtained from the Py-ART config file
                alarm: Bool
                    If true an alarm is sent
                tol_abs: float
                    Margin of tolerance from the reference value. If the
                    current value is above this margin an alarm is sent. If
                    the margin is not specified it is not possible to send any
                    alarm
                tol_trend: float
                    Margin of tolerance from the reference value. If the
                    trend of the last X events is above this margin an alarm
                    is sent. If the margin is not specified it is not possible
                    to send any alarm
                plot_until_year_end: Bool
                    If true will always set the xmax of the timeseries to the end of the current year
                nevents_min: int
                    Minimum number of events with sufficient points to send an
                    alarm related to the trend. If not specified it is not
                    possible to send any alarm
                sender: str
                    The mail of the alarm sender. If not specified it is not
                    possible to send any alarm
                receiver_list: list of str
                    The list of emails of the people that will receive the
                    alarm.. If not specified it is not possible to send any
                    alarm
        'PPI_HISTOGRAM': Plots a histogram of data at a particular
            elevation angle.
            User defined parameters:
                anglenr: int
                    The elevation angle number to plot
        'SAVEVOL': Saves the monitoring data in a C/F radar file. The data
            field contains histograms of data for each pair of azimuth and
            elevation angles
        'VOL_HISTOGRAM': Plots a histogram of data collected from all the
            radar volume.
            User defined parameters:
                write_data: bool
                    If true the resultant histogram is also saved in a csv
                    file. Default True.
        'VOL_TS': Computes statistics of the gathered data and writes them in
            a csv file and plots a time series of those statistics.
            User defined parameters:
                quantiles: list of 3 floats
                    the quantiles to compute. Default 25., 50., 75.
                ref_value: float
                    The reference value. Default 0
                sort_by_date: Bool
                    If true when reading the csv file containing the
                    statistics the data is sorted by date. Default False
                rewrite: Bool
                    If true the csv file containing the statistics is
                    rewritten
                add_data_in_fname: Bool
                    If true and the data used is cumulative the year is
                    written in the csv file name and the plot file name
                npoints_min: int
                    Minimum number of points to use the data point in the
                    plotting and to send an alarm. Default 0
                vmin, vmax: float or None
                    Limits of the Y-axis (data value). If None the limits
                    are obtained from the Py-ART config file
                alarm: Bool
                    If true an alarm is sent
                tol_abs: float
                    Margin of tolerance from the reference value. If the
                    current value is above this margin an alarm is sent. If
                    the margin is not specified it is not possible to send any
                    alarm
                tol_trend: float
                    Margin of tolerance from the reference value. If the
                    trend of the last X events is above this margin an alarm
                    is sent. If the margin is not specified it is not possible
                    to send any alarm
                plot_until_year_end: Bool
                    If true will always set the xmax of the timeseries to the end of the current year
                nevents_min: int
                    Minimum number of events with sufficient points to send an
                    alarm related to the trend. If not specified it is not
                    possible to send any alarm
                sender: str
                    The mail of the alarm sender. If not specified it is not
                    possible to send any alarm
                receiver_list: list of str
                    The list of emails of the people that will receive the
                    alarm.. If not specified it is not possible to send any
                    alarm

    Parameters
    ----------
    dataset : dictionary
        dictionary containing a histogram object and some metadata

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    filename : str
        the name of the file created. None otherwise

    """

    # check the type of dataset required
    hist_type = prdcfg.get('hist_type', 'cumulative')

    if dataset['hist_type'] != hist_type:
        return None

    hist_obj = dataset['hist_obj']

    dssavedir = prdcfg['dsname']
    if 'dssavename' in prdcfg:
        dssavedir = prdcfg['dssavename']

    if prdcfg['type'] == 'VOL_HISTOGRAM':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in hist_obj.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        write_data = prdcfg.get('write_data', 1)

        timeformat = '%Y%m%d'
        titl = (
            pyart.graph.common.generate_radar_time_begin(
                hist_obj).strftime('%Y-%m-%d') + '\n' +
            get_field_name(hist_obj.fields[field_name], field_name))
        if hist_type == 'instant':
            timeformat = '%Y%m%d%H%M%S'
            titl = (
                pyart.graph.common.generate_radar_time_begin(
                    hist_obj).isoformat() + 'Z' + '\n' +
                get_field_name(hist_obj.fields[field_name], field_name))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=dataset['timeinfo'])

        fname_list = make_filename(
            'histogram', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=dataset['timeinfo'], timeformat=timeformat)

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        labelx = get_colobar_label(hist_obj.fields[field_name], field_name)

        bin_centers = hist_obj.range['data']
        hist = np.sum(hist_obj.fields[field_name]['data'], axis=0)
        plot_histogram2(
            bin_centers, hist, fname_list, labelx=labelx,
            labely='Number of Samples', titl=titl)

        print('----- save to '+' '.join(fname_list))

        if write_data:
            fname = savedir+make_filename(
                'histogram', prdcfg['dstype'], prdcfg['voltype'],
                ['csv'], timeinfo=dataset['timeinfo'],
                timeformat=timeformat)[0]

            step = bin_centers[1]-bin_centers[0]
            bin_edges = np.append(
                bin_centers-step/2., bin_centers[-1]+step/2.)
            write_histogram(
                bin_edges, hist, fname, datatype=prdcfg['voltype'], step=step)
            print('----- save to '+fname)

            return fname

        return fname_list

    if prdcfg['type'] == 'PPI_HISTOGRAM':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in hist_obj.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        el_vec = np.sort(hist_obj.fixed_angle['data'])
        el = el_vec[prdcfg['anglenr']]
        ind_el = np.where(hist_obj.fixed_angle['data'] == el)[0][0]

        timeformat = '%Y%m%d'
        titl = (
            '{:.1f}'.format(el)+' Deg. ' +
            pyart.graph.common.generate_radar_time_begin(
                hist_obj).strftime('%Y-%m-%d') + '\n' +
            get_field_name(hist_obj.fields[field_name], field_name))
        if hist_type == 'instant':
            timeformat = '%Y%m%d%H%M%S'
            titl = (
                '{:.1f}'.format(el)+' Deg. ' +
                pyart.graph.common.generate_radar_time_begin(
                    hist_obj).isoformat() + 'Z' + '\n' +
                get_field_name(hist_obj.fields[field_name], field_name))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=dataset['timeinfo'])

        fname_list = make_filename(
            'ppi', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='el'+'{:.1f}'.format(el),
            timeinfo=dataset['timeinfo'], timeformat=timeformat)

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        labelx = get_colobar_label(hist_obj.fields[field_name], field_name)

        sweep_start = hist_obj.sweep_start_ray_index['data'][ind_el]
        sweep_end = hist_obj.sweep_end_ray_index['data'][ind_el]
        values = hist_obj.fields[field_name]['data'][sweep_start:sweep_end, :]
        plot_histogram2(
            hist_obj.range['data'], np.sum(values, axis=0),
            fname_list, labelx=labelx, labely='Number of Samples',
            titl=titl)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'ANGULAR_DENSITY':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in hist_obj.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        el_vec = np.sort(hist_obj.fixed_angle['data'])
        el = el_vec[prdcfg['anglenr']]
        ind_el = np.where(hist_obj.fixed_angle['data'] == el)[0][0]

        timeformat = '%Y%m%d'
        if hist_type == 'instant':
            timeformat = '%Y%m%d%H%M%S'

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=dataset['timeinfo'])

        fname_list = make_filename(
            'ppi', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='el'+'{:.1f}'.format(el),
            timeinfo=dataset['timeinfo'], timeformat=timeformat)

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        quantiles = prdcfg.get('quantiles', np.array([25., 50., 75.]))
        ref_value = prdcfg.get('ref_value', 0.)
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        plot_density(
            hist_obj, hist_type, field_name, ind_el, prdcfg, fname_list,
            quantiles=quantiles, ref_value=ref_value, vmin=vmin, vmax=vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'VOL_TS':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in hist_obj.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # put time info in file path and name
        csvtimeinfo_path = None
        csvtimeinfo_file = None
        timeformat = None
        if hist_type == 'instant':
            csvtimeinfo_path = dataset['timeinfo']
            csvtimeinfo_file = dataset['timeinfo']
            timeformat = '%Y%m%d'
        if prdcfg.get('add_date_in_fname', False):
            csvtimeinfo_file = dataset['timeinfo']
            timeformat = '%Y'

        quantiles = prdcfg.get('quantiles', np.array([25., 50., 75.]))
        ref_value = prdcfg.get('ref_value', 0.)
        sort_by_date = prdcfg.get('sort_by_date', False)
        rewrite = prdcfg.get('rewrite', False)
        plot_until_year_end = prdcfg.get('plot_until_year_end', False)        

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=csvtimeinfo_path)

        csvfname = make_filename(
            'ts', prdcfg['dstype'], prdcfg['voltype'], ['csv'],
            timeinfo=csvtimeinfo_file, timeformat=timeformat,
            runinfo=prdcfg['runinfo'])[0]

        csvfname = savedir+csvfname

        quantiles, values = compute_quantiles_from_hist(
            hist_obj.range['data'],
            np.ma.sum(hist_obj.fields[field_name]['data'], axis=0),
            quantiles=quantiles)

        start_time = pyart.graph.common.generate_radar_time_begin(hist_obj)
        np_t = np.ma.sum(hist_obj.fields[field_name]['data'], dtype=int)
        if np.ma.getmaskarray(np_t):
            np_t = 0

        write_monitoring_ts(
            start_time, np_t, values, quantiles, prdcfg['voltype'],
            csvfname)
        print('saved CSV file: '+csvfname)

        date, np_t_vec, cquant_vec, lquant_vec, hquant_vec = (
            read_monitoring_ts(csvfname, sort_by_date=sort_by_date))

        if date is None:
            warn(
                'Unable to plot time series. No valid data')
            return None

        if rewrite:
            val_vec = np.ma.asarray(
                [lquant_vec, cquant_vec, hquant_vec]).T
            write_monitoring_ts(
                date, np_t_vec, val_vec, quantiles, prdcfg['voltype'],
                csvfname, rewrite=True)

        figtimeinfo = None
        titldate = ''
        if hist_type == 'instant':
            figtimeinfo = date[0]
            titldate = date[0].strftime('%Y-%m-%d')
        else:
            titldate = (date[0].strftime('%Y%m%d')+'-' +
                        date[-1].strftime('%Y%m%d'))
            if prdcfg.get('add_date_in_fname', False):
                figtimeinfo = date[0]
                timeformat = '%Y'

        figfname_list = make_filename(
            'ts', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=figtimeinfo, timeformat=timeformat,
            runinfo=prdcfg['runinfo'])

        for i, figfname in enumerate(figfname_list):
            figfname_list[i] = savedir+figfname

        titl = (prdcfg['runinfo']+' Monitoring '+titldate)

        labely = generate_field_name_str(prdcfg['voltype'])

        np_min = prdcfg.get('npoints_min', 0)
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        plot_monitoring_ts(
            date, np_t_vec, cquant_vec, lquant_vec, hquant_vec, field_name,
            figfname_list, ref_value=ref_value, vmin=vmin, vmax=vmax,
            np_min=np_min, labelx='Time UTC', labely=labely, titl=titl,
            plot_until_year_end = plot_until_year_end)
        print('----- save to '+' '.join(figfname_list))

        # generate alarms if needed
        alarm = prdcfg.get('alarm', False)
        if not alarm:
            return figfname_list

        if 'tol_abs' not in prdcfg:
            warn('unable to send alarm. Missing tolerance on target')
            return None

        if 'tol_trend' not in prdcfg:
            warn('unable to send alarm. Missing tolerance in trend')
            return None

        if 'nevents_min' not in prdcfg:
            warn('unable to send alarm. ' +
                 'Missing minimum number of events to compute trend')
            return None

        if 'sender' not in prdcfg:
            warn('unable to send alarm. Missing email sender')
            return None
        if 'receiver_list' not in prdcfg:
            warn('unable to send alarm. Missing email receivers')
            return None

        tol_abs = prdcfg['tol_abs']
        tol_trend = prdcfg['tol_trend']
        nevents_min = prdcfg['nevents_min']
        sender = prdcfg['sender']
        receiver_list = prdcfg['receiver_list']

        np_last = np_t_vec[-1]
        value_last = cquant_vec[-1]

        if np_last < np_min:
            warn('No valid data on day '+date[-1].strftime('%d-%m-%Y'))
            return None

        # check if absolute value exceeded
        abs_exceeded = False
        if ((value_last > ref_value+tol_abs) or
                (value_last < ref_value-tol_abs)):
            warn('Value '+str(value_last)+' exceeds target '+str(ref_value) +
                 ' +/- '+str(tol_abs))
            abs_exceeded = True

        # compute trend and check if last value exceeds it
        mask = np.ma.getmaskarray(cquant_vec)
        ind = np.where(np.logical_and(
            np.logical_not(mask), np_t_vec >= np_min))[0]
        nvalid = len(ind)
        if nvalid <= nevents_min:
            warn('Not enough points to compute reliable trend')
            np_trend = 0
            value_trend = np.ma.masked
        else:
            np_trend_vec = np_t_vec[ind][-(nevents_min+1):-1]
            data_trend_vec = cquant_vec[ind][-(nevents_min+1):-1]

            np_trend = np.sum(np_trend_vec)
            value_trend = np.sum(data_trend_vec*np_trend_vec)/np_trend

        trend_exceeded = False
        if np_trend > 0:
            if ((value_last > value_trend+tol_trend) or
                    (value_last < value_trend-tol_trend)):
                warn('Value '+str(value_last)+'exceeds trend ' +
                     str(value_trend)+' +/- '+str(tol_trend))
                trend_exceeded = True

        if abs_exceeded is False and trend_exceeded is False:
            return None

        alarm_dir = savedir+'/alarms/'
        if not os.path.isdir(alarm_dir):
            os.makedirs(alarm_dir)
        alarm_fname = make_filename(
            'alarm', prdcfg['dstype'], prdcfg['voltype'], ['txt'],
            timeinfo=start_time, timeformat='%Y%m%d')[0]
        alarm_fname = alarm_dir+alarm_fname

        field_dict = pyart.config.get_metadata(field_name)
        param_name = get_field_name(field_dict, field_name)
        param_name_unit = param_name+' ['+field_dict['units']+']'

        write_alarm_msg(
            prdcfg['RadarName'][0], param_name_unit, start_time, ref_value,
            tol_abs, np_trend, value_trend, tol_trend, nevents_min, np_last,
            value_last, alarm_fname)

        print('----- saved monitoring alarm to '+alarm_fname)

        subject = ('NO REPLY: '+param_name+' monitoring alarm for radar ' +
                   prdcfg['RadarName'][0]+' on day ' +
                   start_time.strftime('%d-%m-%Y'))
        send_msg(sender, receiver_list, subject, alarm_fname)

        return alarm_fname

    if prdcfg['type'] == 'CUMUL_VOL_TS':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in hist_obj.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # put time info in file path and name
        csvtimeinfo_path = dataset['timeinfo']
        csvtimeinfo_file = dataset['timeinfo']
        timeformat = '%Y%m%d'

        quantiles = prdcfg.get('quantiles', np.array([25., 50., 75.]))
        ref_value = prdcfg.get('ref_value', 0.)
        sort_by_date = prdcfg.get('sort_by_date', False)
        rewrite = prdcfg.get('rewrite', False)
        plot_until_year_end = prdcfg.get('plot_until_year_end', False)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prod_ref'], timeinfo=csvtimeinfo_path)

        csvfname = make_filename(
            'ts', prdcfg['dstype'], prdcfg['voltype'], ['csv'],
            timeinfo=csvtimeinfo_file, timeformat=timeformat,
            runinfo=prdcfg['runinfo'])[0]

        csvfname = savedir+csvfname

        date, np_t_vec, cquant_vec, lquant_vec, hquant_vec = (
            read_monitoring_ts(csvfname))

        if date is None:
            warn(
                'Unable to plot time series. No valid data')
            return None

        cquant = np.ma.average(cquant_vec, weights=np_t_vec)
        lquant = np.ma.average(lquant_vec, weights=np_t_vec)
        hquant = np.ma.average(hquant_vec, weights=np_t_vec)
        values = np.ma.asarray([lquant, cquant, hquant])
        start_time = date[0]
        np_t = np.ma.sum(np_t_vec, dtype=int)
        if np.ma.getmaskarray(np_t):
            np_t = 0

        csvtimeinfo_path = None
        csvtimeinfo_file = None
        timeformat = None
        if prdcfg.get('add_date_in_fname', False):
            csvtimeinfo_file = dataset['timeinfo']
            timeformat = '%Y'

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=csvtimeinfo_path)

        csvfname = make_filename(
            'ts', prdcfg['dstype'], prdcfg['voltype'], ['csv'],
            timeinfo=csvtimeinfo_file, timeformat=timeformat,
            runinfo=prdcfg['runinfo'])[0]

        csvfname = savedir+csvfname

        write_monitoring_ts(
            start_time, np_t, values, quantiles, prdcfg['voltype'], csvfname)
        print('saved CSV file: '+csvfname)

        date, np_t_vec, cquant_vec, lquant_vec, hquant_vec = (
            read_monitoring_ts(csvfname, sort_by_date=sort_by_date))

        if date is None:
            warn(
                'Unable to plot time series. No valid data')
            return None

        if rewrite:
            val_vec = np.ma.asarray(
                [lquant_vec, cquant_vec, hquant_vec]).T
            write_monitoring_ts(
                date, np_t_vec, val_vec, quantiles, prdcfg['voltype'],
                csvfname, rewrite=True)

        figtimeinfo = None
        titldate = ''
        titldate = (date[0].strftime('%Y%m%d')+'-' +
                    date[-1].strftime('%Y%m%d'))
        if prdcfg.get('add_date_in_fname', False):
            figtimeinfo = date[0]
            timeformat = '%Y'

        figfname_list = make_filename(
            'ts', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=figtimeinfo, timeformat=timeformat,
            runinfo=prdcfg['runinfo'])

        for i, figfname in enumerate(figfname_list):
            figfname_list[i] = savedir+figfname

        titl = (prdcfg['runinfo']+' Monitoring '+titldate)

        labely = generate_field_name_str(prdcfg['voltype'])

        np_min = prdcfg.get('npoints_min', 0)
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        plot_monitoring_ts(
            date, np_t_vec, cquant_vec, lquant_vec, hquant_vec, field_name,
            figfname_list, ref_value=ref_value, vmin=vmin, vmax=vmax,
            np_min=np_min, labelx='Time UTC', labely=labely, titl=titl,
            plot_until_year_end = plot_until_year_end)
        print('----- save to '+' '.join(figfname_list))

        # generate alarms if needed
        alarm = prdcfg.get('alarm', 0)
        if not alarm:
            return figfname_list

        if 'tol_abs' not in prdcfg:
            warn('unable to send alarm. Missing tolerance on target')
            return None

        if 'tol_trend' not in prdcfg:
            warn('unable to send alarm. Missing tolerance in trend')
            return None

        if 'nevents_min' not in prdcfg:
            warn('unable to send alarm. ' +
                 'Missing minimum number of events to compute trend')
            return None

        if 'sender' not in prdcfg:
            warn('unable to send alarm. Missing email sender')
            return None
        if 'receiver_list' not in prdcfg:
            warn('unable to send alarm. Missing email receivers')
            return None

        tol_abs = prdcfg['tol_abs']
        tol_trend = prdcfg['tol_trend']
        nevents_min = prdcfg['nevents_min']
        sender = prdcfg['sender']
        receiver_list = prdcfg['receiver_list']

        np_last = np_t_vec[-1]
        value_last = cquant_vec[-1]

        if np_last < np_min:
            warn('No valid data on day '+date[-1].strftime('%d-%m-%Y'))
            return None

        # check if absolute value exceeded
        abs_exceeded = False
        if ((value_last > ref_value+tol_abs) or
                (value_last < ref_value-tol_abs)):
            warn('Value '+str(value_last)+' exceeds target '+str(ref_value) +
                 ' +/- '+str(tol_abs))
            abs_exceeded = True

        # compute trend and check if last value exceeds it
        mask = np.ma.getmaskarray(cquant_vec)
        ind = np.where(np.logical_and(
            np.logical_not(mask), np_t_vec >= np_min))[0]
        nvalid = len(ind)
        if nvalid <= nevents_min:
            warn('Not enough points to compute reliable trend')
            np_trend = 0
            value_trend = np.ma.masked
        else:
            np_trend_vec = np_t_vec[ind][-(nevents_min+1):-1]
            data_trend_vec = cquant_vec[ind][-(nevents_min+1):-1]

            np_trend = np.sum(np_trend_vec)
            value_trend = np.sum(data_trend_vec*np_trend_vec)/np_trend

        trend_exceeded = False
        if np_trend > 0:
            if ((value_last > value_trend+tol_trend) or
                    (value_last < value_trend-tol_trend)):
                warn('Value '+str(value_last)+'exceeds trend ' +
                     str(value_trend)+' +/- '+str(tol_trend))
                trend_exceeded = True

        if abs_exceeded is False and trend_exceeded is False:
            return None

        alarm_dir = savedir+'/alarms/'
        if not os.path.isdir(alarm_dir):
            os.makedirs(alarm_dir)
        alarm_fname = make_filename(
            'alarm', prdcfg['dstype'], prdcfg['voltype'], ['txt'],
            timeinfo=start_time, timeformat='%Y%m%d')[0]
        alarm_fname = alarm_dir+alarm_fname

        field_dict = pyart.config.get_metadata(field_name)
        param_name = get_field_name(field_dict, field_name)
        param_name_unit = param_name+' ['+field_dict['units']+']'

        write_alarm_msg(
            prdcfg['RadarName'][0], param_name_unit, start_time, ref_value,
            tol_abs, np_trend, value_trend, tol_trend, nevents_min, np_last,
            value_last, alarm_fname)

        print('----- saved monitoring alarm to '+alarm_fname)

        subject = ('NO REPLY: '+param_name+' monitoring alarm for radar ' +
                   prdcfg['RadarName'][0]+' on day ' +
                   start_time.strftime('%d-%m-%Y'))
        send_msg(sender, receiver_list, subject, alarm_fname)

        return alarm_fname

    if prdcfg['type'] == 'SAVEVOL':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in hist_obj.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        new_dataset = deepcopy(hist_obj)
        new_dataset.fields = dict()
        new_dataset.add_field(field_name, hist_obj.fields[field_name])

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=dataset['timeinfo'])

        fname = make_filename(
            'savevol', prdcfg['dstype'], prdcfg['voltype'], ['nc'],
            timeinfo=dataset['timeinfo'])[0]

        fname = savedir+fname

        pyart.io.cfradial.write_cfradial(fname, new_dataset)
        print('saved file: '+fname)

        return fname

    warn(' Unsupported product type: ' + prdcfg['type'])
    return None
