"""
pyrad.prod.process_timeseries_products
======================================

Functions for obtaining Pyrad products from a time series datasets

.. autosummary::
    :toctree: generated/

    generate_timeseries_products

"""

from copy import deepcopy
from warnings import warn

import numpy as np
from netCDF4 import num2date

from ..io.io_aux import get_save_dir, make_filename, get_fieldname_pyart
from ..io.io_aux import generate_field_name_str

from ..io.read_data_sensor import get_sensor_data
from ..io.read_data_other import read_timeseries

from ..io.write_data import write_ts_polar_data, write_ts_cum
from ..io.write_data import write_ts_grid_data, write_multiple_points

from ..graph.plots_timeseries import plot_timeseries, plot_timeseries_comp
from ..graph.plots_vol import plot_cappi, plot_traj
from ..graph.plots import plot_scatter_comp

from ..util.radar_utils import rainfall_accumulation


def generate_timeseries_products(dataset, prdcfg):
    """
    Generates time series products. Accepted product types:
        'COMPARE_CUMULATIVE_POINT': Plots in the same graph 2 time series of
            data accumulation (tipically rainfall rate). One time series is
            a point measurement of radar data while the other is from a
            co-located instrument (rain gauge or disdrometer)
            User defined parameters:
                dpi: int
                    The pixel density of the plot. Default 72
                vmin, vmax: float
                    The limits of the Y-axis. If none they will be obtained
                    from the Py-ART config file.
                sensor: str
                    The sensor type. Can be 'rgage' or 'disdro'
                sensorid: str
                    The sensor ID.
                location: str
                    A string identifying the location of the disdrometer
                freq: float
                    The frequency used to retrieve the polarimetric variables
                    of a disdrometer
                ele: float
                    The elevation angle used to retrieve the polarimetric
                    variables of a disdrometer
                ScanPeriod: float
                    The scaning period of the radar in seconds. This parameter
                    is defined in the 'loc' config file
        'COMPARE_POINT': Plots in the same graph 2 time series of
            data . One time series is a point measurement of radar data while
            the other is from a co-located instrument (rain gauge or
            disdrometer)
            User defined parameters:
                dpi: int
                    The pixel density of the plot. Default 72
                vmin, vmax: float
                    The limits of the Y-axis. If none they will be obtained
                    from the Py-ART config file.
                sensor: str
                    The sensor type. Can be 'rgage' or 'disdro'
                sensorid: str
                    The sensor ID.
                location: str
                    A string identifying the location of the disdrometer
                freq: float
                    The frequency used to retrieve the polarimetric variables
                    of a disdrometer
                ele: float
                    The elevation angle used to retrieve the polarimetric
                    variables of a disdrometer
        'COMPARE_TIME_AVG': Creates a scatter plot of average radar data
            versus average sensor data.
            User defined parameters:
                dpi: int
                    The pixel density of the plot. Default 72
                sensor: str
                    The sensor type. Can be 'rgage' or 'disdro'
                sensorid: str
                    The sensor ID.
                location: str
                    A string identifying the location of the disdrometer
                freq: float
                    The frequency used to retrieve the polarimetric variables
                    of a disdrometer
                ele: float
                    The elevation angle used to retrieve the polarimetric
                    variables of a disdrometer
                cum_time: float
                    Data accumulation time [s]. Default 3600.
                base_time: float
                    Starting moment of the accumulation [s from midnight].
                    Default 0.
        'PLOT_AND_WRITE': Writes and plots a trajectory time series.
            User defined parameters:
                ymin, ymax: float
                    The minimum and maximum value of the Y-axis. If none it
                    will be obtained from the Py-ART config file.
        'PLOT_AND_WRITE_POINT': Plots and writes a time series of radar data
            at a particular point
            User defined parameters:
                dpi: int
                    The pixel density of the plot. Default 72
                vmin, vmax: float
                    The limits of the Y-axis. If none they will be obtained
                    from the Py-ART config file.
        'PLOT_CUMULATIVE_POINT': Plots a time series of radar data
            accumulation at a particular point.
            User defined parameters:
                dpi: int
                    The pixel density of the plot. Default 72
                vmin, vmax: float
                    The limits of the Y-axis. If none they will be obtained
                    from the Py-ART config file.
                ScanPeriod: float
                    The scaning period of the radar in seconds. This parameter
                    is defined in the 'loc' config file
        'PLOT_HIST': plots and writes a histogram of all the data gathered
            during the trajectory processing
            User defined parameters:
                step: float or None
                    The quantization step of the data. If None it will be
                    obtained from the Py-ART config file
        'TRAJ_CAPPI_IMAGE': Creates a CAPPI image with the trajectory position
            overplot on it.
            User defined parameters:
                color_ref: str
                    The meaning of the color code with which the trajectory is
                    plotted. Can be 'None', 'altitude' (the absolute
                    altitude), 'rel_altitude' (altitude relative to the CAPPI
                    altitude), 'time' (trajectory time respect of the start of
                    the radar scan leading to the CAPPI)
                altitude: float
                    The CAPPI altitude [m]
                wfunc: str
                    Function used in the gridding of the radar data. The
                    function types are defined in pyart.map.grid_from_radars.
                    Default 'NEAREST'
                res: float
                    The CAPPI resolution [m]. Default 500.

    Parameters
    ----------
    dataset : dictionary
        radar object

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    no return

    """
    dssavedir = prdcfg['dsname']
    if 'dssavename' in prdcfg:
        dssavedir = prdcfg['dssavename']

    prdsavedir = prdcfg['prdname']
    if 'prdsavedir' in prdcfg:
        prdsavedir = prdcfg['prdsavedir']

    if prdcfg['type'] == 'WRITE_MULTIPLE_POINTS':

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir, prdsavedir,
            timeinfo=dataset['ref_time'])

        csvfname = make_filename(
            'POI', prdcfg['dstype'], dataset['datatype'], ['csv'],
            timeinfo=dataset['ref_time'])[0]

        csvfname = savedir+csvfname
        write_multiple_points(dataset, csvfname)
        print('saved CSV file: '+csvfname)

        return csvfname

    if prdcfg['type'] == 'PLOT_AND_WRITE_POINT':
        set_time_info = prdcfg.get('set_time_info', True)
        if 'antenna_coordinates_az_el_r' in dataset:
            az = '{:.1f}'.format(
                dataset['antenna_coordinates_az_el_r'][0])
            el = '{:.1f}'.format(
                dataset['antenna_coordinates_az_el_r'][1])
            r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
            gateinfo = ('az'+az+'r'+r+'el'+el)
        else:
            lon = '{:.3f}'.format(
                dataset['point_coordinates_WGS84_lon_lat_alt'][0])
            lat = '{:.3f}'.format(
                dataset['point_coordinates_WGS84_lon_lat_alt'][1])
            alt = '{:.1f}'.format(
                dataset['point_coordinates_WGS84_lon_lat_alt'][2])
            gateinfo = ('lon'+lon+'lat'+lat+'alt'+alt)

        timeformat = None
        timeinfo = None
        if set_time_info:
            timeformat = '%Y%m%d'
            timeinfo = dataset['ref_time']

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir, prdsavedir,
            timeinfo=timeinfo)

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], ['csv'],
            prdcfginfo=gateinfo, timeinfo=timeinfo,
            timeformat=timeformat)[0]

        csvfname = savedir+csvfname

        if not dataset['final']:
            if 'antenna_coordinates_az_el_r' in dataset:
                write_ts_polar_data(dataset, csvfname)
            else:
                write_ts_grid_data(dataset, csvfname)
            print('saved CSV file: '+csvfname)

        dpi = prdcfg.get('dpi', 72)
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)
        plot_only_final = prdcfg.get('plot_only_final', False)

        if plot_only_final and not dataset['final']:
            return None
        if not plot_only_final and dataset['final']:
            return None

        date, value = read_timeseries(csvfname)
        if date is None:
            warn(
                'Unable to plot time series. No valid data')
            return None

        timeinfo_fig = None
        if set_time_info:
            timeinfo_fig = date[0]

        figfname_list = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=timeinfo_fig, timeformat=timeformat)

        for i, figfname in enumerate(figfname_list):
            figfname_list[i] = savedir+figfname

        if 'antenna_coordinates_az_el_r' in dataset:
            label1 = 'Radar (az, el, r): ('+az+', '+el+', '+r+')'
        else:
            label1 = 'Grid (lon, lat, alt): ('+lon+', '+lat+', '+alt+')'
        titl = ('Time Series '+date[0].strftime('%Y-%m-%d'))

        labely = generate_field_name_str(dataset['datatype'])

        plot_timeseries(
            date, [value], figfname_list, labelx='Time UTC',
            labely=labely, labels=[label1], title=titl, dpi=dpi,
            ymin=vmin, ymax=vmax)
        print('----- save to '+' '.join(figfname_list))

        return figfname_list

    if prdcfg['type'] == 'PLOT_CUMULATIVE_POINT':
        dpi = prdcfg.get('dpi', 72)
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)
        plot_only_final = prdcfg.get('plot_only_final', False)
        set_time_info = prdcfg.get('set_time_info', True)

        if plot_only_final and not dataset['final']:
            return None
        if not plot_only_final and dataset['final']:
            return None

        if 'antenna_coordinates_az_el_r' in dataset:
            az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
            el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
            r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
            gateinfo = ('az'+az+'r'+r+'el'+el)
        else:
            lon = '{:.3f}'.format(
                dataset['point_coordinates_WGS84_lon_lat_alt'][0])
            lat = '{:.3f}'.format(
                dataset['point_coordinates_WGS84_lon_lat_alt'][1])
            alt = '{:.1f}'.format(
                dataset['point_coordinates_WGS84_lon_lat_alt'][2])
            gateinfo = ('lon'+lon+'lat'+lat+'alt'+alt)

        timeformat = None
        timeinfo = None
        if set_time_info:
            timeformat = '%Y%m%d'
            timeinfo = dataset['ref_time']

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdid'], timeinfo=timeinfo)

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], ['csv'],
            prdcfginfo=gateinfo, timeinfo=timeinfo,
            timeformat=timeformat)[0]

        csvfname = savedir+csvfname

        date, value = read_timeseries(csvfname)

        if date is None:
            warn(
                'Unable to plot accumulation time series. No valid data')
            return None

        timeinfo_fig = None
        if set_time_info:
            timeinfo_fig = date[0]

        figfname_list = make_filename(
            'ts_cum', prdcfg['dstype'], dataset['datatype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=timeinfo_fig, timeformat=timeformat)

        for i, figfname in enumerate(figfname_list):
            figfname_list[i] = savedir+figfname

        if 'antenna_coordinates_az_el_r' in dataset:
            label1 = 'Radar (az, el, r): ('+az+', '+el+', '+r+')'
        else:
            label1 = 'Grid (lon, lat, alt): ('+lon+', '+lat+', '+alt+')'
        titl = ('Time Series Acc. '+date[0].strftime('%Y-%m-%d'))

        labely = 'Radar estimated rainfall accumulation (mm)'

        plot_timeseries(
            date, [value], figfname_list, labelx='Time UTC',
            labely=labely, labels=[label1], title=titl,
            period=prdcfg['ScanPeriod']*60.,
            ymin=vmin, ymax=vmax, dpi=dpi)
        print('----- save to '+' '.join(figfname_list))

        return figfname_list

    if prdcfg['type'] == 'COMPARE_POINT':
        dpi = prdcfg.get('dpi', 72)
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)
        plot_only_final = prdcfg.get('plot_only_final', False)
        set_time_info = prdcfg.get('set_time_info', True)

        if plot_only_final and not dataset['final']:
            return None
        if not plot_only_final and dataset['final']:
            return None

        if 'antenna_coordinates_az_el_r' in dataset:
            az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
            el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
            r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
            gateinfo = ('az'+az+'r'+r+'el'+el)
        else:
            lon = '{:.3f}'.format(
                dataset['point_coordinates_WGS84_lon_lat_alt'][0])
            lat = '{:.3f}'.format(
                dataset['point_coordinates_WGS84_lon_lat_alt'][1])
            alt = '{:.1f}'.format(
                dataset['point_coordinates_WGS84_lon_lat_alt'][2])
            gateinfo = ('lon'+lon+'lat'+lat+'alt'+alt)

        timeformat = None
        timeinfo = None
        if set_time_info:
            timeformat = '%Y%m%d'
            timeinfo = dataset['ref_time']

        savedir_ts = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdid'], timeinfo=timeinfo)

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], ['csv'],
            prdcfginfo=gateinfo, timeinfo=timeinfo,
            timeformat=timeformat)[0]

        csvfname = savedir_ts+csvfname

        radardate, radarvalue = read_timeseries(csvfname)
        if radardate is None:
            warn(
                'Unable to plot sensor comparison at point of interest. ' +
                'No valid radar data')
            return None

        sensordate, sensorvalue, sensortype, _ = get_sensor_data(
            radardate[0], dataset['datatype'], prdcfg)
        if sensordate is None:
            warn(
                'Unable to plot sensor comparison at point of interest. ' +
                'No valid sensor data')
            return None

        timeinfo_fig = None
        if set_time_info:
            timeinfo_fig = radardate[0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=timeinfo_fig)

        figfname_list = make_filename(
            'ts_comp', prdcfg['dstype'], dataset['datatype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=timeinfo_fig, timeformat=timeformat)

        for i, figfname in enumerate(figfname_list):
            figfname_list[i] = savedir+figfname

        if 'antenna_coordinates_az_el_r' in dataset:
            label1 = 'Radar (az, el, r): ('+az+', '+el+', '+r+')'
        else:
            label1 = 'Grid (lon, lat, alt): ('+lon+', '+lat+', '+alt+')'
        label2 = sensortype+' '+prdcfg['sensorid']
        titl = 'Time Series Comp. '+radardate[0].strftime('%Y-%m-%d')
        labely = generate_field_name_str(dataset['datatype'])

        plot_timeseries_comp(
            radardate, radarvalue, sensordate, sensorvalue, figfname_list,
            labelx='Time UTC', labely=labely, label1=label1, label2=label2,
            titl=titl, ymin=vmin, ymax=vmax, dpi=dpi)
        print('----- save to '+' '.join(figfname_list))

        return figfname_list

    if prdcfg['type'] == 'COMPARE_CUMULATIVE_POINT':
        dpi = prdcfg.get('dpi', 72)
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)
        plot_only_final = prdcfg.get('plot_only_final', False)
        set_time_info = prdcfg.get('set_time_info', True)

        if plot_only_final and not dataset['final']:
            return None
        if not plot_only_final and dataset['final']:
            return None

        if 'antenna_coordinates_az_el_r' in dataset:
            az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
            el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
            r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
            gateinfo = ('az'+az+'r'+r+'el'+el)
        else:
            lon = '{:.3f}'.format(
                dataset['point_coordinates_WGS84_lon_lat_alt'][0])
            lat = '{:.3f}'.format(
                dataset['point_coordinates_WGS84_lon_lat_alt'][1])
            alt = '{:.1f}'.format(
                dataset['point_coordinates_WGS84_lon_lat_alt'][2])
            gateinfo = ('lon'+lon+'lat'+lat+'alt'+alt)

        timeformat = None
        timeinfo = None
        if set_time_info:
            timeformat = '%Y%m%d'
            timeinfo = dataset['ref_time']

        savedir_ts = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdid'], timeinfo=timeinfo)

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], ['csv'],
            prdcfginfo=gateinfo, timeinfo=timeinfo,
            timeformat=timeformat)[0]

        csvfname = savedir_ts+csvfname

        radardate, radarvalue = read_timeseries(csvfname)
        if radardate is None:
            warn(
                'Unable to plot sensor comparison at point of interest. ' +
                'No valid radar data')
            return None

        sensordate, sensorvalue, sensortype, period2 = get_sensor_data(
            radardate[0], dataset['datatype'], prdcfg)
        if sensordate is None:
            warn(
                'Unable to plot sensor comparison at point of interest. ' +
                'No valid sensor data')
            return None

        timeinfo_fig = None
        if set_time_info:
            timeinfo_fig = radardate[0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=timeinfo_fig)

        figfname_list = make_filename(
            'ts_cumcomp', prdcfg['dstype'], dataset['datatype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=timeinfo_fig, timeformat=timeformat)

        for i, figfname in enumerate(figfname_list):
            figfname_list[i] = savedir+figfname

        if 'antenna_coordinates_az_el_r' in dataset:
            label1 = 'Radar (az, el, r): ('+az+', '+el+', '+r+')'
        else:
            label1 = 'Grid (lon, lat, alt): ('+lon+', '+lat+', '+alt+')'
        label2 = sensortype+' '+prdcfg['sensorid']
        titl = ('Time Series Acc. Comp. ' +
                radardate[0].strftime('%Y-%m-%d'))
        labely = 'Rainfall accumulation (mm)'

        plot_timeseries_comp(
            radardate, radarvalue, sensordate, sensorvalue,
            figfname_list, labelx='Time UTC', labely=labely,
            label1=label1, label2=label2, titl=titl,
            period1=prdcfg['ScanPeriod']*60., period2=period2,
            ymin=vmin, ymax=vmax, dpi=dpi)
        print('----- save to '+' '.join(figfname_list))

        return figfname_list

    if prdcfg['type'] == 'COMPARE_TIME_AVG':
        dpi = prdcfg.get('dpi', 72)
        plot_only_final = prdcfg.get('plot_only_final', False)

        if plot_only_final and not dataset['final']:
            return None
        if not plot_only_final and dataset['final']:
            return None

        az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
        el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
        r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
        gateinfo = ('az'+az+'r'+r+'el'+el)

        savedir_ts = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdid'], timeinfo=dataset['time'])

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], ['csv'],
            prdcfginfo=gateinfo, timeinfo=dataset['time'],
            timeformat='%Y%m%d')[0]

        radardate, radarvalue = read_timeseries(savedir_ts+csvfname)
        if radardate is None:
            warn(
                'Unable to compared time averaged data at POI. ' +
                'No valid radar data')
            return None

        sensordate, sensorvalue, sensortype, period2 = get_sensor_data(
            radardate[0], dataset['datatype'], prdcfg)
        if sensordate is None:
            warn(
                'Unable to compared time averaged data at POI. ' +
                'No valid sensor data')
            return None

        cum_time = prdcfg.get('cum_time', 3600)
        base_time = prdcfg.get('base_time', 0)

        sensordate_cum, sensorvalue_cum, np_sensor_cum = rainfall_accumulation(
            sensordate, sensorvalue, cum_time=cum_time, base_time=base_time,
            dropnan=False)

        radardate_cum, radarvalue_cum, np_radar_cum = rainfall_accumulation(
            radardate, radarvalue, cum_time=cum_time, base_time=base_time,
            dropnan=False)

        # find common time stamps
        ind = np.where(np.in1d(radardate_cum, sensordate_cum))[0]
        if ind.size == 0:
            warn('No sensor data for radar data time stamps')
        radardate_cum2 = radardate_cum[ind]
        radarvalue_cum2 = radarvalue_cum[ind]
        np_radar_cum2 = np_radar_cum[ind]

        ind = np.where(np.in1d(sensordate_cum, radardate_cum2))[0]
        sensorvalue_cum2 = sensorvalue_cum[ind]
        np_sensor_cum2 = np_sensor_cum[ind]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=radardate[0])

        fname = make_filename(
            str(cum_time)+'s_acc_ts_comp', prdcfg['dstype'],
            dataset['datatype'], ['csv'], prdcfginfo=gateinfo,
            timeinfo=radardate[0], timeformat='%Y%m%d')[0]

        fname = savedir+fname

        new_dataset = deepcopy(dataset)
        new_dataset.update({
            'time': radardate_cum2,
            'sensor_value': sensorvalue_cum2,
            'np_sensor': np_sensor_cum2,
            'radar_value': radarvalue_cum2,
            'np_radar': np_radar_cum2,
            'cum_time': cum_time})
        new_dataset.update(prdcfg)

        write_ts_cum(new_dataset, fname)

        print('saved CSV file: '+fname)

        figfname_list = make_filename(
            str(cum_time)+'s_acc_ts_comp', prdcfg['dstype'],
            dataset['datatype'], prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=radardate[0], timeformat='%Y%m%d')

        for i, figfname in enumerate(figfname_list):
            figfname_list[i] = savedir+figfname

        labelx = sensortype+' '+prdcfg['sensorid']+' (mm)'
        labely = 'Radar (az, el, r): ('+az+', '+el+', '+r+') (mm)'
        titl = (str(cum_time)+'s Acc. Comp. ' +
                radardate_cum[0].strftime('%Y-%m-%d'))

        plot_scatter_comp(
            sensorvalue_cum2, radarvalue_cum2, figfname_list, labelx=labelx,
            labely=labely, titl=titl, axis='equal', dpi=dpi)

        print('----- save to '+' '.join(figfname_list))

        return figfname_list

    # ================================================================
    if prdcfg['type'] == 'PLOT_AND_WRITE':
        if not dataset['final']:
            return None

        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['ts_dict']:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        ts = dataset['ts_dict'][field_name]
        timeinfo = ts.time_vector[0]

        savedir = get_save_dir(prdcfg['basepath'], prdcfg['procname'],
                               dssavedir, prdsavedir,
                               timeinfo=timeinfo)

        dstype_str = prdcfg['dstype'].lower().replace('_', '')
        fname = make_filename('ts', dstype_str, ts.datatype,
                              ['csv'],
                              prdcfginfo=None, timeinfo=timeinfo,
                              timeformat='%Y%m%d%H%M%S',
                              runinfo=prdcfg['runinfo'])

        ts.write(savedir + fname[0])

        fname = make_filename('ts', dstype_str, ts.datatype,
                              prdcfg['imgformat'],
                              prdcfginfo=None, timeinfo=timeinfo,
                              timeformat='%Y%m%d%H%M%S',
                              runinfo=prdcfg['runinfo'])

        ymin = prdcfg.get('ymin', None)
        ymax = prdcfg.get('ymax', None)

        ts.plot(savedir + fname[0], ymin=ymin, ymax=ymax)

        return None

    if prdcfg['type'] == 'PLOT_HIST':
        if not dataset['final']:
            return None

        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['ts_dict']:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        ts = dataset['ts_dict'][field_name]
        timeinfo = ts.time_vector[0]

        savedir = get_save_dir(prdcfg['basepath'], prdcfg['procname'],
                               dssavedir, prdsavedir,
                               timeinfo=timeinfo)

        dstype_str = prdcfg['dstype'].lower().replace('_', '')

        fname = make_filename('hist', dstype_str, ts.datatype,
                              prdcfg['imgformat'],
                              prdcfginfo=None, timeinfo=timeinfo,
                              timeformat='%Y%m%d%H%M%S',
                              runinfo=prdcfg['runinfo'])

        step = prdcfg.get('step', None)

        ts.plot_hist(savedir + fname[0], step=step)

        return None

    if prdcfg['type'] == 'TRAJ_CAPPI_IMAGE':
        if dataset['final']:
            return None

        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        savedir = get_save_dir(prdcfg['basepath'], prdcfg['procname'],
                               dssavedir, prdsavedir,
                               timeinfo=prdcfg['timeinfo'])

        color_ref = prdcfg.get('color_ref', 'None')
        if color_ref == 'altitude':
            prdtype = 'cappi_alt'
        elif color_ref == 'rel_altitude':
            prdtype = 'cappi_rel_alt'
        elif color_ref == 'time':
            prdtype = 'cappi_time'
        else:
            prdtype = 'cappi'

        fname_list = make_filename(
            prdtype, prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            prdcfginfo='alt'+'{:.1f}'.format(prdcfg['altitude']),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        fig, ax = plot_cappi(
            dataset['radar'], field_name, prdcfg['altitude'], prdcfg,
            fname_list, save_fig=False)

        # start time of radar object
        t_start = dataset['radar'].time['data'].min()
        dt_start = num2date(t_start, dataset['radar'].time['units'],
                            dataset['radar'].time['calendar'])

        fname_list = plot_traj(
            dataset['rng_traj'], dataset['azi_traj'], dataset['ele_traj'],
            dataset['time_traj'], prdcfg, fname_list,
            rad_alt=dataset['radar'].altitude['data'], rad_tstart=dt_start,
            ax=ax, fig=fig, save_fig=True)

        print('----- saved to '+' '.join(fname_list))

        return None

    # ================================================================
    warn(' Unsupported product type: ' + prdcfg['type'])
    return None
