"""
pyrad.io.write_data
====================

Functions for writing pyrad output data

.. autosummary::
    :toctree: generated/

    write_centroids
    write_proc_periods
    write_fixed_angle
    write_ts_lightning
    send_msg
    write_alarm_msg
    write_last_state
    write_smn
   *write_timeseries_point
    write_trt_info
    write_trt_cell_data
    write_trt_thundertracking_data
    write_trt_cell_scores
    write_trt_cell_lightning
    write_trt_rpc
    write_rhi_profile
    write_field_coverage
    write_cdf
    write_histogram
    write_quantiles
    write_ts_polar_data
    write_ts_grid_data
    write_ts_ml
    write_ts_stats
    write_ts_cum
    write_monitoring_ts
    write_excess_gates
    write_intercomp_scores_ts
    write_colocated_gates
    write_colocated_data
    write_colocated_data_time_avg
    write_sun_hits
    write_sun_retrieval

"""

from __future__ import print_function
import glob
import csv
import os

from warnings import warn
import smtplib
from email.message import EmailMessage
import fcntl
import errno
import time
from datetime import datetime as dt

import numpy as np

from pyart.config import get_fillvalue

from .io_aux import generate_field_name_str


def write_centroids(fname, centroids_dict, var_names):
    """
    writes the centroids of a semi-supervised hydrometeor classification

    Parameters
    ----------
    centroids_dict : dict
        dictionary containing the centroids
    var_names : tuple
        List of variables

    Returns
    -------
    fname : str
        the name of the file containing the content

    """
    with open(fname, 'w', newline='') as csvfile:
        csvfile.write("# Centroids file\n")
        csvfile.write("# Header lines with comments are preceded by '#'\n")
        csvfile.write("#\n")

        field_names = np.append('hydro_type', var_names)

        writer = csv.DictWriter(csvfile, field_names)
        writer.writeheader()

        for hydro_type in centroids_dict.keys():
            dict_row = {'hydro_type': hydro_type}
            for ivar, var_name in enumerate(var_names):
                dict_row.update({var_name: centroids_dict[hydro_type][ivar]})
            writer.writerow(dict_row)
        csvfile.close()

    return fname


def write_proc_periods(start_times, end_times, fname):
    """
    writes an output file containing start and stop times of periods to
    process

    Parameters
    ----------
    start_times, end_times : datetime object
        The starting and ending times of the periods
    fname : str
        The name of the file where to write

    Returns
    -------
    fname : str
        the name of the file containing the content

    """
    with open(fname, 'w', newline='') as csvfile:
        field_names = ['start_time', 'end_time']
        writer = csv.DictWriter(csvfile, field_names)
        writer.writeheader()

        for start_time, end_time in zip(start_times, end_times):
            dict_row = {
                'start_time': start_time.strftime('%Y%m%d%H%M%S'),
                'end_time': end_time.strftime('%Y%m%d%H%M%S')
            }
            writer.writerow(dict_row)
        csvfile.close()

    return fname


def write_fixed_angle(time_data, fixed_angle, rad_lat, rad_lon, rad_alt,
                      fname):
    """
    writes an output file with the fixed angle data

    Parameters
    ----------
    time_data : datetime object
        The scan time
    fixed_angle : float
        The first fixed angle in the scan
    rad_lat, rad_lon, rad_alt : float
        Latitude, longitude [deg] and altitude [m MSL] of the radar
    fname : str
        The name of the file where to write

    Returns
    -------
    fname : str
        the name of the file containing the content

    """
    filelist = glob.glob(fname)
    if not filelist:
        with open(fname, 'w', newline='') as csvfile:
            fieldnames = [
                'date_time', 'rad_lat', 'rad_lon', 'rad_alt', 'fixed_angle']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            writer.writerow({
                'date_time': time_data.strftime('%Y%m%d%H%M%S'),
                'rad_lat': rad_lat,
                'rad_lon': rad_lon,
                'rad_alt': rad_alt,
                'fixed_angle': fixed_angle})

            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = [
                'date_time', 'rad_lat', 'rad_lon', 'rad_alt', 'fixed_angle']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writerow({
                'date_time': time_data.strftime('%Y%m%d%H%M%S'),
                'rad_lat': rad_lat,
                'rad_lon': rad_lon,
                'rad_alt': rad_alt,
                'fixed_angle': fixed_angle})

            csvfile.close()


def write_ts_lightning(flashnr, time_data, time_in_flash, lat, lon, alt, dBm,
                       vals_list, fname, pol_vals_labels):
    """
    writes the LMA sources data and the value of the colocated polarimetric
    variables

    Parameters
    ----------
    flashnr : int
        flash number
    time_data : datetime object
        flash source time
    time_in_flash : float
        seconds since start of flash
    lat, lon, alt : float
        latitude, longitude [deg] and altitude [m MSL] of the flash source
    dBm : float
        flash power
    vals_list : list of arrays
        List containing the data for each polarimetric variable
    fname : str
        the name of the file containing the content
    pol_values_labels : list of strings
        List containing strings identifying each polarimetric variable

    Returns
    -------
    fname : str
        the name of the file containing the content

    """
    with open(fname, 'w', newline='') as csvfile:
        vals_list_aux = []
        for j, label in enumerate(pol_vals_labels):
            vals_list_aux.append(
                vals_list[j].filled(fill_value=get_fillvalue()))

        csvfile.write("# Weather radar timeseries data file\n")
        csvfile.write("# Project: MALSplus\n")
        if time_data.size > 0:
            csvfile.write("# Start : %s UTC\n" %
                          time_data[0].strftime("%Y-%m-%d %H:%M:%S"))
            csvfile.write("# End   : %s UTC\n" %
                          time_data[-1].strftime("%Y-%m-%d %H:%M:%S"))
        csvfile.write("# Header lines with comments are preceded by '#'\n")
        csvfile.write("#\n")

        field_names = [
            'flashnr', 'time_data', 'time_in_flash', 'lat', 'lon', 'alt',
            'dBm']
        field_names.extend(pol_vals_labels)

        writer = csv.DictWriter(csvfile, field_names)
        writer.writeheader()

        if flashnr.size == 0.:
            warn('No data to write in file '+fname)
            csvfile.close()
            return fname

        for i, flash in enumerate(flashnr):
            dict_row = {
                'flashnr': int(flash),
                'time_data': time_data[i].strftime('%Y-%m-%d %H:%M:%S.%f'),
                'time_in_flash': time_in_flash[i],
                'lat': lat[i],
                'lon': lon[i],
                'alt': alt[i],
                'dBm': dBm[i],
            }
            for j, label in enumerate(pol_vals_labels):
                dict_row.update(
                    {label: vals_list_aux[j][i]})

            writer.writerow(dict_row)
        csvfile.close()

    return fname


def send_msg(sender, receiver_list, subject, fname):
    """
    sends the content of a text file by email

    Parameters
    ----------
    sender : str
        the email address of the sender
    receiver_list : list of string
        list with the email addresses of the receiver
    subject : str
        the subject of the email
    fname : str
        name of the file containing the content of the email message

    Returns
    -------
    fname : str
        the name of the file containing the content

    """
    # Create the container email message.
    msg = EmailMessage()
    msg['Subject'] = subject
    msg['From'] = sender
    msg['To'] = ', '.join(receiver_list)

    # Open the plain text file whose name is in fname for reading.
    with open(fname) as fp:
        # Create a text/plain message
        msg.set_content(fp.read())

    # Send the message via our own SMTP server.
    with smtplib.SMTP('localhost') as s:
        s.send_message(msg)

    return fname


def write_alarm_msg(radar_name, param_name_unit, date_last, target, tol_abs,
                    np_trend, value_trend, tol_trend, nevents, np_last,
                    value_last, fname):
    """
    writes an alarm file

    Parameters
    ----------
    radar_name : str
        Name of the radar being controlled
    param_name_unit : str
        Parameter and units
    date_last : datetime object
        date of the current event
    target, tol_abs : float
        Target value and tolerance
    np_trend : int
        Total number of points in trend
    value_trend, tol_trend : float
        Trend value and tolerance
    nevents: int
        Number of events in trend
    np_last : int
        Number of points in the current event
    value_last : float
        Value of the current event
    fname : str
        Name of file where to store the alarm information

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    with open(fname, 'w', newline='') as txtfile:
        txtfile.write(
            'Weather radar polarimetric parameters monitoring alarm\n')
        if radar_name is not None:
            txtfile.write('Radar name: '+radar_name+'\n')
        txtfile.write('Parameter [Unit]: '+param_name_unit+'\n')
        txtfile.write('Date : '+date_last.strftime('%Y%m%d')+'\n')
        txtfile.write('Target value: '+str(target)+' +/- '+str(tol_abs)+'\n')
        if np_trend > 0:
            txtfile.write(
                'Number of points in trend: '+str(np_trend)+' in ' +
                str(nevents)+' events\n')
            txtfile.write(
                'Trend value: '+str(value_trend)+' (tolerance: +/- ' +
                str(tol_trend)+')\n')
        else:
            txtfile.write('Number of points in trend: NA\n')
            txtfile.write('Trend value: NA\n')
        txtfile.write('Number of points: '+str(np_last)+'\n')
        txtfile.write('Value: '+str(value_last)+'\n')

        txtfile.close()

    return fname


def write_last_state(datetime_last, fname):
    """
    writes SwissMetNet data in format datetime,avg_value, std_value

    Parameters
    ----------
    datetime_last : datetime object
        date and time of the last state
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    try:
        with open(fname, 'w', newline='') as txtfile:
            txtfile.write(datetime_last.strftime('%Y%m%d%H%M%S'))
            txtfile.close()

            return fname
    except EnvironmentError:
        warn('Unable to write on file '+fname)
        return None


def write_smn(datetime_vec, value_avg_vec, value_std_vec, fname):
    """
    writes SwissMetNet data in format datetime,avg_value, std_value

    Parameters
    ----------
    datetime_vec : datetime array
        array containing the measurement time
    value_avg_vec : float array
        array containing the average value
    value_std_vec : float array
        array containing the standard deviation
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    with open(fname, 'w', newline='') as csvfile:
        fieldnames = ['datetime', 'avg', 'std']
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for i, value_avg in enumerate(value_avg_vec):
            writer.writerow({
                'datetime': datetime_vec[i].strftime('%Y%m%d%H%M%S'),
                'avg': value_avg,
                'std': value_std_vec[i]})

        csvfile.close()

    return fname


def write_timeseries_point(fname, data, dstype, text, timeformat=None, timeinfo=None):
    """
    Write one timesample of a time series to a file

    Parameters
    ----------
    data      : time series structure with the fields:
      time    : Time sample in Julian calendar format.
      label[] : Header of each value column
      value[] : Data value of each column
    dataType  : Type of the data
    cfginfo   : Information string describing the kind of data
                in the text file. This string is part of the
                generated file name.
    text      : Text with description of the data in the
                file. Must be a string array. Each string
                is written to a new line an the comment symbol
                is set at the beginning.

    Keywords
    --------
    timeformat: If set, the date time column is set according to the
                time format number. If not set, the time format
                is 'YYYY-MM-DD, <seconds of day>'
    timeinfo  : Used for the time stamp of the output file. If not set,
                data.time is used to generate it.

    """

    datatime = data['time'].strftime("%Y-%m-%d %H:%M:%S")
    datavalue = data['value']
    datatypename = dstype
    unit = data['unit']

    print("----- Write timeseries ", fname)

    if os.path.isfile(fname) == False:
        #File does not exist. Open it and fill in header info.
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write("# Weather radar timeseries data file\n")
            csvfile.write("# Project: MALSplus\n")
            csvfile.write("# Data/Unit : " +  datatypename + " [" + unit + "]\n")
            csvfile.write("# Start : " + datatime + " UTC\n")
            csvfile.write("# Header lines with comments are preceded by '#'\n")
            for line in text:
                csvfile.write("# " + line +'\n')
            csvfile.write("#\n")

            if timeformat is None:
                label_str = "# Date [YYYY-MM-DD hh:mm:ss]"
                tformat = "%Y-%m-%d %H:%M:%S"
                time_str_old = dt.strptime(datatime, tformat)
                time_str = time_str_old.strftime(tformat)
            else:
                label_str = "# Date ["+ timeformat +"]"
                tformat = timeformat
                time_str = datatime.strftime(tformat)

            for label in data['label']:
                label_str = label_str + ", " + label
            csvfile.write(label_str + '\n')

            for value in data['value']:
                time_str = time_str + ", " + ('%.4f'% value)
            csvfile.write(time_str + '\n')

    else:
        if timeformat is None:
            tformat = "%Y-%m-%d %H:%M:%S"
            time_str_old = dt.strptime(datatime, tformat)
            time_str = time_str_old.strftime(tformat)
        else:
            tformat = timeformat
            time_str = datatime.strftime(tformat)

        with open(fname, 'a', newline='') as csvfile:

            for value in data['value']:
                time_str = time_str + ", " + ('%.4f'% value)
            csvfile.write(time_str + '\n')

    csvfile.close()

    return fname


def write_trt_info(ids, max_rank, nscans, time_start, time_end, fname):
    """
    writes TRT info of the thundertracking

    Parameters
    ----------
    ids, max_rank, nscans, time_start, time_end: array
        the cell parameters
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    with open(fname, 'w', newline='') as csvfile:
        fieldnames = [
            'id', 'max_rank', 'nscans_Xband', 'time_start', 'time_end']
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for i, id_cell in enumerate(ids):
            writer.writerow({
                'id': id_cell,
                'max_rank': max_rank[i],
                'nscans_Xband': nscans[i],
                'time_start': time_start[i].strftime('%Y%m%d%H%M'),
                'time_end': time_end[i].strftime('%Y%m%d%H%M')
            })

        csvfile.close()

    return fname


def write_trt_cell_data(
        traj_ID, yyyymmddHHMM, lon, lat, ell_L, ell_S, ell_or, area,
        vel_x, vel_y, det, RANKr, CG_n, CG_p, CG, CG_percent_p, ET45,
        ET45m, ET15, ET15m, VIL, maxH, maxHm, POH, RANK, Dvel_x,
        Dvel_y, cell_contour, fname):
    """
    writes TRT cell data

    Parameters
    ----------
    traj_ID, yyyymmddHHMM, lon, lat, ell_L, ell_S, ell_or, area,
    vel_x, vel_y, det, RANKr, CG_n, CG_p, CG, CG_percent_p, ET45,
    ET45m, ET15, ET15m, VIL, maxH, maxHm, POH, RANK, Dvel_x,
    Dvel_y, cell_contour:
        the cell parameters
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    try:
        with open(fname, 'w', newline='') as csvfile:
            fieldnames = [
                'traj_ID', 'yyyymmddHHMM', 'lon', 'lat', 'ell_L', 'ell_S',
                'ell_or', 'area', 'vel_x', 'vel_y', 'det', 'RANKr', 'CG-',
                'CG+', 'CG', '%CG+', 'ET45', 'ET45m', 'ET15', 'ET15m',
                'VIL', 'maxH', 'maxHm', 'POH', 'RANK', 'Dvel_x', 'Dvel_y',
                'cell_contour_lon-lat']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            for i, traj_ID_el in enumerate(traj_ID):
                cell_contour_aux = cell_contour[i]
                npoints_contour = len(cell_contour_aux['lon'])
                cell_contour_arr = np.empty(2*npoints_contour, dtype=float)
                cell_contour_arr[0:-1:2] = cell_contour_aux['lon']
                cell_contour_arr[1::2] = cell_contour_aux['lat']
                cell_contour_str = str(cell_contour_arr[0])
                for j in range(1, 2*npoints_contour):
                    cell_contour_str += ' '+str(cell_contour_arr[j])

                writer.writerow({
                    'traj_ID': traj_ID_el,
                    'yyyymmddHHMM': yyyymmddHHMM[i].strftime('%Y%m%d%H%M'),
                    'lon': lon[i],
                    'lat': lat[i],
                    'ell_L': ell_L[i],
                    'ell_S': ell_S[i],
                    'ell_or': ell_or[i],
                    'area': area[i],
                    'vel_x': vel_x[i],
                    'vel_y': vel_y[i],
                    'det': det[i],
                    'RANKr': RANKr[i],
                    'CG-': CG_n[i],
                    'CG+': CG_p[i],
                    'CG': CG[i],
                    '%CG+': CG_percent_p[i],
                    'ET45': ET45[i],
                    'ET45m': ET45m[i],
                    'ET15': ET15[i],
                    'ET15m': ET15m[i],
                    'VIL': VIL[i],
                    'maxH': maxH[i],
                    'maxHm': maxHm[i],
                    'POH': POH[i],
                    'RANK': RANK[i],
                    'Dvel_x': Dvel_x[i],
                    'Dvel_y': Dvel_y[i],
                    'cell_contour_lon-lat': cell_contour_str
                })

            csvfile.close()

            return fname
    except EnvironmentError:
        warn('Unable to write on file '+fname)
        return None


def write_trt_thundertracking_data(
        traj_ID, scan_ordered_time, scan_time, azi, rng, yyyymmddHHMM, lon,
        lat, ell_L, ell_S, ell_or, area, vel_x, vel_y, det, RANKr, CG_n, CG_p,
        CG, CG_percent_p, ET45, ET45m, ET15, ET15m, VIL, maxH, maxHm, POH,
        RANK, Dvel_x, Dvel_y, cell_contour, fname):
    """
    writes TRT cell data of the thundertracking scan

    Parameters
    ----------
    traj_ID, scan_ordered_time, scan_time, azi, rng, yyyymmddHHMM, lon, lat,
    ell_L, ell_S, ell_or, area, vel_x, vel_y, det, RANKr, CG_n, CG_p, CG,
    CG_percent_p, ET45, ET45m, ET15, ET15m, VIL, maxH, maxHm, POH, RANK,
    Dvel_x, Dvel_y, cell_contour:
        the cell parameters
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    try:
        with open(fname, 'w', newline='') as csvfile:
            fieldnames = [
                'traj_ID', 'scan_ordered_time', 'scan_time', 'azi', 'rng',
                'yyyymmddHHMM', 'lon', 'lat', 'ell_L', 'ell_S', 'ell_or',
                'area', 'vel_x', 'vel_y', 'det', 'RANKr', 'CG-', 'CG+', 'CG',
                '%CG+', 'ET45', 'ET45m', 'ET15', 'ET15m', 'VIL', 'maxH',
                'maxHm', 'POH', 'RANK', 'Dvel_x', 'Dvel_y',
                'cell_contour_lon-lat']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            for i, traj_ID_el in enumerate(traj_ID):
                cell_contour_aux = cell_contour[i]
                npoints_contour = len(cell_contour_aux['lon'])
                cell_contour_arr = np.empty(2*npoints_contour, dtype=float)
                cell_contour_arr[0:-1:2] = cell_contour_aux['lon']
                cell_contour_arr[1::2] = cell_contour_aux['lat']
                cell_contour_str = str(cell_contour_arr[0])
                for j in range(1, 2*npoints_contour):
                    cell_contour_str += ' '+str(cell_contour_arr[j])

                if np.ma.is_masked(scan_time[i]):
                    scan_time_str = get_fillvalue()
                else:
                    scan_time_str = scan_time[i].strftime('%Y%m%d%H%M%S')
                writer.writerow({
                    'traj_ID': traj_ID_el,
                    'scan_ordered_time': scan_ordered_time[i].strftime(
                        '%Y%m%d%H%M%S.%f'),
                    'scan_time': scan_time_str,
                    'azi': azi[i],
                    'rng': rng[i],
                    'yyyymmddHHMM': yyyymmddHHMM[i].strftime('%Y%m%d%H%M'),
                    'lon': lon[i],
                    'lat': lat[i],
                    'ell_L': ell_L[i],
                    'ell_S': ell_S[i],
                    'ell_or': ell_or[i],
                    'area': area[i],
                    'vel_x': vel_x[i],
                    'vel_y': vel_y[i],
                    'det': det[i],
                    'RANKr': RANKr[i],
                    'CG-': CG_n[i],
                    'CG+': CG_p[i],
                    'CG': CG[i],
                    '%CG+': CG_percent_p[i],
                    'ET45': ET45[i],
                    'ET45m': ET45m[i],
                    'ET15': ET15[i],
                    'ET15m': ET15m[i],
                    'VIL': VIL[i],
                    'maxH': maxH[i],
                    'maxHm': maxHm[i],
                    'POH': POH[i],
                    'RANK': RANK[i],
                    'Dvel_x': Dvel_x[i],
                    'Dvel_y': Dvel_y[i],
                    'cell_contour_lon-lat': cell_contour_str
                })

            csvfile.close()

            return fname
    except EnvironmentError:
        warn('Unable to write on file '+fname)
        return None


def write_trt_cell_scores(
        traj_ID, flash_density_max_time, flash_density_max_rank,
        nflashes_max_list, area_flash_max_list, flash_density_max,
        rank_max_time, rank_max, fname):
    """
    writes TRT cells scores

    Parameters
    ----------
    traj_ID : array of ints
        The ID of the cells
    flash_density_max_time : array of date times
        The time at which the maximum flash density was reached for each cell
    flash_density_max_rank : array of floats
        The rank when the maximum flash density was reached for each cell
    nflashes_max_list : array of ints
        the number of flashes when the max flash density was reached
    area_flash_max_list : array of floats
        The area when the max flash density was reached
    flash_density_max : array of floats
        The maximum flash density for each cell
    rank_max_time : array of datetime
        the time at wich the maximum rank of each cell was reached
    rank_max : array of float
        the rank when the maximum rank of each cell was reached
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    with open(fname, 'w', newline='') as csvfile:
        fieldnames = [
            'traj ID', 'max flash density time',
            'max flash density rank', 'max flash density flashes',
            'max flash density area', 'max flash density', 'max rank time',
            'max rank']
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for i, traj_ID_el in enumerate(traj_ID):
            writer.writerow({
                'traj ID': traj_ID_el,
                'max flash density time': flash_density_max_time[i].strftime(
                    '%Y-%m-%d %H:%M:%S'),
                'max flash density rank': flash_density_max_rank[i],
                'max flash density flashes': nflashes_max_list[i],
                'max flash density area': area_flash_max_list[i],
                'max flash density': flash_density_max[i],
                'max rank time': rank_max_time[i].strftime(
                    '%Y-%m-%d %H:%M:%S'),
                'max rank': rank_max[i],
            })

        csvfile.close()

    return fname


def write_trt_cell_lightning(
        cell_ID, cell_time, lon, lat, area, rank, nflash, flash_density,
        fname, timeformat='%Y%m%d%H%M'):
    """
    writes the lightning data for each TRT cell

    Parameters
    ----------
    cell_ID : array of ints
        the cell ID
    cell_time : array of datetime
        the time step
    lon, lat : array of floats
        the latitude and longitude of the center of the cell
    area : array of floats
        the area of the cell
    rank : array of floats
        the rank of the cell
    nflash : array of ints
        the number of flashes/sources within the cell
    flash_density : array of floats
        the flash/source density
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    nflash = nflash.filled(fill_value=get_fillvalue())
    flash_density = flash_density.filled(fill_value=get_fillvalue())
    with open(fname, 'w', newline='') as csvfile:
        fieldnames = [
            'traj_ID', 'yyyymmddHHMM', 'lon', 'lat', 'area', 'RANKr',
            'nflashes', 'flash_dens']
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for i, traj_ID_el in enumerate(cell_ID):
            writer.writerow({
                'traj_ID': traj_ID_el,
                'yyyymmddHHMM': cell_time[i].strftime(timeformat),
                'lon': lon[i],
                'lat': lat[i],
                'area': area[i],
                'RANKr': rank[i],
                'nflashes': nflash[i],
                'flash_dens': flash_density[i],
            })

        csvfile.close()

    return fname


def write_trt_rpc(cell_ID, cell_time, lon, lat, area, rank, hmin, hmax, freq,
                  fname, timeformat='%Y%m%d%H%M'):
    """
    writes the rimed particles column data for a TRT cell

    Parameters
    ----------
    cell_ID : array of ints
        the cell ID
    cell_time : array of datetime
        the time step
    lon, lat : array of floats
        the latitude and longitude of the center of the cell
    area : array of floats
        the area of the cell
    rank : array of floats
        the rank of the cell
    hmin, hmax : array of floats
        Minimum and maximum altitude of the rimed particle column
    freq : array of floats
        Frequency of the species constituting the rime particle column within
        the limits of it
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    hmin = hmin.filled(fill_value=get_fillvalue())
    hmax = hmax.filled(fill_value=get_fillvalue())
    freq = freq.filled(fill_value=get_fillvalue())
    with open(fname, 'w', newline='') as csvfile:
        fieldnames = [
            'traj_ID', 'yyyymmddHHMM', 'lon', 'lat', 'area', 'RANKr',
            'hmin', 'hmax', 'freq']
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for i, traj_ID_el in enumerate(cell_ID):
            writer.writerow({
                'traj_ID': traj_ID_el,
                'yyyymmddHHMM': cell_time[i].strftime(timeformat),
                'lon': lon[i],
                'lat': lat[i],
                'area': area[i],
                'RANKr': rank[i],
                'hmin': hmin[i],
                'hmax': hmax[i],
                'freq': freq[i]
            })

        csvfile.close()

    return fname


def write_rhi_profile(hvec, data, nvalid_vec, labels, fname, datatype=None,
                      timeinfo=None, sector=None):
    """
    writes the values of an RHI profile in a text file

    Parameters
    ----------
    hvec : float array
        array containing the alitude in m MSL
    data : list of float array
        the quantities at each altitude
    nvalid_vec : int array
        number of valid data points used to compute the quantiles
    labels : list of strings
        label specifying the quantitites in data
    fname : str
        file name where to store the data
    datatype : str
        the data type
    timeinfo : datetime object
        time of the rhi profile
    sector : dict
        dictionary specying the sector limits

    Returns
    -------
    fname : str
        the name of the file where data has been written

    """
    data_aux = list()
    for data_points in data:
        data_aux.append(data_points.filled(fill_value=get_fillvalue()))
    with open(fname, 'w', newline='') as csvfile:
        csvfile.write('# Height Profile\n')
        csvfile.write('# ==============\n')
        csvfile.write('#\n')
        if datatype is None:
            csvfile.write('# Datatype (Unit) : Not specified\n')
        else:
            csvfile.write('# Datatype (Unit) : '+datatype+'\n')
        csvfile.write('# Fill value : '+str(get_fillvalue())+'\n')
        if timeinfo is None:
            csvfile.write('# Time     : Not specified\n')
        else:
            csvfile.write('# Time     : ' +
                          timeinfo.strftime('%Y-%m-%d %H:%M:%S')+'\n')
        if sector is None:
            csvfile.write('# Sector specification: None\n')
        else:
            csvfile.write('# Sector specification:\n')
            csvfile.write('#   Azimuth         : ' + str(sector['az']) +
                          ' deg\n')
            csvfile.write('#   Range start     : ' +
                          str(sector['rmin'])+' m\n')
            csvfile.write('#   Range stop      : ' +
                          str(sector['rmax'])+' m\n')

        fieldnames = ['Altitude [m MSL]']
        fieldnames.extend(labels)
        fieldnames.append('N valid')
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for j, height in enumerate(hvec):
            data_dict = {fieldnames[0]: height}
            for i in range(len(labels)):
                data_dict.update({fieldnames[i+1]: data_aux[i][j]})
            data_dict.update({fieldnames[-1]: nvalid_vec[j]})
            writer.writerow(data_dict)

        csvfile.close()

    return fname


def write_field_coverage(quantiles, values, ele_start, ele_stop, azi_start,
                         azi_stop, threshold, nvalid_min, datatype, timeinfo,
                         fname):
    """
    writes the quantiles of the coverage on a particular sector

    Parameters
    ----------
    quantiles : datetime array
        array containing the quantiles computed
    values : float array
        quantile value
    ele_start, ele_stop, azi_start, azi_stop : float
        The limits of the sector
    threshold : float
        The minimum value to consider the data valid
    nvalid_min : int
        the minimum number of points to consider that there are values in a
        ray
    datatype : str
        data type and units
    timeinfo : datetime object
        the time stamp of the data
    fname : str
        name of the file where to write the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    with open(fname, 'w', newline='') as csvfile:
        csvfile.write('# Quantiles of the field coverage\n')
        csvfile.write('# Datatype (Unit) : '+datatype+'\n')
        csvfile.write('# Time : '+timeinfo.strftime('%Y-%m-%d %H:%M:%S')+'\n')
        csvfile.write('# Sector specification:\n')
        csvfile.write('#   Azimuth start   : '+str(azi_start)+' deg\n')
        csvfile.write('#   Azimuth stop    : '+str(azi_stop)+' deg\n')
        csvfile.write('#   Elevation start : '+str(ele_start)+' deg\n')
        csvfile.write('#   Elevation stop  : '+str(ele_stop)+' deg\n')
        csvfile.write('# Threshold : '+str(threshold)+'\n')
        csvfile.write('# Minimum number of valid gates per ray : ' +
                      str(nvalid_min)+'\n')
        fieldnames = ['Quantile [%]', 'Rain extension [m]']
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for i, quant in enumerate(quantiles):
            writer.writerow({
                'Quantile [%]': quant,
                'Rain extension [m]': values[i]})

        csvfile.close()

    return fname


def write_cdf(quantiles, values, ntot, nnan, nclut, nblocked, nprec_filter,
              noutliers, ncdf, fname, use_nans=False, nan_value=0.,
              filterprec=[], vismin=None, sector=None, datatype=None,
              timeinfo=None):
    """
    writes a cumulative distribution function

    Parameters
    ----------
    quantiles : datetime array
        array containing the measurement time
    values : float array
        array containing the average value
    fname : float array
        array containing the standard deviation
    sector : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    hydrotype_list = ['NC', 'DS', 'CR', 'LR', 'GR', 'RN', 'VI', 'WS', 'MH',
                      'IH/HDG']
    with open(fname, 'w', newline='') as txtfile:
        txtfile.write('Statistical analysis\n')
        txtfile.write('====================\n\n')
        if datatype is None:
            txtfile.write('Datatype (Unit) : Not specified\n')
        else:
            txtfile.write('Datatype (Unit) : '+datatype+'\n')
        if timeinfo is None:
            txtfile.write('Time     : Not specified\n')
        else:
            txtfile.write('Time     : ' +
                          timeinfo.strftime('%Y-%m-%d %H:%M:%S')+'\n')
        if sector is None:
            txtfile.write('Sector specification: None\n')
        else:
            txtfile.write('Sector specification:\n')
            if sector['rmin'] is None:
                txtfile.write('  Range start     : Not specified\n')
            else:
                txtfile.write('  Range start     : ' +
                              str(sector['rmin'])+' m\n')

            if sector['rmax'] is None:
                txtfile.write('  Range stop      : Not specified\n')
            else:
                txtfile.write('  Range stop      : ' +
                              str(sector['rmax'])+' m\n')

            if sector['azmin'] is None:
                txtfile.write('  Azimuth start   : Not specified\n')
            else:
                txtfile.write('  Azimuth start   : ' +
                              str(sector['azmin'])+' deg\n')

            if sector['azmax'] is None:
                txtfile.write('  Azimuth stop    : Not specified\n')
            else:
                txtfile.write('  Azimuth stop    : ' +
                              str(sector['azmax'])+' deg\n')

            if sector['elmin'] is None:
                txtfile.write('  Elevation start : Not specified\n')
            else:
                txtfile.write('  Elevation start : ' +
                              str(sector['elmin'])+' deg\n')

            if sector['elmax'] is None:
                txtfile.write('  Elevation stop  : Not specified\n')
            else:
                txtfile.write('  Elevation stop  : ' +
                              str(sector['elmax'])+' deg\n')

            if sector['hmin'] is None:
                txtfile.write('  Height start    : Not specified\n')
            else:
                txtfile.write('  Height start    : ' +
                              str(sector['hmin'])+' m\n')

            if sector['hmax'] is None:
                txtfile.write('  Height stop     : Not specified\n')
            else:
                txtfile.write('  Height stop     : ' +
                              str(sector['hmax'])+' m\n')
            txtfile.write('')
        txtfile.write('Total number of gates in sector      : ' +
                      str(ntot)+'\n')
        txtfile.write('Number of gates with no value (NaNs) : ' +
                      str(nnan)+'\n')
        if use_nans:
            txtfile.write('  NaNs are set to          : '+str(nan_value)+'\n')
        else:
            txtfile.write('  NaNs are ignored!\n')
        if nclut == -1:
            txtfile.write('Clutter contaminated gates           : ' +
                          'Not checked\n')
        else:
            txtfile.write('Clutter contaminated gates           : ' +
                          str(nclut)+'\n')
        if nblocked == -1:
            txtfile.write('Blocked gates                        : ' +
                          'Not checked\n')
        else:
            txtfile.write('Blocked gates (vismin = '+str(int(vismin)) +
                          ') : '+str(nblocked)+'\n')
        if nprec_filter == -1:
            txtfile.write('Filtered precipitation gates         : None\n')
        else:
            txtfile.write('Filtered precipitation gates         : ' +
                          str(nprec_filter)+'\n')
            txtfile.write('  precipitation types filtered: ')
            for ind_hydro in filterprec:
                txtfile.write(hydrotype_list[ind_hydro]+' ')
            txtfile.write('\n')
        txtfile.write('Number of outliers                   : ' +
                      str(noutliers)+'\n')
        txtfile.write('Number of gates used for histogram   : '+str(ncdf) +
                      ' ('+str(int(ncdf*100./ntot))+'%)\n\n')
        txtfile.write('Quantiles\n')
        txtfile.write('=========\n')
        for i, value in enumerate(values):
            txtfile.write('  Quantile_'+str(int(quantiles[i]))+' = ' +
                          '{:5.1f}'.format(value)+'\n')

        txtfile.close()

    return fname


def write_histogram(bin_edges, values, fname, datatype='undefined',
                    step=0):
    """
    writes a histogram

    Parameters
    ----------
    bin_edges : float array
        array containing the histogram bin edges
    values : int array
        array containing the number of points in each bin
    fname : str
        file name
    datatype :str
        The data type
    step : str
        The bin step

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    with open(fname, 'w', newline='') as csvfile:
        datatype_str = 'undefined'
        if datatype != 'undefined':
            datatype_str = generate_field_name_str(datatype)
        csvfile.write(
            '# Weather radar data histogram file\n' +
            '# Comment lines are preceded by "#"\n' +
            '# Description:\n' +
            '# Histogram of weather radar data.\n' +
            '# Data: '+datatype_str+'\n' +
            '# Step: '+str(step)+'\n'
            '#\n')
        fieldnames = ['bin_edge_left', 'bin_edge_right', 'value']
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for i, val in enumerate(values):
            writer.writerow({
                'bin_edge_left': bin_edges[i],
                'bin_edge_right': bin_edges[i+1],
                'value': val})
        csvfile.close()

    return fname


def write_quantiles(quantiles, values, fname, datatype='undefined'):
    """
    writes quantiles

    Parameters
    ----------
    quantiles : float array
        array containing the quantiles to write
    values : float array
        array containing the value of each quantile
    fname : str
        file name
    datatype :str
        The data type

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    values_aux = values.filled(fill_value=get_fillvalue())
    with open(fname, 'w', newline='') as csvfile:
        csvfile.write(
            '# Weather radar data histogram file\n' +
            '# Comment lines are preceded by "#"\n' +
            '# Description:\n' +
            '# Histogram of weather radar data.\n' +
            '# Data: '+generate_field_name_str(datatype)+'\n' +
            '#\n')
        fieldnames = ['quantile', 'value']
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for i, quant in enumerate(quantiles):
            writer.writerow({
                'quantile': quant,
                'value': values_aux[i]})
        csvfile.close()

    return fname


def write_ts_polar_data(dataset, fname):
    """
    writes time series of data

    Parameters
    ----------
    dataset : dict
        dictionary containing the time series parameters

    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    filelist = glob.glob(fname)
    nsamples = len(dataset['used_antenna_coordinates_az_el_r'][0])
    if not filelist:
        with open(fname, 'w', newline='') as csvfile:
            if nsamples > 1:
                start_time = dataset['time'][0]
            else:
                start_time = dataset['time']
            csvfile.write('# Weather radar timeseries data file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('# Description: \n')
            csvfile.write('# Time series of a weather radar data over a ' +
                          'fixed location.\n')
            csvfile.write(
                '# Location [lon, lat, alt]: ' +
                str(dataset['point_coordinates_WGS84_lon_lat_alt']) + '\n')
            csvfile.write(
                '# Nominal antenna coordinates used [az, el, r]: ' +
                str(dataset['antenna_coordinates_az_el_r'])+'\n')
            csvfile.write(
                '# Data: '+generate_field_name_str(dataset['datatype'])+'\n')
            csvfile.write('# Fill Value: '+str(get_fillvalue())+'\n')
            csvfile.write(
                '# Start: ' +
                start_time.strftime('%Y-%m-%d %H:%M:%S UTC')+'\n')
            csvfile.write('#\n')

            fieldnames = ['date', 'az', 'el', 'r', 'value']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()

            if nsamples == 1:
                writer.writerow({
                    'date': dataset['time'],
                    'az': dataset['used_antenna_coordinates_az_el_r'][0][0],
                    'el': dataset['used_antenna_coordinates_az_el_r'][1][0],
                    'r': dataset['used_antenna_coordinates_az_el_r'][2][0],
                    'value': dataset['value']})
            else:
                for i in range(nsamples):
                    writer.writerow({
                        'date': dataset['time'][i],
                        'az':
                            dataset['used_antenna_coordinates_az_el_r'][0][i],
                        'el':
                            dataset['used_antenna_coordinates_az_el_r'][1][i],
                        'r':
                            dataset['used_antenna_coordinates_az_el_r'][2][i],
                        'value': dataset['value'][i]})
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = ['date', 'az', 'el', 'r', 'value']
            writer = csv.DictWriter(csvfile, fieldnames)
            if nsamples == 1:
                writer.writerow({
                    'date': dataset['time'],
                    'az': dataset['used_antenna_coordinates_az_el_r'][0][0],
                    'el': dataset['used_antenna_coordinates_az_el_r'][1][0],
                    'r': dataset['used_antenna_coordinates_az_el_r'][2][0],
                    'value': dataset['value']})
            else:
                for i in range(nsamples):
                    writer.writerow({
                        'date': dataset['time'][i],
                        'az':
                            dataset['used_antenna_coordinates_az_el_r'][0][i],
                        'el':
                            dataset['used_antenna_coordinates_az_el_r'][1][i],
                        'r':
                            dataset['used_antenna_coordinates_az_el_r'][2][i],
                        'value': dataset['value'][i]})
            csvfile.close()

    return fname


def write_ts_grid_data(dataset, fname):
    """
    writes time series of data

    Parameters
    ----------
    dataset : dict
        dictionary containing the time series parameters

    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    filelist = glob.glob(fname)
    if not filelist:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write('# Gridded data timeseries data file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('# Description: \n')
            csvfile.write('# Time series of a gridded data over a ' +
                          'fixed location.\n')
            csvfile.write(
                '# Nominal location [lon, lat, alt]: ' +
                str(dataset['point_coordinates_WGS84_lon_lat_alt']) + '\n')
            csvfile.write(
                '# Grid points used [iz, iy, ix]: ' +
                str(dataset['grid_points_iz_iy_ix'])+'\n')
            csvfile.write(
                '# Data: '+generate_field_name_str(dataset['datatype'])+'\n')
            csvfile.write('# Fill Value: '+str(get_fillvalue())+'\n')
            csvfile.write(
                '# Start: ' +
                dataset['time'].strftime('%Y-%m-%d %H:%M:%S UTC')+'\n')
            csvfile.write('#\n')

            fieldnames = ['date', 'lon', 'lat', 'alt', 'value']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            writer.writerow(
                {'date': dataset['time'].strftime('%Y-%m-%d %H:%M:%S.%f'),
                 'lon': dataset['used_coordinates_WGS84_lon_lat_alt'][0],
                 'lat': dataset['used_coordinates_WGS84_lon_lat_alt'][1],
                 'alt': dataset['used_coordinates_WGS84_lon_lat_alt'][2],
                 'value': dataset['value']})
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = ['date', 'lon', 'lat', 'alt', 'value']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writerow(
                {'date': dataset['time'].strftime('%Y-%m-%d %H:%M:%S.%f'),
                 'lon': dataset['used_coordinates_WGS84_lon_lat_alt'][0],
                 'lat': dataset['used_coordinates_WGS84_lon_lat_alt'][1],
                 'alt': dataset['used_coordinates_WGS84_lon_lat_alt'][2],
                 'value': dataset['value']})
            csvfile.close()

    return fname


def write_ts_ml(dt_ml, ml_top_avg, ml_top_std, thick_avg, thick_std,
                nrays_valid, nrays_total, fname):
    """
    writes time series of melting layer data

    Parameters
    ----------
    dt_ml : date time array
        array of time steps
    ml_top_avg, ml_top_std: float arrays
        the average and the standard deviation of the melting layer top height
    thick_avg, thick_std: float arrays
        the average and the standard deviation of the metling layer thickness
    nrays_valid, nrays_total: int arrays
        the number of rays where melting layer has been identified and the
        total number of arrays in the scan
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    filelist = glob.glob(fname)
    if not filelist:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write(
                '# Weather radar detected melting layer data file\n' +
                '# Comment lines are preceded by "#"\n' +
                '# Description: \n' +
                '# Time series of melting layer data detected by weather radar.\n' +
                '# Fill Value: '+str(get_fillvalue())+'\n' +
                '# Start: '+dt_ml.strftime('%Y-%m-%d %H:%M:%S UTC')+'\n' +
                '#\n')

            fieldnames = [
                'date-time [UTC]', 'mean ml top height [m MSL]',
                'std ml top height [m MSL]', 'mean ml thickness [m]',
                'std ml thickness [m]', 'N valid rays', 'rays total']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            writer.writerow(
                {'date-time [UTC]': dt_ml.strftime('%Y-%m-%d %H:%M:%S'),
                 'mean ml top height [m MSL]': ml_top_avg.filled(
                     fill_value=get_fillvalue()),
                 'std ml top height [m MSL]': ml_top_std.filled(
                     fill_value=get_fillvalue()),
                 'mean ml thickness [m]': thick_avg.filled(
                     fill_value=get_fillvalue()),
                 'std ml thickness [m]': thick_std.filled(
                     fill_value=get_fillvalue()),
                 'N valid rays': nrays_valid,
                 'rays total': nrays_total})
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = [
                'date-time [UTC]', 'mean ml top height [m MSL]',
                'std ml top height [m MSL]', 'mean ml thickness [m]',
                'std ml thickness [m]', 'N valid rays', 'rays total']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writerow(
                {'date-time [UTC]': dt_ml.strftime('%Y-%m-%d %H:%M:%S'),
                 'mean ml top height [m MSL]': ml_top_avg.filled(
                     fill_value=get_fillvalue()),
                 'std ml top height [m MSL]': ml_top_std.filled(
                     fill_value=get_fillvalue()),
                 'mean ml thickness [m]': thick_avg.filled(
                     fill_value=get_fillvalue()),
                 'std ml thickness [m]': thick_std.filled(
                     fill_value=get_fillvalue()),
                 'N valid rays': nrays_valid,
                 'rays total': nrays_total})
            csvfile.close()

    return fname


def write_ts_stats(dt, value, fname, stat='mean'):
    """
    writes time series of statistics

    Parameters
    ----------
    dt : date time array
        array of time steps
    value: float arrays
        the average and the standard deviation of the melting layer top height
    fname : str
        file name where to store the data
    stat : str
        Statistic that is written

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    filelist = glob.glob(fname)
    if not filelist:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write(
                '# Statistics\n' +
                '# Comment lines are preceded by "#"\n' +
                '# Description: \n' +
                '# Time series of '+stat+'.\n' +
                '# Fill Value: '+str(get_fillvalue())+'\n' +
                '# Start: '+dt.strftime('%Y-%m-%d %H:%M:%S UTC')+'\n' +
                '#\n')

            fieldnames = ['date-time [UTC]', 'value']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            writer.writerow(
                {'date-time [UTC]': dt.strftime('%Y-%m-%d %H:%M:%S'),
                 'value': value.filled(
                     fill_value=get_fillvalue())[0]})
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = ['date-time [UTC]', 'value']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writerow(
                {'date-time [UTC]': dt.strftime('%Y-%m-%d %H:%M:%S'),
                 'value': value.filled(
                     fill_value=get_fillvalue())[0]})
            csvfile.close()

    return fname


def write_ts_cum(dataset, fname):
    """
    writes time series accumulation of data

    Parameters
    ----------
    dataset : dict
        dictionary containing the time series parameters

    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    radar_value = dataset['radar_value'].filled(fill_value=get_fillvalue())
    sensor_value = dataset['sensor_value'].filled(fill_value=get_fillvalue())
    np_radar = dataset['np_radar']
    np_sensor = dataset['np_sensor']

    with open(fname, 'w', newline='') as csvfile:
        csvfile.write('# Precipitation accumulation data file\n')
        csvfile.write('# Comment lines are preceded by "#"\n')
        csvfile.write('# Description: \n')
        csvfile.write('# Time series of a precipitation accumulation of ' +
                      'weather radar data and another sensor over a ' +
                      'fixed location.\n')
        csvfile.write(
            '# Location [lon, lat, alt]: ' +
            str(dataset['point_coordinates_WGS84_lon_lat_alt']) + '\n')
        csvfile.write(
            '# Nominal antenna coordinates used [az, el, r]: ' +
            str(dataset['antenna_coordinates_az_el_r'])+'\n')
        csvfile.write('# sensor type: ' + dataset['sensor']+'\n')
        csvfile.write('# sensor ID: ' + dataset['sensorid']+'\n')
        csvfile.write('# Data: Precipitation accumulation over ' +
                      str(dataset['cum_time'])+' s\n')
        csvfile.write('# Fill Value: '+str(get_fillvalue())+'\n')
        csvfile.write(
            '# Start: ' +
            dataset['time'][0].strftime('%Y-%m-%d %H:%M:%S UTC')+'\n')
        csvfile.write('#\n')

        fieldnames = ['date', 'np_radar', 'radar_value', 'np_sensor',
                      'sensor_value']
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for i, rad_val in enumerate(radar_value):
            writer.writerow({
                'date': dataset['time'][i],
                'np_radar': np_radar[i],
                'radar_value': rad_val,
                'np_sensor': np_sensor[i],
                'sensor_value': sensor_value[i]})

        csvfile.close()

    return fname


def write_monitoring_ts(start_time, np_t, values, quantiles, datatype, fname,
                        rewrite=False):
    """
    writes time series of data

    Parameters
    ----------
    start_time : datetime object or array of date time objects
        the time of the monitoring
    np_t : int or array of ints
        the total number of points
    values: float array with 3 elements of array of arrays
        the values at certain quantiles
    quantiles: float array with 3 elements
        the quantiles computed
    datatype : str
        The data type
    fname : str
        file name where to store the data
    rewrite : bool
        if True a new file is created

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    nvalues = np.size(start_time)
    if nvalues == 1:
        start_time_aux = np.asarray([start_time])
        np_t_aux = np.asarray([np_t])
        values_aux = np.asarray([values.filled(fill_value=get_fillvalue())])
    else:
        start_time_aux = np.asarray(start_time)
        values_aux = values.filled(fill_value=get_fillvalue())
        np_t_aux = np_t

    if rewrite:
        file_exists = False
    else:
        filelist = glob.glob(fname)
        if not filelist:
            file_exists = False
        else:
            file_exists = True

    if not file_exists:
        with open(fname, 'w', newline='') as csvfile:
            while True:
                try:
                    fcntl.flock(csvfile, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except OSError as e:
                    if e.errno != errno.EAGAIN:
                        raise
                    else:
                        time.sleep(0.1)

            csvfile.write(
                '# Weather radar monitoring timeseries data file\n' +
                '# Comment lines are preceded by "#"\n' +
                '# Description: \n' +
                '# Time series of a monitoring of weather radar data.\n' +
                '# Quantiles: '+str(quantiles[1])+', '+str(quantiles[0])+', ' +
                str(quantiles[2])+' percent.\n' +
                '# Data: '+generate_field_name_str(datatype)+'\n' +
                '# Fill Value: '+str(get_fillvalue())+'\n' +
                '# Start: '+start_time_aux[0].strftime(
                    '%Y-%m-%d %H:%M:%S UTC')+'\n' +
                '#\n')

            fieldnames = ['date', 'NP', 'central_quantile', 'low_quantile',
                          'high_quantile']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            for i, np_t_el in enumerate(np_t_aux):
                writer.writerow({
                    'date': start_time_aux[i].strftime('%Y%m%d%H%M%S'),
                    'NP': np_t_el,
                    'central_quantile': values_aux[i, 1],
                    'low_quantile': values_aux[i, 0],
                    'high_quantile': values_aux[i, 2]})

            fcntl.flock(csvfile, fcntl.LOCK_UN)
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            while True:
                try:
                    fcntl.flock(csvfile, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except OSError as e:
                    if e.errno != errno.EAGAIN:
                        raise
                    else:
                        time.sleep(0.1)

            fieldnames = ['date', 'NP', 'central_quantile', 'low_quantile',
                          'high_quantile']
            writer = csv.DictWriter(csvfile, fieldnames)
            for i, np_t_el in enumerate(np_t_aux):
                writer.writerow({
                    'date': start_time_aux[i].strftime('%Y%m%d%H%M%S'),
                    'NP': np_t_el,
                    'central_quantile': values_aux[i, 1],
                    'low_quantile': values_aux[i, 0],
                    'high_quantile': values_aux[i, 2]})
            fcntl.flock(csvfile, fcntl.LOCK_UN)
            csvfile.close()
    return fname


def write_excess_gates(excess_dict, fname):
    """
    Writes the position and values of gates that have a frequency of
    occurrence higher than a particular threshold

    Parameters
    ----------
    excess_dict : dict
        dictionary containing the gates parameters
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    ngates = len(excess_dict['ray_ind'])
    with open(fname, 'w', newline='') as csvfile:
        csvfile.write('# Gates exceeding '+str(excess_dict['quant_min']) +
                      ' percentile data file\n')
        csvfile.write('# Comment lines are preceded by "#"\n')
        csvfile.write(
            '# Data collection start time: ' +
            excess_dict['starttime'].strftime('%Y-%m-%d %H:%M:%S UTC')+'\n')
        csvfile.write(
            '# Data collection end time: ' +
            excess_dict['endtime'].strftime('%Y-%m-%d %H:%M:%S UTC')+'\n')
        csvfile.write('# Number of gates in file: '+str(ngates)+'\n')
        csvfile.write('#\n')

        fieldnames = [
            'ray_ind', 'rng_ind', 'ele', 'azi', 'rng', 'nsamples',
            'occurrence', 'freq_occu']
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for i, ray_ind in enumerate(excess_dict['ray_ind']):
            writer.writerow({
                'ray_ind': ray_ind,
                'rng_ind': excess_dict['rng_ind'][i],
                'ele': excess_dict['ele'][i],
                'azi': excess_dict['azi'][i],
                'rng': excess_dict['rng'][i],
                'nsamples': excess_dict['nsamples'][i],
                'occurrence': excess_dict['occurrence'][i],
                'freq_occu': excess_dict['freq_occu'][i]})

        csvfile.close()

    return fname


def write_intercomp_scores_ts(start_time, stats, field_name, fname,
                              rad1_name='RADAR001', rad2_name='RADAR002',
                              rewrite=False):
    """
    writes time series of radar intercomparison scores

    Parameters
    ----------
    start_time : datetime object or array of date time objects
        the time of the intercomparison
    stats : dict
        dictionary containing the statistics
    field_name : str
        The name of the field
    fname : str
        file name where to store the data
    rad1_name, rad2_name : str
        Name of the radars intercompared
    rewrite : bool
        if True a new file is created

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    nvalues = np.size(start_time)

    meanbias = stats['meanbias'].filled(fill_value=get_fillvalue())
    medianbias = stats['medianbias'].filled(fill_value=get_fillvalue())
    quant25bias = stats['quant25bias'].filled(fill_value=get_fillvalue())
    quant75bias = stats['quant75bias'].filled(fill_value=get_fillvalue())
    modebias = stats['modebias'].filled(fill_value=get_fillvalue())
    corr = stats['corr'].filled(fill_value=get_fillvalue())
    slope = stats['slope'].filled(fill_value=get_fillvalue())
    intercep = stats['intercep'].filled(fill_value=get_fillvalue())
    intercep_slope_1 = stats['intercep_slope_1'].filled(
        fill_value=get_fillvalue())

    if nvalues == 1:
        start_time_aux = np.asarray([start_time])
        meanbias = np.asarray([meanbias])
        medianbias = np.asarray([medianbias])
        quant25bias = np.asarray([quant25bias])
        quant75bias = np.asarray([quant75bias])
        modebias = np.asarray([modebias])
        corr = np.asarray([corr])
        slope = np.asarray([slope])
        intercep = np.asarray([intercep])
        intercep_slope_1 = np.asarray([intercep_slope_1])
        np_t = np.asarray([stats['npoints']])
    else:
        start_time_aux = np.asarray(start_time)
        np_t = stats['npoints']

    if rewrite:
        file_exists = False
    else:
        filelist = glob.glob(fname)
        if not filelist:
            file_exists = False
        else:
            file_exists = True

    if not file_exists:
        with open(fname, 'w', newline='') as csvfile:
            while True:
                try:
                    fcntl.flock(csvfile, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except OSError as e:
                    if e.errno != errno.EAGAIN:
                        raise
                    else:
                        time.sleep(0.1)

            csvfile.write(
                '# Weather radar intercomparison scores timeseries file\n' +
                '# Comment lines are preceded by "#"\n' +
                '# Description: \n' +
                '# Time series of the intercomparison between two radars.\n' +
                '# Radar 1: '+rad1_name+'\n' +
                '# Radar 2: '+rad2_name+'\n' +
                '# Field name: '+field_name+'\n' +
                '# Fill Value: '+str(get_fillvalue())+'\n' +
                '# Start: '+start_time_aux[0].strftime(
                    '%Y-%m-%d %H:%M:%S UTC')+'\n' +
                '#\n')

            fieldnames = ['date', 'NP', 'mean_bias', 'median_bias',
                          'quant25_bias', 'quant75_bias', 'mode_bias', 'corr',
                          'slope_of_linear_regression',
                          'intercep_of_linear_regression',
                          'intercep_of_linear_regression_of_slope_1']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()

            for i, dt in enumerate(start_time_aux):
                writer.writerow({
                    'date': dt.strftime('%Y%m%d%H%M%S'),
                    'NP': np_t[i],
                    'mean_bias': meanbias[i],
                    'median_bias': medianbias[i],
                    'quant25_bias': quant25bias[i],
                    'quant75_bias': quant75bias[i],
                    'mode_bias': modebias[i],
                    'corr': corr[i],
                    'slope_of_linear_regression': slope[i],
                    'intercep_of_linear_regression': intercep[i],
                    'intercep_of_linear_regression_of_slope_1': (
                        intercep_slope_1[i])
                    })

            fcntl.flock(csvfile, fcntl.LOCK_UN)
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            while True:
                try:
                    fcntl.flock(csvfile, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except OSError as e:
                    if e.errno != errno.EAGAIN:
                        raise
                    else:
                        time.sleep(0.1)

            fieldnames = ['date', 'NP', 'mean_bias', 'median_bias',
                          'quant25_bias', 'quant75_bias', 'mode_bias', 'corr',
                          'slope_of_linear_regression',
                          'intercep_of_linear_regression',
                          'intercep_of_linear_regression_of_slope_1']
            writer = csv.DictWriter(csvfile, fieldnames)
            for i, dt in enumerate(start_time_aux):
                writer.writerow({
                    'date': dt.strftime('%Y%m%d%H%M%S'),
                    'NP': np_t[i],
                    'mean_bias': meanbias[i],
                    'median_bias': medianbias[i],
                    'quant25_bias': quant25bias[i],
                    'quant75_bias': quant75bias[i],
                    'mode_bias': modebias[i],
                    'corr': corr[i],
                    'slope_of_linear_regression': slope[i],
                    'intercep_of_linear_regression': intercep[i],
                    'intercep_of_linear_regression_of_slope_1': (
                        intercep_slope_1[i])
                    })
            fcntl.flock(csvfile, fcntl.LOCK_UN)
            csvfile.close()

    return fname


def write_colocated_gates(coloc_gates, fname):
    """
    Writes the position of gates colocated with two radars

    Parameters
    ----------
    coloc_gates : dict
        dictionary containing the colocated gates parameters
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    with open(fname, 'w', newline='') as csvfile:
        csvfile.write('# Colocated radar gates data file\n')
        csvfile.write('# Comment lines are preceded by "#"\n')
        csvfile.write('#\n')

        fieldnames = [
            'rad1_ray_ind', 'rad1_rng_ind', 'rad1_ele', 'rad1_azi', 'rad1_rng',
            'rad2_ray_ind', 'rad2_rng_ind', 'rad2_ele', 'rad2_azi', 'rad2_rng']
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for i, rad1_ray_ind in enumerate(coloc_gates['rad1_ray_ind']):
            writer.writerow({
                'rad1_ray_ind': rad1_ray_ind,
                'rad1_rng_ind': coloc_gates['rad1_rng_ind'][i],
                'rad1_ele': coloc_gates['rad1_ele'][i],
                'rad1_azi': coloc_gates['rad1_azi'][i],
                'rad1_rng': coloc_gates['rad1_rng'][i],
                'rad2_ray_ind': coloc_gates['rad2_ray_ind'][i],
                'rad2_rng_ind': coloc_gates['rad2_rng_ind'][i],
                'rad2_ele': coloc_gates['rad2_ele'][i],
                'rad2_azi': coloc_gates['rad2_azi'][i],
                'rad2_rng': coloc_gates['rad2_rng'][i]})

        csvfile.close()

    return fname


def write_colocated_data(coloc_data, fname):
    """
    Writes the data of gates colocated with two radars

    Parameters
    ----------
    coloc_data : dict
        dictionary containing the colocated data parameters
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    filelist = glob.glob(fname)
    if not filelist:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write('# Colocated radar gates data file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('#\n')

            fieldnames = [
                'rad1_time', 'rad1_ray_ind', 'rad1_rng_ind', 'rad1_ele',
                'rad1_azi', 'rad1_rng', 'rad1_val', 'rad2_time',
                'rad2_ray_ind', 'rad2_rng_ind', 'rad2_ele', 'rad2_azi',
                'rad2_rng', 'rad2_val']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            for i, rad1_time in enumerate(coloc_data['rad1_time']):
                writer.writerow({
                    'rad1_time': rad1_time.strftime('%Y%m%d%H%M%S'),
                    'rad1_ray_ind': coloc_data['rad1_ray_ind'][i],
                    'rad1_rng_ind': coloc_data['rad1_rng_ind'][i],
                    'rad1_ele': coloc_data['rad1_ele'][i],
                    'rad1_azi': coloc_data['rad1_azi'][i],
                    'rad1_rng': coloc_data['rad1_rng'][i],
                    'rad1_val': coloc_data['rad1_val'][i],
                    'rad2_time': (
                        coloc_data['rad2_time'][i].strftime('%Y%m%d%H%M%S')),
                    'rad2_ray_ind': coloc_data['rad2_ray_ind'][i],
                    'rad2_rng_ind': coloc_data['rad2_rng_ind'][i],
                    'rad2_ele': coloc_data['rad2_ele'][i],
                    'rad2_azi': coloc_data['rad2_azi'][i],
                    'rad2_rng': coloc_data['rad2_rng'][i],
                    'rad2_val': coloc_data['rad2_val'][i]})

            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = [
                'rad1_time', 'rad1_ray_ind', 'rad1_rng_ind', 'rad1_ele',
                'rad1_azi', 'rad1_rng', 'rad1_val', 'rad2_time',
                'rad2_ray_ind', 'rad2_rng_ind', 'rad2_ele', 'rad2_azi',
                'rad2_rng', 'rad2_val']
            writer = csv.DictWriter(csvfile, fieldnames)
            for i, rad1_time in enumerate(coloc_data['rad1_time']):
                writer.writerow({
                    'rad1_time': rad1_time.strftime('%Y%m%d%H%M%S'),
                    'rad1_ray_ind': coloc_data['rad1_ray_ind'][i],
                    'rad1_rng_ind': coloc_data['rad1_rng_ind'][i],
                    'rad1_ele': coloc_data['rad1_ele'][i],
                    'rad1_azi': coloc_data['rad1_azi'][i],
                    'rad1_rng': coloc_data['rad1_rng'][i],
                    'rad1_val': coloc_data['rad1_val'][i],
                    'rad2_time': (
                        coloc_data['rad2_time'][i].strftime('%Y%m%d%H%M%S')),
                    'rad2_ray_ind': coloc_data['rad2_ray_ind'][i],
                    'rad2_rng_ind': coloc_data['rad2_rng_ind'][i],
                    'rad2_ele': coloc_data['rad2_ele'][i],
                    'rad2_azi': coloc_data['rad2_azi'][i],
                    'rad2_rng': coloc_data['rad2_rng'][i],
                    'rad2_val': coloc_data['rad2_val'][i]})
            csvfile.close()

    return fname


def write_colocated_data_time_avg(coloc_data, fname):
    """
    Writes the time averaged data of gates colocated with two radars

    Parameters
    ----------
    coloc_data : dict
        dictionary containing the colocated data parameters
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    filelist = glob.glob(fname)
    if not filelist:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write('# Colocated radar gates data file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('#\n')

            fieldnames = [
                'rad1_time', 'rad1_ray_ind', 'rad1_rng_ind', 'rad1_ele',
                'rad1_azi', 'rad1_rng', 'rad1_dBZavg', 'rad1_PhiDPavg',
                'rad1_Flagavg', 'rad2_time', 'rad2_ray_ind', 'rad2_rng_ind',
                'rad2_ele', 'rad2_azi', 'rad2_rng', 'rad2_dBZavg',
                'rad2_PhiDPavg', 'rad2_Flagavg']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            for i, rad1_time in enumerate(coloc_data['rad1_time']):
                writer.writerow({
                    'rad1_time': rad1_time.strftime('%Y%m%d%H%M%S'),
                    'rad1_ray_ind': coloc_data['rad1_ray_ind'][i],
                    'rad1_rng_ind': coloc_data['rad1_rng_ind'][i],
                    'rad1_ele': coloc_data['rad1_ele'][i],
                    'rad1_azi': coloc_data['rad1_azi'][i],
                    'rad1_rng': coloc_data['rad1_rng'][i],
                    'rad1_dBZavg': coloc_data['rad1_dBZavg'][i],
                    'rad1_PhiDPavg': coloc_data['rad1_PhiDPavg'][i],
                    'rad1_Flagavg': coloc_data['rad1_Flagavg'][i],
                    'rad2_time': (
                        coloc_data['rad2_time'][i].strftime('%Y%m%d%H%M%S')),
                    'rad2_ray_ind': coloc_data['rad2_ray_ind'][i],
                    'rad2_rng_ind': coloc_data['rad2_rng_ind'][i],
                    'rad2_ele': coloc_data['rad2_ele'][i],
                    'rad2_azi': coloc_data['rad2_azi'][i],
                    'rad2_rng': coloc_data['rad2_rng'][i],
                    'rad2_dBZavg': coloc_data['rad2_dBZavg'][i],
                    'rad2_PhiDPavg': coloc_data['rad2_PhiDPavg'][i],
                    'rad2_Flagavg': coloc_data['rad2_Flagavg'][i]})

            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = [
                'rad1_time', 'rad1_ray_ind', 'rad1_rng_ind', 'rad1_ele',
                'rad1_azi', 'rad1_rng', 'rad1_dBZavg', 'rad1_PhiDPavg',
                'rad1_Flagavg', 'rad2_time', 'rad2_ray_ind', 'rad2_rng_ind',
                'rad2_ele', 'rad2_azi', 'rad2_rng', 'rad2_dBZavg',
                'rad2_PhiDPavg', 'rad2_Flagavg']
            writer = csv.DictWriter(csvfile, fieldnames)
            for i, rad1_time in enumerate(coloc_data['rad1_time']):
                writer.writerow({
                    'rad1_time': rad1_time.strftime('%Y%m%d%H%M%S'),
                    'rad1_ray_ind': coloc_data['rad1_ray_ind'][i],
                    'rad1_rng_ind': coloc_data['rad1_rng_ind'][i],
                    'rad1_ele': coloc_data['rad1_ele'][i],
                    'rad1_azi': coloc_data['rad1_azi'][i],
                    'rad1_rng': coloc_data['rad1_rng'][i],
                    'rad1_dBZavg': coloc_data['rad1_dBZavg'][i],
                    'rad1_PhiDPavg': coloc_data['rad1_PhiDPavg'][i],
                    'rad1_Flagavg': coloc_data['rad1_Flagavg'][i],
                    'rad2_time': (
                        coloc_data['rad2_time'][i].strftime('%Y%m%d%H%M%S')),
                    'rad2_ray_ind': coloc_data['rad2_ray_ind'][i],
                    'rad2_rng_ind': coloc_data['rad2_rng_ind'][i],
                    'rad2_ele': coloc_data['rad2_ele'][i],
                    'rad2_azi': coloc_data['rad2_azi'][i],
                    'rad2_rng': coloc_data['rad2_rng'][i],
                    'rad2_dBZavg': coloc_data['rad2_dBZavg'][i],
                    'rad2_PhiDPavg': coloc_data['rad2_PhiDPavg'][i],
                    'rad2_Flagavg': coloc_data['rad2_Flagavg'][i]})
            csvfile.close()

    return fname


def write_sun_hits(sun_hits, fname):
    """
    Writes sun hits data.

    Parameters
    ----------
    sun_hits : dict
        dictionary containing the sun hits parameters

    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    dBm_sun_hit = sun_hits['dBm_sun_hit'].filled(fill_value=get_fillvalue())
    std_dBm_sun_hit = sun_hits['std(dBm_sun_hit)'].filled(
        fill_value=get_fillvalue())
    dBmv_sun_hit = sun_hits['dBmv_sun_hit'].filled(fill_value=get_fillvalue())
    std_dBmv_sun_hit = sun_hits['std(dBmv_sun_hit)'].filled(
        fill_value=get_fillvalue())
    zdr_sun_hit = sun_hits['ZDR_sun_hit'].filled(fill_value=get_fillvalue())
    std_zdr_sun_hit = sun_hits['std(ZDR_sun_hit)'].filled(
        fill_value=get_fillvalue())

    filelist = glob.glob(fname)
    if not filelist:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write('# Weather radar sun hits data file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('# Fill Value: '+str(get_fillvalue())+'\n')
            csvfile.write('#\n')

            fieldnames = [
                'time', 'ray', 'NPrng',
                'rad_el', 'rad_az', 'sun_el', 'sun_az',
                'dBm_sun_hit', 'std(dBm_sun_hit)', 'NPh', 'NPhval',
                'dBmv_sun_hit', 'std(dBmv_sun_hit)', 'NPv', 'NPvval',
                'ZDR_sun_hit', 'std(ZDR_sun_hit)', 'NPzdr', 'NPzdrval']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            for i, sh_time in enumerate(sun_hits['time']):
                writer.writerow({
                    'time': sh_time.strftime('%Y-%m-%d %H:%M:%S.%f'),
                    'ray': sun_hits['ray'][i],
                    'NPrng': sun_hits['NPrng'][i],
                    'rad_el': sun_hits['rad_el'][i],
                    'rad_az': sun_hits['rad_az'][i],
                    'sun_el': sun_hits['sun_el'][i],
                    'sun_az': sun_hits['sun_az'][i],
                    'dBm_sun_hit': dBm_sun_hit[i],
                    'std(dBm_sun_hit)': std_dBm_sun_hit[i],
                    'NPh': sun_hits['NPh'][i],
                    'NPhval': sun_hits['NPhval'][i],
                    'dBmv_sun_hit': dBmv_sun_hit[i],
                    'std(dBmv_sun_hit)': std_dBmv_sun_hit[i],
                    'NPv': sun_hits['NPv'][i],
                    'NPvval': sun_hits['NPvval'][i],
                    'ZDR_sun_hit': zdr_sun_hit[i],
                    'std(ZDR_sun_hit)': std_zdr_sun_hit[i],
                    'NPzdr': sun_hits['NPzdr'][i],
                    'NPzdrval': sun_hits['NPzdrval'][i]})
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = [
                'time', 'ray', 'NPrng',
                'rad_el', 'rad_az', 'sun_el', 'sun_az',
                'dBm_sun_hit', 'std(dBm_sun_hit)', 'NPh', 'NPhval',
                'dBmv_sun_hit', 'std(dBmv_sun_hit)', 'NPv', 'NPvval',
                'ZDR_sun_hit', 'std(ZDR_sun_hit)', 'NPzdr', 'NPzdrval']
            writer = csv.DictWriter(csvfile, fieldnames)
            for i, sh_time in enumerate(sun_hits['time']):
                writer.writerow({
                    'time': sh_time.strftime('%Y-%m-%d %H:%M:%S.%f'),
                    'ray': sun_hits['ray'][i],
                    'NPrng': sun_hits['NPrng'][i],
                    'rad_el': sun_hits['rad_el'][i],
                    'rad_az': sun_hits['rad_az'][i],
                    'sun_el': sun_hits['sun_el'][i],
                    'sun_az': sun_hits['sun_az'][i],
                    'dBm_sun_hit': dBm_sun_hit[i],
                    'std(dBm_sun_hit)': std_dBm_sun_hit[i],
                    'NPh': sun_hits['NPh'][i],
                    'NPhval': sun_hits['NPhval'][i],
                    'dBmv_sun_hit': dBmv_sun_hit[i],
                    'std(dBmv_sun_hit)': std_dBmv_sun_hit[i],
                    'NPv': sun_hits['NPv'][i],
                    'NPvval': sun_hits['NPvval'][i],
                    'ZDR_sun_hit': zdr_sun_hit[i],
                    'std(ZDR_sun_hit)': std_zdr_sun_hit[i],
                    'NPzdr': sun_hits['NPzdr'][i],
                    'NPzdrval': sun_hits['NPzdrval'][i]})
            csvfile.close()

    return fname


def write_sun_retrieval(sun_retrieval, fname):
    """
    Writes sun retrieval data.

    Parameters
    ----------
    sun_retrieval : dict
        dictionary containing the sun retrieval parameters

    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    first_hit_time = sun_retrieval['first_hit_time'].strftime('%Y%m%d%H%M%S')
    last_hit_time = sun_retrieval['last_hit_time'].strftime('%Y%m%d%H%M%S')
    el_width_h = sun_retrieval['el_width_h'].filled(fill_value=get_fillvalue())
    az_width_h = sun_retrieval['az_width_h'].filled(fill_value=get_fillvalue())
    el_bias_h = sun_retrieval['el_bias_h'].filled(fill_value=get_fillvalue())
    az_bias_h = sun_retrieval['az_bias_h'].filled(fill_value=get_fillvalue())
    dBm_sun_est = sun_retrieval['dBm_sun_est'].filled(
        fill_value=get_fillvalue())
    std_dBm_sun_est = sun_retrieval['std(dBm_sun_est)'].filled(
        fill_value=get_fillvalue())
    sf_h = sun_retrieval['sf_h'].filled(fill_value=get_fillvalue())

    el_width_v = sun_retrieval['el_width_v'].filled(fill_value=get_fillvalue())
    az_width_v = sun_retrieval['az_width_v'].filled(fill_value=get_fillvalue())
    el_bias_v = sun_retrieval['el_bias_v'].filled(fill_value=get_fillvalue())
    az_bias_v = sun_retrieval['az_bias_v'].filled(fill_value=get_fillvalue())
    dBmv_sun_est = sun_retrieval['dBmv_sun_est'].filled(
        fill_value=get_fillvalue())
    std_dBmv_sun_est = sun_retrieval['std(dBmv_sun_est)'].filled(
        fill_value=get_fillvalue())
    sf_v = sun_retrieval['sf_v'].filled(fill_value=get_fillvalue())

    zdr_sun_est = sun_retrieval['ZDR_sun_est'].filled(
        fill_value=get_fillvalue())
    std_zdr_sun_est = sun_retrieval['std(ZDR_sun_est)'].filled(
        fill_value=get_fillvalue())
    sf_ref = sun_retrieval['sf_ref'].filled(
        fill_value=get_fillvalue())
    ref_time = 'None'
    if sun_retrieval['ref_time'] is not None:
        ref_time = sun_retrieval['ref_time'].strftime('%Y%m%d%H%M%S')

    filelist = glob.glob(fname)
    if not filelist:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write('# Weather radar sun retrievals data file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('# Fill Value: '+str(get_fillvalue())+'\n')
            csvfile.write('#\n')

            fieldnames = [
                'first_hit_time', 'last_hit_time',
                'nhits_h', 'el_width_h', 'az_width_h',
                'el_bias_h', 'az_bias_h', 'dBm_sun_est', 'std(dBm_sun_est)',
                'sf_h',
                'nhits_v', 'el_width_v', 'az_width_v',
                'el_bias_v', 'az_bias_v', 'dBmv_sun_est', 'std(dBmv_sun_est)',
                'sf_v',
                'nhits_zdr', 'ZDR_sun_est', 'std(ZDR_sun_est)', 'sf_ref',
                'ref_time']

            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            writer.writerow(
                {'first_hit_time': first_hit_time,
                 'last_hit_time': last_hit_time,
                 'nhits_h': sun_retrieval['nhits_h'],
                 'el_width_h': el_width_h,
                 'az_width_h': az_width_h,
                 'el_bias_h': el_bias_h,
                 'az_bias_h': az_bias_h,
                 'dBm_sun_est': dBm_sun_est,
                 'std(dBm_sun_est)': std_dBm_sun_est,
                 'sf_h': sf_h,
                 'nhits_v': sun_retrieval['nhits_v'],
                 'el_width_v': el_width_v,
                 'az_width_v': az_width_v,
                 'el_bias_v': el_bias_v,
                 'az_bias_v': az_bias_v,
                 'dBmv_sun_est': dBmv_sun_est,
                 'std(dBmv_sun_est)': std_dBmv_sun_est,
                 'sf_v': sf_v,
                 'nhits_zdr': sun_retrieval['nhits_zdr'],
                 'ZDR_sun_est': zdr_sun_est,
                 'std(ZDR_sun_est)': std_zdr_sun_est,
                 'sf_ref': sf_ref,
                 'ref_time': ref_time})
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = [
                'first_hit_time', 'last_hit_time',
                'nhits_h', 'el_width_h', 'az_width_h',
                'el_bias_h', 'az_bias_h', 'dBm_sun_est', 'std(dBm_sun_est)',
                'sf_h',
                'nhits_v', 'el_width_v', 'az_width_v',
                'el_bias_v', 'az_bias_v', 'dBmv_sun_est', 'std(dBmv_sun_est)',
                'sf_v',
                'nhits_zdr', 'ZDR_sun_est', 'std(ZDR_sun_est)', 'sf_ref',
                'ref_time']

            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writerow(
                {'first_hit_time': first_hit_time,
                 'last_hit_time': last_hit_time,
                 'nhits_h': sun_retrieval['nhits_h'],
                 'el_width_h': el_width_h,
                 'az_width_h': az_width_h,
                 'el_bias_h': el_bias_h,
                 'az_bias_h': az_bias_h,
                 'dBm_sun_est': dBm_sun_est,
                 'std(dBm_sun_est)': std_dBm_sun_est,
                 'sf_h': sf_h,
                 'nhits_v': sun_retrieval['nhits_v'],
                 'el_width_v': el_width_v,
                 'az_width_v': az_width_v,
                 'el_bias_v': el_bias_v,
                 'az_bias_v': az_bias_v,
                 'dBmv_sun_est': dBmv_sun_est,
                 'std(dBmv_sun_est)': std_dBmv_sun_est,
                 'sf_v': sf_v,
                 'nhits_zdr': sun_retrieval['nhits_zdr'],
                 'ZDR_sun_est': zdr_sun_est,
                 'std(ZDR_sun_est)': std_zdr_sun_est,
                 'sf_ref': sf_ref,
                 'ref_time': ref_time})
            csvfile.close()

    return fname
