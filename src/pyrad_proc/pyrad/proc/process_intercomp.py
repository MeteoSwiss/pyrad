"""
pyrad.proc.process_intercomp
============================

Functions used in the inter-comparison between radars

.. autosummary::
    :toctree: generated/

    process_time_stats
    process_time_stats2
    process_time_avg
    process_weighted_time_avg
    process_time_avg_flag
    process_colocated_gates
    process_intercomp
    process_intercomp_time_avg
    process_fields_diff
    process_intercomp_fields

"""

from copy import deepcopy
from warnings import warn
import datetime
import numpy as np
import scipy
from netCDF4 import num2date

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..io.io_aux import get_save_dir, make_filename
from ..io.read_data_other import read_colocated_gates, read_colocated_data
from ..io.read_data_other import read_colocated_data_time_avg
from ..io.read_data_radar import interpol_field

from ..util.radar_utils import time_avg_range, get_range_bins_to_avg
from ..util.radar_utils import find_colocated_indexes


def process_time_stats(procstatus, dscfg, radar_list=None):
    """
    computes the temporal statistics of a field

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        period : float. Dataset keyword
            the period to average [s]. If -1 the statistics are going to be
            performed over the entire data. Default 3600.
        start_average : float. Dataset keyword
            when to start the average [s from midnight UTC]. Default 0.
        lin_trans: int. Dataset keyword
            If 1 apply linear transformation before averaging
        use_nan : bool. Dataset keyword
            If true non valid data will be used
        nan_value : float. Dataset keyword
            The value of the non valid data. Default 0
        stat: string. Dataset keyword
            Statistic to compute: Can be mean, std, cov, min, max. Default
            mean
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_name = get_fieldname_pyart(datatype)
        break
    ind_rad = int(radarnr[5:8])-1

    start_average = dscfg.get('start_average', 0.)
    period = dscfg.get('period', 3600.)
    lin_trans = dscfg.get('lin_trans', 0)
    use_nan = dscfg.get('use_nan', 0)
    nan_value = dscfg.get('nan_value', 0.)
    stat = dscfg.get('stat', 'mean')

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if field_name not in radar.fields:
            warn(field_name+' not available.')
            return None, None

        # Prepare auxiliary radar
        field = deepcopy(radar.fields[field_name])
        if stat in ('mean', 'std', 'cov'):
            if lin_trans:
                field['data'] = np.ma.power(10., 0.1*field['data'])

            if use_nan:
                field['data'] = np.ma.asarray(field['data'].filled(nan_value))

            if stat in ('std', 'cov'):
                sum2_dict = pyart.config.get_metadata('sum_squared')
                sum2_dict['data'] = field['data']*field['data']
        else:
            if use_nan:
                field['data'] = np.ma.asarray(field['data'].filled(nan_value))

        npoints_dict = pyart.config.get_metadata('number_of_samples')
        npoints_dict['data'] = np.ma.asarray(
            np.logical_not(np.ma.getmaskarray(field['data'])), dtype=int)

        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()
        radar_aux.add_field(field_name, field)
        radar_aux.add_field('number_of_samples', npoints_dict)

        if stat in ('std', 'cov'):
            radar_aux.add_field('sum_squared', sum2_dict)

        # first volume: initialize start and end time of averaging
        if dscfg['initialized'] == 0:
            avg_par = dict()
            if period != -1:
                date_00 = dscfg['timeinfo'].replace(
                    hour=0, minute=0, second=0, microsecond=0)

                avg_par.update(
                    {'starttime': date_00+datetime.timedelta(
                        seconds=start_average)})
                avg_par.update(
                    {'endtime': avg_par['starttime']+datetime.timedelta(
                        seconds=period)})
            else:
                avg_par.update({'starttime': dscfg['timeinfo']})
                avg_par.update({'endtime': dscfg['timeinfo']})

            avg_par.update({'timeinfo': dscfg['timeinfo']})
            dscfg['global_data'] = avg_par
            dscfg['initialized'] = 1

        if dscfg['initialized'] == 0:
            return None, None

        dscfg['global_data']['timeinfo'] = dscfg['timeinfo']
        # no radar object in global data: create it
        if 'radar_out' not in dscfg['global_data']:
            if period != -1:
                # get start and stop times of new radar object
                (dscfg['global_data']['starttime'],
                 dscfg['global_data']['endtime']) = (
                     time_avg_range(
                         dscfg['timeinfo'], dscfg['global_data']['starttime'],
                         dscfg['global_data']['endtime'], period))

                # check if volume time older than starttime
                if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
                    dscfg['global_data'].update({'radar_out': radar_aux})
            else:
                dscfg['global_data'].update({'radar_out': radar_aux})

            return None, None

        # still accumulating: add field to global field
        if (period == -1 or
                dscfg['timeinfo'] < dscfg['global_data']['endtime']):

            if period == -1:
                dscfg['global_data']['endtime'] = dscfg['timeinfo']

            field_interp = interpol_field(
                dscfg['global_data']['radar_out'], radar_aux, field_name)
            npoints_interp = interpol_field(
                dscfg['global_data']['radar_out'], radar_aux,
                'number_of_samples')

            if use_nan:
                field_interp['data'] = np.ma.asarray(
                    field_interp['data'].filled(nan_value))
                dscfg['global_data']['radar_out'].fields[
                    'number_of_samples']['data'] += np.ma.asarray(
                        npoints_interp['data'].filled(fill_value=1),
                        dtype=int)
            else:
                dscfg['global_data']['radar_out'].fields[
                    'number_of_samples']['data'] += np.ma.asarray(
                        npoints_interp['data'].filled(fill_value=0),
                        dtype=int)

            if stat in ('mean', 'std', 'cov'):
                masked_sum = np.ma.getmaskarray(
                    dscfg['global_data']['radar_out'].fields[
                        field_name]['data'])
                valid_sum = np.logical_and(
                    np.logical_not(masked_sum),
                    np.logical_not(np.ma.getmaskarray(field_interp['data'])))

                dscfg['global_data']['radar_out'].fields[
                    field_name]['data'][masked_sum] = (
                        field_interp['data'][masked_sum])

                dscfg['global_data']['radar_out'].fields[
                    field_name]['data'][valid_sum] += (
                        field_interp['data'][valid_sum])

                if stat in ('cov', 'std'):
                    dscfg['global_data']['radar_out'].fields[
                        'sum_squared']['data'][masked_sum] = (
                            field_interp['data'][masked_sum] *
                            field_interp['data'][masked_sum])

                    dscfg['global_data']['radar_out'].fields[
                        'sum_squared']['data'][valid_sum] += (
                            field_interp['data'][valid_sum] *
                            field_interp['data'][valid_sum])

            elif stat == 'max':
                dscfg['global_data']['radar_out'].fields[
                    field_name]['data'] = np.maximum(
                        dscfg['global_data']['radar_out'].fields[
                            field_name]['data'].filled(fill_value=-1.e300),
                        field_interp['data'].filled(fill_value=-1.e300))

                dscfg['global_data']['radar_out'].fields[
                    field_name]['data'] = np.ma.masked_values(
                        dscfg['global_data']['radar_out'].fields[
                            field_name]['data'], -1.e300)
            elif stat == 'min':
                dscfg['global_data']['radar_out'].fields[
                    field_name]['data'] = np.minimum(
                        dscfg['global_data']['radar_out'].fields[
                            field_name]['data'].filled(fill_value=1.e300),
                        field_interp['data'].filled(fill_value=1.e300))

                dscfg['global_data']['radar_out'].fields[
                    field_name]['data'] = np.ma.masked_values(
                        dscfg['global_data']['radar_out'].fields[
                            field_name]['data'], 1.e300)

            return None, None

        # we have reached the end of the accumulation period: do the averaging
        # and start a new object (only reachable if period != -1)
        if stat in ('mean', 'std', 'cov'):
            field_mean = (
                dscfg['global_data']['radar_out'].fields[field_name]['data'] /
                dscfg['global_data']['radar_out'].fields[
                    'number_of_samples']['data'])

            if stat == 'mean':
                if lin_trans:
                    dscfg['global_data']['radar_out'].fields[
                        field_name]['data'] = 10.*np.ma.log10(field_mean)
                else:
                    dscfg['global_data']['radar_out'].fields[
                        field_name]['data'] = field_mean
            elif stat in ('std', 'cov'):
                field_std = np.ma.sqrt(
                    dscfg['global_data']['radar_out'].fields[
                        'sum_squared']['data'] /
                    dscfg['global_data']['radar_out'].fields[
                        'number_of_samples']['data']-field_mean*field_mean)

                if stat == 'std':
                    if lin_trans:
                        dscfg['global_data']['radar_out'].fields[
                            field_name]['data'] = 10.*np.ma.log10(field_std)
                    else:
                        dscfg['global_data']['radar_out'].fields[
                            field_name]['data'] = field_std
                else:
                    if lin_trans:
                        dscfg['global_data']['radar_out'].fields[
                            field_name]['data'] = 10.*np.ma.log10(
                                field_std/field_mean)
                    else:
                        dscfg['global_data']['radar_out'].fields[
                            field_name]['data'] = field_std/field_mean

        new_dataset = {
            'radar_out': deepcopy(dscfg['global_data']['radar_out']),
            'timeinfo': dscfg['global_data']['endtime']}

        dscfg['global_data']['starttime'] += datetime.timedelta(
            seconds=period)
        dscfg['global_data']['endtime'] += datetime.timedelta(seconds=period)

        # remove old radar object from global_data dictionary
        dscfg['global_data'].pop('radar_out', None)

        # get start and stop times of new radar object
        dscfg['global_data']['starttime'], dscfg['global_data']['endtime'] = (
            time_avg_range(
                dscfg['timeinfo'], dscfg['global_data']['starttime'],
                dscfg['global_data']['endtime'], period))

        # check if volume time older than starttime
        if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
            dscfg['global_data'].update({'radar_out': radar_aux})

        return new_dataset, ind_rad

    # no more files to process if there is global data pack it up
    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None, None
        if 'radar_out' not in dscfg['global_data']:
            return None, None

        if stat in ('mean', 'std', 'cov'):
            field_mean = (
                dscfg['global_data']['radar_out'].fields[field_name]['data'] /
                dscfg['global_data']['radar_out'].fields[
                    'number_of_samples']['data'])

            if stat == 'mean':
                if lin_trans:
                    dscfg['global_data']['radar_out'].fields[
                        field_name]['data'] = 10.*np.ma.log10(field_mean)
                else:
                    dscfg['global_data']['radar_out'].fields[
                        field_name]['data'] = field_mean

            elif stat in ('std', 'cov'):
                field_std = np.ma.sqrt(
                    dscfg['global_data']['radar_out'].fields[
                        'sum_squared']['data'] /
                    dscfg['global_data']['radar_out'].fields[
                        'number_of_samples']['data']-field_mean*field_mean)
                if stat == 'std':
                    if lin_trans:
                        dscfg['global_data']['radar_out'].fields[
                            field_name]['data'] = 10.*np.ma.log10(field_std)
                    else:
                        dscfg['global_data']['radar_out'].fields[
                            field_name]['data'] = field_std
                else:
                    dscfg['global_data']['radar_out'].fields[
                        field_name]['data'] = field_std/field_mean

        new_dataset = {
            'radar_out': deepcopy(dscfg['global_data']['radar_out']),
            'timeinfo': dscfg['global_data']['endtime']}

        return new_dataset, ind_rad


def process_time_stats2(procstatus, dscfg, radar_list=None):
    """
    computes the temporal mean of a field

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        period : float. Dataset keyword
            the period to average [s]. If -1 the statistics are going to be
            performed over the entire data. Default 3600.
        start_average : float. Dataset keyword
            when to start the average [s from midnight UTC]. Default 0.
        stat: string. Dataset keyword
            Statistic to compute: Can be median, mode, percentileXX
        use_nan : bool. Dataset keyword
            If true non valid data will be used
        nan_value : float. Dataset keyword
            The value of the non valid data. Default 0
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_name = get_fieldname_pyart(datatype)
        break
    ind_rad = int(radarnr[5:8])-1

    start_average = dscfg.get('start_average', 0.)
    period = dscfg.get('period', 3600.)
    use_nan = dscfg.get('use_nan', 0)
    nan_value = dscfg.get('nan_value', 0.)
    stat = dscfg.get('stat', 'median')
    if 'percentile' in stat:
        percentile = float(stat.replace('percentile', ''))

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if field_name not in radar.fields:
            warn(field_name+' not available.')
            return None, None

        # prepare auxiliary radar
        field = deepcopy(radar.fields[field_name])
        if use_nan:
            field['data'] = np.ma.asarray(field['data'].filled(nan_value))
        npoints_dict = pyart.config.get_metadata('number_of_samples')
        npoints_dict['data'] = np.ma.asarray(
            np.logical_not(np.ma.getmaskarray(field['data'])), dtype=int)

        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()
        radar_aux.add_field(field_name, field)
        radar_aux.add_field('number_of_samples', npoints_dict)

        # first volume: initialize start and end time of averaging
        if dscfg['initialized'] == 0:
            avg_par = dict()
            if period != -1:
                date_00 = dscfg['timeinfo'].replace(
                    hour=0, minute=0, second=0, microsecond=0)

                avg_par.update(
                    {'starttime': date_00+datetime.timedelta(
                        seconds=start_average)})
                avg_par.update(
                    {'endtime': avg_par['starttime']+datetime.timedelta(
                        seconds=period)})
            else:
                avg_par.update({'starttime': dscfg['timeinfo']})
                avg_par.update({'endtime': dscfg['timeinfo']})

            avg_par.update({'timeinfo': dscfg['timeinfo']})
            dscfg['global_data'] = avg_par
            dscfg['initialized'] = 1

        if dscfg['initialized'] == 0:
            return None, None

        dscfg['global_data']['timeinfo'] = dscfg['timeinfo']
        # no radar object in global data: create it
        if 'radar_out' not in dscfg['global_data']:
            if period != -1:
                # get start and stop times of new radar object
                (dscfg['global_data']['starttime'],
                 dscfg['global_data']['endtime']) = (
                     time_avg_range(
                         dscfg['timeinfo'], dscfg['global_data']['starttime'],
                         dscfg['global_data']['endtime'], period))

                # check if volume time older than starttime
                if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
                    dscfg['global_data'].update({'radar_out': radar_aux})
                    dscfg['global_data'].update(
                        {'field_data': np.atleast_3d(
                            radar_aux.fields[field_name]['data'])})
            else:
                dscfg['global_data'].update({'radar_out': radar_aux})
                dscfg['global_data'].update(
                    {'field_data': np.atleast_3d(
                        radar_aux.fields[field_name]['data'])})

            return None, None

        # still accumulating: add field to global field
        if (period == -1 or
                dscfg['timeinfo'] < dscfg['global_data']['endtime']):

            if period == -1:
                dscfg['global_data']['endtime'] = dscfg['timeinfo']

            field_interp = interpol_field(
                dscfg['global_data']['radar_out'], radar_aux, field_name)
            npoints_interp = interpol_field(
                dscfg['global_data']['radar_out'], radar_aux,
                'number_of_samples')

            if use_nan:
                field_interp['data'] = np.ma.asarray(
                    field_interp['data'].filled(nan_value))
                dscfg['global_data']['radar_out'].fields[
                    'number_of_samples']['data'] += np.ma.asarray(
                        npoints_interp['data'].filled(fill_value=1),
                        dtype=int)
            else:
                dscfg['global_data']['radar_out'].fields[
                    'number_of_samples']['data'] += np.ma.asarray(
                        npoints_interp['data'].filled(fill_value=0),
                        dtype=int)

            dscfg['global_data']['field_data'] = np.ma.append(
                dscfg['global_data']['field_data'],
                np.atleast_3d(field_interp['data']), axis=2)

            return None, None

        # we have reached the end of the accumulation period: do the averaging
        # and start a new object (only reachable if period != -1)
        if stat == 'median':
            dscfg['global_data']['radar_out'].fields[
                field_name]['data'] = np.ma.median(
                    dscfg['global_data']['field_data'], axis=2)
        elif stat == 'mode':
            mode_data, _ = scipy.stats.mode(
                dscfg['global_data']['field_data'].filled(fill_value=np.nan),
                axis=2, nan_policy='omit')
            dscfg['global_data']['radar_out'].fields[field_name]['data'] = (
                np.ma.masked_invalid(np.squeeze(mode_data, axis=2)))
        elif 'percentile' in stat:
            percent_data = np.nanpercentile(
                dscfg['global_data']['field_data'].filled(fill_value=np.nan),
                percentile, axis=2)
            dscfg['global_data']['radar_out'].fields[field_name]['data'] = (
                np.ma.masked_invalid(percent_data))

        new_dataset = {
            'radar_out': deepcopy(dscfg['global_data']['radar_out']),
            'timeinfo': dscfg['global_data']['endtime']}

        dscfg['global_data']['starttime'] += datetime.timedelta(
            seconds=period)
        dscfg['global_data']['endtime'] += datetime.timedelta(seconds=period)

        # remove old radar object from global_data dictionary
        dscfg['global_data'].pop('radar_out', None)

        # get start and stop times of new radar object
        dscfg['global_data']['starttime'], dscfg['global_data']['endtime'] = (
            time_avg_range(
                dscfg['timeinfo'], dscfg['global_data']['starttime'],
                dscfg['global_data']['endtime'], period))

        # check if volume time older than starttime
        if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
            dscfg['global_data'].update({'radar_out': radar_aux})

        return new_dataset, ind_rad

    # no more files to process if there is global data pack it up
    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None, None
        if 'radar_out' not in dscfg['global_data']:
            return None, None

        if stat == 'median':
            dscfg['global_data']['radar_out'].fields[field_name]['data'] = (
                np.ma.median(dscfg['global_data']['field_data'], axis=2))
        elif stat == 'mode':
            mode_data, _ = scipy.stats.mode(
                dscfg['global_data']['field_data'].filled(fill_value=np.nan),
                axis=2, nan_policy='omit')
            dscfg['global_data']['radar_out'].fields[field_name]['data'] = (
                np.ma.masked_invalid(np.squeeze(mode_data, axis=2)))
        elif 'percentile' in stat:
            percent_data = np.nanpercentile(
                dscfg['global_data']['field_data'].filled(fill_value=np.nan),
                percentile, axis=2)
            dscfg['global_data']['radar_out'].fields[field_name]['data'] = (
                np.ma.masked_invalid(percent_data))

        new_dataset = {
            'radar_out': deepcopy(dscfg['global_data']['radar_out']),
            'timeinfo': dscfg['global_data']['endtime']}

        return new_dataset, ind_rad


def process_time_avg(procstatus, dscfg, radar_list=None):
    """
    computes the temporal mean of a field

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        period : float. Dataset keyword
            the period to average [s]. Default 3600.
        start_average : float. Dataset keyword
            when to start the average [s from midnight UTC]. Default 0.
        lin_trans: int. Dataset keyword
            If 1 apply linear transformation before averaging
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_name = get_fieldname_pyart(datatype)
        break
    ind_rad = int(radarnr[5:8])-1

    lin_trans = dscfg.get('lin_trans', 0)

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if field_name not in radar.fields:
            warn(field_name+' not available.')
            return None, None

        period = dscfg.get('period', 3600.)

        field = deepcopy(radar.fields[field_name])
        if lin_trans:
            field['data'] = np.ma.power(10., 0.1*field['data'])

        field['data'] = field['data'].filled(fill_value=0.)
        field['data'] = np.ma.asarray(field['data'])

        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()
        radar_aux.add_field(field_name, field)
        npoints_dict = pyart.config.get_metadata('number_of_samples')
        npoints_dict['data'] = np.ma.ones(
            (radar.nrays, radar.ngates), dtype=int)
        radar_aux.add_field('number_of_samples', npoints_dict)

        # first volume: initialize start and end time of averaging
        if dscfg['initialized'] == 0:
            start_average = dscfg.get('start_average', 0.)

            date_00 = dscfg['timeinfo'].replace(
                hour=0, minute=0, second=0, microsecond=0)

            avg_par = dict()
            avg_par.update(
                {'starttime': date_00+datetime.timedelta(
                    seconds=start_average)})
            avg_par.update(
                {'endtime': avg_par['starttime']+datetime.timedelta(
                    seconds=period)})
            avg_par.update({'timeinfo': dscfg['timeinfo']})
            dscfg['global_data'] = avg_par
            dscfg['initialized'] = 1

        if dscfg['initialized'] == 0:
            return None, None

        dscfg['global_data']['timeinfo'] = dscfg['timeinfo']
        # no radar object in global data: create it
        if 'radar_out' not in dscfg['global_data']:
            # get start and stop times of new radar object
            (dscfg['global_data']['starttime'],
             dscfg['global_data']['endtime']) = (
                 time_avg_range(
                     dscfg['timeinfo'], dscfg['global_data']['starttime'],
                     dscfg['global_data']['endtime'], period))

            # check if volume time older than starttime
            if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
                dscfg['global_data'].update({'radar_out': radar_aux})

            return None, None

        # still accumulating: add field to global field
        if dscfg['timeinfo'] < dscfg['global_data']['endtime']:
            field_interp = interpol_field(
                dscfg['global_data']['radar_out'], radar_aux, field_name)
            npoints_interp = interpol_field(
                dscfg['global_data']['radar_out'], radar_aux,
                'number_of_samples')
            dscfg['global_data']['radar_out'].fields[field_name]['data'] += (
                field_interp['data'].filled(fill_value=0))
            dscfg['global_data']['radar_out'].fields[
                'number_of_samples']['data'] += (
                    npoints_interp['data'].filled(fill_value=0)).astype('int')

            return None, None

        # we have reached the end of the accumulation period: do the averaging
        # and start a new object
        dscfg['global_data']['radar_out'].fields[field_name]['data'] /= (
            dscfg['global_data']['radar_out'].fields[
                'number_of_samples']['data'])
        if lin_trans:
            dscfg['global_data']['radar_out'].fields[field_name]['data'] = (
                10.*np.ma.log10(
                    dscfg['global_data']['radar_out'].fields[
                        field_name]['data']))

        new_dataset = {
            'radar_out': deepcopy(dscfg['global_data']['radar_out']),
            'timeinfo': dscfg['global_data']['endtime']}

        dscfg['global_data']['starttime'] += datetime.timedelta(
            seconds=period)
        dscfg['global_data']['endtime'] += datetime.timedelta(seconds=period)

        # remove old radar object from global_data dictionary
        dscfg['global_data'].pop('radar_out', None)

        # get start and stop times of new radar object
        dscfg['global_data']['starttime'], dscfg['global_data']['endtime'] = (
            time_avg_range(
                dscfg['timeinfo'], dscfg['global_data']['starttime'],
                dscfg['global_data']['endtime'], period))

        # check if volume time older than starttime
        if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
            dscfg['global_data'].update({'radar_out': radar_aux})

        return new_dataset, ind_rad

    # no more files to process if there is global data pack it up
    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None, None
        if 'radar_out' not in dscfg['global_data']:
            return None, None

        (dscfg['global_data']['radar_out'].fields[field_name][
            'data']) /= (
                dscfg['global_data']['radar_out'].fields[
                    'number_of_samples']['data'])
        if lin_trans:
            dscfg['global_data']['radar_out'].fields[field_name]['data'] = (
                10.*np.ma.log10(
                    dscfg['global_data']['radar_out'].fields[
                        field_name]['data']))

        new_dataset = {
            'radar_out': deepcopy(dscfg['global_data']['radar_out']),
            'timeinfo': dscfg['global_data']['endtime']}

        return new_dataset, ind_rad


def process_weighted_time_avg(procstatus, dscfg, radar_list=None):
    """
    computes the temporal mean of a field weighted by the reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        period : float. Dataset keyword
            the period to average [s]. Default 3600.
        start_average : float. Dataset keyword
            when to start the average [s from midnight UTC]. Default 0.
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : Radar
        radar object
    ind_rad : int
        radar index

    """
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ('dBZ', 'dBZc', 'dBuZ', 'dBZv', 'dBZvc', 'dBuZv'):
            refl_name = get_fieldname_pyart(datatype)
        else:
            field_name = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]
        if field_name not in radar.fields or refl_name not in radar.fields:
            warn('Unable to compute weighted average. Missing data')
            return None, None

        period = dscfg.get('period', 3600.)

        field = deepcopy(radar.fields[field_name])
        field['data'] = field['data'].filled(fill_value=0.)
        field['data'] = np.ma.asarray(field['data'])

        refl_field = deepcopy(radar.fields[refl_name])
        refl_field['data'] = np.ma.power(10., 0.1*refl_field['data'])
        refl_field['data'] = refl_field['data'].filled(fill_value=0.)
        refl_field['data'] = np.ma.asarray(refl_field['data'])

        field['data'] *= refl_field['data']

        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()
        radar_aux.add_field(field_name, field)
        radar_aux.add_field(refl_name, refl_field)

        # first volume: initialize start and end time of averaging
        if dscfg['initialized'] == 0:
            start_average = dscfg.get('start_average', 0.)

            date_00 = dscfg['timeinfo'].replace(
                hour=0, minute=0, second=0, microsecond=0)

            avg_par = dict()
            avg_par.update(
                {'starttime': date_00+datetime.timedelta(
                    seconds=start_average)})
            avg_par.update(
                {'endtime': avg_par['starttime']+datetime.timedelta(
                    seconds=period)})
            avg_par.update({'timeinfo': dscfg['timeinfo']})
            dscfg['global_data'] = avg_par
            dscfg['initialized'] = 1

        if dscfg['initialized'] == 0:
            return None, None

        dscfg['global_data']['timeinfo'] = dscfg['timeinfo']
        # no radar object in global data: create it
        if 'radar_out' not in dscfg['global_data']:
            # get start and stop times of new radar object
            (dscfg['global_data']['starttime'],
             dscfg['global_data']['endtime']) = (
                 time_avg_range(
                     dscfg['timeinfo'], dscfg['global_data']['starttime'],
                     dscfg['global_data']['endtime'], period))

            # check if volume time older than starttime
            if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
                dscfg['global_data'].update({'radar_out': radar_aux})

            return None, None

        # still accumulating: add field to global field
        if dscfg['timeinfo'] < dscfg['global_data']['endtime']:
            field_interp = interpol_field(
                dscfg['global_data']['radar_out'], radar_aux, field_name)
            dscfg['global_data']['radar_out'].fields[field_name]['data'] += (
                field_interp['data'].filled(fill_value=0))

            refl_interp = interpol_field(
                dscfg['global_data']['radar_out'], radar_aux, refl_name)
            dscfg['global_data']['radar_out'].fields[refl_name]['data'] += (
                refl_interp['data'].filled(fill_value=0))

            return None, None

        # we have reached the end of the accumulation period: do the averaging
        # and start a new object
        dscfg['global_data']['radar_out'].fields[field_name]['data'] /= (
            dscfg['global_data']['radar_out'].fields[refl_name]['data'])

        new_dataset = {
            'radar_out': deepcopy(dscfg['global_data']['radar_out']),
            'timeinfo': dscfg['global_data']['endtime']}

        dscfg['global_data']['starttime'] += datetime.timedelta(
            seconds=period)
        dscfg['global_data']['endtime'] += datetime.timedelta(seconds=period)

        # remove old radar object from global_data dictionary
        dscfg['global_data'].pop('radar_out', None)

        # get start and stop times of new radar object
        dscfg['global_data']['starttime'], dscfg['global_data']['endtime'] = (
            time_avg_range(
                dscfg['timeinfo'], dscfg['global_data']['starttime'],
                dscfg['global_data']['endtime'], period))

        # check if volume time older than starttime
        if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
            dscfg['global_data'].update({'radar_out': radar_aux})

        return new_dataset, ind_rad

    # no more files to process if there is global data pack it up
    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None, None
        if 'radar_out' not in dscfg['global_data']:
            return None, None

        dscfg['global_data']['radar_out'].fields[field_name]['data'] /= (
            dscfg['global_data']['radar_out'].fields[refl_name]['data'])

        new_dataset = {
            'radar_out': deepcopy(dscfg['global_data']['radar_out']),
            'timeinfo': dscfg['global_data']['endtime']}

        return new_dataset, ind_rad


def process_time_avg_flag(procstatus, dscfg, radar_list=None):
    """
    computes a flag field describing the conditions of the data used while
    averaging

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        period : float. Dataset keyword
            the period to average [s]. Default 3600.
        start_average : float. Dataset keyword
            when to start the average [s from midnight UTC]. Default 0.
        phidpmax: float. Dataset keyword
            maximum PhiDP
        beamwidth : float. Dataset keyword
            the antenna beamwidth [deg]. If None that of the keys
            radar_beam_width_h or radar_beam_width_v in attribute
            instrument_parameters of the radar object will be used. If the key
            or the attribute are not present the beamwidth will be set to None
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : Radar
        radar object
    ind_rad : int
        radar index

    """
    temp_name = None
    hydro_name = None
    iso0_name = None
    echo_name = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ('PhiDP', 'PhiDPc'):
            phidp_name = get_fieldname_pyart(datatype)
        elif datatype == 'echoID':
            echo_name = get_fieldname_pyart(datatype)
        elif datatype == 'hydro':
            hydro_name = get_fieldname_pyart(datatype)
        elif datatype == 'TEMP':
            temp_name = get_fieldname_pyart(datatype)
        elif datatype == 'H_ISO0':
            iso0_name = 'height_over_iso0'

    ind_rad = int(radarnr[5:8])-1

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        phidpmax = dscfg.get('phidpmax', 60.)
        period = dscfg.get('period', 3600.)

        time_avg_flag = pyart.config.get_metadata('time_avg_flag')
        time_avg_flag['data'] = np.ma.zeros(
            (radar.nrays, radar.ngates), dtype=int)

        if phidp_name not in radar.fields:
            warn('Missing PhiDP data')
            time_avg_flag['data'] += 1
        else:
            phidp_field = radar.fields[phidp_name]
            time_avg_flag['data'][phidp_field['data'] > phidpmax] += 1

        if echo_name is not None:
            if echo_name not in radar.fields:
                warn('Missing echo ID data')
                time_avg_flag['data'] += 100
            else:
                echo_field = radar.fields[echo_name]
                time_avg_flag['data'][echo_field['data'] == 2] += 100

        if hydro_name is not None and echo_name is not None:
            if ((hydro_name not in radar.fields) or
                    (echo_name not in radar.fields)):
                warn('Missing hydrometeor classification data')
                time_avg_flag['data'] += 10000
            else:
                hydro_field = radar.fields[hydro_name]
                # check where is no rain
                is_not_rain = np.logical_and(
                    hydro_field['data'] != 4, hydro_field['data'] != 6)
                # where is no rain should be precip
                is_not_rain = np.logical_and(
                    is_not_rain, echo_field['data'] == 3)
                time_avg_flag['data'][is_not_rain] += 10000
        elif temp_name is not None:
            if temp_name not in radar.fields:
                warn('Missing temperature data')
                time_avg_flag['data'] += 10000
            else:
                beamwidth = dscfg.get('beamwidth', None)
                if beamwidth is None:
                    if radar.instrument_parameters is not None:
                        if ('radar_beam_width_h' in
                                radar.instrument_parameters):
                            beamwidth = radar.instrument_parameters[
                                'radar_beam_width_h']['data'][0]
                        elif ('radar_beam_width_v' in
                              radar.instrument_parameters):
                            beamwidth = radar.instrument_parameters[
                                'radar_beam_width_v']['data'][0]
                if beamwidth is None:
                    warn('Antenna beam width unknown.')

                mask_fzl, _ = pyart.correct.get_mask_fzl(
                    radar, fzl=None, doc=None, min_temp=0., max_h_iso0=0.,
                    thickness=700., beamwidth=beamwidth,
                    temp_field=temp_name, iso0_field=iso0_name,
                    temp_ref='temperature')
                time_avg_flag['data'][mask_fzl] += 10000
        elif iso0_name is not None:
            if iso0_name not in radar.fields:
                warn('Missing height relative to iso0 data')
                time_avg_flag['data'] += 10000
            else:
                beamwidth = dscfg.get('beamwidth', None)
                if beamwidth is None:
                    if radar.instrument_parameters is not None:
                        if ('radar_beam_width_h' in
                                radar.instrument_parameters):
                            beamwidth = radar.instrument_parameters[
                                'radar_beam_width_h']['data'][0]
                        elif ('radar_beam_width_v' in
                              radar.instrument_parameters):
                            beamwidth = radar.instrument_parameters[
                                'radar_beam_width_v']['data'][0]
                if beamwidth is None:
                    warn('Antenna beam width unknown.')

                mask_fzl, _ = pyart.correct.get_mask_fzl(
                    radar, fzl=None, doc=None, min_temp=0., max_h_iso0=0.,
                    thickness=700., beamwidth=beamwidth,
                    temp_field=temp_name, iso0_field=iso0_name,
                    temp_ref='height_over_iso0')
                time_avg_flag['data'][mask_fzl] += 10000

        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()
        radar_aux.add_field('time_avg_flag', time_avg_flag)

        # first volume: initialize start and end time of averaging
        if dscfg['initialized'] == 0:
            start_average = dscfg.get('start_average', 0.)

            date_00 = dscfg['timeinfo'].replace(
                hour=0, minute=0, second=0, microsecond=0)

            avg_par = dict()
            avg_par.update(
                {'starttime': date_00+datetime.timedelta(
                    seconds=start_average)})
            avg_par.update(
                {'endtime': avg_par['starttime']+datetime.timedelta(
                    seconds=period)})
            avg_par.update({'timeinfo': dscfg['timeinfo']})
            dscfg['global_data'] = avg_par
            dscfg['initialized'] = 1

        if dscfg['initialized'] == 0:
            return None, None

        dscfg['global_data']['timeinfo'] = dscfg['timeinfo']
        # no radar object in global data: create it
        if 'radar_out' not in dscfg['global_data']:
            # get start and stop times of new radar object
            (dscfg['global_data']['starttime'],
             dscfg['global_data']['endtime']) = (
                 time_avg_range(
                     dscfg['timeinfo'], dscfg['global_data']['starttime'],
                     dscfg['global_data']['endtime'], period))

            # check if volume time older than starttime
            if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
                dscfg['global_data'].update({'radar_out': radar_aux})

            return None, None

        # still accumulating: add field to global field
        if dscfg['timeinfo'] < dscfg['global_data']['endtime']:
            flag_interp = interpol_field(
                dscfg['global_data']['radar_out'], radar_aux, 'time_avg_flag')
            dscfg['global_data']['radar_out'].fields[
                'time_avg_flag']['data'] += (
                    flag_interp['data'].filled(fill_value=0)).astype(int)

            return None, None

        # we have reached the end of the accumulation: start a new object
        new_dataset = {
            'radar_out': deepcopy(dscfg['global_data']['radar_out']),
            'timeinfo': dscfg['global_data']['endtime']}

        dscfg['global_data']['starttime'] += datetime.timedelta(
            seconds=period)
        dscfg['global_data']['endtime'] += datetime.timedelta(seconds=period)

        # remove old radar object from global_data dictionary
        dscfg['global_data'].pop('radar_out', None)

        # get start and stop times of new radar object
        dscfg['global_data']['starttime'], dscfg['global_data']['endtime'] = (
            time_avg_range(
                dscfg['timeinfo'], dscfg['global_data']['starttime'],
                dscfg['global_data']['endtime'], period))

        # check if volume time older than starttime
        if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
            dscfg['global_data'].update({'radar_out': radar_aux})

        return new_dataset, ind_rad

    # no more files to process if there is global data pack it up
    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None, None
        if 'radar_out' not in dscfg['global_data']:
            return None, None

        new_dataset = {
            'radar_out': deepcopy(dscfg['global_data']['radar_out']),
            'timeinfo': dscfg['global_data']['endtime']}

        return new_dataset, ind_rad


def process_colocated_gates(procstatus, dscfg, radar_list=None):
    """
    Find colocated gates within two radars

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        h_tol : float. Dataset keyword
            Tolerance in altitude difference between radar gates [m].
            Default 100.
        latlon_tol : float. Dataset keyword
            Tolerance in latitude and longitude position between radar gates
            [deg]. Default 0.0005
        vol_d_tol : float. Dataset keyword
            Tolerance in pulse volume diameter [m]. Default 100.
        vismin : float. Dataset keyword
            Minimum visibility [percent]. Default None.
        hmin : float. Dataset keyword
            Minimum altitude [m MSL]. Default None.
        hmax : float. Dataset keyword
            Maximum altitude [m MSL]. Default None.
        rmin : float. Dataset keyword
            Minimum range [m]. Default None.
        rmax : float. Dataset keyword
            Maximum range [m]. Default None.
        elmin : float. Dataset keyword
            Minimum elevation angle [deg]. Default None.
        elmax : float. Dataset keyword
            Maximum elevation angle [deg]. Default None.
        azrad1min : float. Dataset keyword
            Minimum azimuth angle [deg] for radar 1. Default None.
        azrad1max : float. Dataset keyword
            Maximum azimuth angle [deg] for radar 1. Default None.
        azrad2min : float. Dataset keyword
            Minimum azimuth angle [deg] for radar 2. Default None.
        azrad2max : float. Dataset keyword
            Maximum azimuth angle [deg] for radar 2. Default None.

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : radar object
        radar object containing the flag field
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # check how many radars are there
    radarnr_dict = dict()
    ind_radar_list = set()
    for datatypedescr in dscfg['datatype']:
        radarnr = datatypedescr.split(':')[0]
        radarnr_dict.update({radarnr: []})
        ind_radar_list.add(int(radarnr[5:8])-1)

    ind_radar_list = list(ind_radar_list)

    if (len(radarnr_dict) != 2) or (len(radar_list) < 2):
        warn('Intercomparison requires data from two different radars')
        return None, None

    # create the list of data types for each radar
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if radarnr in radarnr_dict:
            radarnr_dict[radarnr].append(get_fieldname_pyart(datatype))

    radar1 = radar_list[ind_radar_list[0]]
    radar2 = radar_list[ind_radar_list[1]]

    if radar1 is None or radar2 is None:
        warn('Unable to inter-compare radars. Missing radar')

    if 'instrument_name' in radar1.metadata:
        print('Radar 1: '+radar1.metadata['instrument_name'])
    if 'instrument_name' in radar2.metadata:
        print('Radar 2: '+radar2.metadata['instrument_name'])

    coloc_gates_field = 'colocated_gates'

    h_tol = dscfg.get('h_tol', 100.)
    latlon_tol = dscfg.get('latlon_tol', 0.0005)
    vol_d_tol = dscfg.get('vol_d_tol', 100.)
    vismin = dscfg.get('vismin', None)
    hmin = dscfg.get('hmin', None)
    hmax = dscfg.get('hmax', None)
    rmin = dscfg.get('rmin', None)
    rmax = dscfg.get('rmax', None)
    elmin = dscfg.get('elmin', None)
    elmax = dscfg.get('elmax', None)
    azrad1min = dscfg.get('azrad1min', None)
    azrad1max = dscfg.get('azrad1max', None)
    azrad2min = dscfg.get('azrad2min', None)
    azrad2max = dscfg.get('azrad2max', None)

    visib_field = None
    if 'visibility' in radarnr_dict['RADAR'+'{:03d}'.format(
            ind_radar_list[0]+1)]:
        visib_field = 'visibility'
    if vismin is not None and visib_field is None:
        warn('Unable to filter data according to visibility. ' +
             'Visibility field for RADAR'+'{:03d}'.format(
                 ind_radar_list[0]+1)+' not available')

    gate_coloc_rad1_dict = pyart.util.intersection(
        radar1, radar2,
        h_tol=h_tol, latlon_tol=latlon_tol, vol_d_tol=vol_d_tol,
        vismin=vismin, hmin=hmin, hmax=hmax, rmin=rmin, rmax=rmax,
        elmin=elmin, elmax=elmax, azmin=azrad1min, azmax=azrad1max,
        visib_field=visib_field, intersec_field=coloc_gates_field)

    visib_field = None
    if 'visibility' in radarnr_dict['RADAR'+'{:03d}'.format(
            ind_radar_list[1]+1)]:
        visib_field = 'visibility'

    if vismin is not None and visib_field is None:
        warn('Unable to filter data according to visibility. ' +
             'Visibility field for RADAR'+'{:03d}'.format(
                 ind_radar_list[1]+1)+' not available')

    gate_coloc_rad2_dict = pyart.util.intersection(
        radar2, radar1,
        h_tol=h_tol, latlon_tol=latlon_tol, vol_d_tol=vol_d_tol,
        vismin=vismin, hmin=hmin, hmax=hmax, rmin=rmin, rmax=rmax,
        elmin=elmin, elmax=elmax, azmin=azrad2min, azmax=azrad2max,
        visib_field=visib_field, intersec_field=coloc_gates_field)

    new_rad1 = deepcopy(radar1)
    new_rad1.fields = dict()
    new_rad1.add_field('colocated_gates', gate_coloc_rad1_dict)

    new_rad2 = deepcopy(radar2)
    new_rad2.fields = dict()
    new_rad2.add_field('colocated_gates', gate_coloc_rad2_dict)

    coloc_rad1_dict, new_rad1.fields['colocated_gates'] = (
        pyart.util.colocated_gates(
            new_rad1, new_rad2, h_tol=h_tol,
            latlon_tol=latlon_tol, coloc_gates_field=coloc_gates_field))

    coloc_rad2_dict, new_rad2.fields['colocated_gates'] = (
        pyart.util.colocated_gates(
            new_rad2, new_rad1, h_tol=h_tol,
            latlon_tol=latlon_tol, coloc_gates_field=coloc_gates_field))

    # prepare output
    rad1_dict = {
        'coloc_dict': coloc_rad1_dict,
        'radar_out': new_rad1}

    rad2_dict = {
        'coloc_dict': coloc_rad2_dict,
        'radar_out': new_rad2}

    new_dataset = {
        'RADAR'+'{:03d}'.format(ind_radar_list[0]+1): rad1_dict,
        'RADAR'+'{:03d}'.format(ind_radar_list[1]+1): rad2_dict}

    return new_dataset, ind_radar_list


def process_intercomp(procstatus, dscfg, radar_list=None):
    """
    intercomparison between two radars at co-located gates. The variables
    compared must be of the same type.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types        
        colocgatespath : string.
            base path to the file containing the coordinates of the co-located
            gates
        coloc_radars_name : string. Dataset keyword
            string identifying the radar names
        rays_are_indexed : bool. Dataset keyword
            If True it is considered that the rays are indexed and that the
            data can be selected simply looking at the ray number.
            Default false
        azi_tol : float. Dataset keyword
            azimuth tolerance between the two radars. Default 0.5 deg
        ele_tol : float. Dataset keyword
            elevation tolerance between the two radars. Default 0.5 deg
        rng_tol : float. Dataset keyword
            range tolerance between the two radars. Default 50 m
        coloc_data_dir : string. Dataset keyword
            name of the directory containing the csv file with colocated data
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing a dictionary with intercomparison data and the
        key "final" which contains a boolean that is true when all volumes
        have been processed
    ind_rad : int
        radar index

    """
    # check how many radars have to be compared and which datatype to use
    ind_radar_list = []
    field_name_list = []
    datatype_list = []
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_name = get_fieldname_pyart(datatype)
        ind_radar_list.append(int(radarnr[5:8])-1)
        datatype_list.append(datatype)
        field_name_list.append(field_name)
    ind_radar_list = np.array(ind_radar_list)
    datatype_list = np.array(datatype_list)
    field_name_list = np.array(field_name_list)
    
    if (ind_radar_list.size != 2) or (np.unique(ind_radar_list).size < 2):
        warn('Intercomparison requires data from two different radars')
        return None, None
    if np.unique(field_name_list).size > 1:
        warn('Intercomparison must be performed on the same variable')
        return None, None
        
    if procstatus == 0:
        savedir = dscfg['colocgatespath']+dscfg['coloc_radars_name']+'/'

        prdtype = 'info'
        if 'prdtype' in dscfg:
            prdtype = dscfg['prdtype']

        fname = make_filename(
            prdtype, 'COLOCATED_GATES', dscfg['coloc_radars_name'], ['csv'],
            timeinfo=None)[0]

        (rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi, rad1_rng,
         rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi, rad2_rng) = (
             read_colocated_gates(savedir+fname))

        if rad1_ele is None:
            raise ValueError('Unable to intercompare radars. ' +
                             'Missing colocated gates file')

        dscfg['global_data'] = {
            'rad1_ray_ind': rad1_ray_ind,
            'rad1_rng_ind': rad1_rng_ind,
            'rad1_ele': rad1_ele,
            'rad1_azi': rad1_azi,
            'rad1_rng': rad1_rng,
            'rad2_ray_ind': rad2_ray_ind,
            'rad2_rng_ind': rad2_rng_ind,
            'rad2_ele': rad2_ele,
            'rad2_azi': rad2_azi,
            'rad2_rng': rad2_rng}

        return None, None

    if procstatus == 1:
        radar1 = radar_list[ind_radar_list[0]]
        radar2 = radar_list[ind_radar_list[1]]

        if radar1 is None or radar2 is None:
            warn('Unable to inter-compare radars. Missing radar')
            return None, None

        if ((field_name not in radar1.fields) or
                (field_name not in radar2.fields)):
            warn('Unable to get values of field '+field_name +
                 ' at colocated range bins. ' +
                 'Field missing in one of the radars')
            return None, None

        if not dscfg['initialized']:
            dscfg['global_data'].update({'timeinfo': dscfg['timeinfo']})
            dscfg['global_data'].update(
                {'rad1_name': dscfg['RadarName'][ind_radar_list[0]]})
            dscfg['global_data'].update(
                {'rad2_name': dscfg['RadarName'][ind_radar_list[1]]})
            dscfg['initialized'] = 1

        rad1_field = radar1.fields[field_name]['data']
        rad2_field = radar2.fields[field_name]['data']

        intercomp_dict = {
            'rad1_time': [],
            'rad1_ray_ind': [],
            'rad1_rng_ind': [],
            'rad1_ele': [],
            'rad1_azi': [],
            'rad1_rng': [],
            'rad1_val': [],
            'rad2_time': [],
            'rad2_ray_ind': [],
            'rad2_rng_ind': [],
            'rad2_ele': [],
            'rad2_azi': [],
            'rad2_rng': [],
            'rad2_val': []}

        # determine if radar data has to be averaged
        avg_rad1, avg_rad2, avg_rad_lim = get_range_bins_to_avg(
            radar1.range['data'], radar2.range['data'])

        # rays are indexed to regular grid
        rays_are_indexed = dscfg.get('rays_are_indexed', False)
        if not rays_are_indexed:
            azi_tol = dscfg.get('azi_tol', 0.5)
            ele_tol = dscfg.get('ele_tol', 0.5)
            rng_tol = dscfg.get('rng_tol', 50.)

            rad1_ray_ind, rad1_rng_ind, rad2_ray_ind, rad2_rng_ind = (
                find_colocated_indexes(
                    radar1, radar2, dscfg['global_data']['rad1_ele'],
                    dscfg['global_data']['rad1_azi'],
                    dscfg['global_data']['rad1_rng'],
                    dscfg['global_data']['rad2_ele'],
                    dscfg['global_data']['rad2_azi'],
                    dscfg['global_data']['rad2_rng'], ele_tol=ele_tol,
                    azi_tol=azi_tol, rng_tol=rng_tol))
        else:
            rad1_ray_ind = deepcopy(dscfg['global_data']['rad1_ray_ind'])
            rad1_rng_ind = deepcopy(dscfg['global_data']['rad1_rng_ind'])
            rad2_ray_ind = deepcopy(dscfg['global_data']['rad2_ray_ind'])
            rad2_rng_ind = deepcopy(dscfg['global_data']['rad2_rng_ind'])

        # keep only indices of valid gates
        val1_vec = rad1_field[rad1_ray_ind, rad1_rng_ind]
        val2_vec = rad2_field[rad1_ray_ind, rad1_rng_ind]

        mask_val1 = np.ma.getmaskarray(val1_vec)
        mask_val2 = np.ma.getmaskarray(val2_vec)

        isvalid = np.logical_not(np.logical_or(mask_val1, mask_val2))

        rad1_ray_ind = rad1_ray_ind[isvalid]
        rad1_rng_ind = rad1_rng_ind[isvalid]
        rad2_ray_ind = rad2_ray_ind[isvalid]
        rad2_rng_ind = rad2_rng_ind[isvalid]

        # if averaging required loop over valid gates and average
        if avg_rad1:
            ngates_valid = len(rad1_ray_ind)
            val1_vec = np.ma.masked_all(ngates_valid, dtype=float)
            is_valid_avg = np.zeros(ngates_valid, dtype=bool)
            for i in range(ngates_valid):
                if rad1_rng_ind[i]+avg_rad_lim[1] >= radar1.ngates:
                    continue
                if rad1_rng_ind[i]+avg_rad_lim[0] < 0:
                    continue
                ind_rng = list(range(rad1_rng_ind[i]+avg_rad_lim[0],
                                     rad1_rng_ind[i]+avg_rad_lim[1]+1))

                if np.any(np.ma.getmaskarray(
                        rad1_field[rad1_ray_ind[i], ind_rng])):
                    continue

                val1_vec[i] = np.ma.asarray(np.ma.mean(
                    rad1_field[rad1_ray_ind[i], ind_rng]))

                is_valid_avg[i] = True

            rad1_ray_ind = rad1_ray_ind[is_valid_avg]
            rad1_rng_ind = rad1_rng_ind[is_valid_avg]
            rad2_ray_ind = rad2_ray_ind[is_valid_avg]
            rad2_rng_ind = rad2_rng_ind[is_valid_avg]

            val1_vec = val1_vec[is_valid_avg]
            val2_vec = rad2_field[rad2_ray_ind, rad2_rng_ind]

        elif avg_rad2:
            ngates_valid = len(rad2_ray_ind)
            val2_vec = np.ma.masked_all(ngates_valid, dtype=float)
            is_valid_avg = np.zeros(ngates_valid, dtype=bool)
            for i in range(ngates_valid):
                if rad2_rng_ind[i]+avg_rad_lim[1] >= radar2.ngates:
                    continue
                if rad2_rng_ind[i]+avg_rad_lim[0] < 0:
                    continue
                ind_rng = list(range(rad2_rng_ind[i]+avg_rad_lim[0],
                                     rad2_rng_ind[i]+avg_rad_lim[1]+1))

                if np.any(np.ma.getmaskarray(
                        rad2_field[rad2_ray_ind[i], ind_rng])):
                    continue

                val2_vec[i] = np.ma.asarray(np.ma.mean(
                    rad2_field[rad2_ray_ind[i], ind_rng]))

                is_valid_avg[i] = True

            rad1_ray_ind = rad1_ray_ind[is_valid_avg]
            rad1_rng_ind = rad1_rng_ind[is_valid_avg]
            rad2_ray_ind = rad2_ray_ind[is_valid_avg]
            rad2_rng_ind = rad2_rng_ind[is_valid_avg]

            val2_vec = val2_vec[is_valid_avg]
            val1_vec = rad1_field[rad1_ray_ind, rad1_rng_ind]
        else:
            val1_vec = val1_vec[isvalid]
            val2_vec = val2_vec[isvalid]

        intercomp_dict['rad1_time'] = num2date(
            radar1.time['data'][rad1_ray_ind], radar1.time['units'],
            radar1.time['calendar'])
        intercomp_dict['rad1_ray_ind'] = rad1_ray_ind
        intercomp_dict['rad1_rng_ind'] = rad1_rng_ind
        intercomp_dict['rad1_ele'] = radar1.elevation['data'][rad1_ray_ind]
        intercomp_dict['rad1_azi'] = radar1.azimuth['data'][rad1_ray_ind]
        intercomp_dict['rad1_rng'] = radar1.range['data'][rad1_rng_ind]
        intercomp_dict['rad1_val'] = val1_vec

        intercomp_dict['rad2_time'] = num2date(
            radar2.time['data'][rad2_ray_ind], radar2.time['units'],
            radar2.time['calendar'])
        intercomp_dict['rad2_ray_ind'] = rad2_ray_ind
        intercomp_dict['rad2_rng_ind'] = rad2_rng_ind
        intercomp_dict['rad2_ele'] = radar2.elevation['data'][rad2_ray_ind]
        intercomp_dict['rad2_azi'] = radar2.azimuth['data'][rad2_ray_ind]
        intercomp_dict['rad2_rng'] = radar2.range['data'][rad2_rng_ind]
        intercomp_dict['rad2_val'] = val2_vec

        new_dataset = {'intercomp_dict': intercomp_dict,
                       'timeinfo': dscfg['global_data']['timeinfo'],
                       'final': False}
        return new_dataset, None

    if procstatus == 2:
        savedir = get_save_dir(
            dscfg['basepath'], dscfg['procname'], dscfg['dsname'],
            dscfg['coloc_data_dir'],
            timeinfo=dscfg['global_data']['timeinfo'], create_dir=False)

        fname = make_filename(
            'colocated_data', dscfg['type'], datatype, ['csv'],
            timeinfo=dscfg['global_data']['timeinfo'], timeformat='%Y%m%d')

        fname = savedir+fname[0]

        coloc_data = read_colocated_data(fname)

        intercomp_dict = {
            'rad1_name': dscfg['global_data']['rad1_name'],
            'rad1_time': coloc_data[0],
            'rad1_ray_ind': coloc_data[1],
            'rad1_rng_ind': coloc_data[2],
            'rad1_ele': coloc_data[3],
            'rad1_azi': coloc_data[4],
            'rad1_rng': coloc_data[5],
            'rad1_val': coloc_data[6],
            'rad2_name': dscfg['global_data']['rad2_name'],
            'rad2_time': coloc_data[7],
            'rad2_ray_ind': coloc_data[8],
            'rad2_rng_ind': coloc_data[9],
            'rad2_ele': coloc_data[10],
            'rad2_azi': coloc_data[11],
            'rad2_rng': coloc_data[12],
            'rad2_val': coloc_data[13]}

        new_dataset = {'intercomp_dict': intercomp_dict,
                       'timeinfo': dscfg['global_data']['timeinfo'],
                       'final': True}

        return new_dataset, None


def process_intercomp_time_avg(procstatus, dscfg, radar_list=None):
    """
    intercomparison between the average reflectivity of two radars

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        colocgatespath : string.
            base path to the file containing the coordinates of the co-located
            gates
        coloc_data_dir : string. Dataset keyword
            name of the directory containing the csv file with colocated data
        coloc_radars_name : string. Dataset keyword
            string identifying the radar names
        rays_are_indexed : bool. Dataset keyword
            If True it is considered that the rays are indexed and that the
            data can be selected simply looking at the ray number.
            Default false
        azi_tol : float. Dataset keyword
            azimuth tolerance between the two radars. Default 0.5 deg
        ele_tol : float. Dataset keyword
            elevation tolerance between the two radars. Default 0.5 deg
        rng_tol : float. Dataset keyword
            range tolerance between the two radars. Default 50 m
        clt_max : int. Dataset keyword
            maximum number of samples that can be clutter contaminated.
            Default 100 i.e. all
        phi_excess_max : int. Dataset keyword
            maximum number of samples that can have excess instantaneous
            PhiDP. Default 100 i.e. all
        non_rain_max : int. Dataset keyword
            maximum number of samples that can be no rain. Default 100 i.e. all
        phi_avg_max : float. Dataset keyword
            maximum average PhiDP allowed. Default 600 deg i.e. any

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing a dictionary with intercomparison data and the
        key "final" which contains a boolean that is true when all volumes
        have been processed
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        savedir = dscfg['colocgatespath']+dscfg['coloc_radars_name']+'/'

        prdtype = 'info'
        if 'prdtype' in dscfg:
            prdtype = dscfg['prdtype']

        fname = make_filename(
            prdtype, 'COLOCATED_GATES', dscfg['coloc_radars_name'], ['csv'],
            timeinfo=None)[0]

        (rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi, rad1_rng,
         rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi, rad2_rng) = (
             read_colocated_gates(savedir+fname))

        if rad1_ele is None:
            raise ValueError('Unable to intercompare radars. ' +
                             'Missing colocated gates file')

        dscfg['global_data'] = {
            'rad1_ray_ind': rad1_ray_ind,
            'rad1_rng_ind': rad1_rng_ind,
            'rad1_ele': rad1_ele,
            'rad1_azi': rad1_azi,
            'rad1_rng': rad1_rng,
            'rad2_ray_ind': rad2_ray_ind,
            'rad2_rng_ind': rad2_rng_ind,
            'rad2_ele': rad2_ele,
            'rad2_azi': rad2_azi,
            'rad2_rng': rad2_rng}

        return None, None

    if procstatus == 1:
        # check how many radars are there
        ind_radar_list = set()
        for datatypedescr in dscfg['datatype']:
            radarnr = datatypedescr.split(':')[0]
            ind_radar_list.add(int(radarnr[5:8])-1)

        ind_radar_list = list(ind_radar_list)

        if (len(ind_radar_list) != 2) or (len(radar_list) < 2):
            warn('Intercomparison requires data from two different radars')
            return None, None

        radarnr_list = ['RADAR'+'{:03d}'.format(ind_radar_list[0]+1),
                        'RADAR'+'{:03d}'.format(ind_radar_list[1]+1)]

        # get field names
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            if radarnr == radarnr_list[0]:
                if (datatype in (
                        'dBZ', 'dBZc', 'dBuZ', 'dBZv', 'dBZvc', 'dBuZv')):
                    rad1_refl_field = get_fieldname_pyart(datatype)
                elif datatype in ('PhiDP', 'PhiDPc'):
                    rad1_phidp_field = get_fieldname_pyart(datatype)
                elif datatype == 'time_avg_flag':
                    rad1_flag_field = get_fieldname_pyart(datatype)
            elif radarnr == radarnr_list[1]:
                if (datatype in (
                        'dBZ', 'dBZc', 'dBuZ', 'dBZv', 'dBZvc', 'dBuZv')):
                    rad2_refl_field = get_fieldname_pyart(datatype)
                elif datatype in ('PhiDP', 'PhiDPc'):
                    rad2_phidp_field = get_fieldname_pyart(datatype)
                elif datatype == 'time_avg_flag':
                    rad2_flag_field = get_fieldname_pyart(datatype)

        radar1 = radar_list[ind_radar_list[0]]
        radar2 = radar_list[ind_radar_list[1]]

        if radar1 is None or radar2 is None:
            warn('Unable to inter-compare radars. Missing radar')
            return None, None

        if ((rad1_refl_field not in radar1.fields) or
                (rad1_phidp_field not in radar1.fields) or
                (rad1_flag_field not in radar1.fields) or
                (rad2_refl_field not in radar2.fields) or
                (rad2_phidp_field not in radar2.fields) or
                (rad2_flag_field not in radar2.fields)):
            warn('Unable to compare radar time avg fields. ' +
                 'Fields missing')
            return None, None

        if not dscfg['initialized']:
            dscfg['global_data'].update({'timeinfo': dscfg['timeinfo']})
            dscfg['global_data'].update(
                {'rad1_name': dscfg['RadarName'][ind_radar_list[0]]})
            dscfg['global_data'].update(
                {'rad2_name': dscfg['RadarName'][ind_radar_list[1]]})
            dscfg['initialized'] = 1

        refl1 = radar1.fields[rad1_refl_field]['data']
        refl2 = radar2.fields[rad2_refl_field]['data']

        phidp1 = radar1.fields[rad1_phidp_field]['data']
        phidp2 = radar2.fields[rad2_phidp_field]['data']

        flag1 = radar1.fields[rad1_flag_field]['data']
        flag2 = radar2.fields[rad2_flag_field]['data']

        intercomp_dict = {
            'rad1_time': [],
            'rad1_ray_ind': [],
            'rad1_rng_ind': [],
            'rad1_ele': [],
            'rad1_azi': [],
            'rad1_rng': [],
            'rad1_dBZavg': [],
            'rad1_PhiDPavg': [],
            'rad1_Flagavg': [],
            'rad2_time': [],
            'rad2_ray_ind': [],
            'rad2_rng_ind': [],
            'rad2_ele': [],
            'rad2_azi': [],
            'rad2_rng': [],
            'rad2_dBZavg': [],
            'rad2_PhiDPavg': [],
            'rad2_Flagavg': []}

        # determine if radar data has to be averaged
        avg_rad1, avg_rad2, avg_rad_lim = get_range_bins_to_avg(
            radar1.range['data'], radar2.range['data'])

        # rays are indexed to regular grid
        rays_are_indexed = dscfg.get('rays_are_indexed', False)
        # get current radars gates indices
        if not rays_are_indexed:
            azi_tol = dscfg.get('azi_tol', 0.5)
            ele_tol = dscfg.get('ele_tol', 0.5)
            rng_tol = dscfg.get('rng_tol', 50.)

            rad1_ray_ind, rad1_rng_ind, rad2_ray_ind, rad2_rng_ind = (
                find_colocated_indexes(
                    radar1, radar2, dscfg['global_data']['rad1_ele'],
                    dscfg['global_data']['rad1_azi'],
                    dscfg['global_data']['rad1_rng'],
                    dscfg['global_data']['rad2_ele'],
                    dscfg['global_data']['rad2_azi'],
                    dscfg['global_data']['rad2_rng'], ele_tol=ele_tol,
                    azi_tol=azi_tol, rng_tol=rng_tol))
        else:
            rad1_ray_ind = deepcopy(dscfg['global_data']['rad1_ray_ind'])
            rad1_rng_ind = deepcopy(dscfg['global_data']['rad1_rng_ind'])
            rad2_ray_ind = deepcopy(dscfg['global_data']['rad2_ray_ind'])
            rad2_rng_ind = deepcopy(dscfg['global_data']['rad2_rng_ind'])

        # keep only indices and data of valid gates
        refl1_vec = refl1[rad1_ray_ind, rad1_rng_ind]
        phidp1_vec = phidp1[rad1_ray_ind, rad1_rng_ind]
        flag1_vec = flag1[rad1_ray_ind, rad1_rng_ind]

        refl2_vec = refl2[rad2_ray_ind, rad2_rng_ind]
        phidp2_vec = phidp2[rad2_ray_ind, rad2_rng_ind]
        flag2_vec = flag2[rad2_ray_ind, rad2_rng_ind]

        mask_refl1 = np.ma.getmaskarray(refl1_vec)
        mask_phidp1 = np.ma.getmaskarray(phidp1_vec)
        mask_refl2 = np.ma.getmaskarray(refl2_vec)
        mask_phidp2 = np.ma.getmaskarray(phidp2_vec)

        isvalid = np.logical_not(
            np.logical_or(np.logical_or(mask_refl1, mask_refl2),
                          np.logical_or(mask_phidp1, mask_phidp2)))

        rad1_ray_ind = rad1_ray_ind[isvalid]
        rad1_rng_ind = rad1_rng_ind[isvalid]
        rad2_ray_ind = rad2_ray_ind[isvalid]
        rad2_rng_ind = rad2_rng_ind[isvalid]

        # if averaging required loop over valid gates and average
        # only if all gates valid
        if avg_rad1:
            ngates_valid = len(rad1_ray_ind)
            refl1_vec = np.ma.masked_all(ngates_valid, dtype=float)
            phidp1_vec = np.ma.masked_all(ngates_valid, dtype=float)
            flag1_vec = np.ma.masked_all(ngates_valid, dtype=int)
            is_valid_avg = np.zeros(ngates_valid, dtype=bool)
            for i in range(ngates_valid):
                if rad1_rng_ind[i]+avg_rad_lim[1] >= radar1.ngates:
                    continue
                if rad1_rng_ind[i]+avg_rad_lim[0] < 0:
                    continue
                ind_rng = list(range(rad1_rng_ind[i]+avg_rad_lim[0],
                                     rad1_rng_ind[i]+avg_rad_lim[1]+1))

                if np.any(np.ma.getmaskarray(
                        refl1[rad1_ray_ind[i], ind_rng])):
                    continue
                if np.any(np.ma.getmaskarray(
                        phidp1[rad1_ray_ind[i], ind_rng])):
                    continue

                refl1_vec[i] = np.ma.asarray(np.ma.mean(
                    refl1[rad1_ray_ind[i], ind_rng]))
                phidp1_vec[i] = np.ma.asarray(np.ma.mean(
                    phidp1[rad1_ray_ind[i], ind_rng]))

                rad1_flag = flag1[rad1_ray_ind[i], ind_rng]

                rad1_excess_phi = rad1_flag % 100
                rad1_clt = ((rad1_flag-rad1_excess_phi) % 10000) / 100
                rad1_prec = (
                    ((rad1_flag-rad1_clt*100-rad1_excess_phi) % 1000000) /
                    10000)

                flag1_vec[i] = int(
                    10000*np.max(rad1_prec)+100*np.max(rad1_clt) +
                    np.max(rad1_excess_phi))
                is_valid_avg[i] = True

            rad1_ray_ind = rad1_ray_ind[is_valid_avg]
            rad1_rng_ind = rad1_rng_ind[is_valid_avg]
            rad2_ray_ind = rad2_ray_ind[is_valid_avg]
            rad2_rng_ind = rad2_rng_ind[is_valid_avg]

            refl1_vec = refl1_vec[is_valid_avg]
            phidp1_vec = phidp1_vec[is_valid_avg]
            flag1_vec = flag1_vec[is_valid_avg]

            refl2_vec = refl2[rad2_ray_ind, rad2_rng_ind]
            phidp2_vec = phidp2[rad2_ray_ind, rad2_rng_ind]
            flag2_vec = flag2[rad2_ray_ind, rad2_rng_ind]

        elif avg_rad2:
            ngates_valid = len(rad2_ray_ind)
            refl2_vec = np.ma.masked_all(ngates_valid, dtype=float)
            phidp2_vec = np.ma.masked_all(ngates_valid, dtype=float)
            flag2_vec = np.ma.masked_all(ngates_valid, dtype=int)
            is_valid_avg = np.zeros(ngates_valid, dtype=bool)
            for i in range(ngates_valid):
                if rad2_rng_ind[i]+avg_rad_lim[1] >= radar2.ngates:
                    continue
                if rad2_rng_ind[i]+avg_rad_lim[0] < 0:
                    continue
                ind_rng = list(range(rad2_rng_ind[i]+avg_rad_lim[0],
                                     rad2_rng_ind[i]+avg_rad_lim[1]+1))

                if np.any(np.ma.getmaskarray(
                        refl2[rad2_ray_ind[i], ind_rng])):
                    continue
                if np.any(np.ma.getmaskarray(
                        phidp2[rad2_ray_ind[i], ind_rng])):
                    continue

                refl2_vec[i] = np.ma.asarray(np.ma.mean(
                    refl2[rad2_ray_ind[i], ind_rng]))
                phidp2_vec[i] = np.ma.asarray(np.ma.mean(
                    phidp2[rad2_ray_ind[i], ind_rng]))

                rad2_flag = flag2[rad2_ray_ind[i], ind_rng]

                rad2_excess_phi = rad2_flag % 100
                rad2_clt = ((rad2_flag-rad2_excess_phi) % 10000) / 100
                rad2_prec = (
                    ((rad2_flag-rad2_clt*100-rad2_excess_phi) % 1000000) /
                    10000)

                flag2_vec[i] = int(
                    10000*np.max(rad2_prec)+100*np.max(rad2_clt) +
                    np.max(rad2_excess_phi))
                is_valid_avg[i] = True

            rad1_ray_ind = rad1_ray_ind[is_valid_avg]
            rad1_rng_ind = rad1_rng_ind[is_valid_avg]
            rad2_ray_ind = rad2_ray_ind[is_valid_avg]
            rad2_rng_ind = rad2_rng_ind[is_valid_avg]

            refl2_vec = refl2_vec[is_valid_avg]
            phidp2_vec = phidp2_vec[is_valid_avg]
            flag2_vec = flag2_vec[is_valid_avg]

            refl1_vec = refl1[rad1_ray_ind, rad1_rng_ind]
            phidp1_vec = phidp1[rad1_ray_ind, rad1_rng_ind]
            flag1_vec = flag1[rad1_ray_ind, rad1_rng_ind]
        else:
            refl1_vec = refl1_vec[isvalid]
            phidp1_vec = phidp1_vec[isvalid]
            flag1_vec = flag1_vec[isvalid]
            refl2_vec = refl2_vec[isvalid]
            phidp2_vec = phidp2_vec[isvalid]
            flag2_vec = flag2_vec[isvalid]

        intercomp_dict['rad1_time'] = np.empty(
            len(rad1_ray_ind), dtype=datetime.datetime)
        intercomp_dict['rad1_time'][:] = dscfg['global_data']['timeinfo']
        intercomp_dict['rad1_ray_ind'] = rad1_ray_ind
        intercomp_dict['rad1_rng_ind'] = rad1_rng_ind
        intercomp_dict['rad1_ele'] = radar1.elevation['data'][rad1_ray_ind]
        intercomp_dict['rad1_azi'] = radar1.azimuth['data'][rad1_ray_ind]
        intercomp_dict['rad1_rng'] = radar1.range['data'][rad1_rng_ind]
        intercomp_dict['rad1_dBZavg'] = refl1_vec
        intercomp_dict['rad1_PhiDPavg'] = phidp1_vec
        intercomp_dict['rad1_Flagavg'] = flag1_vec

        intercomp_dict['rad2_time'] = deepcopy(intercomp_dict['rad1_time'])
        intercomp_dict['rad2_ray_ind'] = rad2_ray_ind
        intercomp_dict['rad2_rng_ind'] = rad2_rng_ind
        intercomp_dict['rad2_ele'] = radar2.elevation['data'][rad2_ray_ind]
        intercomp_dict['rad2_azi'] = radar2.azimuth['data'][rad2_ray_ind]
        intercomp_dict['rad2_rng'] = radar2.range['data'][rad2_rng_ind]
        intercomp_dict['rad2_dBZavg'] = refl2_vec
        intercomp_dict['rad2_PhiDPavg'] = phidp2_vec
        intercomp_dict['rad2_Flagavg'] = flag2_vec

        new_dataset = {'intercomp_dict': intercomp_dict,
                       'timeinfo': dscfg['global_data']['timeinfo'],
                       'final': False}
        return new_dataset, None

    if procstatus == 2:
        # get field name
        refl_type = None
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            if datatype in ('dBZ', 'dBZc', 'dBuZ', 'dBZv', 'dBZvc', 'dBuZv'):
                refl_type = datatype
                break

        if refl_type is None:
            warn('Unknown reflectivity type')
            return None, None

        savedir = get_save_dir(
            dscfg['basepath'], dscfg['procname'], dscfg['dsname'],
            dscfg['coloc_data_dir'],
            timeinfo=dscfg['global_data']['timeinfo'], create_dir=False)

        fname = make_filename(
            'colocated_data', dscfg['type'], refl_type, ['csv'],
            timeinfo=dscfg['global_data']['timeinfo'], timeformat='%Y%m%d')

        fname = savedir+fname[0]

        (rad1_time, rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi, rad1_rng,
         rad1_dBZ, rad1_phi, rad1_flag, rad2_time, rad2_ray_ind, rad2_rng_ind,
         rad2_ele, rad2_azi, rad2_rng, rad2_dBZ, rad2_phi, rad2_flag) = (
             read_colocated_data_time_avg(fname))

        rad1_excess_phi = (rad1_flag % 100).astype(int)
        rad2_excess_phi = (rad2_flag % 100).astype(int)

        rad1_clt = (((rad1_flag-rad1_excess_phi) % 10000) / 100).astype(int)
        rad2_clt = (((rad2_flag-rad2_excess_phi) % 10000) / 100).astype(int)

        rad1_non_rain = (
            ((rad1_flag-rad1_clt*100-rad1_excess_phi) % 1000000) /
            10000).astype(int)
        rad2_non_rain = (
            ((rad2_flag-rad2_clt*100-rad2_excess_phi) % 1000000) /
            10000).astype(int)

        clt_max = dscfg.get('clt_max', 100)
        phi_excess_max = dscfg.get('phi_excess_max', 100)
        non_rain_max = dscfg.get('non_rain_max', 100)
        phi_avg_max = dscfg.get('phi_avg_max', 600.)

        # filter out invalid data
        ind_val = np.where(
            np.logical_and.reduce((
                rad1_clt <= clt_max, rad2_clt <= clt_max,
                rad1_excess_phi <= phi_excess_max,
                rad2_excess_phi <= phi_excess_max,
                rad1_non_rain <= non_rain_max, rad2_non_rain <= non_rain_max,
                rad1_phi <= phi_avg_max, rad2_phi <= phi_avg_max)))[0]

        intercomp_dict = {
            'rad1_name': dscfg['global_data']['rad1_name'],
            'rad1_time': rad1_time[ind_val],
            'rad1_ray_ind': rad1_ray_ind[ind_val],
            'rad1_rng_ind': rad1_rng_ind[ind_val],
            'rad1_ele': rad1_ele[ind_val],
            'rad1_azi': rad1_azi[ind_val],
            'rad1_rng': rad1_rng[ind_val],
            'rad1_val': rad1_dBZ[ind_val],
            'rad2_name': dscfg['global_data']['rad2_name'],
            'rad2_time': rad2_time[ind_val],
            'rad2_ray_ind': rad1_ray_ind[ind_val],
            'rad2_rng_ind': rad1_rng_ind[ind_val],
            'rad2_ele': rad2_ele[ind_val],
            'rad2_azi': rad2_azi[ind_val],
            'rad2_rng': rad2_rng[ind_val],
            'rad2_val': rad2_dBZ[ind_val]}

        new_dataset = {'intercomp_dict': intercomp_dict,
                       'timeinfo': dscfg['global_data']['timeinfo'],
                       'final': True}

        return new_dataset, None


def process_fields_diff(procstatus, dscfg, radar_list=None):
    """
    Computes the field difference between RADAR001 and radar002,
    i.e. RADAR001-RADAR002. Assumes both radars have the same geometry

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing a radar object containing the field differences
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # check how many radars are there
    radarnr_dict = dict()
    ind_radar_list = set()
    for datatypedescr in dscfg['datatype']:
        radarnr = datatypedescr.split(':')[0]
        radarnr_dict.update({radarnr: []})
        ind_radar_list.add(int(radarnr[5:8])-1)

    ind_radar_list = list(ind_radar_list)

    if (len(radarnr_dict) != 2) or (len(radar_list) < 2):
        warn('Intercomparison requires data from two different radars')
        return None, None

    # create the list of data types for each radar
    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
    field_name_1 = get_fieldname_pyart(datatype)

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][1])
    field_name_2 = get_fieldname_pyart(datatype)

    radar1 = radar_list[ind_radar_list[0]]
    radar2 = radar_list[ind_radar_list[1]]

    if radar1 is None or radar2 is None:
        warn('Unable to inter-compare radars. Missing radar')
        return None, None

    if ((field_name_1 not in radar1.fields) or
            (field_name_2 not in radar2.fields)):
        warn('Unable to compare fields '+field_name_1+'and '+field_name_2 +
             '. Field missing in one of the radars')
        return None, None

    field_diff = pyart.config.get_metadata('fields_difference')
    field_diff['data'] = (
        radar1.fields[field_name_1]['data'] -
        radar2.fields[field_name_2]['data'])
    field_diff['long_name'] = field_name_1+' - '+field_name_2

    rad_diff = deepcopy(radar1)
    rad_diff.fields = dict()
    rad_diff.add_field('fields_difference', field_diff)

    new_dataset = {'radar_out': rad_diff}

    return new_dataset, None


def process_intercomp_fields(procstatus, dscfg, radar_list=None):
    """
    intercomparison between two radars

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing a dictionary with intercomparison data
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # check how many radars are there
    radarnr_dict = dict()
    ind_radar_list = set()
    for datatypedescr in dscfg['datatype']:
        radarnr = datatypedescr.split(':')[0]
        radarnr_dict.update({radarnr: []})
        ind_radar_list.add(int(radarnr[5:8])-1)

    ind_radar_list = list(ind_radar_list)

    if (len(radarnr_dict) != 2) or (len(radar_list) < 2):
        warn('Intercomparison requires data from two different radars')
        return None, None

    # create the list of data types for each radar
    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
    field_name_1 = get_fieldname_pyart(datatype)

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][1])
    field_name_2 = get_fieldname_pyart(datatype)

    radar1 = radar_list[ind_radar_list[0]]
    radar2 = radar_list[ind_radar_list[1]]

    if radar1 is None or radar2 is None:
        warn('Unable to inter-compare radars. Missing radar')
        return None, None

    if ((field_name_1 not in radar1.fields) or
            (field_name_2 not in radar2.fields)):
        warn('Unable to compare fields '+field_name_1+' and '+field_name_2 +
             '. Field missing in one of the radars')
        return None, None

    data1 = deepcopy(radar1.fields[field_name_1]['data'])
    data2 = deepcopy(radar2.fields[field_name_2]['data'])

    mask1 = np.ma.getmaskarray(data1)
    mask2 = np.ma.getmaskarray(data2)

    data1[mask2] = np.ma.masked
    data2[mask1] = np.ma.masked

    intercomp_dict = {
        'rad1_name': dscfg['RadarName'][ind_radar_list[0]],
        'rad1_val': data1.compressed(),
        'rad2_name': dscfg['RadarName'][ind_radar_list[1]],
        'rad2_val': data2.compressed()}

    new_dataset = {'intercomp_dict': intercomp_dict,
                   'timeinfo': dscfg['timeinfo'],
                   'final': False}

    return new_dataset, None
