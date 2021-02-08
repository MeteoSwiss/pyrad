"""
pyrad.proc.process_calib
===========================

Functions for monitoring data quality and correct bias and noise effects

.. autosummary::
    :toctree: generated/

    process_correct_bias
    process_correct_noise_rhohv
    process_gc_monitoring
    process_occurrence
    process_time_avg_std
    process_occurrence_period
    process_sun_hits
   *process_sunscan

"""

from copy import deepcopy
from warnings import warn
import numpy as np
import sys

from scipy.constants import c as c_speed
from netCDF4 import num2date
from datetime import datetime as dt

try:
    import pysolar
    _PYSOLAR_AVAILABLE = True
except ImportError:
    _PYSOLAR_AVAILABLE = False

import pyart
from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..io.read_data_sun import read_sun_hits_multiple_days, read_solar_flux
from ..io.read_data_other import read_excess_gates
from ..io.read_data_radar import interpol_field

from ..util.radar_utils import get_closest_solar_flux, get_histogram_bins
from ..util.radar_utils import find_ray_index, find_rng_index


def process_correct_bias(procstatus, dscfg, radar_list=None):
    """
    Corrects a bias on the data

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type to correct for bias
        bias : float. Dataset keyword
            The bias to be corrected [dB]. Default 0
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        break
    field_name = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if field_name not in radar.fields:
        warn('Unable to correct for bias field ' + field_name +
             '. Field not available')
        return None, None

    bias = dscfg.get('bias', 0.)
    corrected_field = pyart.correct.correct_bias(
        radar, bias=bias, field_name=field_name)

    if field_name.startswith('corrected_'):
        new_field_name = field_name
    else:
        new_field_name = 'corrected_'+field_name

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(new_field_name, corrected_field)

    return new_dataset, ind_rad


def process_correct_noise_rhohv(procstatus, dscfg, radar_list=None):
    """
    identifies echoes as 0: No data, 1: Noise, 2: Clutter,
    3: Precipitation

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The data types used in the correction
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'uRhoHV':
            urhohv = 'uncorrected_cross_correlation_ratio'
        if datatype == 'SNRh':
            snr = 'signal_to_noise_ratio_hh'
        if datatype == 'ZDR':
            zdr = 'differential_reflectivity'
        if datatype == 'ZDRc':
            zdr = 'corrected_differential_reflectivity'
        if datatype == 'Nh':
            nh = 'noisedBZ_hh'
        if datatype == 'Nv':
            nv = 'noisedBZ_vv'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if ((urhohv not in radar.fields) or
            (snr not in radar.fields) or
            (zdr not in radar.fields) or
            (nh not in radar.fields) or
            (nv not in radar.fields)):
        warn('Unable to correct RhoHV field for noise. Missing fields')
        return None, None

    rhohv = pyart.correct.correct_noise_rhohv(
        radar, urhohv_field=urhohv, snr_field=snr, zdr_field=zdr,
        nh_field=nh, nv_field=nv, rhohv_field='cross_correlation_ratio')

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field('cross_correlation_ratio', rhohv)

    return new_dataset, ind_rad


def process_gc_monitoring(procstatus, dscfg, radar_list=None):
    """
    computes ground clutter monitoring statistics

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        excessgatespath : str. Config keyword
            The path to the gates in excess of quantile location
        excessgates_fname : str. Dataset keyword
            The name of the gates in excess of quantile file
        datatype : list of string. Dataset keyword
            The input data types
        step : float. Dataset keyword
            The width of the histogram bin. Default is None. In that case the
            default step in function get_histogram_bins is used
        regular_grid : Boolean. Dataset keyword
            Whether the radar has a Boolean grid or not. Default False
        val_min : Float. Dataset keyword
            Minimum value to consider that the gate has signal. Default None
        filter_prec : str. Dataset keyword
            Give which type of volume should be filtered. None, no filtering;
            keep_wet, keep wet volumes; keep_dry, keep dry volumes.
        rmax_prec : float. Dataset keyword
            Maximum range to consider when looking for wet gates [m]
        percent_prec_max : float. Dataset keyword
            Maxim percentage of wet gates to consider the volume dry
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : Radar
        radar object containing histogram data
    ind_rad : int
        radar index

    """
    echoid_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'echoID':
            echoid_field = get_fieldname_pyart(datatype)
        else:
            field_name = get_fieldname_pyart(datatype)
    ind_rad = int(radarnr[5:8])-1

    if procstatus == 0:
        savedir = dscfg['excessgatespath']
        fname = dscfg['excessgates_fname']
        ray_ind, rng_ind, ele, azi, rng, nsamples, occurrence, freq_occu = (
            read_excess_gates(savedir+fname))

        dscfg['global_data'] = {
            'ray_ind': ray_ind,
            'rng_ind': rng_ind,
            'ele': ele,
            'azi': azi,
            'rng': rng,
            'nsamples': nsamples,
            'occurrence': occurrence,
            'freq_occu': freq_occu}

        return None, None

    if dscfg['global_data']['ray_ind'] is None:
        warn('Unable to get statistics of clutter')
        return None, None

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = deepcopy(radar_list[ind_rad])

        if field_name not in radar.fields:
            warn(field_name+' not available.')
            return None, None

        # filter out low values
        val_min = dscfg.get('val_min', None)
        mask = np.ma.getmaskarray(radar.fields[field_name]['data'])
        if val_min is not None:
            mask = np.logical_or(
                mask, radar.fields[field_name]['data'] < val_min)

        field = deepcopy(radar.fields[field_name]['data'])
        field[mask] = np.ma.masked

        # filter wet or dry volumes
        filter_prec = dscfg.get('filter_prec', 'None')
        if filter_prec in ('keep_wet', 'keep_dry'):
            if echoid_field not in radar.fields:
                warn('Unable to determine if there is precipitation ' +
                     'close to the radar. Missing echoID field.')
                return None, None

            # Put invalid values to noise
            echoid = deepcopy(radar.fields[echoid_field]['data'])
            echoid[mask] = 1

            rmax_prec = dscfg.get('rmax_prec', 0.)
            percent_prec_max = dscfg.get('percent_prec_max', 10.)

            ngates = radar.ngates
            if rmax_prec > 0.:
                ngates = len(
                    radar.range['data'][radar.range['data'] < rmax_prec])
            ngates_total = ngates*radar.nrays

            prec_field = echoid[:, :ngates]
            ngates_prec = np.size(prec_field[prec_field == 3])

            percent_prec = ngates_prec/ngates_total*100.
            warn('Percent gates with precipitation: '+str(percent_prec)+'\n')
            if percent_prec > percent_prec_max:
                if filter_prec == 'keep_dry':
                    warn('Radar volume is precipitation contaminated.\n' +
                         'Maximum percentage allowed: '+str(percent_prec_max))
                    return None, None
            else:
                if filter_prec == 'keep_wet':
                    warn('Radar volume has not enough precipitation.\n' +
                         'Minimum percentage required: ' +
                         str(percent_prec_max))
                    return None, None

        step = dscfg.get('step', None)
        bin_edges = get_histogram_bins(field_name, step=step)
        nbins = len(bin_edges)-1
        step = bin_edges[1]-bin_edges[0]
        bin_centers = bin_edges[:-1]+step/2.

        # create histogram object from radar object
        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()
        radar_aux.range['data'] = bin_centers
        radar_aux.ngates = nbins
        radar_aux.nrays = 1

        field_dict = pyart.config.get_metadata(field_name)
        field_dict['data'] = np.ma.zeros((1, nbins), dtype=int)

        # rays are indexed to regular grid
        regular_grid = dscfg.get('regular_grid', False)
        if regular_grid:
            ray_ind = dscfg['global_data']['ray_ind']
            rng_ind = dscfg['global_data']['rng_ind']
            field = field[ray_ind, rng_ind].compressed()
        else:
            azi_tol = dscfg.get('azi_tol', 0.5)
            ele_tol = dscfg.get('ele_tol', 0.5)
            rng_tol = dscfg.get('rng_tol', 50.)

            # get indexes of gates close to target
            ngc = np.size(dscfg['global_data']['ray_ind'])

            ray_ind = np.ma.masked_all(ngc, dtype=int)
            rng_ind = np.ma.masked_all(ngc, dtype=int)
            for i in range(ngc):
                ind_ray_rad = find_ray_index(
                    radar.elevation['data'], radar.azimuth['data'],
                    dscfg['global_data']['ele'][i],
                    dscfg['global_data']['azi'][i],
                    ele_tol=ele_tol, azi_tol=azi_tol)
                if ind_ray_rad is None:
                    continue
                ind_rng_rad = find_rng_index(
                    radar.range['data'], dscfg['global_data']['rng'][i],
                    rng_tol=rng_tol)
                if ind_rng_rad is None:
                    continue
                ray_ind[i] = ind_ray_rad
                rng_ind[i] = ind_rng_rad
            ray_ind = ray_ind.compressed()
            rng_ind = rng_ind.compressed()
            field = field[ray_ind, rng_ind].compressed()

        # put gates with values off limits to limit
        # and compute histogram
        field[field < bin_centers[0]] = bin_centers[0]
        field[field > bin_centers[-1]] = bin_centers[-1]

        field_dict['data'][0, :], bin_edges = np.histogram(
            field, bins=bin_edges)
        radar_aux.add_field(field_name, field_dict)
        start_time = pyart.graph.common.generate_radar_time_begin(radar_aux)

        # Put histogram in Memory or add to existing histogram
        if dscfg['initialized'] == 0:
            dscfg['global_data'].update({
                'hist_obj': radar_aux,
                'timeinfo': start_time})
            dscfg['initialized'] = 1
        else:
            dscfg['global_data']['hist_obj'].fields[field_name]['data'] += (
                field_dict['data'].filled(fill_value=0)).astype('int64')

        #    dscfg['global_data']['timeinfo'] = dscfg['timeinfo']

        dataset = dict()
        dataset.update({'hist_obj': radar_aux})
        dataset.update({'hist_type': 'instant'})
        dataset.update({'timeinfo': start_time})

        return dataset, ind_rad

    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None, None

        dataset = dict()
        dataset.update({'hist_obj': dscfg['global_data']['hist_obj']})
        dataset.update({'hist_type': 'cumulative'})
        dataset.update({'timeinfo': dscfg['global_data']['timeinfo']})

        return dataset, ind_rad


def process_occurrence(procstatus, dscfg, radar_list=None):
    """
    computes the frequency of occurrence of data. It looks only for gates
    where data is present.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        regular_grid : Boolean. Dataset keyword
            Whether the radar has a Boolean grid or not. Default False
        rmin, rmax : float. Dataset keyword
            minimum and maximum ranges where the computation takes place. If
            -1 the whole range is considered. Default is -1
        val_min : Float. Dataset keyword
            Minimum value to consider that the gate has signal. Default None
        filter_prec : str. Dataset keyword
            Give which type of volume should be filtered. None, no filtering;
            keep_wet, keep wet volumes; keep_dry, keep dry volumes.
        rmax_prec : float. Dataset keyword
            Maximum range to consider when looking for wet gates [m]
        percent_prec_max : float. Dataset keyword
            Maxim percentage of wet gates to consider the volume dry
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """
    echoid_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'echoID':
            echoid_field = get_fieldname_pyart(datatype)
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

        if field_name not in radar.fields:
            warn(field_name+' not available.')
            return None, None

        # filter out low values
        val_min = dscfg.get('val_min', None)
        mask = np.ma.getmaskarray(radar.fields[field_name]['data'])
        if val_min is not None:
            mask = np.logical_or(
                mask, radar.fields[field_name]['data'] < val_min)

        filter_prec = dscfg.get('filter_prec', 'None')
        if filter_prec in ('keep_wet', 'keep_dry'):
            if echoid_field not in radar.fields:
                warn('Unable to determine if there is precipitation ' +
                     'close to the radar. Missing echoID field.')
                return None, None

            # Put invalid values to noise
            echoid = deepcopy(radar.fields[echoid_field]['data'])
            echoid[mask] = 1

            rmax_prec = dscfg.get('rmax_prec', 0.)
            percent_prec_max = dscfg.get('percent_prec_max', 10.)

            ngates = radar.ngates
            if rmax_prec > 0.:
                ngates = len(
                    radar.range['data'][radar.range['data'] < rmax_prec])
            ngates_total = ngates*radar.nrays

            prec_field = echoid[:, :ngates]
            ngates_prec = np.size(prec_field[prec_field == 3])

            percent_prec = ngates_prec/ngates_total*100.
            warn('Percent gates with precipitation: '+str(percent_prec)+'\n')
            if percent_prec > percent_prec_max:
                if filter_prec == 'keep_dry':
                    warn('Radar volume is precipitation contaminated.\n' +
                         'Maximum percentage allowed: '+str(percent_prec_max))
                    return None, None
            else:
                if filter_prec == 'keep_wet':
                    warn('Radar volume has not enough precipitation.\n' +
                         'Minimum percentage required: ' +
                         str(percent_prec_max))
                    return None, None

        # prepare field number of samples and occurrence
        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()

        npoints_dict = pyart.config.get_metadata('number_of_samples')
        npoints_dict['data'] = np.ma.ones(
            (radar.nrays, radar.ngates), dtype=int)
        radar_aux.add_field('number_of_samples', npoints_dict)

        occu_dict = pyart.config.get_metadata('occurrence')
        occu_dict['data'] = np.ma.zeros(
            (radar.nrays, radar.ngates), dtype=int)
        occu_dict['data'][np.logical_not(mask)] = 1

        # filter out out of range data
        rmin = dscfg.get('rmin', -1.)
        rmax = dscfg.get('rmax', -1.)
        if rmin >= 0.:
            ind_min = np.where(radar_aux.range['data'] < rmin)[0]
            if ind_min:
                ind_min = ind_min[-1]
                occu_dict['data'][:, 0:ind_min+1] = 0
        if rmax >= 0.:
            ind_max = np.where(radar_aux.range['data'] > rmax)[0]
            if ind_max:
                ind_max = ind_max[0]
                occu_dict['data'][:, ind_max:radar_aux.ngates] = 0

        radar_aux.add_field('occurrence', occu_dict)

        # first volume: initialize radar object
        if dscfg['initialized'] == 0:
            new_dataset = {
                'radar_out': radar_aux,
                'starttime': dscfg['timeinfo'],
                'endtime': dscfg['timeinfo'],
                'occu_final': False}

            dscfg['global_data'] = new_dataset
            dscfg['initialized'] = 1

            return new_dataset, ind_rad

        # accumulate data
        regular_grid = False
        if 'regular_grid' in dscfg:
            regular_grid = dscfg['regular_grid']

        if not regular_grid:
            occu_interp = interpol_field(
                dscfg['global_data']['radar_out'], radar_aux, 'occurrence')
            npoints_interp = interpol_field(
                dscfg['global_data']['radar_out'], radar_aux,
                'number_of_samples')
        else:
            if radar_aux.nrays != dscfg['global_data']['radar_out'].nrays:
                warn('Unable to accumulate radar object. ' +
                     'Number of rays of current radar different from ' +
                     'reference. nrays current: '+str(radar_aux.nrays) +
                     ' nrays ref: ' +
                     str(dscfg['global_data']['radar_out'].nrays))
                return None, None
            occu_interp = radar_aux.fields['occurrence']
            npoints_interp = radar_aux.fields['number_of_samples']

        dscfg['global_data']['radar_out'].fields['occurrence']['data'] += (
            np.ma.asarray(
                occu_interp['data'].filled(fill_value=0)).astype('int'))
        dscfg['global_data']['radar_out'].fields['number_of_samples'][
            'data'] += (np.ma.asarray(
                npoints_interp['data'].filled(fill_value=0)).astype('int'))
        dscfg['global_data']['endtime'] = dscfg['timeinfo']

        new_dataset = {
            'radar_out': dscfg['global_data']['radar_out'],
            'starttime': dscfg['global_data']['starttime'],
            'endtime': dscfg['global_data']['endtime'],
            'occu_final': False}

        return new_dataset, ind_rad

    # no more files to process. Compute frequency of occurrence
    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None, None
        if 'radar_out' not in dscfg['global_data']:
            return None, None

        radar = dscfg['global_data']['radar_out']

        freq_occu_dict = pyart.config.get_metadata('frequency_of_occurrence')
        freq_occu_dict['data'] = (100.*radar.fields['occurrence']['data'] /
                                  radar.fields['number_of_samples']['data'])

        radar.add_field('frequency_of_occurrence', freq_occu_dict)

        new_dataset = {
            'radar_out': dscfg['global_data']['radar_out'],
            'starttime': dscfg['global_data']['starttime'],
            'endtime': dscfg['global_data']['endtime'],
            'occu_final': True}

        return new_dataset, ind_rad


def process_time_avg_std(procstatus, dscfg, radar_list=None):
    """
    computes the average and standard deviation of data. It looks only for
    gates where data is present.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        regular_grid : Boolean. Dataset keyword
            Whether the radar has a Boolean grid or not. Default False
        rmin, rmax : float. Dataset keyword
            minimum and maximum ranges where the computation takes place. If
            -1 the whole range is considered. Default is -1
        val_min : Float. Dataset keyword
            Minimum reflectivity value to consider that the gate has signal.
            Default None
        filter_prec : str. Dataset keyword
            Give which type of volume should be filtered. None, no filtering;
            keep_wet, keep wet volumes; keep_dry, keep dry volumes.
        rmax_prec : float. Dataset keyword
            Maximum range to consider when looking for wet gates [m]
        percent_prec_max : float. Dataset keyword
            Maxim percentage of wet gates to consider the volume dry
        lin_trans : Boolean. Dataset keyword
            If True the data will be transformed into linear units. Default
            False
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """
    echoid_field = None
    refl_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'echoID':
            echoid_field = get_fieldname_pyart(datatype)
        elif (datatype in ('dBZ', 'dBZc', 'dBZv', 'dBZvc', 'dBuZ', 'dBuZc') and
              refl_field is None):
            refl_field = get_fieldname_pyart(datatype)
        else:
            field_name = get_fieldname_pyart(datatype)
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

        # filter out low reflectivity values
        val_min = dscfg.get('val_min', None)
        mask = np.ma.getmaskarray(radar.fields[field_name]['data'])
        if val_min is not None and refl_field is not None:
            mask = np.logical_or(
                mask, radar.fields[refl_field]['data'] < val_min)

        filter_prec = dscfg.get('filter_prec', 'None')
        if filter_prec in ('keep_wet', 'keep_dry'):
            if echoid_field not in radar.fields:
                warn('Unable to determine if there is precipitation ' +
                     'close to the radar. Missing echoID field.')
                return None, None

            # Put invalid values to noise
            echoid = deepcopy(radar.fields[echoid_field]['data'])
            echoid[mask] = 1

            rmax_prec = dscfg.get('rmax_prec', 0.)
            percent_prec_max = dscfg.get('percent_prec_max', 10.)

            ngates = radar.ngates
            if rmax_prec > 0.:
                ngates = len(
                    radar.range['data'][radar.range['data'] < rmax_prec])
            ngates_total = ngates*radar.nrays

            prec_field = echoid[:, :ngates]
            ngates_prec = np.size(prec_field[prec_field == 3])

            percent_prec = ngates_prec/ngates_total*100.
            warn('Percent gates with precipitation: '+str(percent_prec)+'\n')
            if percent_prec > percent_prec_max:
                if filter_prec == 'keep_dry':
                    warn('Radar volume is precipitation contaminated.\n' +
                         'Maximum percentage allowed: '+str(percent_prec_max))
                    return None, None
            else:
                if filter_prec == 'keep_wet':
                    warn('Radar volume has not enough precipitation.\n' +
                         'Minimum percentage required: ' +
                         str(percent_prec_max))
                    return None, None

        # filter out out of range data
        rmin = dscfg.get('rmin', -1.)
        rmax = dscfg.get('rmax', -1.)
        if rmin >= 0.:
            ind_min = np.where(radar.range['data'] < rmin)[0]
            if ind_min:
                ind_min = ind_min[-1]
                mask[:, 0:ind_min+1] = 1
        if rmax >= 0.:
            ind_max = np.where(radar.range['data'] > rmax)[0]
            if ind_max:
                ind_max = ind_max[0]
                mask[:, ind_max:radar.ngates] = 1

        # prepare field number of samples and values sum
        field = deepcopy(radar.fields[field_name]['data'])
        if lin_trans:
            field = np.ma.power(10., 0.1*field)

        field = np.ma.masked_where(mask, field)
        field = np.ma.asarray(field)

        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()

        sum_dict = pyart.config.get_metadata('sum')
        sum_dict['data'] = field
        radar_aux.add_field('sum', sum_dict)

        sum2_dict = pyart.config.get_metadata('sum_squared')
        sum2_dict['data'] = field*field
        radar_aux.add_field('sum_squared', sum2_dict)

        npoints_dict = pyart.config.get_metadata('number_of_samples')
        npoints_dict['data'] = np.ma.asarray(np.logical_not(mask), dtype=int)
        radar_aux.add_field('number_of_samples', npoints_dict)

        # first volume: initialize radar object
        if dscfg['initialized'] == 0:
            new_dataset = {
                'radar_out': radar_aux,
                'starttime': dscfg['timeinfo'],
                'endtime': dscfg['timeinfo'],
                'occu_final': False}

            dscfg['global_data'] = new_dataset
            dscfg['initialized'] = 1

            return new_dataset, ind_rad

        # accumulate data
        regular_grid = False
        if 'regular_grid' in dscfg:
            regular_grid = dscfg['regular_grid']

        if not regular_grid:
            sum_interp = interpol_field(
                dscfg['global_data']['radar_out'], radar_aux, 'sum')
            sum2_interp = interpol_field(
                dscfg['global_data']['radar_out'], radar_aux, 'sum_squared')
            npoints_interp = interpol_field(
                dscfg['global_data']['radar_out'], radar_aux,
                'number_of_samples')
        else:
            if radar_aux.nrays != dscfg['global_data']['radar_out'].nrays:
                warn('Unable to accumulate radar object. ' +
                     'Number of rays of current radar different from ' +
                     'reference. nrays current: '+str(radar_aux.nrays) +
                     ' nrays ref: ' +
                     str(dscfg['global_data']['radar_out'].nrays))
                return None, None
            sum_interp = radar_aux.fields['sum']
            sum2_interp = radar_aux.fields['sum_squared']
            npoints_interp = radar_aux.fields['number_of_samples']

        valid = np.logical_not(np.ma.getmaskarray(sum_interp))

        dscfg['global_data']['radar_out'].fields['sum']['data'][valid] += (
            np.ma.asarray(sum_interp['data'][valid]))
        dscfg['global_data']['radar_out'].fields['sum_squared'][
            'data'][valid] += np.ma.asarray(sum2_interp['data'][valid])
        dscfg['global_data']['radar_out'].fields['number_of_samples'][
            'data'][valid] += np.ma.asarray(npoints_interp['data'][valid])
        dscfg['global_data']['endtime'] = dscfg['timeinfo']

        new_dataset = {
            'radar_out': dscfg['global_data']['radar_out'],
            'starttime': dscfg['global_data']['starttime'],
            'endtime': dscfg['global_data']['endtime'],
            'occu_final': False}

        return new_dataset, ind_rad

    # no more files to process. Compute mean and standard deviation
    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None, None
        if 'radar_out' not in dscfg['global_data']:
            return None, None

        field_mean = (
            dscfg['global_data']['radar_out'].fields['sum']['data'] /
            dscfg['global_data']['radar_out'].fields[
                'number_of_samples']['data'])
        field_std = np.ma.sqrt(
            dscfg['global_data']['radar_out'].fields['sum_squared']['data'] /
            dscfg['global_data']['radar_out'].fields[
                'number_of_samples']['data']-field_mean*field_mean)
        if lin_trans:
            field_mean = 10.*np.ma.log10(field_mean)
            field_std = 10.*np.ma.log10(field_std)

        radar = dscfg['global_data']['radar_out']

        mean_dict = pyart.config.get_metadata(field_name)
        mean_dict['data'] = field_mean
        radar.add_field(field_name, mean_dict)

        std_dict = pyart.config.get_metadata('standard_deviation')
        std_dict['data'] = field_std
        radar.add_field('standard_deviation', std_dict)

        new_dataset = {
            'radar_out': radar,
            'starttime': dscfg['global_data']['starttime'],
            'endtime': dscfg['global_data']['endtime'],
            'occu_final': True}

        return new_dataset, ind_rad


def process_occurrence_period(procstatus, dscfg, radar_list=None):
    """
    computes the frequency of occurrence over a long period of time by adding
    together shorter periods

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        regular_grid : Boolean. Dataset keyword
            Whether the radar has a Boolean grid or not. Default False
        rmin, rmax : float. Dataset keyword
            minimum and maximum ranges where the computation takes place. If
            -1 the whole range is considered. Default is -1
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
        if datatype == 'occurrence':
            occu_field = get_fieldname_pyart(datatype)
        elif datatype == 'nsamples':
            nsamples_field = get_fieldname_pyart(datatype)
    ind_rad = int(radarnr[5:8])-1

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if ((occu_field not in radar.fields) or
                (nsamples_field not in radar.fields)):
            warn('Unable to compute frequency of occurrence. Missing data')
            return None, None

        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()
        radar_aux.add_field('occurrence', radar.fields['occurrence'])
        radar_aux.add_field(
            'number_of_samples', radar.fields['number_of_samples'])

        # filter out out of range data
        rmin = dscfg.get('rmin', -1.)
        rmax = dscfg.get('rmax', -1.)
        if rmin >= 0.:
            ind_min = np.where(radar_aux.range['data'] < rmin)[0]
            if ind_min:
                ind_min = ind_min[-1]
                radar_aux.fields['occurrence']['data'][:, 0:ind_min+1] = 0
        if rmax >= 0.:
            ind_max = np.where(radar_aux.range['data'] > rmax)[0]
            if ind_max:
                ind_max = ind_max[0]
                radar_aux.fields['occurrence']['data'][
                    :, ind_max:radar_aux.ngates] = 0

        # first volume: initialize radar object
        if dscfg['initialized'] == 0:
            new_dataset = {
                'radar_out': radar_aux,
                'starttime': dscfg['timeinfo'],
                'endtime': dscfg['timeinfo'],
                'occu_final': False}

            dscfg['global_data'] = new_dataset
            dscfg['initialized'] = 1

            return new_dataset, ind_rad

        # accumulate data
        regular_grid = False
        if 'regular_grid' in dscfg:
            regular_grid = dscfg['regular_grid']

        if not regular_grid:
            occu_interp = interpol_field(
                dscfg['global_data']['radar_out'], radar_aux, 'occurrence')
            npoints_interp = interpol_field(
                dscfg['global_data']['radar_out'], radar_aux,
                'number_of_samples')
        else:
            if radar_aux.nrays != dscfg['global_data']['radar_out'].nrays:
                warn('Unable to accumulate radar object. ' +
                     'Number of rays of current radar different from ' +
                     'reference. nrays current: '+str(radar_aux.nrays) +
                     ' nrays ref: ' +
                     str(dscfg['global_data']['radar_out'].nrays))
                return None, None
            occu_interp = radar_aux.fields['occurrence']
            npoints_interp = radar_aux.fields['number_of_samples']

        dscfg['global_data']['radar_out'].fields['occurrence']['data'] += (
            np.ma.asarray(
                occu_interp['data'].filled(fill_value=0)).astype('int'))
        dscfg['global_data']['radar_out'].fields['number_of_samples'][
            'data'] += np.ma.asarray(
                npoints_interp['data'].filled(fill_value=0)).astype('int')
        dscfg['global_data']['endtime'] = dscfg['timeinfo']

        new_dataset = {
            'radar_out': dscfg['global_data']['radar_out'],
            'starttime': dscfg['global_data']['starttime'],
            'endtime': dscfg['global_data']['endtime'],
            'occu_final': False}

        return new_dataset, ind_rad

    # no more files to process. Compute frequency of occurrence
    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None, None
        if 'radar_out' not in dscfg['global_data']:
            return None, None

        radar = dscfg['global_data']['radar_out']

        freq_occu_dict = pyart.config.get_metadata('frequency_of_occurrence')
        freq_occu_dict['data'] = (100.*radar.fields['occurrence']['data'] /
                                  radar.fields['number_of_samples']['data'])

        radar.add_field('frequency_of_occurrence', freq_occu_dict)

        new_dataset = {
            'radar_out': dscfg['global_data']['radar_out'],
            'starttime': dscfg['global_data']['starttime'],
            'endtime': dscfg['global_data']['endtime'],
            'occu_final': True}

        return new_dataset, ind_rad


def process_sun_hits(procstatus, dscfg, radar_list=None):
    """
    monitoring of the radar using sun hits

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        delev_max : float. Dataset keyword
            maximum elevation distance from nominal radar elevation where to
            look for a sun hit signal [deg]. Default 1.5
        dazim_max : float. Dataset keyword
            maximum azimuth distance from nominal radar elevation where to
            look for a sun hit signal [deg]. Default 1.5
        elmin : float. Dataset keyword
            minimum radar elevation where to look for sun hits [deg].
            Default 1.
        attg : float. Dataset keyword
            gaseous attenuation. Default None
        sun_position : string. Datset keyword
            The function to compute the sun position to use. Can be 'MF' or
            'pysolar'
        sun_hit_method : str. Dataset keyword
            Method used to estimate the power of the sun hit. Can be HS
            (Hildebrand and Sekhon 1974) or Ivic (Ivic 2013)
        rmin : float. Dataset keyword
            minimum range where to look for a sun hit signal [m]. Used in HS
            method. Default 50000.
        hmin : float. Dataset keyword
            minimum altitude where to look for a sun hit signal [m MSL].
            Default 10000. The actual range from which a sun hit signal will
            be search will be the minimum between rmin and the range from
            which the altitude is higher than hmin. Used in HS method. Default
            10000.
        nbins_min : int. Dataset keyword.
            minimum number of range bins that have to contain signal to
            consider the ray a potential sun hit. Default 20 for HS and 8000
            for Ivic.
        npulses_ray : int
            Default number of pulses used in the computation of the ray. If the
            number of pulses is not in radar.instrument_parameters this will be
            used instead. Used in Ivic method. Default 30
        iterations: int
            number of iterations in step 7 of Ivic method. Default 10.
        max_std_pwr : float. Dataset keyword
            maximum standard deviation of the signal power to consider the
            data a sun hit [dB]. Default 2. Used in HS method
        max_std_zdr : float. Dataset keyword
            maximum standard deviation of the ZDR to consider the
            data a sun hit [dB]. Default 2.
        az_width_co : float. Dataset keyword
            co-polar antenna azimuth width (convoluted with sun width) [deg].
            Default None
        el_width_co : float. Dataset keyword
            co-polar antenna elevation width (convoluted with sun width)
            [deg]. Default None
        az_width_cross : float. Dataset keyword
            cross-polar antenna azimuth width (convoluted with sun width)
            [deg]. Default None
        el_width_cross : float. Dataset keyword
            cross-polar antenna elevation width (convoluted with sun width)
            [deg]. Default None
        ndays : int. Dataset keyword
            number of days used in sun retrieval. Default 1
        coeff_band : float. Dataset keyword
            multiplicate coefficient to transform pulse width into receiver
            bandwidth
        frequency : float. Dataset keyword
            the radar frequency [Hz]. If None that of the key
            frequency in attribute instrument_parameters of the radar
            object will be used. If the key or the attribute are not present
            frequency dependent parameters will not be computed
        beamwidth : float. Dataset keyword
            the antenna beamwidth [deg]. If None that of the keys
            radar_beam_width_h or radar_beam_width_v in attribute
            instrument_parameters of the radar object will be used. If the key
            or the attribute are not present the beamwidth dependent
            parameters will not be computed
        pulse_width : float. Dataset keyword
            the pulse width [s]. If None that of the key
            pulse_width in attribute instrument_parameters of the radar
            object will be used. If the key or the attribute are not present
            the pulse width dependent parameters will not be computed
        ray_angle_res : float. Dataset keyword
            the ray angle resolution [deg]. If None that of the key
            ray_angle_res in attribute instrument_parameters of the radar
            object will be used. If the key or the attribute are not present
            the ray angle resolution parameters will not be computed
        AntennaGainH, AntennaGainV : float. Dataset keyword
            the horizontal (vertical) polarization antenna gain [dB].
            If None that of the attribute instrument_parameters of the radar
            object will be used. If the key or the attribute are not present
            the ray angle resolution parameters will not be computed
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    sun_hits_dict : dict
        dictionary containing a radar object, a sun_hits dict and a
        sun_retrieval dictionary
    ind_rad : int
        radar index

    """


    if procstatus == 0:
        return None, None

    if procstatus == 1:
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            if datatype == 'dBm':
                pwrh_field = 'signal_power_hh'
            if datatype == 'dBmv':
                pwrv_field = 'signal_power_vv'
            if datatype == 'ZDRu':
                zdr_field = 'unfiltered_differential_reflectivity'
            if datatype == 'ZDRuc':
                zdr_field = 'corrected_unfiltered_differential_reflectivity'
            if datatype == 'ZDR':
                zdr_field = 'differential_reflectivity'

        ind_rad = int(radarnr[5:8])-1
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if ((pwrh_field not in radar.fields) or
                (pwrv_field not in radar.fields) or
                (zdr_field not in radar.fields)):
            warn('Unable to get sun hits. Missing data')
            return None, None

        # initialize dataset
        if dscfg['initialized'] == 0:
            radar_par = dict()

            freq = dscfg.get('frequency', None)
            if freq is None:
                if (radar.instrument_parameters is not None and
                        'frequency' in radar.instrument_parameters):
                    freq = radar.instrument_parameters['frequency']['data'][0]
            if freq is None:
                warn('Radar frequency unknown.')
            else:
                radar_par.update({'wavelen': 3e8/freq})

            beamwidth = dscfg.get('beamwidth', None)
            if beamwidth is None:
                if radar.instrument_parameters is not None:
                    if 'radar_beam_width_h' in radar.instrument_parameters:
                        beamwidth = radar.instrument_parameters[
                            'radar_beam_width_h']['data'][0]
                    elif 'radar_beam_width_v' in radar.instrument_parameters:
                        beamwidth = radar.instrument_parameters[
                            'radar_beam_width_v']['data'][0]
            if beamwidth is None:
                warn('Antenna beam width unknown.')
            else:
                radar_par.update({'beamwidth': beamwidth})

            pulse_width = dscfg.get('pulse_width', None)
            if pulse_width is None:
                if (radar.instrument_parameters is not None and
                        'pulse_width' in radar.instrument_parameters):
                    pulse_width = radar.instrument_parameters['pulse_width'][
                        'data'][0]
            if pulse_width is None:
                warn('Pulse width unknown.')
            else:
                radar_par.update({'pulse_width': pulse_width})

            ray_angle_res = dscfg.get('ray_angle_res', None)
            if ray_angle_res is None:
                if radar.ray_angle_res is not None:
                    ray_angle_res = radar.ray_angle_res['data'][0]
            if ray_angle_res is None:
                warn('Angular resolution unknown.')
            else:
                radar_par.update({'angle_step': ray_angle_res})

            antenna_gain_h = dscfg.get('AntennaGainH', None)
            if antenna_gain_h is None:
                if (radar.instrument_parameters is not None and
                        'radar_antenna_gain_h' in radar.instrument_parameters):
                    antenna_gain_h = (
                        radar.instrument_parameters['radar_antenna_gain_h'][
                            'data'][0])
            if antenna_gain_h is None:
                warn('Horizontal antenna gain unknown.')
            else:
                radar_par.update({'antenna_gain_h': antenna_gain_h})

            antenna_gain_v = dscfg.get('AntennaGainV', None)
            if antenna_gain_v is None:
                if (radar.instrument_parameters is not None and
                        'radar_antenna_gain_v' in radar.instrument_parameters):
                    antenna_gain_v = (
                        radar.instrument_parameters['radar_antenna_gain_v'][
                            'data'][0])
            if antenna_gain_v is None:
                warn('Vertical antenna gain unknown.')
            else:
                radar_par.update({'antenna_gain_v': antenna_gain_v})

            radar_par.update({'timeinfo': dscfg['timeinfo']})
            dscfg['global_data'] = radar_par
            dscfg['initialized'] = 1

        dscfg['global_data']['timeinfo'] = dscfg['timeinfo']

        # user values
        delev_max = dscfg.get('delev_max', 1.5)
        dazim_max = dscfg.get('dazim_max', 1.5)
        elmin = dscfg.get('elmin', 1.)
        attg = dscfg.get('attg', None)
        max_std_zdr = dscfg.get('max_std_zdr', 2.)
        sun_hit_method = dscfg.get('sun_hit_method', 'HS')
        sun_position = dscfg.get('sun_position', 'MF')

        if sun_hit_method == 'HS':
            rmin = dscfg.get('rmin', 50000.)
            hmin = dscfg.get('hmin', 10000.)
            nbins_min = dscfg.get('nbins_min', 20)
            max_std_pwr = dscfg.get('max_std_pwr', 2.)

            sun_hits, new_radar = pyart.correct.get_sun_hits(
                radar, delev_max=delev_max, dazim_max=dazim_max, elmin=elmin,
                rmin=rmin, hmin=hmin, nbins_min=nbins_min,
                max_std_pwr=max_std_pwr, max_std_zdr=max_std_zdr,
                attg=attg, sun_position=sun_position, pwrh_field=pwrh_field,
                pwrv_field=pwrv_field, zdr_field=zdr_field)
        else:
            npulses_ray = dscfg.get('npulses_ray', 30)

            nbins_min = dscfg.get('nbins_min', 800)
            iterations = dscfg.get('iterations', 10)

            sun_hits, new_radar = pyart.correct.get_sun_hits_ivic(
                radar, delev_max=delev_max, dazim_max=dazim_max, elmin=elmin,
                npulses_ray=npulses_ray, nbins_min=nbins_min,
                iterations=iterations, attg=attg, max_std_zdr=max_std_zdr,
                sun_position=sun_position, pwrh_field=pwrh_field,
                pwrv_field=pwrv_field, zdr_field=zdr_field)


        if sun_hits is None:
            return None, None

        sun_hits_dataset = dict()
        sun_hits_dataset.update({'sun_hits': sun_hits})
        sun_hits_dataset.update({'radar_out': new_radar})
        sun_hits_dataset.update({'timeinfo': dscfg['timeinfo']})

        return sun_hits_dataset, ind_rad

    if procstatus == 2:
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            break

        ind_rad = int(radarnr[5:8])-1

        # user values
        az_width_co = dscfg.get('az_width_co', None)
        el_width_co = dscfg.get('el_width_co', None)
        az_width_cross = dscfg.get('az_width_cross', None)
        el_width_cross = dscfg.get('el_width_cross', None)
        nfiles = dscfg.get('ndays', 1)

        sun_hits = read_sun_hits_multiple_days(
            dscfg, dscfg['global_data']['timeinfo'], nfiles=nfiles)

        if sun_hits[0] is None:
            return None, None

        sun_pwr_h = sun_hits[7]
        sun_pwr_v = sun_hits[11]

        # get DRAO reference
        sf_ref = np.ma.asarray(np.ma.masked)
        ref_time = None
        if 'wavelen' in dscfg['global_data']:

            flx_dt, flx_val = read_solar_flux(
                dscfg['solarfluxpath']+'fluxtable.txt')

            if flx_dt is not None:
                flx_dt_closest, flx_val_closest = get_closest_solar_flux(
                    sun_hits[0], flx_dt, flx_val)

                # flux at radar wavelength
                sf_radar = pyart.correct.solar_flux_lookup(
                    flx_val_closest, dscfg['global_data']['wavelen'])

                sf_ref = np.ma.asarray(sf_radar[-1])
                ref_time = flx_dt_closest[-1]

                # scaling of the power to account for solar flux variations.
                # The last sun hit is the reference. The scale factor is in dB
                scale_factor = -10.*np.log10(sf_radar/sf_ref)
                sun_pwr_h += scale_factor
                sun_pwr_v += scale_factor
            else:
                warn('Unable to compute solar power reference. ' +
                     'Missing DRAO data')
        else:
            warn('Unable to compute solar power reference. ' +
                 'Missing radar wavelength')

        sun_retrieval_h = pyart.correct.sun_retrieval(
            sun_hits[4], sun_hits[6], sun_hits[3], sun_hits[5],
            sun_pwr_h, sun_hits[8],
            az_width_co=az_width_co, el_width_co=el_width_co,
            az_width_cross=az_width_cross, el_width_cross=el_width_cross,
            is_zdr=False)

        sun_retrieval_v = pyart.correct.sun_retrieval(
            sun_hits[4], sun_hits[6], sun_hits[3], sun_hits[5],
            sun_pwr_v, sun_hits[12],
            az_width_co=az_width_co, el_width_co=el_width_co,
            az_width_cross=az_width_cross, el_width_cross=el_width_cross,
            is_zdr=False)

        sun_retrieval_zdr = pyart.correct.sun_retrieval(
            sun_hits[4], sun_hits[6], sun_hits[3], sun_hits[5],
            sun_hits[15], sun_hits[16],
            az_width_co=az_width_co, el_width_co=el_width_co,
            az_width_cross=az_width_cross, el_width_cross=el_width_cross,
            is_zdr=True)

        sun_retrieval_dict = {
            'first_hit_time': sun_hits[0][0],
            'last_hit_time': sun_hits[0][-1],
            'dBm_sun_est': np.ma.asarray(np.ma.masked),
            'std(dBm_sun_est)': np.ma.asarray(np.ma.masked),
            'sf_h': np.ma.asarray(np.ma.masked),
            'az_bias_h': np.ma.asarray(np.ma.masked),
            'el_bias_h': np.ma.asarray(np.ma.masked),
            'az_width_h': np.ma.asarray(np.ma.masked),
            'el_width_h': np.ma.asarray(np.ma.masked),
            'nhits_h': 0,
            'par_h': None,
            'dBmv_sun_est': np.ma.asarray(np.ma.masked),
            'std(dBmv_sun_est)': np.ma.asarray(np.ma.masked),
            'sf_v': np.ma.asarray(np.ma.masked),
            'az_bias_v': np.ma.asarray(np.ma.masked),
            'el_bias_v': np.ma.asarray(np.ma.masked),
            'az_width_v': np.ma.asarray(np.ma.masked),
            'el_width_v': np.ma.asarray(np.ma.masked),
            'nhits_v': 0,
            'par_v': None,
            'ZDR_sun_est': np.ma.asarray(np.ma.masked),
            'std(ZDR_sun_est)': np.ma.asarray(np.ma.masked),
            'az_bias_zdr': np.ma.asarray(np.ma.masked),
            'el_bias_zdr': np.ma.asarray(np.ma.masked),
            'nhits_zdr': 0,
            'par_zdr': None,
            'sf_ref': np.ma.asarray(sf_ref),
            'ref_time': ref_time,
            'lant': np.ma.asarray(np.ma.masked)}

        if sun_retrieval_h is not None:
            # correct for scanning losses and the polarization of the antenna
            if (('angle_step' in dscfg['global_data']) and
                    ('beamwidth' in dscfg['global_data'])):
                lant = pyart.correct.scanning_losses(
                    dscfg['global_data']['angle_step'],
                    dscfg['global_data']['beamwidth'])

            else:
                warn('Unable to estimate scanning losses. ' +
                     'Missing radar parameters. ' +
                     'Antenna losses will be neglected')
                lant = 0.
            ptoa_h = sun_retrieval_h[0]+lant+3.

            # compute observed solar flux
            if (('pulse_width' in dscfg['global_data']) and
                    ('wavelen' in dscfg['global_data']) and
                    ('antenna_gain_h' in dscfg['global_data'])):
                sf_h = pyart.correct.ptoa_to_sf(
                    ptoa_h, dscfg['global_data']['pulse_width'],
                    dscfg['global_data']['wavelen'],
                    dscfg['global_data']['antenna_gain_h'])
            else:
                warn('Unable to estimate observed solar flux. ' +
                     'Missing radar parameters')
                sf_h = np.ma.asarray(np.ma.masked)

            sun_retrieval_dict['dBm_sun_est'] = np.ma.asarray(ptoa_h)
            sun_retrieval_dict['std(dBm_sun_est)'] = np.ma.asarray(
                sun_retrieval_h[1])
            sun_retrieval_dict['sf_h'] = np.ma.asarray(sf_h)
            sun_retrieval_dict['az_bias_h'] = np.ma.asarray(
                sun_retrieval_h[2])
            sun_retrieval_dict['el_bias_h'] = np.ma.asarray(
                sun_retrieval_h[3])
            sun_retrieval_dict['az_width_h'] = np.ma.asarray(
                sun_retrieval_h[4])
            sun_retrieval_dict['el_width_h'] = np.ma.asarray(
                sun_retrieval_h[5])
            sun_retrieval_dict['nhits_h'] = np.asarray(sun_retrieval_h[6])
            sun_retrieval_dict['par_h'] = np.ma.asarray(sun_retrieval_h[7])
            sun_retrieval_dict['lant'] = np.ma.asarray(lant)
        if sun_retrieval_v is not None:
            # correct for scanning losses and the polarization of the antenna
            if (('angle_step' in dscfg['global_data']) and
                    ('beamwidth' in dscfg['global_data'])):
                lant = pyart.correct.scanning_losses(
                    dscfg['global_data']['angle_step'],
                    dscfg['global_data']['beamwidth'])
            else:
                lant = 0.
                warn('Unable to estimate scanning losses. ' +
                     'Missing radar parameters. ' +
                     'Antenna losses will be neglected')
            ptoa_v = sun_retrieval_v[0]+lant+3.

            # compute observed solar flux
            if (('pulse_width' in dscfg['global_data']) and
                    ('wavelen' in dscfg['global_data']) and
                    ('antenna_gain_v' in dscfg['global_data'])):
                sf_v = pyart.correct.ptoa_to_sf(
                    ptoa_v, dscfg['global_data']['pulse_width'],
                    dscfg['global_data']['wavelen'],
                    dscfg['global_data']['antenna_gain_v'])
            else:
                warn('Unable to estimate observed solar flux. ' +
                     'Missing radar parameters')
                sf_v = np.ma.asarray(np.ma.masked)

            sun_retrieval_dict['dBmv_sun_est'] = np.ma.asarray(ptoa_v)
            sun_retrieval_dict['std(dBmv_sun_est)'] = np.ma.asarray(
                sun_retrieval_v[1])
            sun_retrieval_dict['sf_v'] = np.ma.asarray(sf_v)
            sun_retrieval_dict['az_bias_v'] = np.ma.asarray(
                sun_retrieval_v[2])
            sun_retrieval_dict['el_bias_v'] = np.ma.asarray(
                sun_retrieval_v[3])
            sun_retrieval_dict['az_width_v'] = np.ma.asarray(
                sun_retrieval_v[4])
            sun_retrieval_dict['el_width_v'] = np.ma.asarray(
                sun_retrieval_v[5])
            sun_retrieval_dict['nhits_v'] = np.ma.asarray(sun_retrieval_v[6])
            sun_retrieval_dict['par_v'] = np.ma.asarray(sun_retrieval_v[7])
            sun_retrieval_dict['lant'] = np.ma.asarray(lant)
        if sun_retrieval_zdr is not None:
            sun_retrieval_dict['ZDR_sun_est'] = np.ma.asarray(
                sun_retrieval_zdr[0])
            sun_retrieval_dict['std(ZDR_sun_est)'] = np.ma.asarray(
                sun_retrieval_zdr[1])
            sun_retrieval_dict['az_bias_zdr'] = np.ma.asarray(
                sun_retrieval_zdr[2])
            sun_retrieval_dict['el_bias_dzr'] = np.ma.asarray(
                sun_retrieval_zdr[3])
            sun_retrieval_dict['nhits_zdr'] = np.asarray(sun_retrieval_zdr[6])
            sun_retrieval_dict['par_zdr'] = np.asarray(sun_retrieval_zdr[7])

        sun_hits_dict = {
            'time': sun_hits[0],
            'ray': sun_hits[1],
            'NPrng': sun_hits[2],
            'rad_el': sun_hits[3],
            'rad_az': sun_hits[4],
            'sun_el': sun_hits[5],
            'sun_az': sun_hits[6],
            'dBm_sun_hit': sun_hits[7],
            'std(dBm_sun_hit)': sun_hits[8],
            'NPh': sun_hits[9],
            'NPhval': sun_hits[10],
            'dBmv_sun_hit': sun_hits[11],
            'std(dBmv_sun_hit)': sun_hits[12],
            'NPv': sun_hits[13],
            'NPvval': sun_hits[14],
            'ZDR_sun_hit': sun_hits[15],
            'std(ZDR_sun_hit)': sun_hits[16],
            'NPzdr': sun_hits[17],
            'NPzdrval': sun_hits[18]}

        sun_hits_dataset = dict()
        sun_hits_dataset.update({'sun_hits_final': sun_hits_dict})
        sun_hits_dataset.update({'sun_retrieval': sun_retrieval_dict})
        sun_hits_dataset.update(
            {'timeinfo': dscfg['global_data']['timeinfo']})

        return sun_hits_dataset, ind_rad


def process_sunscan(procstatus, dscfg, radar_list=None):
    """
    Processing of automatic sun scans for monitoring purposes of the radar system.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        delev_max : float. Dataset keyword
            maximum elevation distance from nominal radar elevation where to
            look for a sun hit signal [deg]. Default 1.5
        dazim_max : float. Dataset keyword
            maximum azimuth distance from nominal radar elevation where to
            look for a sun hit signal [deg]. Default 1.5
        elmin : float. Dataset keyword
            minimum radar elevation where to look for sun hits [deg].
            Default 1.
        attg : float. Dataset keyword
            gaseous attenuation. Default None
        sun_position : string. Datset keyword
            The function to compute the sun position to use. Can be 'MF' or
            'pysolar'
        sun_hit_method : str. Dataset keyword
            Method used to estimate the power of the sun hit. Should be PSR. HS
            (Hildebrand and Sekhon 1974) or Ivic (Ivic 2013) are implemented but not tested.
        n_noise_bins : int. Dataset keyword
            Number of bins to use for noise estimation
        noise_threshold : float. Dataset keyword
            Distance over the noise level in [dBm]
        min_num_samples : int. Dataset keyword
            Minimal number of samples above the noise level
        max_fit_stddev : float. Dataset keyword
            Maximal allowed standard deviation for a valid sun fit [dBm]
        do_second_noise_est : string ('Yes' or 'No'). Dataset keyword
            Used to trigger a second noise estimation based on the first fit
            Requires another dataset keyword: n_indfar_bins
        n_indfar_bins : int. Dataset keyword
            Number of samples most remote from the sun center
        az_width_co : float. Dataset keyword
            co-polar antenna azimuth width (convoluted with sun width) [deg].
            Default None
        el_width_co : float. Dataset keyword
            co-polar antenna elevation width (convoluted with sun width)
            [deg]. Default None
        az_width_cross : float. Dataset keyword
            cross-polar antenna azimuth width (convoluted with sun width)
            [deg]. Default None
        el_width_cross : float. Dataset keyword
            cross-polar antenna elevation width (convoluted with sun width)
            [deg]. Default None
        rmin : float. Dataset keyword
            minimum range where to look for a sun hit signal [m]. Used in HS
            method. Default 50000.
        hmin : float. Dataset keyword
            minimum altitude where to look for a sun hit signal [m MSL].
            Default 10000. The actual range from which a sun hit signal will
            be search will be the minimum between rmin and the range from
            which the altitude is higher than hmin. Used in HS method. Default
            10000.
        nbins_min : int. Dataset keyword.
            minimum number of range bins that have to contain signal to
            consider the ray a potential sun hit. Default 20 for HS and 8000
            for Ivic.
        npulses_ray : int
            Default number of pulses used in the computation of the ray. If the
            number of pulses is not in radar.instrument_parameters this will be
            used instead. Used in Ivic method. Default 30
        flat_reg_wlen : int
            Length of the flat region window [m]. Used in Ivic method. Default
            8000.
        iterations: int
            number of iterations in step 7 of Ivic method. Default 10.
        max_std_pwr : float. Dataset keyword
            maximum standard deviation of the signal power to consider the
            data a sun hit [dB]. Default 2. Used in HS method
        max_std_zdr : float. Dataset keyword
            maximum standard deviation of the ZDR to consider the
            data a sun hit [dB]. Default 2.
        ndays : int. Dataset keyword
            number of days used in sun retrieval. Default 1
        coeff_band : float. Dataset keyword
            multiplicate coefficient to transform pulse width into receiver
            bandwidth
        frequency : float. Dataset keyword
            the radar frequency [Hz]. If None that of the key
            frequency in attribute instrument_parameters of the radar
            object will be used. If the key or the attribute are not present
            frequency dependent parameters will not be computed
        beamwidth : float. Dataset keyword
            the antenna beamwidth [deg]. If None that of the keys
            radar_beam_width_h or radar_beam_width_v in attribute
            instrument_parameters of the radar object will be used. If the key
            or the attribute are not present the beamwidth dependent
            parameters will not be computed
        pulse_width : float. Dataset keyword
            the pulse width [s]. If None that of the key
            pulse_width in attribute instrument_parameters of the radar
            object will be used. If the key or the attribute are not present
            the pulse width dependent parameters will not be computed
        ray_angle_res : float. Dataset keyword
            the ray angle resolution [deg]. If None that of the key
            ray_angle_res in attribute instrument_parameters of the radar
            object will be used. If the key or the attribute are not present
            the ray angle resolution parameters will not be computed
        AntennaGainH, AntennaGainV : float. Dataset keyword
            the horizontal (vertical) polarization antenna gain [dB].
            If None that of the attribute instrument_parameters of the radar
            object will be used. If the key or the attribute are not present
            the ray angle resolution parameters will not be computed
    radar_list : list of Radar objects
        Optional. list of radar objects



    Returns
    -------
    sunscan_dataset : dict
        dictionary containing a radar object, a sun_hits dict, a
        sun_retrieval dictionary, field_name and timeinfo
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    pwrh_field = None
    pwrv_field = None
    zdr_field = None

    sun_hit_method = dscfg.get('sun_hit_method', 'PSR')
    sun_position = dscfg.get('sun_position', 'MF')
    n_noise_bins = dscfg.get('n_noise_bins', 8)
    noise_threshold = dscfg.get('noise_threshold', 1.5)
    min_num_samples = dscfg.get('min_num_samples', 60)
    max_fit_stddev = dscfg.get('max_fit_stddev', 0.8)
    do_second_noise_est = dscfg.get('do_second_noise_est', 'Yes')
    n_indfar_bins = dscfg.get('n_indfar_bins', 10)

    az_width_co = dscfg.get('az_width_co', None)
    el_width_co = dscfg.get('el_width_co', None)
    az_width_cross = dscfg.get('az_width_cross', None)
    el_width_cross = dscfg.get('el_width_cross', None)

    delev_max = dscfg.get('delev_max', 3.0)
    dazim_max = dscfg.get('dazim_max', 3.0)
    elmin = dscfg.get('elmin', 5.)
    attg = dscfg.get('attg', None)

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if sun_hit_method == "PSR":
            if datatype == 'NdBmh':
                pwrh_field = get_fieldname_pyart(datatype)
            elif datatype == 'NdBmv':
                pwrv_field = get_fieldname_pyart(datatype)
            elif datatype == 'ZDR':
                zdr_field = get_fieldname_pyart(datatype)
            else:
                warn('ERROR: No valid datatype')
        else:
            if datatype == 'dBm':
                pwrh_field = 'signal_power_hh'
            if datatype == 'dBmv':
                pwrv_field = 'signal_power_vv'
            if datatype == 'ZDRu':
                zdr_field = 'unfiltered_differential_reflectivity'
            if datatype == 'ZDRuc':
                zdr_field = 'corrected_unfiltered_differential_reflectivity'
            if datatype == 'ZDR':
                zdr_field = 'differential_reflectivity'

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None

    radar = radar_list[ind_rad]

    if ((pwrh_field not in radar.fields) and
            (pwrv_field not in radar.fields) and
            (zdr_field not in radar.fields)):
        warn('Unable to process sunscan. Missing data')
        return None, None

    # initialize dataset
    if dscfg['initialized'] == 0:
        radar_par = dict()

        freq = dscfg.get('frequency', None)
        if freq is None:
            if (radar.instrument_parameters is not None and
                    'frequency' in radar.instrument_parameters):
                freq = radar.instrument_parameters['frequency']['data'][0]
        if freq is None:
            warn('Radar frequency unknown.')
        else:
            radar_par.update({'wavelen': 3e8/freq})

        beamwidth = dscfg.get('beamwidth', None)
        if beamwidth is None:
            if radar.instrument_parameters is not None:
                if 'radar_beam_width_h' in radar.instrument_parameters:
                    beamwidth = radar.instrument_parameters[
                        'radar_beam_width_h']['data'][0]
                elif 'radar_beam_width_v' in radar.instrument_parameters:
                    beamwidth = radar.instrument_parameters[
                        'radar_beam_width_v']['data'][0]
        if beamwidth is None:
            warn('Antenna beam width unknown.')
        else:
            radar_par.update({'beamwidth': beamwidth})

        pulse_width = dscfg.get('pulse_width', None)
        if pulse_width is None:
            if (radar.instrument_parameters is not None and
                    'pulse_width' in radar.instrument_parameters):
                pulse_width = radar.instrument_parameters['pulse_width'][
                    'data'][0]
        if pulse_width is None:
            warn('Pulse width unknown.')
        else:
            radar_par.update({'pulse_width': pulse_width})

        ray_angle_res = dscfg.get('ray_angle_res', None)
        if ray_angle_res is None:
            if radar.ray_angle_res is not None:
                ray_angle_res = radar.ray_angle_res['data'][0]
        if ray_angle_res is None:
            warn('Angular resolution unknown.')
        else:
            radar_par.update({'angle_step': ray_angle_res})

        antenna_gain_h = dscfg.get('AntennaGainH', None)
        if antenna_gain_h is None:
            if (radar.instrument_parameters is not None and
                    'radar_antenna_gain_h' in radar.instrument_parameters):
                antenna_gain_h = (
                    radar.instrument_parameters['radar_antenna_gain_h'][
                        'data'][0])
        if antenna_gain_h is None:
            warn('Horizontal antenna gain unknown.')
        else:
            radar_par.update({'antenna_gain_h': antenna_gain_h})

        antenna_gain_v = dscfg.get('AntennaGainV', None)
        if antenna_gain_v is None:
            if (radar.instrument_parameters is not None and
                    'radar_antenna_gain_v' in radar.instrument_parameters):
                antenna_gain_v = (
                    radar.instrument_parameters['radar_antenna_gain_v'][
                        'data'][0])
        if antenna_gain_v is None:
            warn('Vertical antenna gain unknown.')
        else:
            radar_par.update({'antenna_gain_v': antenna_gain_v})

        radar_par.update({'timeinfo': dscfg['timeinfo']})
        dscfg['global_data'] = radar_par
        dscfg['initialized'] = 1

    dscfg['global_data']['timeinfo'] = dscfg['timeinfo']

    # get time for sunscan
    time = num2date(radar.time['data'], radar.time['units'],
                    radar.time['calendar'])
    center_index = int(len(time) / 2.)
    sunscan_time = time[center_index]

    if sun_hit_method == 'PSR':
        if datatype == 'NdBmh':
            sunvol = radar.fields[pwrh_field]['data'][:, 0]
        elif datatype == 'NdBmv':
            sunvol = radar.fields[pwrv_field]['data'][:, 0]
        else:
            warn('ERROR: No valid datatype')

        indsorted = np.ma.argsort(sunvol)
        valsorted = sunvol[indsorted]
        valmax = np.ma.max(valsorted)
        noise_level1 = np.mean(valsorted[0:n_noise_bins-1])
        noise_thres = noise_level1 + noise_threshold
        indval = np.ma.where(sunvol > noise_thres)
        nval = np.ma.size(indval)

        if nval < min_num_samples:
            warn('WARNING: Sun scan processing failed! Not enough valid samples')
            return None, None

        sun_hits = pyart.correct.get_sun_hits_psr(
            radar, delev_max=delev_max, dazim_max=dazim_max, elmin=elmin,
            noise_thres=noise_thres, attg=attg, sun_position=sun_position,
            pwrh_field=pwrh_field, pwrv_field=pwrv_field)

    elif sun_hit_method == 'HS':
        rmin = dscfg.get('rmin', 50000.)
        hmin = dscfg.get('hmin', 10000.)
        nbins_min = dscfg.get('nbins_min', 20)
        max_std_zdr = dscfg.get('max_std_zdr', 2.)
        max_std_pwr = dscfg.get('max_std_pwr', 2.)

        sun_hits, new_radar = pyart.correct.get_sun_hits(
            radar, delev_max=delev_max, dazim_max=dazim_max, elmin=elmin,
            rmin=rmin, hmin=hmin, nbins_min=nbins_min,
            max_std_pwr=max_std_pwr, max_std_zdr=max_std_zdr,
            attg=attg, sun_position=sun_position, pwrh_field=pwrh_field,
            pwrv_field=pwrv_field, zdr_field=zdr_field)

    elif sun_hit_method == 'IVIC':
        npulses_ray = dscfg.get('npulses_ray', 30)
        flat_reg_wlen_rng = dscfg.get('flat_reg_wlen', 8000.)
        nbins_min = dscfg.get('nbins_min', 800)
        iterations = dscfg.get('iterations', 10)
        max_std_zdr = dscfg.get('max_std_zdr', 2.)

        r_res = radar.range['data'][1]-radar.range['data'][0]
        flat_reg_wlen = int(flat_reg_wlen_rng/r_res)

        sun_hits, new_radar = pyart.correct.get_sun_hits_ivic(
            radar, delev_max=delev_max, dazim_max=dazim_max, elmin=elmin,
            npulses_ray=npulses_ray, flat_reg_wlen=flat_reg_wlen,
            nbins_min=nbins_min, iterations=iterations, attg=attg,
            max_std_zdr=max_std_zdr, sun_position=sun_position,
            pwrh_field=pwrh_field, pwrv_field=pwrv_field,
            zdr_field=zdr_field)
    else:
        warn('Warning: No valid sun_hit_method specified.')

    if sun_hits is None:
        return None, None

    sun_pwr_h = sun_hits['dBm_sun_hit']
    sun_pwr_v = sun_hits['dBmv_sun_hit']

    #Substraction of noise
    vals_lin_h = 10. ** (sun_pwr_h / 10.)
    vals_lin_v = 10. ** (sun_pwr_v / 10.)
    noise_lin = 10. ** (noise_level1 / 10.)
    vals_nonoise_lin_h = vals_lin_h - noise_lin
    vals_nonoise_lin_v = vals_lin_v - noise_lin

    indnotval_h = np.ma.where(vals_nonoise_lin_h <= 1e-37)
    indnotval_v = np.ma.where(vals_nonoise_lin_v <= 1e-37)

    if np.ma.size(indnotval_h) > 0:
        warn('WARNING: SunScan: Too small linear values!')
        vals_nonoise_lin_h[indnotval_h] = np.ma.masked
    if np.ma.size(indnotval_v) > 0:
        warn('WARNING: SunScan: Too small linear values!')
        vals_nonoise_lin_v[indnotval_v] = np.ma.masked

    if pwrh_field is not None:
        vals_nonoise_db = 10. * np.ma.log10(vals_nonoise_lin_h)

    if pwrv_field is not None:
        vals_nonoise_db = 10. * np.ma.log10(vals_nonoise_lin_v)

    valmax_nonoise = np.ma.max(vals_nonoise_db)

    # get DRAO reference
    sf_ref = np.ma.asarray(np.ma.masked)
    ref_time = None
    if 'wavelen' in dscfg['global_data']:

        flx_dt, flx_val = read_solar_flux(
            dscfg['solarfluxpath']+'fluxtable.txt')

        if flx_dt is not None:
            flx_dt_closest, flx_val_closest = get_closest_solar_flux(
                sun_hits['time'], flx_dt, flx_val)

            # flux at radar wavelength
            sf_radar = pyart.correct.solar_flux_lookup(
                flx_val_closest, dscfg['global_data']['wavelen'])

            sf_ref = np.ma.asarray(sf_radar[-1])
            ref_time = flx_dt_closest[-1]

            # scaling of the power to account for solar flux variations.
            # The last sun hit is the reference. The scale factor is in dB
            scale_factor = -10.*np.log10(sf_radar/sf_ref)
            if sun_hit_method == 'PSR':
                vals_nonoise_db += scale_factor
            else:
                sun_pwr_h += scale_factor
                sun_pwr_v += scale_factor
        else:
            warn('Unable to compute solar power reference. ' +
                 'Missing DRAO data')
    else:
        warn('Unable to compute solar power reference. ' +
             'Missing radar wavelength')

    sun_retrieval_h = None
    sun_retrieval_v = None

    if sun_hit_method == 'PSR':
        if pwrh_field is not None:
            sun_retrieval_h = pyart.correct.sun_retrieval(
                sun_hits['rad_az'], sun_hits['sun_az'], sun_hits['rad_el'], sun_hits['sun_el'],
                vals_nonoise_db, sun_hits['std(dBm_sun_hit)'],
                az_width_co=az_width_co, el_width_co=el_width_co,
                az_width_cross=az_width_cross, el_width_cross=el_width_cross,
                is_zdr=False)

            azoff = np.ma.asarray(sun_retrieval_h[2])
            eloff = np.ma.asarray(sun_retrieval_h[3])
            par = np.ma.asarray(sun_retrieval_h[7])
            fit_stddev = np.ma.asarray(sun_retrieval_h[1])
            if fit_stddev >= max_fit_stddev:
                return None, None

        if pwrv_field is not None:
            sun_retrieval_v = pyart.correct.sun_retrieval(
                sun_hits['rad_az'], sun_hits['sun_az'], sun_hits['rad_el'], sun_hits['sun_el'],
                vals_nonoise_db, sun_hits['std(dBmv_sun_hit)'],
                az_width_co=az_width_co, el_width_co=el_width_co,
                az_width_cross=az_width_cross, el_width_cross=el_width_cross,
                is_zdr=False)

            azoff = np.ma.asarray(sun_retrieval_v[2])
            eloff = np.ma.asarray(sun_retrieval_v[3])
            par = np.ma.asarray(sun_retrieval_v[7])
            fit_stddev = np.ma.asarray(sun_retrieval_v[1])
            if fit_stddev >= max_fit_stddev:
                return None, None

        sundist = np.zeros_like(sunvol)
        sunpwrmat = np.zeros_like(sunvol)
        sunvalmat = np.zeros_like(sunvol)
        sunpos_el = np.zeros_like(sunvol)
        sunpos_az = np.zeros_like(sunvol)

        # get time at each ray
        time = num2date(radar.time['data'], radar.time['units'],
                        radar.time['calendar'])

        for ray in range(radar.nrays):
            if _PYSOLAR_AVAILABLE and sun_position == 'pysolar':
                elev_sun, azim_sun = pyart.correct.sun_position_pysolar(
                    time[ray], radar.latitude['data'][0],
                    radar.longitude['data'][0])
            else:
                elev_sun, azim_sun = pyart.correct.sun_position_mfr(
                    time[ray], radar.latitude['data'][0],
                    radar.longitude['data'][0], refraction=True)

            #azshift?
            delev = np.ma.abs(radar.elevation['data'][ray]-elev_sun)
            dazim = np.ma.abs(
                (radar.azimuth['data'][ray]-azim_sun) *
                np.ma.cos(elev_sun*np.pi/180.))
            if dazim > 360.:
                dazim -= 360.

            sundist[ray] = np.sqrt((dazim - azoff)**2 + (delev - eloff)**2)
            sunpwrmat[ray] = (par[0] + par[1] * dazim + par[2] * delev +
                              par[3] * dazim**2 + par[4] * delev**2)
            sunvalmat[ray] = sunvol[ray]
            sunpos_el[ray] = elev_sun
            sunpos_az[ray] = azim_sun

        #Second noise estimation: removal of sunpower influence
        if do_second_noise_est == 'Yes':
            #Find the samples most remote from the sun center
            inddistsorted = np.argsort(sundist)[::-1]
            indfar = inddistsorted[0:n_indfar_bins-1]
            noise_far_lin = 10. ** (sunvalmat[indfar] / 10.) - 10. ** (sunpwrmat[indfar] / 10.)
            noise_level2 = 10. * np.log10(np.mean(noise_far_lin))
            noise_thres2 = noise_level2 + noise_threshold

            sun_hits = pyart.correct.get_sun_hits_psr(
                radar, delev_max=delev_max, dazim_max=dazim_max, elmin=elmin,
                noise_thres=noise_thres2, attg=attg, sun_position=sun_position,
                pwrh_field=pwrh_field, pwrv_field=pwrv_field)

            if sun_hits is None:
                return None

            sun_pwr_h = sun_hits['dBm_sun_hit']
            sun_pwr_v = sun_hits['dBmv_sun_hit']

            #Substraction of noise
            vals_lin_h = 10. ** (sun_pwr_h / 10.)
            vals_lin_v = 10. ** (sun_pwr_v / 10.)
            noise_lin = 10. ** (noise_level2 / 10.)
            vals_nonoise_lin_h = vals_lin_h - noise_lin
            vals_nonoise_lin_v = vals_lin_v - noise_lin

            indnotval_h = np.ma.where(vals_nonoise_lin_h <= 1e-37)
            indnotval_v = np.ma.where(vals_nonoise_lin_v <= 1e-37)

            if np.ma.size(indnotval_h) > 0:
                warn('WARNING: SunScan: Too small linear values!')
                vals_nonoise_lin_h[indnotval_h] = np.ma.masked
            if np.ma.size(indnotval_v) > 0:
                warn('WARNING: SunScan: Too small linear values!')
                vals_nonoise_lin_v[indnotval_v] = np.ma.masked

            if pwrh_field is not None:
                vals_nonoise_db = 10. * np.ma.log10(vals_nonoise_lin_h)

            if pwrv_field is not None:
                vals_nonoise_db = 10. * np.ma.log10(vals_nonoise_lin_v)

            valmax_nonoise = np.ma.max(vals_nonoise_db)

            # get DRAO reference
            sf_ref = np.ma.asarray(np.ma.masked)
            ref_time = None
            if 'wavelen' in dscfg['global_data']:

                flx_dt, flx_val = read_solar_flux(
                    dscfg['solarfluxpath']+'fluxtable.txt')

                if flx_dt is not None:
                    flx_dt_closest, flx_val_closest = get_closest_solar_flux(
                        sun_hits['time'], flx_dt, flx_val)

                    # flux at radar wavelength
                    sf_radar = pyart.correct.solar_flux_lookup(
                        flx_val_closest, dscfg['global_data']['wavelen'])

                    sf_ref = np.ma.asarray(sf_radar[-1])
                    ref_time = flx_dt_closest[-1]

                    # scaling of the power to account for solar flux variations.
                    # The last sun hit is the reference. The scale factor is in dB
                    scale_factor = -10.*np.log10(sf_radar/sf_ref)
                    vals_nonoise_db += scale_factor
                else:
                    warn('Unable to compute solar power reference. ' +
                         'Missing DRAO data')
            else:
                warn('Unable to compute solar power reference. ' +
                     'Missing radar wavelength')

            sun_retrieval_h = None
            sun_retrieval_v = None

            if pwrh_field is not None:
                sun_retrieval_h = pyart.correct.sun_retrieval(
                    sun_hits['rad_az'], sun_hits['sun_az'], sun_hits['rad_el'], sun_hits['sun_el'],
                    vals_nonoise_db, sun_hits['std(dBm_sun_hit)'],
                    az_width_co=az_width_co, el_width_co=el_width_co,
                    az_width_cross=az_width_cross, el_width_cross=el_width_cross,
                    is_zdr=False)

            if pwrv_field is not None:
                sun_retrieval_v = pyart.correct.sun_retrieval(
                    sun_hits['rad_az'], sun_hits['sun_az'], sun_hits['rad_el'], sun_hits['sun_el'],
                    vals_nonoise_db, sun_hits['std(dBmv_sun_hit)'],
                    az_width_co=az_width_co, el_width_co=el_width_co,
                    az_width_cross=az_width_cross, el_width_cross=el_width_cross,
                    is_zdr=False)

    else:
        sun_retrieval_h = pyart.correct.sun_retrieval(
            sun_hits[4], sun_hits[6], sun_hits[3], sun_hits[5],
            sun_pwr_h, sun_hits[8],
            az_width_co=az_width_co, el_width_co=el_width_co,
            az_width_cross=az_width_cross, el_width_cross=el_width_cross,
            is_zdr=False)

        sun_retrieval_v = pyart.correct.sun_retrieval(
            sun_hits[4], sun_hits[6], sun_hits[3], sun_hits[5],
            sun_pwr_v, sun_hits[12],
            az_width_co=az_width_co, el_width_co=el_width_co,
            az_width_cross=az_width_cross, el_width_cross=el_width_cross,
            is_zdr=False)

        sun_retrieval_zdr = pyart.correct.sun_retrieval(
            sun_hits[4], sun_hits[6], sun_hits[3], sun_hits[5],
            sun_hits[15], sun_hits[16],
            az_width_co=az_width_co, el_width_co=el_width_co,
            az_width_cross=az_width_cross, el_width_cross=el_width_cross,
            is_zdr=True)

    sun_retrieval_dict = {
        'first_hit_time': sun_hits['time'][0],
        'sunscan_time': sunscan_time,
        'last_hit_time': sun_hits['time'][-1],
        'sunpos_el': sunpos_el[center_index],
        'sunpos_az': sunpos_az[center_index],
        'sun_maxpwr_noise': np.ma.asarray(np.ma.masked),
        'sun_maxpwr_nonoise': np.ma.asarray(np.ma.masked),
        'noise_pwr': np.ma.asarray(np.ma.masked),
        'dBm_sun_est': np.ma.asarray(np.ma.masked),
        'dBm_sun_est_toa': np.ma.asarray(np.ma.masked),
        'std(dBm_sun_est)': np.ma.asarray(np.ma.masked),
        'sf_h': np.ma.asarray(np.ma.masked),
        'az_bias_h': np.ma.asarray(np.ma.masked),
        'el_bias_h': np.ma.asarray(np.ma.masked),
        'az_width_h': np.ma.asarray(np.ma.masked),
        'el_width_h': np.ma.asarray(np.ma.masked),
        'nhits_h': 0,
        'par_h': None,
        'dBmv_sun_est': np.ma.asarray(np.ma.masked),
        'dBmv_sun_est_toa': np.ma.asarray(np.ma.masked),
        'std(dBmv_sun_est)': np.ma.asarray(np.ma.masked),
        'sf_v': np.ma.asarray(np.ma.masked),
        'az_bias_v': np.ma.asarray(np.ma.masked),
        'el_bias_v': np.ma.asarray(np.ma.masked),
        'az_width_v': np.ma.asarray(np.ma.masked),
        'el_width_v': np.ma.asarray(np.ma.masked),
        'nhits_v': 0,
        'par_v': None,
        'ZDR_sun_est': np.ma.asarray(np.ma.masked),
        'std(ZDR_sun_est)': np.ma.asarray(np.ma.masked),
        'az_bias_zdr': np.ma.asarray(np.ma.masked),
        'el_bias_zdr': np.ma.asarray(np.ma.masked),
        'nhits_zdr': 0,
        'par_zdr': None,
        'sf_ref': np.ma.asarray(sf_ref),
        'ref_time': ref_time,
        'lant': np.ma.asarray(np.ma.masked)}

    if sun_retrieval_h is not None:
        # correct for scanning losses and the aneotenna polarization
        if (('angle_step' in dscfg['global_data']) and
                ('beamwidth' in dscfg['global_data'])):
            lant = pyart.correct.scanning_losses(
                dscfg['global_data']['angle_step'],
                dscfg['global_data']['beamwidth'])
        else:
            warn('Unable to estimate scanning losses. ' +
                 'Missing radar parameters. ' +
                 'Antenna losses will be neglected')
            lant = 0.
        ptoa_h = sun_retrieval_h[0]+lant+3.

        # compute observed solar flux
        if (('pulse_width' in dscfg['global_data']) and
                ('wavelen' in dscfg['global_data']) and
                ('antenna_gain_h' in dscfg['global_data'])):
            sf_h = pyart.correct.ptoa_to_sf(
                ptoa_h, dscfg['global_data']['pulse_width'],
                dscfg['global_data']['wavelen'],
                dscfg['global_data']['antenna_gain_h'])
        else:
            warn('Unable to estimate observed solar flux. ' +
                 'Missing radar parameters')
            sf_h = np.ma.asarray(np.ma.masked)

        sun_retrieval_dict['sun_maxpwr_noise'] = np.ma.asarray(valmax)
        sun_retrieval_dict['sun_maxpwr_nonoise'] = np.ma.asarray(valmax_nonoise)
        if do_second_noise_est == 'Yes':
            sun_retrieval_dict['noise_pwr'] = np.ma.asarray(noise_level2)
        else:
            sun_retrieval_dict['noise_pwr'] = np.ma.asarray(noise_level1)
        sun_retrieval_dict['dBm_sun_est'] = np.ma.asarray(sun_retrieval_h[0])
        sun_retrieval_dict['dBm_sun_est_toa'] = np.ma.asarray(ptoa_h)
        sun_retrieval_dict['std(dBm_sun_est)'] = np.ma.asarray(
            sun_retrieval_h[1])
        sun_retrieval_dict['sf_h'] = np.ma.asarray(sf_h)
        sun_retrieval_dict['az_bias_h'] = np.ma.asarray(
            sun_retrieval_h[2])
        sun_retrieval_dict['el_bias_h'] = np.ma.asarray(
            sun_retrieval_h[3])
        sun_retrieval_dict['az_width_h'] = np.ma.asarray(
            sun_retrieval_h[4])
        sun_retrieval_dict['el_width_h'] = np.ma.asarray(
            sun_retrieval_h[5])
        sun_retrieval_dict['nhits_h'] = np.asarray(sun_retrieval_h[6])
        sun_retrieval_dict['par_h'] = np.ma.asarray(sun_retrieval_h[7])
        sun_retrieval_dict['lant'] = np.ma.asarray(lant)

    if sun_retrieval_v is not None:
        # correct for scanning losses and the polarization of the antenna
        if (('angle_step' in dscfg['global_data']) and
                ('beamwidth' in dscfg['global_data'])):
            lant = pyart.correct.scanning_losses(
                dscfg['global_data']['angle_step'],
                dscfg['global_data']['beamwidth'])
        else:
            lant = 0.
            warn('Unable to estimate scanning losses. ' +
                 'Missing radar parameters. ' +
                 'Antenna losses will be neglected')
        ptoa_v = sun_retrieval_v[0]+lant+3.

        # compute observed solar flux
        if (('pulse_width' in dscfg['global_data']) and
                ('wavelen' in dscfg['global_data']) and
                ('antenna_gain_v' in dscfg['global_data'])):
            sf_v = pyart.correct.ptoa_to_sf(
                ptoa_v, dscfg['global_data']['pulse_width'],
                dscfg['global_data']['wavelen'],
                dscfg['global_data']['antenna_gain_v'])
        else:
            warn('Unable to estimate observed solar flux. ' +
                 'Missing radar parameters')
            sf_v = np.ma.asarray(np.ma.masked)

        sun_retrieval_dict['sun_maxpwr_noise'] = np.ma.asarray(valmax)
        sun_retrieval_dict['sun_maxpwr_nonoise'] = np.ma.asarray(valmax_nonoise)
        if do_second_noise_est == 'Yes':
            sun_retrieval_dict['noise_pwr'] = np.ma.asarray(noise_level2)
        else:
            sun_retrieval_dict['noise_pwr'] = np.ma.asarray(noise_level1)
        sun_retrieval_dict['dBmv_sun_est'] = np.ma.asarray(sun_retrieval_v[0])
        sun_retrieval_dict['dBmv_sun_est_toa'] = np.ma.asarray(ptoa_v)
        sun_retrieval_dict['std(dBmv_sun_est)'] = np.ma.asarray(
            sun_retrieval_v[1])
        sun_retrieval_dict['sf_v'] = np.ma.asarray(sf_v)
        sun_retrieval_dict['az_bias_v'] = np.ma.asarray(
            sun_retrieval_v[2])
        sun_retrieval_dict['el_bias_v'] = np.ma.asarray(
            sun_retrieval_v[3])
        sun_retrieval_dict['az_width_v'] = np.ma.asarray(
            sun_retrieval_v[4])
        sun_retrieval_dict['el_width_v'] = np.ma.asarray(
            sun_retrieval_v[5])
        sun_retrieval_dict['nhits_v'] = np.ma.asarray(sun_retrieval_v[6])
        sun_retrieval_dict['par_v'] = np.ma.asarray(sun_retrieval_v[7])
        sun_retrieval_dict['lant'] = np.ma.asarray(lant)

    sun_hits_dict = {
        'time': sun_hits['time'],
        'ray': sun_hits['ray'],
        'NPrng': sun_hits['NPrng'],
        'nhits': sun_hits['nhits'],
        'rad_el': sun_hits['rad_el'],
        'rad_az': sun_hits['rad_az'],
        'sun_el': sun_hits['sun_el'],
        'sun_az': sun_hits['sun_az']}

    sunscan_dataset = dict()
    sunscan_dataset.update({'sun_hits': sun_hits_dict})
    sunscan_dataset.update({'sun_retrieval': sun_retrieval_dict})
    sunscan_dataset.update({'field_name': get_fieldname_pyart(datatype)})
    sunscan_dataset.update({'radar_out': radar})
    sunscan_dataset.update({'timeinfo': dscfg['global_data']['timeinfo']})

    return sunscan_dataset, ind_rad
