"""
pyrad.proc.process_retrieve
===========================

Functions for retrieving new moments and products

.. autosummary::
    :toctree: generated/

    process_ccor
    process_signal_power
    process_rcs_pr
    process_rcs
    process_vol_refl
    process_snr
    process_radial_noise_hs
    process_radial_noise_ivic
    process_l
    process_cdr
    process_vpr
    process_rainrate
    process_rainfall_accumulation
    process_bird_density


"""

import datetime
from copy import deepcopy
from warnings import warn

import numpy as np

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..io.read_data_radar import interpol_field

from ..util.radar_utils import time_avg_range


def process_ccor(procstatus, dscfg, radar_list=None):
    """
    Computes the Clutter Correction Ratio, i.e. the ratio between the
    signal without Doppler filtering and the signal with Doppler filtering

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
        dictionary containing the output
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ('dBZ', 'dBZv'):
            filt_field = get_fieldname_pyart(datatype)
        elif datatype in ('dBuZ', 'dBuZv'):
            unfilt_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if filt_field not in radar.fields or unfilt_field not in radar.fields:
        warn('Unable to compute CCOR. Missing fields')
        return None, None

    ccor_field = 'clutter_correction_ratio_hh'
    if 'vv' in filt_field:
        ccor_field = 'clutter_correction_ratio_vv'

    ccor = pyart.retrieve.compute_ccor(
        radar, filt_field=filt_field, unfilt_field=unfilt_field,
        ccor_field=ccor_field)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(ccor_field, ccor)

    return new_dataset, ind_rad


def process_signal_power(procstatus, dscfg, radar_list=None):
    """
    Computes the signal power in dBm

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        mflossh, mflossv : float. Dataset keyword
            The matching filter losses of the horizontal (vertical) channel
            [dB]. If None it will be obtained from the attribute
            radar_calibration of the radar object. Defaults to 0
        radconsth, radconstv : float. Dataset keyword
            The horizontal (vertical) channel radar constant. If None it will
            be obtained from the attribute radar_calibration of the radar
            object
        lrxh, lrxv : float. Global keyword
            The horizontal (vertical) receiver losses from the antenna feed to
            the reference point. [dB] positive value. Default 0
        lradomeh, lradomev : float. Global keyword
            The 1-way dry radome horizontal (vertical) channel losses.
            [dB] positive value. Default 0.
        attg : float. Dataset keyword
            The gas attenuation [dB/km]. If none it will be obtained from the
            attribute radar_calibration of the radar object or assigned
            according to the radar frequency. Defaults to 0.
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
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'dBuZ':
            refl_field = 'unfiltered_reflectivity'
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'dBuZc':
            refl_field = 'corrected_unfiltered_reflectivity'
        if datatype == 'dBZv':
            refl_field = 'reflectivity_vv'
        if datatype == 'dBuZv':
            refl_field = 'unfiltered_reflectivity_vv'
        if datatype == 'dBuZvc':
            refl_field = 'corrected_unfiltered_reflectivity_vv'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if refl_field not in radar.fields:
        warn('Unable to obtain signal power. Missing field '+refl_field)
        return None, None

    lrx = 0.
    lradome = 0.
    if refl_field.endswith('_vv'):
        pwr_field = 'signal_power_vv'

        lmf = dscfg.get('mflossv', None)
        radconst = dscfg.get('radconstv', None)
        if 'lrxv' in dscfg:
            lrx = dscfg['lrxv'][ind_rad]
        if 'lradomev' in dscfg:
            lradome = dscfg['lradomev'][ind_rad]
    else:
        pwr_field = 'signal_power_hh'

        lmf = dscfg.get('mflossh', None)
        radconst = dscfg.get('radconsth', None)
        if 'lrxh' in dscfg:
            lrx = dscfg['lrxh'][ind_rad]
        if 'lradomeh' in dscfg:
            lradome = dscfg['lradomeh'][ind_rad]

    attg = dscfg.get('attg', None)

    s_pwr = pyart.retrieve.compute_signal_power(
        radar, lmf=lmf, attg=attg, radconst=radconst, lrx=lrx,
        lradome=lradome, refl_field=refl_field, pwr_field=pwr_field)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(pwr_field, s_pwr)

    return new_dataset, ind_rad


def process_rcs_pr(procstatus, dscfg, radar_list=None):
    """
    Computes the radar cross-section (assuming a point target) from radar
    reflectivity by first computing the received power and then the RCS from
    it.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        AntennaGainH, AntennaGainV : float. Dataset keyword
            The horizontal (vertical) polarization antenna gain [dB]. If None
            it will be obtained from the attribute instrument_parameters of
            the radar object
        txpwrh, txpwrv : float. Dataset keyword
            The transmitted power of the horizontal (vertical) channel [dBm].
            If None it will be obtained from the attribute radar_calibration
            of the radar object
        mflossh, mflossv : float. Dataset keyword
            The matching filter losses of the horizontal (vertical) channel
            [dB]. If None it will be obtained from the attribute
            radar_calibration of the radar object. Defaults to 0
        radconsth, radconstv : float. Dataset keyword
            The horizontal (vertical) channel radar constant. If None it will
            be obtained from the attribute radar_calibration of the radar
            object
        lrxh, lrxv : float. Global keyword
            The horizontal (vertical) receiver losses from the antenna feed to
            the reference point. [dB] positive value. Default 0
        ltxh, ltxv : float. Global keyword
            The horizontal (vertical) transmitter losses from the output of the
            high power amplifier to the antenna feed. [dB] positive value.
            Default 0
        lradomeh, lradomev : float. Global keyword
            The 1-way dry radome horizontal (vertical) channel losses.
            [dB] positive value. Default 0.
        attg : float. Dataset keyword
            The gas attenuation [dB/km]. If none it will be obtained from the
            attribute radar_calibration of the radar object or assigned
            according to the radar frequency. Defaults to 0.
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
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'dBuZ':
            refl_field = 'unfiltered_reflectivity'
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'dBuZc':
            refl_field = 'corrected_unfiltered_reflectivity'
        if datatype == 'dBZv':
            refl_field = 'reflectivity_vv'
        if datatype == 'dBuZv':
            refl_field = 'unfiltered_reflectivity_vv'
        if datatype == 'dBuZvc':
            refl_field = 'corrected_unfiltered_reflectivity_vv'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if refl_field not in radar.fields:
        warn('Unable to obtain RCS. Missing field '+refl_field)
        return None, None

    lrx = 0.
    lradome = 0.
    ltx = 0.
    if refl_field.endswith('_vv'):
        rcs_field = 'radar_cross_section_vv'

        lmf = dscfg.get('mflossv', None)
        radconst = dscfg.get('radconstv', None)
        antenna_gain = dscfg.get('AntennaGainH', None)
        tx_pwr = dscfg.get('txpwrv', None)
        if 'lrxv' in dscfg:
            lrx = dscfg['lrxv'][ind_rad]
        if 'ltxv' in dscfg:
            ltx = dscfg['ltxv'][ind_rad]
        if 'lradomev' in dscfg:
            lradome = dscfg['lradomev'][ind_rad]
    else:
        rcs_field = 'radar_cross_section_hh'

        lmf = dscfg.get('mflossh', None)
        radconst = dscfg.get('radconsth', None)
        antenna_gain = dscfg.get('AntennaGainV', None)
        tx_pwr = dscfg.get('txpwrh', None)
        if 'lrxh' in dscfg:
            lrx = dscfg['lrxh'][ind_rad]
        if 'ltxh' in dscfg:
            ltx = dscfg['ltxh'][ind_rad]
        if 'lradomeh' in dscfg:
            lradome = dscfg['lradomeh'][ind_rad]

    attg = dscfg.get('attg', None)

    rcs_dict = pyart.retrieve.compute_rcs_from_pr(
        radar, lmf=lmf, attg=attg, radconst=radconst, tx_pwr=tx_pwr,
        antenna_gain=antenna_gain, lrx=lrx, ltx=ltx,
        lradome=lradome, refl_field=refl_field, rcs_field=rcs_field)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(rcs_field, rcs_dict)

    return new_dataset, ind_rad


def process_rcs(procstatus, dscfg, radar_list=None):
    """
    Computes the radar cross-section (assuming a point target) from radar
    reflectivity.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        kw2 : float. Dataset keyowrd
            The water constant
        pulse_width : float. Dataset keyowrd
            The pulse width [s]
        beamwidthv : float. Global keyword
            The vertical polarization antenna beamwidth [deg]. Used if input
            is vertical reflectivity
        beamwidthh : float. Global keyword
            The horizontal polarization antenna beamwidth [deg]. Used if input
            is horizontal reflectivity
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
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'dBuZ':
            refl_field = 'unfiltered_reflectivity'
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'dBuZc':
            refl_field = 'corrected_unfiltered_reflectivity'
        if datatype == 'dBZv':
            refl_field = 'reflectivity_vv'
        if datatype == 'dBuZv':
            refl_field = 'unfiltered_reflectivity_vv'
        if datatype == 'dBuZvc':
            refl_field = 'corrected_unfiltered_reflectivity_vv'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if refl_field not in radar.fields:
        warn('Unable to obtain RCS. Missing field '+refl_field)
        return None, None

    if refl_field.endswith('_vv'):
        rcs_field = 'radar_cross_section_vv'

        beamwidth = dscfg.get('beamwidthv', None)
    else:
        rcs_field = 'radar_cross_section_hh'

        beamwidth = dscfg.get('beamwidthh', None)

    pulse_width = dscfg.get('PulseWidth', None)
    kw2 = dscfg.get('kw2', 0.93)

    rcs_dict = pyart.retrieve.compute_rcs(
        radar, kw2=kw2, pulse_width=pulse_width, beamwidth=beamwidth,
        refl_field=refl_field, rcs_field=rcs_field)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(rcs_field, rcs_dict)

    return new_dataset, ind_rad


def process_vol_refl(procstatus, dscfg, radar_list=None):
    """
    Computes the volumetric reflectivity in 10log10(cm^2 km^-3)

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        freq : float. Dataset keyword
            The radar frequency
        kw : float. Dataset keyword
            The water constant
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
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'dBuZ':
            refl_field = 'unfiltered_reflectivity'
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'dBuZc':
            refl_field = 'corrected_unfiltered_reflectivity'
        if datatype == 'dBZv':
            refl_field = 'reflectivity_vv'
        if datatype == 'dBuZv':
            refl_field = 'unfiltered_reflectivity_vv'
        if datatype == 'dBuZvc':
            refl_field = 'corrected_unfiltered_reflectivity_vv'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if refl_field not in radar.fields:
        warn('Unable to obtain signal power. Missing field '+refl_field)
        return None, None

    if refl_field.endswith('_vv'):
        vol_refl_field = 'volumetric_reflectivity_vv'
    else:
        vol_refl_field = 'volumetric_reflectivity'

    freq = dscfg.get('freq', None)
    kw = dscfg.get('kw', None)

    vol_refl_dict = pyart.retrieve.compute_vol_refl(
        radar, kw=kw, freq=freq, refl_field=refl_field,
        vol_refl_field=vol_refl_field)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(vol_refl_field, vol_refl_dict)

    return new_dataset, ind_rad


def process_snr(procstatus, dscfg, radar_list=None):
    """
    Computes SNR

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The input data type
        output_type : string. Dataset keyword
            The output data type. Either SNRh or SNRv
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
        if datatype == 'dBZ':
            refl = 'reflectivity'
        if datatype == 'dBuZ':
            refl = 'unfiltered_reflectivity'
        if datatype == 'dBZv':
            refl = 'reflectivity_vv'
        if datatype == 'dBuZv':
            refl = 'unfiltered_reflectivity_vv'
        if datatype == 'Nh':
            noise = 'noisedBZ_hh'
        if datatype == 'Nv':
            noise = 'noisedBZ_vv'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if (refl not in radar.fields) or (noise not in radar.fields):
        warn('Unable to compute SNR. Missing data')
        return None, None

    output_type = dscfg.get('output_type', 'SNRh')
    if output_type == 'SNRh':
        snr_field = 'signal_to_noise_ratio_hh'
    else:
        snr_field = 'signal_to_noise_ratio_vv'

    snr = pyart.retrieve.compute_snr(
        radar, refl_field=refl, noise_field=noise,
        snr_field=snr_field)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(snr_field, snr)

    return new_dataset, ind_rad


def process_radial_noise_hs(procstatus, dscfg, radar_list=None):
    """
    Computes the radial noise from the signal power using the
    Hildebrand and Sekhon 1974 method

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The input data type
        rmin : float. Dataset keyword
            The minimum range from which to start the computation
        nbins_min : int. Dataset keyword
            The minimum number of noisy gates to consider the estimation valid
        max_std_pwr : float. Dataset keyword
            The maximum standard deviation of the noise power to consider the
            estimation valid
        get_noise_pos : bool. Dataset keyword
            If True a field flagging the position of the noisy gets will be
            returned
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

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
    pwr_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if pwr_field not in radar.fields:
        warn('Unable to compute radial noise. Missing signal power')
        return None, None

    # user values
    rmin = dscfg.get('rmin', 10000.)
    nbins_min = dscfg.get('nbins_min', 50)
    max_std_pwr = dscfg.get('max_std_pwr', 2.)
    get_noise_pos = dscfg.get('get_noise_pos', False)

    if 'hh' in pwr_field:
        noise_field = 'noisedBm_hh'
        noise_pos_field = 'noise_pos_h'
    else:
        noise_field = 'noisedBm_vv'
        noise_pos_field = 'noise_pos_v'

    ind_rmin = np.where(radar.range['data'] >= rmin)[0]
    if ind_rmin.size == 0:
        warn('No data at ranges further than rmin '+str(rmin)+'  m.')
        return None, None
    ind_rmin = ind_rmin[0]

    noise, noise_pos = pyart.retrieve.compute_radial_noise_hs(
        radar, ind_rmin=ind_rmin, nbins_min=nbins_min, max_std_pwr=max_std_pwr,
        pwr_field=pwr_field, noise_field=noise_field,
        get_noise_pos=get_noise_pos)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(noise_field, noise)
    if noise_pos is not None:
        new_dataset['radar_out'].add_field(noise_pos_field, noise_pos)

    return new_dataset, ind_rad


def process_radial_noise_ivic(procstatus, dscfg, radar_list=None):
    """
    Computes the radial noise from the signal power using the
    Ivic 2013 method

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The input data type
        npulses_ray : int
            Default number of pulses used in the computation of the ray. If
            the number of pulses is not in radar.instrument_parameters this
            will be used instead. Default 30
        ngates_min: int
            minimum number of gates with noise to consider the retrieval
            valid. Default 800
        iterations: int
            number of iterations in step 7. Default 10.
        get_noise_pos : bool
            If true an additional field with gates containing noise according
            to the algorithm is produced
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

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
    pwr_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if pwr_field not in radar.fields:
        warn('Unable to compute radial noise. Missing signal power')
        return None, None

    # user values
    npulses_ray = dscfg.get('npulses_ray', 30)
    ngates_min = dscfg.get('ngates_min', 800)
    iterations = dscfg.get('iterations', 10)
    get_noise_pos = dscfg.get('get_noise_pos', False)

    if 'hh' in pwr_field:
        noise_field = 'noisedBm_hh'
        noise_pos_field = 'noise_pos_h'
    else:
        noise_field = 'noisedBm_vv'
        noise_pos_field = 'noise_pos_v'

    noise, noise_pos = pyart.retrieve.compute_radial_noise_ivic(
        radar, npulses_ray=npulses_ray, ngates_min=ngates_min,
        iterations=iterations, pwr_field=pwr_field, noise_field=noise_field,
        get_noise_pos=get_noise_pos)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(noise_field, noise)
    if noise_pos is not None:
        new_dataset['radar_out'].add_field(noise_pos_field, noise_pos)

    return new_dataset, ind_rad


def process_l(procstatus, dscfg, radar_list=None):
    """
    Computes L parameter

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The input data type
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

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'])
    rhohv = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if rhohv not in radar.fields:
        print('Unable to compute L. Missing RhoHV field')
        return None, None

    l = pyart.retrieve.compute_l(
        radar, rhohv_field=rhohv,
        l_field='logarithmic_cross_correlation_ratio')

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(
        'logarithmic_cross_correlation_ratio', l)

    return new_dataset, ind_rad


def process_cdr(procstatus, dscfg, radar_list=None):
    """
    Computes Circular Depolarization Ratio

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The input data type
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
        if datatype == 'RhoHV':
            rhohv = 'cross_correlation_ratio'
        if datatype == 'uRhoHV':
            rhohv = 'uncorrected_cross_correlation_ratio'
        if datatype == 'RhoHVu':
            rhohv = 'unfiltered_cross_correlation_ratio'
        if datatype == 'ZDR':
            zdr = 'differential_reflectivity'
        if datatype == 'ZDRc':
            zdr = 'corrected_differential_reflectivity'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if ((rhohv not in radar.fields) or
            (zdr not in radar.fields)):
        warn('Unable to compute CDR field. Missing data')
        return None, None

    cdr = pyart.retrieve.compute_cdr(
        radar, rhohv_field=rhohv, zdr_field=zdr,
        cdr_field='circular_depolarization_ratio')

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field('circular_depolarization_ratio', cdr)

    return new_dataset, ind_rad


def process_vpr(procstatus, dscfg, radar_list=None):
    """
    Computes the vertical profile of reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The input data type
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
        if datatype in ('dBZ', 'dBZc'):
            refl_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if refl_field not in radar.fields:
        warn('ERROR: Unable to compute VPR. Missing data')
        return None, None

    refl_corr = pyart.correct.correct_vpr(radar, refl_field=refl_field)
    # User defined variables
    nvalid_min = dscfg.get('nvalid_min', 20)
    angle_min = dscfg.get('angle_min', 0.)
    angle_max = dscfg.get('angle_max', 4.)
    ml_thickness_min = dscfg.get('ml_thickness_min', 200)
    ml_thickness_max = dscfg.get('ml_thickness_max', 800)
    ml_thickness_step = dscfg.get('ml_thickness_step', 200)
    iso0_max = dscfg.get('iso0_max', 5000.)
    ml_top_diff_max = dscfg.get('ml_top_diff_max', 200)
    ml_top_step = dscfg.get('ml_top_step', 200)
    ml_peak_min = dscfg.get('ml_peak_min', 1.)
    ml_peak_max = dscfg.get('ml_peak_max', 6.)
    ml_peak_step = dscfg.get('ml_peak_step', 1.)
    dr_min = dscfg.get('dr_min', -6.)
    dr_max = dscfg.get('dr_max', -1.5)
    dr_step = dscfg.get('dr_step', 1.5)
    dr_default = dscfg.get('dr_default', -4.5)
    dr_alt = dscfg.get('dr_alt', 800.)
    h_max = dscfg.get('h_max', 6000.)
    h_res = dscfg.get('h_res', 1.)
    max_weight = dscfg.get('max_weight', 9.)
    rmin_obs = dscfg.get('rmin_obs', 5000.)
    rmax_obs = dscfg.get('rmax_obs', 150000.)
    use_ml = dscfg.get('use_ml', False)
    ml_datatypedescr = dscfg.get('ml_datatype', None)
    z_datatypedescr = dscfg.get('z_datatype', None)
    vpr_theo_datatypedescr = dscfg.get('vpr_theo_datatype', None)
    vpr_memory_max = dscfg.get('vpr_memory_max', 0.)
    filter_vpr_memory_max = dscfg.get('filter_vpr_memory_max', 0.)
    weight_mem = dscfg.get('weight_mem', 0.75)

    iso0 = None
    if use_ml:
        if (ml_datatypedescr is None or dscfg['loadbasepath'] is None
                or dscfg['loadname'] is None):
            warn('unable to find files containing melting layer information')
            iso0 = None
        else:
            flist = get_file_list(
                ml_datatypedescr, [dscfg['timeinfo']], [dscfg['timeinfo']],
                dscfg)
            if not flist:
                warn('unable to find files containing'
                     ' melting layer information')
                iso0 = None
            else:
                radar_ml = pyart.io.read_cfradial(flist[0])
                if radar_ml is None:
                    warn('Unable to use retrieved melting layer data')
                    iso0 = None
                else:
                    print('Using file {} with melting layer information'.format(
                        flist[0]))
                    iso0 = np.ma.mean(
                        radar_ml.fields['melting_layer_height']['data'][:, 1])
                    ml_bottom = np.ma.mean(
                        radar_ml.fields['melting_layer_height']['data'][:, 0])
                    ml_thickness = iso0-ml_bottom
                    ml_thickness_min = ml_thickness-ml_thickness_step
                    ml_thickness_max = ml_thickness+ml_thickness_step

    radar_mem_list = None
    if vpr_memory_max > 0:
        if (z_datatypedescr is None or dscfg['loadbasepath'] is None
                or dscfg['loadname'] is None):
            warn('unable to find files with azimuthal reflectivity average')
        else:
            # get previous reflectivity files but do not include current one
            flist = get_file_list(
                z_datatypedescr,
                [dscfg['timeinfo']-datetime.timedelta(minutes=vpr_memory_max)],
                [dscfg['timeinfo']-datetime.timedelta(seconds=1)], dscfg)
            if not flist:
                warn('unable to find files containing reflectivity')
            else:
                radar_mem_list = []
                for fname in flist:
                    radar_vpr = pyart.io.read_cfradial(fname)
                    if radar_vpr is None:
                        warn('Unable to use azimuthal reflectivity average')
                        continue
                    print('Using file {} with azimuthal averaged refl'.format(
                          fname))
                    radar_mem_list.append(radar_vpr)
                if not radar_mem_list:
                    radar_mem_list = None

    vpr_theo_dict_mem = None
    if filter_vpr_memory_max > 0.:
        if (vpr_theo_datatypedescr is None or dscfg['loadbasepath'] is None
                or dscfg['loadname'] is None):
            warn('unable to find files with theoretical VPR')
        else:
            # get previous VPR retrieved files but do not include current one
            flist = get_file_list(
                vpr_theo_datatypedescr,
                [dscfg['timeinfo']
                 - datetime.timedelta(minutes=filter_vpr_memory_max)],
                [dscfg['timeinfo']-datetime.timedelta(seconds=1)], dscfg)
            if not flist:
                warn('unable to find files containing retrieved VPR')
            else:
                height, _, vals = read_rhi_profile(flist[-1], labels=['Znorm'])
                vpr_theo_dict_mem = {
                    'value': vals[:, 0],
                    'altitude': height}
                print('Using file {} with previously retrieved VPR'.format(
                      flist[-1]))

    corr_refl_field = 'corrected_reflectivity'
    corr_field = 'vpr_correction'
    refl_corr, vpr_corr, vpr_theo_dict, radar_rhi = pyart.correct.correct_vpr(
        radar, nvalid_min=nvalid_min, angle_min=angle_min,
        angle_max=angle_max, ml_thickness_min=ml_thickness_min,
        ml_thickness_max=ml_thickness_max,
        ml_thickness_step=ml_thickness_step, iso0_max=iso0_max,
        ml_top_diff_max=ml_top_diff_max, ml_top_step=ml_top_step,
        ml_peak_min=ml_peak_min, ml_peak_max=ml_peak_max,
        ml_peak_step=ml_peak_step, dr_min=dr_min, dr_max=dr_max,
        dr_step=dr_step, dr_default=dr_default, dr_alt=dr_alt, h_max=h_max,
        h_res=h_res, max_weight=max_weight, rmin_obs=rmin_obs,
        rmax_obs=rmax_obs, iso0=iso0, weight_mem=weight_mem,
        vpr_theo_dict_mem=vpr_theo_dict_mem, radar_mem_list=radar_mem_list,
        refl_field=refl_field, corr_refl_field=corr_refl_field,
        corr_field=corr_field, temp_field=temp_field, iso0_field=iso0_field,
        temp_ref=temp_ref)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field('corrected_reflectivity', refl_corr)

    return new_dataset, ind_rad


def process_rainrate(procstatus, dscfg, radar_list=None):
    """
    Estimates rainfall rate from polarimetric moments

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The input data type
        RR_METHOD : string. Dataset keyword
            The rainfall rate estimation method. One of the following:
            Z, ZPoly, KDP, A, ZKDP, ZA, hydro
        alpha, beta : float
            factor and exponent of the R-Var power law R = alpha*Var^Beta.
            Default value depending on RR_METHOD. Z (0.0376, 0.6112),
            KDP (None, None), A (None, None)
        alphaz, betaz : float
            factor and exponent of the R-Z power law R = alpha*Z^Beta.
            Default value (0.0376, 0.6112)
        alphazr, betazr : float
            factor and exponent of the R-Z power law R = alpha*Z^Beta applied
            to rain in method hydro. Default value (0.0376, 0.6112)
        alphazs, betazs : float
            factor and exponent of the R-Z power law R = alpha*Z^Beta applied
            to solid precipitation in method hydro. Default value (0.1, 0.5)
        alphakdp, betakdp : float
            factor and exponent of the R-KDP power law R = alpha*KDP^Beta.
            Default value (None, None)
        alphaa, betaa : float
            factor and exponent of the R-Ah power law R = alpha*Ah^Beta.
            Default value (None, None)
        thresh : float
            In hybrid methods, Rainfall rate threshold at which the retrieval
            method used changes [mm/h]. Default value depending on RR_METHOD.
            ZKDP 10, ZA 10, hydro 10
        mp_factor : float
            Factor by which the Z-R relation is multiplied in the melting layer
            in method hydro. Default 0.6
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

    if 'RR_METHOD' not in dscfg:
        raise Exception(
            "ERROR: Undefined parameter 'RR_METHOD' for dataset '%s'"
            % dscfg['dsname'])

    if dscfg['RR_METHOD'] == 'Z':
        radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
        refl_field = get_fieldname_pyart(datatype)

        ind_rad = int(radarnr[5:8])-1
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if refl_field not in radar.fields:
            warn('ERROR: Unable to compute rainfall rate. Missing data')
            return None, None

        # user defined parameters
        alpha = dscfg.get('alpha', 0.0376)
        beta = dscfg.get('beta', 0.6112)

        rain = pyart.retrieve.est_rain_rate_z(
            radar, alpha=alpha, beta=beta, refl_field=refl_field,
            rr_field=None)

    elif dscfg['RR_METHOD'] == 'ZPoly':
        radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
        refl_field = get_fieldname_pyart(datatype)

        ind_rad = int(radarnr[5:8])-1
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if refl_field not in radar.fields:
            warn('Unable to compute rainfall rate. Missing data')
            return None, None

        rain = pyart.retrieve.est_rain_rate_zpoly(
            radar, refl_field=refl_field, rr_field=None)

    elif dscfg['RR_METHOD'] == 'KDP':
        radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
        kdp_field = get_fieldname_pyart(datatype)

        ind_rad = int(radarnr[5:8])-1
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if kdp_field not in radar.fields:
            warn('Unable to compute rainfall rate. Missing data')
            return None, None

        # user defined parameters
        alpha = dscfg.get('alpha', None)
        beta = dscfg.get('beta', None)

        rain = pyart.retrieve.est_rain_rate_kdp(
            radar, alpha=alpha, beta=beta, kdp_field=kdp_field, rr_field=None)

    elif dscfg['RR_METHOD'] == 'A':
        radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
        a_field = get_fieldname_pyart(datatype)

        ind_rad = int(radarnr[5:8])-1
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if a_field not in radar.fields:
            warn('Unable to compute rainfall rate. Missing data')
            return None, None

        # user defined parameters
        alpha = dscfg.get('alpha', None)
        beta = dscfg.get('beta', None)

        rain = pyart.retrieve.est_rain_rate_a(
            radar, alpha=alpha, beta=beta, a_field=a_field, rr_field=None)

    elif dscfg['RR_METHOD'] == 'ZKDP':
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            if datatype == 'dBZc':
                refl_field = 'corrected_reflectivity'
            if datatype == 'KDPc':
                kdp_field = 'corrected_specific_differential_phase'
            if datatype == 'dBZ':
                refl_field = 'reflectivity'
            if datatype == 'KDP':
                kdp_field = 'specific_differential_phase'

        ind_rad = int(radarnr[5:8])-1
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if ((refl_field not in radar.fields) or
                (kdp_field not in radar.fields)):
            warn('Unable to compute rainfall rate. Missing data')
            return None, None

        # user defined parameters
        alphaz = dscfg.get('alphaz', 0.0376)
        betaz = dscfg.get('betaz', 0.6112)
        alphakdp = dscfg.get('alphakdp', None)
        betakdp = dscfg.get('betakdp', None)
        thresh = dscfg.get('thresh', 10.)

        rain = pyart.retrieve.est_rain_rate_zkdp(
            radar, alphaz=alphaz, betaz=betaz, alphakdp=alphakdp,
            betakdp=betakdp, refl_field=refl_field, kdp_field=kdp_field,
            rr_field=None, master_field=refl_field, thresh=thresh,
            thresh_max=True)

    elif dscfg['RR_METHOD'] == 'ZA':
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            if datatype == 'dBZc':
                refl_field = 'corrected_reflectivity'
            if datatype == 'Ahc':
                a_field = 'corrected_specific_attenuation'
            if datatype == 'dBZ':
                refl_field = 'reflectivity'
            if datatype == 'Ah':
                a_field = 'specific_attenuation'

        ind_rad = int(radarnr[5:8])-1
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if ((refl_field not in radar.fields) or
                (a_field not in radar.fields)):
            warn('Unable to compute rainfall rate. Missing data')
            return None, None

        # user defined parameters
        alphaz = dscfg.get('alphaz', 0.0376)
        betaz = dscfg.get('betaz', 0.6112)
        alphaa = dscfg.get('alphaa', None)
        betaa = dscfg.get('betaa', None)
        thresh = dscfg.get('thresh', 5.)

        rain = pyart.retrieve.est_rain_rate_za(
            radar, alphaz=alphaz, betaz=betaz, alphaa=alphaa, betaa=betaa,
            refl_field=refl_field, a_field=a_field, rr_field=None,
            master_field=refl_field, thresh=thresh, thresh_max=True)

    elif dscfg['RR_METHOD'] == 'hydro':
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            if datatype == 'dBZc':
                refl_field = 'corrected_reflectivity'
            if datatype == 'Ahc':
                a_field = 'corrected_specific_attenuation'
            if datatype == 'hydro':
                hydro_field = 'radar_echo_classification'
            if datatype == 'dBZ':
                refl_field = 'reflectivity'
            if datatype == 'Ah':
                a_field = 'specific_attenuation'

        ind_rad = int(radarnr[5:8])-1
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        # user defined parameters
        alphazr = dscfg.get('alphazr', 0.0376)
        betazr = dscfg.get('betazr', 0.6112)
        alphazs = dscfg.get('alphazs', 0.1)
        betazs = dscfg.get('betazs', 0.5)
        alphaa = dscfg.get('alphaa', None)
        betaa = dscfg.get('betaa', None)
        thresh = dscfg.get('thresh', 5.)
        mp_factor = dscfg.get('mp_factor', 0.6)

        if ((refl_field in radar.fields) and
                (a_field in radar.fields) and
                (hydro_field in radar.fields)):
            rain = pyart.retrieve.est_rain_rate_hydro(
                radar, alphazr=alphazr, betazr=betazr, alphazs=alphazs,
                betazs=betazs, alphaa=alphaa, betaa=betaa,
                mp_factor=mp_factor, refl_field=refl_field, a_field=a_field,
                hydro_field=hydro_field, rr_field=None,
                master_field=refl_field, thresh=thresh, thresh_max=True)
        elif refl_field in radar.fields:
            warn('Unable to compute rainfall rate using hydrometeor ' +
                 'classification. Missing data. ' +
                 'A simple Z-R relation will be used instead')
            rain = pyart.retrieve.est_rain_rate_z(
                radar, alpha=alphazr, beta=betazr, refl_field=refl_field,
                rr_field=None)
        else:
            warn('Unable to compute rainfall rate using hydrometeor ' +
                 'classification. Missing data.')
            return None, None
    else:
        raise Exception(
            "ERROR: Unknown rainfall rate retrieval method " +
            dscfg['RR_METHOD'])

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field('radar_estimated_rain_rate', rain)

    return new_dataset, ind_rad


def process_rainfall_accumulation(procstatus, dscfg, radar_list=None):
    """
    Computes rainfall accumulation fields

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
    if procstatus == 0:
        return None, None

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg['datatype'][0])
    if datatype != 'RR':
        warn('Unable to compute rainfall accumulation. Data type: '+datatype +
             '. Expected RR')
        return None, None
    field_name = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1

    start_average = dscfg.get('start_average', 0.)
    period = dscfg.get('period', 3600.)
    use_nan = dscfg.get('use_nan', 0)
    nan_value = dscfg.get('nan_value', 0.)
    rr_acu_name = 'rainfall_accumulation'

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if field_name not in radar.fields:
            warn('ERROR: Unable to compute rainfall accumulation. '
                 'Missing data')
            return None, None

        # Prepare auxiliary radar
        field_data = deepcopy(radar.fields[field_name]['data'])
        if use_nan:
            field_data = np.ma.asarray(field_data.filled(nan_value))

        field_data *= dscfg['ScanPeriod']*60./3600.

        rr_acu_dict = pyart.config.get_metadata(rr_acu_name)
        rr_acu_dict['data'] = field_data

        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()
        radar_aux.add_field(rr_acu_name, rr_acu_dict)

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
                dscfg['global_data']['radar_out'], radar_aux, rr_acu_name)

            if use_nan:
                field_interp['data'] = np.ma.asarray(
                    field_interp['data'].filled(nan_value))

            masked_sum = np.ma.getmaskarray(
                dscfg['global_data']['radar_out'].fields[rr_acu_name]['data'])
            valid_sum = np.logical_and(
                np.logical_not(masked_sum),
                np.logical_not(np.ma.getmaskarray(field_interp['data'])))

            dscfg['global_data']['radar_out'].fields[rr_acu_name]['data'][
                masked_sum] = field_interp['data'][masked_sum]

            dscfg['global_data']['radar_out'].fields[
                rr_acu_name]['data'][valid_sum] += (
                    field_interp['data'][valid_sum])

            return None, None

        # we have reached the end of the accumulation period: start a new
        # object (only reachable if period != -1)
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


def process_bird_density(procstatus, dscfg, radar_list=None):
    """
    Computes the bird density from the volumetric reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        sigma_bird : float. Dataset keyword
            The bird radar cross section
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
        if datatype == 'eta_h':
            vol_refl_field = 'volumetric_reflectivity'
        if datatype == 'eta_v':
            vol_refl_field = 'volumetric_reflectivity_vv'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if vol_refl_field not in radar.fields:
        warn('Unable to obtain bird density. Missing field '+vol_refl_field)
        return None, None

    sigma_bird = dscfg.get('sigma_bird', 11.)
    bird_density_dict = pyart.retrieve.compute_bird_density(
        radar, sigma_bird=sigma_bird, vol_refl_field=vol_refl_field,
        bird_density_field='bird_density')

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field('bird_density', bird_density_dict)

    return new_dataset, ind_rad
