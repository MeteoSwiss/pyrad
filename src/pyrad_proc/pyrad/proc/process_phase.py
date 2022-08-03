"""
pyrad.proc.process_phase
======================================

Functions for PhiDP and KDP processing and attenuation correction

.. autosummary::
    :toctree: generated/

    process_correct_phidp0
    process_smooth_phidp_single_window
    process_smooth_phidp_double_window
    process_kdp_leastsquare_single_window
    process_kdp_leastsquare_double_window
    process_phidp_kdp_Vulpiani
    process_phidp_kdp_Kalman
    process_phidp_kdp_Maesaka
    process_phidp_kdp_lp
    process_selfconsistency_kdp_phidp
    process_selfconsistency_bias
    process_attenuation

"""

from copy import deepcopy
from warnings import warn

import numpy as np

import pyart

from ..io.io_aux import get_datatype_fields


def process_correct_phidp0(procstatus, dscfg, radar_list=None):
    """
    corrects phidp of the system phase

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        rmin : float. Dataset keyword
            The minimum range where to look for valid data [m]
        rmax : float. Dataset keyword
            The maximum range where to look for valid data [m]
        rcell : float. Dataset keyword
            The length of a continuous cell to consider it valid precip [m]
        Zmin : float. Dataset keyword
            The minimum reflectivity [dBZ]
        Zmax : float. Dataset keyword
            The maximum reflectivity [dBZ]
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
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'PhiDP':
            psidp_field = 'differential_phase'
        if datatype == 'PhiDPc':
            psidp_field = 'corrected_differential_phase'
        if datatype == 'uPhiDP':
            psidp_field = 'uncorrected_differential_phase'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if (refl_field not in radar.fields) or (psidp_field not in radar.fields):
        warn('Unable to correct PhiDP system offset. Missing data')
        return None, None

    ind_rmin = np.where(radar.range['data'] > dscfg['rmin'])[0][0]
    ind_rmax = np.where(radar.range['data'] < dscfg['rmax'])[0][-1]
    r_res = radar.range['data'][1]-radar.range['data'][0]
    min_rcons = int(dscfg['rcell']/r_res)

    if psidp_field.startswith('corrected_'):
        phidp_field = psidp_field
    elif psidp_field.startswith('uncorrected_'):
        phidp_field = psidp_field.replace('uncorrected_', 'corrected_', 1)
    else:
        phidp_field = 'corrected_'+psidp_field

    try:
        phidp = pyart.correct_sys_phase(
            radar, ind_rmin=ind_rmin, ind_rmax=ind_rmax, min_rcons=min_rcons,
            zmin=dscfg['Zmin'], zmax=dscfg['Zmax'], psidp_field=psidp_field,
            refl_field=refl_field, phidp_field=phidp_field)
    except AttributeError as e:
        warn('Could not find correct_sys_phase function...')
        warn('Skipping system phase correction...')
        warn('Please use PyART-MCH to get this functionality')
        phidp = radar.fields[phidp_field]
     
    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(phidp_field, phidp)

    return new_dataset, ind_rad


def process_smooth_phidp_single_window(procstatus, dscfg, radar_list=None):
    """
    corrects phidp of the system phase and smoothes it using one window

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        rmin : float. Dataset keyword
            The minimum range where to look for valid data [m]
        rmax : float. Dataset keyword
            The maximum range where to look for valid data [m]
        rcell : float. Dataset keyword
            The length of a continuous cell to consider it valid precip [m]
        rwind : float. Dataset keyword
            The length of the smoothing window [m]
        Zmin : float. Dataset keyword
            The minimum reflectivity [dBZ]
        Zmax : float. Dataset keyword
            The maximum reflectivity [dBZ]
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
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'PhiDP':
            psidp_field = 'differential_phase'
        if datatype == 'PhiDPc':
            psidp_field = 'corrected_differential_phase'
        if datatype == 'uPhiDP':
            psidp_field = 'uncorrected_differential_phase'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if (refl_field not in radar.fields) or (psidp_field not in radar.fields):
        warn('Unable to smooth PhiDP. Missing data')
        return None, None

    ind_rmin = np.where(radar.range['data'] > dscfg['rmin'])[0][0]
    ind_rmax = np.where(radar.range['data'] < dscfg['rmax'])[0][-1]
    r_res = radar.range['data'][1]-radar.range['data'][0]
    min_rcons = int(dscfg['rcell']/r_res)
    wind_len = int(dscfg['rwind']/r_res)
    min_valid = int(wind_len/2+1)

    if psidp_field.startswith('corrected_'):
        phidp_field = psidp_field
    elif psidp_field.startswith('uncorrected_'):
        phidp_field = psidp_field.replace('uncorrected_', 'corrected_', 1)
    else:
        phidp_field = 'corrected_'+psidp_field

    phidp = pyart.correct.smooth_phidp_single_window(
        radar, ind_rmin=ind_rmin, ind_rmax=ind_rmax, min_rcons=min_rcons,
        zmin=dscfg['Zmin'], zmax=dscfg['Zmax'], wind_len=wind_len,
        min_valid=min_valid, psidp_field=psidp_field, refl_field=refl_field,
        phidp_field=phidp_field)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(phidp_field, phidp)

    return new_dataset, ind_rad


def process_smooth_phidp_double_window(procstatus, dscfg, radar_list=None):
    """
    corrects phidp of the system phase and smoothes it using one window

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        rmin : float. Dataset keyword
            The minimum range where to look for valid data [m]
        rmax : float. Dataset keyword
            The maximum range where to look for valid data [m]
        rcell : float. Dataset keyword
            The length of a continuous cell to consider it valid precip [m]
        rwinds : float. Dataset keyword
            The length of the short smoothing window [m]
        rwindl : float. Dataset keyword
            The length of the long smoothing window [m]
        Zmin : float. Dataset keyword
            The minimum reflectivity [dBZ]
        Zmax : float. Dataset keyword
            The maximum reflectivity [dBZ]
        Zthr : float. Dataset keyword
            The threshold defining wich smoothed data to used [dBZ]
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
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'PhiDP':
            psidp_field = 'differential_phase'
        if datatype == 'PhiDPc':
            psidp_field = 'corrected_differential_phase'
        if datatype == 'uPhiDP':
            psidp_field = 'uncorrected_differential_phase'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if (refl_field not in radar.fields) or (psidp_field not in radar.fields):
        warn('Unable to smooth PhiDP. Missing data')
        return None, None

    ind_rmin = np.where(radar.range['data'] > dscfg['rmin'])[0][0]
    ind_rmax = np.where(radar.range['data'] < dscfg['rmax'])[0][-1]
    r_res = radar.range['data'][1]-radar.range['data'][0]
    min_rcons = int(dscfg['rcell']/r_res)
    swind_len = int(dscfg['rwinds']/r_res)
    smin_valid = int(swind_len/2+1)
    lwind_len = int(dscfg['rwindl']/r_res)
    lmin_valid = int(lwind_len/2+1)

    if psidp_field.startswith('corrected_'):
        phidp_field = psidp_field
    elif psidp_field.startswith('uncorrected_'):
        phidp_field = psidp_field.replace('uncorrected_', 'corrected_', 1)
    else:
        phidp_field = 'corrected_'+psidp_field

    phidp = pyart.correct.smooth_phidp_double_window(
        radar, ind_rmin=ind_rmin, ind_rmax=ind_rmax, min_rcons=min_rcons,
        zmin=dscfg['Zmin'], zmax=dscfg['Zmax'], swind_len=swind_len,
        smin_valid=smin_valid, lwind_len=lwind_len, lmin_valid=lmin_valid,
        zthr=dscfg['Zthr'], psidp_field=psidp_field, refl_field=refl_field,
        phidp_field=phidp_field)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(phidp_field, phidp)

    return new_dataset, ind_rad


def process_phidp_kdp_Maesaka(procstatus, dscfg, radar_list=None):
    """
    Estimates PhiDP and KDP using the method by Maesaka. This method only
    retrieves data in rain (i.e. below the melting layer)

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        rmin : float. Dataset keyword
            The minimum range where to look for valid data [m]
        rmax : float. Dataset keyword
            The maximum range where to look for valid data [m]
        rcell : float. Dataset keyword
            The length of a continuous cell to consider it valid precip [m]
        Zmin : float. Dataset keyword
            The minimum reflectivity [dBZ]
        Zmax : float. Dataset keyword
            The maximum reflectivity [dBZ]
        fzl : float. Dataset keyword
            The freezing level height [m]. Default 2000.
        ml_thickness : float. Dataset keyword
            The melting layer thickness in meters. Default 700.
        beamwidth : float. Dataset keyword
            the antenna beamwidth [deg]. If None that of the keys
            radar_beam_width_h or radar_beam_width_v in attribute
            instrument_parameters of the radar object will be used. If the key
            or the attribute are not present the beamwidth will be set to None
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

    temp_field = None
    iso0_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'PhiDP':
            psidp_field = 'differential_phase'
        if datatype == 'PhiDPc':
            psidp_field = 'corrected_differential_phase'
        if datatype == 'uPhiDP':
            psidp_field = 'uncorrected_differential_phase'
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'TEMP':
            temp_field = 'temperature'
        if datatype == 'H_ISO0':
            iso0_field = 'height_over_iso0'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if ((refl_field not in radar.fields) or
            (psidp_field not in radar.fields)):
        warn('Unable to retrieve PhiDP KDP using the Maesaka approach. ' +
             'Missing data')
        return None, None

    # determine which freezing level reference
    temp_ref = 'temperature'
    if temp_field is None and iso0_field is None:
        warn('Field to obtain the freezing level was not specified. ' +
             'Using fixed freezing level height')
        temp_ref = 'fixed_fzl'
    elif temp_field is not None:
        if temp_field not in radar.fields:
            warn('COSMO temperature field not available. ' +
                 'Using fixed freezing level height')
            temp_ref = 'fixed_fzl'
    elif iso0_field is not None:
        if iso0_field not in radar.fields:
            warn('Height over iso0 field not available. ' +
                 'Using fixed freezing level height')
            temp_ref = 'fixed_fzl'
        else:
            temp_ref = 'height_over_iso0'

    # determine freezing level height if necessary
    fzl = None
    if temp_ref == 'fixed_fzl':
        if 'fzl' in dscfg:
            fzl = dscfg['fzl']
        else:
            fzl = 2000.
            warn('Freezing level height not defined. Using default ' +
                 str(fzl)+' m')

    thickness = 700.
    if 'ml_thickness' in dscfg:
        thickness = dscfg['ml_thickness']

    phidp_field = 'corrected_differential_phase'
    kdp_field = 'corrected_specific_differential_phase'

    radar_aux = deepcopy(radar)

    # correct PhiDP0
    ind_rmin = np.where(radar_aux.range['data'] > dscfg['rmin'])[0][0]
    ind_rmax = np.where(radar_aux.range['data'] < dscfg['rmax'])[0][-1]
    r_res = radar_aux.range['data'][1]-radar_aux.range['data'][0]
    min_rcons = int(dscfg['rcell']/r_res)

    try:
        phidp = pyart.correct.correct_sys_phase(
            radar_aux, ind_rmin=ind_rmin, ind_rmax=ind_rmax, min_rcons=min_rcons,
            zmin=dscfg['Zmin'], zmax=dscfg['Zmax'], psidp_field=psidp_field,
            refl_field=refl_field, phidp_field=phidp_field)
    except AttributeError as e:
        warn('Could not find correct_sys_phase function...')
        warn('Skipping system phase correction...')
        warn('Please use PyART-MCH to get this functionality')
        phidp = radar_aux.fields[phidp_field]
        
    # filter out data in an above the melting layer
    mask = np.ma.getmaskarray(phidp['data'])

    beamwidth = dscfg.get('beamwidth', None)
    if beamwidth is None:
        if radar.instrument_parameters is not None:
            if 'radar_beam_width_h' in radar.instrument_parameters:
                beamwidth = radar.instrument_parameters['radar_beam_width_h'][
                    'data'][0]
            elif 'radar_beam_width_v' in radar.instrument_parameters:
                beamwidth = radar.instrument_parameters['radar_beam_width_v'][
                    'data'][0]
    if beamwidth is None:
        warn('Antenna beam width unknown.')
        
    try:
        mask_fzl, _ = pyart.correct.get_mask_fzl(
            radar_aux, fzl=fzl, doc=15, min_temp=0., max_h_iso0=0.,
            thickness=thickness, beamwidth=beamwidth, temp_field=temp_field,
            iso0_field=iso0_field, temp_ref=temp_ref)
        mask = np.logical_or(mask, mask_fzl)
    except AttributeError as e:
        warn("Masking above FZL is only available in PyART-MCH")
        warn("Please use this library instead of ARM's version")
    
        

    # filter out data with invalid reflectivity
    mask_refl = np.ma.getmaskarray(radar_aux.fields[refl_field]['data'])
    mask = np.logical_or(mask, mask_refl)

    phidp['data'] = np.ma.masked_where(mask, phidp['data'])
    fill_value = pyart.config.get_fillvalue()
    radar_aux.add_field(phidp_field, phidp, replace_existing=True)

    # the return data is not a masked array
    kdp, phidpf, _ = pyart.retrieve.kdp_proc.kdp_maesaka(
        radar_aux, gatefilter=None, method='cg', backscatter=None, Clpf=1.,
        length_scale=None, first_guess=0.01, finite_order='low',
        fill_value=fill_value, psidp_field=phidp_field, kdp_field=kdp_field,
        phidp_field=phidp_field)

    kdp['data'] = np.ma.masked_where(mask, kdp['data'])
    phidpf['data'] = np.ma.masked_where(mask, phidpf['data'])

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar_aux)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(phidp_field, phidpf)
    new_dataset['radar_out'].add_field(kdp_field, kdp)

    return new_dataset, ind_rad


def process_phidp_kdp_lp(procstatus, dscfg, radar_list=None):
    """
    Estimates PhiDP and KDP using a linear programming algorithm.
    This method only retrieves data in rain (i.e. below the melting layer)

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        fzl : float. Dataset keyword
            The freezing level height [m]. Default 2000.
        ml_thickness : float. Dataset keyword
            The melting layer thickness in meters. Default 700.
        beamwidth : float. Dataset keyword
            the antenna beamwidth [deg]. If None that of the keys
            radar_beam_width_h or radar_beam_width_v in attribute
            instrument_parameters of the radar object will be used. If the key
            or the attribute are not present the beamwidth will be set to None
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

    temp_field = None
    iso0_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'PhiDP':
            psidp_field = 'differential_phase'
        if datatype == 'PhiDPc':
            psidp_field = 'corrected_differential_phase'
        if datatype == 'uPhiDP':
            psidp_field = 'uncorrected_differential_phase'
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'RhoHVc':
            rhv_field = 'corrected_cross_correlation_ratio'
        if datatype == 'RhoHV':
            rhv_field = 'cross_correlation_ratio'
        if datatype == 'SNRh':
            snr_field = 'signal_to_noise_ratio_hh'
        if datatype == 'TEMP':
            temp_field = 'temperature'
        if datatype == 'H_ISO0':
            iso0_field = 'height_over_iso0'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if ((refl_field not in radar.fields) or
            (psidp_field not in radar.fields) or
            (rhv_field not in radar.fields) or
            (snr_field not in radar.fields)):
        warn('Unable to retrieve PhiDP and KDP using the LP approach. ' +
             'Missing data')
        return None, None

    # determine which freezing level reference
    temp_ref = 'temperature'
    if temp_field is None and iso0_field is None:
        warn('Field to obtain the freezing level was not specified. ' +
             'Using fixed freezing level height')
        temp_ref = 'fixed_fzl'
    elif temp_field is not None:
        if temp_field not in radar.fields:
            warn('COSMO temperature field not available. ' +
                 'Using fixed freezing level height')
            temp_ref = 'fixed_fzl'
    elif iso0_field is not None:
        if iso0_field not in radar.fields:
            warn('Height over iso0 field not available. ' +
                 'Using fixed freezing level height')
            temp_ref = 'fixed_fzl'
        else:
            temp_ref = 'height_over_iso0'

    # determine freezing level height if necessary
    fzl = None
    if temp_ref == 'fixed_fzl':
        if 'fzl' in dscfg:
            fzl = dscfg['fzl']
        else:
            fzl = 2000.
            warn('Freezing level height not defined. Using default ' +
                 str(fzl)+' m')

    thickness = 700.
    if 'ml_thickness' in dscfg:
        thickness = dscfg['ml_thickness']

    radar_aux = deepcopy(radar)

    # user config
    LP_solver = dscfg.get('LP_solver', 'cvxopt')

    # filter out data in an above the melting layer
    mask = np.ma.getmaskarray(radar_aux.fields[psidp_field]['data'])
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

    mask_fzl, _ = pyart.correct.get_mask_fzl(
        radar_aux, fzl=fzl, doc=15, min_temp=0., max_h_iso0=0.,
        thickness=thickness, beamwidth=beamwidth, temp_field=temp_field,
        iso0_field=iso0_field, temp_ref=temp_ref)
    mask = np.logical_or(mask, mask_fzl)

    # filter out data with invalid reflectivity
    mask_refl = np.ma.getmaskarray(radar_aux.fields[refl_field]['data'])
    mask = np.logical_or(mask, mask_refl)

    radar_aux.fields[psidp_field]['data'] = np.ma.masked_where(
        mask, radar_aux.fields[psidp_field]['data'])

    phidp_field = 'corrected_differential_phase'
    kdp_field = 'corrected_specific_differential_phase'

    phidp, kdp = pyart.correct.phase_proc_lp(
        radar_aux, 0, debug=False, self_const=60000.0,
        low_z=10.0, high_z=53.0, min_phidp=0.01, min_ncp=10.,
        min_rhv=0.6, fzl=4000.0, sys_phase=0.0,
        overide_sys_phase=True, nowrap=None, really_verbose=False,
        LP_solver=LP_solver, refl_field=refl_field, ncp_field=snr_field,
        rhv_field=rhv_field, phidp_field=psidp_field, kdp_field=kdp_field,
        unf_field=phidp_field, window_len=35, proc=1)

    kdp['data'] = np.ma.masked_where(mask, kdp['data'])
    phidp['data'] = np.ma.masked_where(mask, phidp['data'])

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(phidp_field, phidp)
    new_dataset['radar_out'].add_field(kdp_field, kdp)

    return new_dataset, ind_rad


def process_kdp_leastsquare_single_window(procstatus, dscfg, radar_list=None):
    """
    Computes specific differential phase using a piecewise least square method

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        rwind : float. Dataset keyword
            The length of the segment for the least square method [m]
        vectorize : bool. Dataset keyword
            Whether to vectorize the KDP processing. Default false
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
        if datatype == 'PhiDP':
            phidp_field = 'differential_phase'
        if datatype == 'PhiDPc':
            phidp_field = 'corrected_differential_phase'
        if datatype == 'uPhiDP':
            phidp_field = 'uncorrected_differential_phase'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if phidp_field not in radar.fields:
        warn('Unable to retrieve KDP from PhiDP using least square. ' +
             'Missing data')
        return None, None

    r_res = radar.range['data'][1]-radar.range['data'][0]
    wind_len = int(dscfg['rwind']/r_res)
    min_valid = int(wind_len/2+1)
    kdp_field = 'corrected_specific_differential_phase'
    vectorize = dscfg.get('vectorize', False)

    kdp = pyart.retrieve.kdp_leastsquare_single_window(
        radar, wind_len=wind_len, min_valid=min_valid, phidp_field=phidp_field,
        kdp_field=kdp_field, vectorize=vectorize)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(kdp_field, kdp)

    return new_dataset, ind_rad


def process_kdp_leastsquare_double_window(procstatus, dscfg, radar_list=None):
    """
    Computes specific differential phase using a piecewise least square method

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        rwinds : float. Dataset keyword
            The length of the short segment for the least square method [m]
        rwindl : float. Dataset keyword
            The length of the long segment for the least square method [m]
        Zthr : float. Dataset keyword
            The threshold defining which estimated data to use [dBZ]
        vectorize : Bool. Dataset keyword
            Whether to vectorize the KDP processing. Default false
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
        if datatype == 'PhiDP':
            phidp_field = 'differential_phase'
        if datatype == 'PhiDPc':
            phidp_field = 'corrected_differential_phase'
        if datatype == 'uPhiDP':
            phidp_field = 'uncorrected_differential_phase'
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if (phidp_field not in radar.fields) or (refl_field not in radar.fields):
        warn('Unable to retrieve KDP from PhiDP using least square. ' +
             'Missing data')
        return None, None

    r_res = radar.range['data'][1]-radar.range['data'][0]
    swind_len = int(dscfg['rwinds']/r_res)
    smin_valid = int(swind_len/2+1)
    lwind_len = int(dscfg['rwindl']/r_res)
    lmin_valid = int(lwind_len/2+1)
    vectorize = dscfg.get('vectorize', False)

    kdp_field = 'corrected_specific_differential_phase'

    kdp = pyart.retrieve.kdp_leastsquare_double_window(
        radar, swind_len=swind_len, smin_valid=smin_valid,
        lwind_len=lwind_len, lmin_valid=lmin_valid, zthr=dscfg['Zthr'],
        phidp_field=phidp_field, refl_field=refl_field, kdp_field=kdp_field,
        vectorize=vectorize)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(kdp_field, kdp)

    return new_dataset, ind_rad


def process_phidp_kdp_Vulpiani(procstatus, dscfg, radar_list=None):
    """
    Computes specific differential phase and differential phase using the
    method developed by Vulpiani et al.
    The data is assumed to be clutter free and monotonous

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        rwind : float. Dataset keyword
            The length of the segment [m]
        n_iter : int. Dataset keyword
            number of iterations
        interp : boolean. Dataset keyword
            if set non valid values are interpolated using neighbouring valid
            values
        parallel : boolean. Dataset keyword
            if set use parallel computing
        get_phidp : boolean. Datset keyword
            if set the PhiDP computed by integrating the resultant KDP is
            added to the radar field
        frequency : float. Dataset keyword
            the radar frequency [Hz]. If None that of the key
            frequency in attribute instrument_parameters of the radar
            object will be used. If the key or the attribute are not present
            it will be assumed that the radar is C band
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
        if datatype == 'PhiDP':
            phidp_field = 'differential_phase'
        if datatype == 'PhiDPc':
            phidp_field = 'corrected_differential_phase'
        if datatype == 'uPhiDP':
            phidp_field = 'uncorrected_differential_phase'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if phidp_field not in radar.fields:
        warn('Unable to retrieve KDP from PhiDP using least square. ' +
             'Missing data')
        return None, None

    # number of iterations
    n_iter = 3
    if 'n_iter' in dscfg:
        n_iter = dscfg['n_iter']

    # window length (must be even)
    r_res = radar.range['data'][1]-radar.range['data'][0]
    wind_len = int(dscfg['rwind']/r_res)
    if wind_len % 2 == 1:
        wind_len += 1

    # interpolate invalid values?
    interp = 0
    if 'interp' in dscfg:
        interp = dscfg['interp']

    # parallel computing?
    parallel = 1
    if 'parallel' in dscfg:
        parallel = dscfg['parallel']

    # get PhiDP computed from KDP?
    get_phidp = 0
    if 'get_phidp' in dscfg:
        get_phidp = dscfg['get_phidp']

    # get band
    freq = dscfg.get('frequency', None)
    if freq is None:
        if (radar.instrument_parameters is not None and
                'frequency' in radar.instrument_parameters):
            freq = radar.instrument_parameters['frequency']['data'][0]

    band = 'C'
    if freq is None:
        warn('Radar frequency unknown. Default C band will be applied')
    else:
        band = pyart.retrieve.get_freq_band(freq)
        if band is None:
            if freq < 2e9:
                band = 'S'
            elif freq > 12e9:
                band = 'X'
            warn('Radar frequency out of range. Valid bands are S, C or X. ' +
                 band + ' band will be applied')

    kdp_field = 'corrected_specific_differential_phase'
    phidpr_field = 'corrected_differential_phase'

    kdp_dict, phidpr_dict = pyart.retrieve.kdp_vulpiani(
        radar, gatefilter=None, fill_value=None, psidp_field=phidp_field,
        kdp_field=kdp_field, phidp_field=phidpr_field, band=band,
        windsize=wind_len, n_iter=n_iter, interp=interp,
        prefilter_psidp=False, filter_opt=None, parallel=parallel)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(kdp_field, kdp_dict)
    if get_phidp:
        new_dataset['radar_out'].add_field(phidpr_field, phidpr_dict)

    return new_dataset, ind_rad


def process_phidp_kdp_Kalman(procstatus, dscfg, radar_list=None):
    """
    Computes specific differential phase and differential phase using the
    Kalman filter as proposed by Schneebeli et al.
    The data is assumed to be clutter free and continous

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        parallel : boolean. Dataset keyword
            if set use parallel computing
        get_phidp : boolean. Datset keyword
            if set the PhiDP computed by integrating the resultant KDP is
            added to the radar field
        frequency : float. Dataset keyword
            the radar frequency [Hz]. If None that of the key
            frequency in attribute instrument_parameters of the radar
            object will be used. If the key or the attribute are not present
            it will be assumed that the radar is C band
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
        if datatype == 'PhiDP':
            phidp_field = 'differential_phase'
        if datatype == 'PhiDPc':
            phidp_field = 'corrected_differential_phase'
        if datatype == 'uPhiDP':
            phidp_field = 'uncorrected_differential_phase'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if phidp_field not in radar.fields:
        warn('Unable to retrieve KDP from PhiDP using least square. ' +
             'Missing data')
        return None, None

    # User defined options
    parallel = dscfg.get('parallel', 1)
    get_phidp = dscfg.get('get_phidp', 0)

    # get band
    freq = dscfg.get('frequency', None)
    if freq is None:
        if (radar.instrument_parameters is not None and
                'frequency' in radar.instrument_parameters):
            freq = radar.instrument_parameters['frequency']['data'][0]

    band = 'C'
    if freq is None:
        warn('Radar frequency unknown. Default C band will be applied')
    else:
        band = pyart.retrieve.get_freq_band(freq)
        if band is None:
            if freq < 2e9:
                band = 'S'
            elif freq > 12e9:
                band = 'X'
            warn('Radar frequency out of range. Valid bands are S, C or X. ' +
                 band + ' band will be applied')

    kdp_field = 'corrected_specific_differential_phase'
    phidpr_field = 'corrected_differential_phase'

    kdp_dict, _, phidpr_dict = pyart.retrieve.kdp_schneebeli(
        radar, gatefilter=None, fill_value=None, psidp_field=phidp_field,
        kdp_field=kdp_field, phidp_field=phidpr_field, band=band, rcov=0,
        pcov=0, prefilter_psidp=False, filter_opt=None, parallel=parallel)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(kdp_field, kdp_dict)
    if get_phidp:
        new_dataset['radar_out'].add_field(phidpr_field, phidpr_dict)

    return new_dataset, ind_rad


def process_attenuation(procstatus, dscfg, radar_list=None):
    """
    Computes specific attenuation and specific differential attenuation using
    the Z-Phi method and corrects reflectivity and differential reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        ATT_METHOD : float. Dataset keyword
            The attenuation estimation method used. One of the following:
            ZPhi, Philin
        fzl : float. Dataset keyword
            The default freezing level height. It will be used if no
            temperature field name is specified or the temperature field is
            not in the radar object. Default 2000.
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

    temp = None
    iso0 = None
    zdr = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'dBZc':
            refl = 'corrected_reflectivity'
        if datatype == 'PhiDPc':
            phidp = 'corrected_differential_phase'
        if datatype == 'ZDRc':
            zdr = 'corrected_differential_reflectivity'
        if datatype == 'dBZ':
            refl = 'reflectivity'
        if datatype == 'PhiDP':
            phidp = 'differential_phase'
        if datatype == 'ZDR':
            zdr = 'differential_reflectivity'
        if datatype == 'TEMP':
            temp = 'temperature'
        if datatype == 'H_ISO0':
            iso0 = 'height_over_iso0'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if phidp not in radar.fields or refl not in radar.fields:
        warn('Unable to compute attenuation. Missing data')
        return None, None

    # determine which freezing level reference
    temp_ref = 'temperature'
    if temp is None and iso0 is None:
        warn('Field to obtain the freezing level was not specified. ' +
             'Using fixed freezing level height')
        temp_ref = 'fixed_fzl'
    elif temp is not None:
        if temp not in radar.fields:
            warn('COSMO temperature field not available. ' +
                 'Using fixed freezing level height')
            temp_ref = 'fixed_fzl'
    elif iso0 is not None:
        if iso0 not in radar.fields:
            warn('Height over iso0 field not available. ' +
                 'Using fixed freezing level height')
            temp_ref = 'fixed_fzl'
        else:
            temp_ref = 'height_over_iso0'

    # determine freezing level height if necessary
    fzl = None
    if temp_ref == 'fixed_fzl':
        if 'fzl' in dscfg:
            fzl = dscfg['fzl']
        else:
            fzl = 2000.
            warn('Freezing level height not defined. Using default ' +
                 str(fzl)+' m')

    att_method = dscfg.get('ATT_METHOD', 'ZPhi')
    if att_method not in ('ZPhi', 'Philin'):
        raise ValueError(
            'Unknown attenuation correction method. ' +
            'Must be one of the following: [ZPhi, Philin]')

    if att_method == 'ZPhi':
        spec_at, pia, cor_z, spec_diff_at, pida, cor_zdr = (
            pyart.correct.calculate_attenuation_zphi(
                radar, doc=15, fzl=fzl, smooth_window_len=0, a_coef=None,
                beta=None, c=None, d=None, refl_field=refl, phidp_field=phidp,
                zdr_field=zdr, temp_field=temp, iso0_field=iso0,
                spec_at_field=None, pia_field=None, corr_refl_field=None,
                spec_diff_at_field=None, pida_field=None, corr_zdr_field=None,
                temp_ref=temp_ref))
    elif att_method == 'Philin':
        spec_at, pia, cor_z, spec_diff_at, pida, cor_zdr = (
            pyart.correct.calculate_attenuation_philinear(
                radar, doc=15, fzl=fzl, pia_coef=None, pida_coef=None,
                refl_field=refl, phidp_field=phidp, zdr_field=zdr,
                temp_field=temp, iso0_field=iso0, spec_at_field=None,
                pia_field=None, corr_refl_field=None, spec_diff_at_field=None,
                pida_field=None, corr_zdr_field=None, temp_ref=temp_ref))

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()

    new_dataset['radar_out'].add_field('specific_attenuation', spec_at)
    new_dataset['radar_out'].add_field('path_integrated_attenuation', pia)
    new_dataset['radar_out'].add_field('corrected_reflectivity', cor_z)

    if (spec_diff_at is not None) and (cor_zdr is not None):
        new_dataset['radar_out'].add_field(
            'specific_differential_attenuation', spec_diff_at)
        new_dataset['radar_out'].add_field(
            'path_integrated_differential_attenuation', pida)
        new_dataset['radar_out'].add_field(
            'corrected_differential_reflectivity', cor_zdr)
    else:
        warn(
            ' Specific differential attenuation and attenuation ' +
            'corrected differential reflectivity not available')

    return new_dataset, ind_rad
