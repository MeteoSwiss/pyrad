"""
pyrad.proc.process_echoclass
===============================

Functions for echo classification and filtering

.. autosummary::
    :toctree: generated/

    process_echo_id
    process_birds_id
    process_clt_to_echo_id
    process_hydro_mf_to_echo_id
    process_hydro_mf_to_hydro
    process_echo_filter
    process_cdf
    process_filter_snr
    process_filter_vel_diff
    process_filter_visibility
    process_outlier_filter
    process_hydroclass
    process_centroids
    process_melting_layer
    process_zdr_column

"""

from copy import deepcopy
from warnings import warn

import sys
import os

import numpy as np

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..io.read_data_other import read_centroids

try:
    import sklearn_extra
    import sklearn
    _SKLEARN_AVAILABLE = True
except ImportError:
    _SKLEARN_AVAILABLE = False


def process_echo_id(procstatus, dscfg, radar_list=None):
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
        if datatype == 'dBZ':
            refl_field = get_fieldname_pyart(datatype)
        if datatype == 'dBuZ':
            refl_field = get_fieldname_pyart(datatype)
        if datatype == 'ZDR':
            zdr_field = get_fieldname_pyart(datatype)
        if datatype == 'ZDRu':
            zdr_field = get_fieldname_pyart(datatype)
        if datatype == 'RhoHV':
            rhv_field = get_fieldname_pyart(datatype)
        if datatype == 'uPhiDP':
            phi_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if ((refl_field not in radar.fields) or
            (zdr_field not in radar.fields) or
            (rhv_field not in radar.fields) or
            (phi_field not in radar.fields)):
        warn('Unable to create radar_echo_id dataset. Missing data')
        return None, None

    echo_id = np.ma.zeros((radar.nrays, radar.ngates), dtype=np.uint8)+3

    # look for clutter
    gatefilter = pyart.filters.moment_and_texture_based_gate_filter(
        radar, zdr_field=zdr_field, rhv_field=rhv_field, phi_field=phi_field,
        refl_field=refl_field, textzdr_field=None, textrhv_field=None,
        textphi_field=None, textrefl_field=None, wind_size=7,
        max_textphi=20., max_textrhv=0.3, max_textzdr=2.85,
        max_textrefl=8., min_rhv=0.6)

    is_clutter = gatefilter.gate_excluded == 1
    echo_id[is_clutter] = 2

    # look for noise
    is_noise = radar.fields[refl_field]['data'].data == (
        pyart.config.get_fillvalue())
    echo_id[is_noise] = 1

    id_field = pyart.config.get_metadata('radar_echo_id')
    id_field['data'] = echo_id
    id_field.update({'_FillValue': 0})

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field('radar_echo_id', id_field)

    return new_dataset, ind_rad


def process_birds_id(procstatus, dscfg, radar_list=None):
    """
    identifies echoes as 0: No data, 1: Noise, 2: Clutter,
    3: Birds

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
        if datatype == 'dBZ':
            refl_field = get_fieldname_pyart(datatype)
        if datatype == 'dBuZ':
            refl_field = get_fieldname_pyart(datatype)
        if datatype == 'ZDR':
            zdr_field = get_fieldname_pyart(datatype)
        if datatype == 'ZDRu':
            zdr_field = get_fieldname_pyart(datatype)
        if datatype == 'RhoHV':
            rhv_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if ((refl_field not in radar.fields) or
            (zdr_field not in radar.fields) or
            (rhv_field not in radar.fields)):
        warn('Unable to create radar_echo_id dataset. Missing data')
        return None, None

    # user defined parameters
    max_zdr = dscfg.get('max_zdr', 3.)
    max_rhv = dscfg.get('max_rhv', 0.9)
    max_refl = dscfg.get('max_refl', 20.)
    rmin = dscfg.get('rmin', 2000.)
    rmax = dscfg.get('rmax', 25000.)
    elmin = dscfg.get('elmin', 1.5)
    elmax = dscfg.get('elmax', 85.)
    echo_id = np.zeros((radar.nrays, radar.ngates), dtype=np.uint8)+3

    # look for clutter
    gatefilter = pyart.filters.birds_gate_filter(
        radar, zdr_field=zdr_field, rhv_field=rhv_field,
        refl_field=refl_field, max_zdr=max_zdr, max_rhv=max_rhv,
        max_refl=max_refl, rmin=rmin, rmax=rmax, elmin=elmin, elmax=elmax)

    is_clutter = gatefilter.gate_excluded == 1
    echo_id[is_clutter] = 2

    # look for noise
    is_noise = radar.fields[refl_field]['data'].data == (
        pyart.config.get_fillvalue())
    echo_id[is_noise] = 1

    id_field = pyart.config.get_metadata('radar_echo_id')
    id_field['data'] = echo_id

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field('radar_echo_id', id_field)

    return new_dataset, ind_rad


def process_clt_to_echo_id(procstatus, dscfg, radar_list=None):
    """
    Converts clutter exit code from rad4alp into pyrad echo ID

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
        if datatype == 'CLT':
            clt_field = get_fieldname_pyart(datatype)
            break

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if clt_field not in radar.fields:
        warn('rad4alp clutter exit code not present. Unable to obtain echoID')
        return None, None

    echo_id = np.zeros((radar.nrays, radar.ngates), dtype=np.uint8)+3
    clt = radar.fields[clt_field]['data']
    echo_id[clt == 1] = 1
    echo_id[clt >= 100] = 2

    id_field = pyart.config.get_metadata('radar_echo_id')
    id_field['data'] = echo_id

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field('radar_echo_id', id_field)

    return new_dataset, ind_rad


def process_hydro_mf_to_echo_id(procstatus, dscfg, radar_list=None):
    """
    Converts MF hydrometeor classification into pyrad echo ID

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
        if datatype == 'hydroMF':
            clt_field = get_fieldname_pyart(datatype)
            break

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if clt_field not in radar.fields:
        warn('MF hydrometeor classification not present.'
             ' Unable to obtain echoID')
        return None, None

    echo_id = np.zeros((radar.nrays, radar.ngates), dtype=np.uint8)+1
    clt = radar.fields[clt_field]['data']
    echo_id[clt >= 8] = 3  # precip
    echo_id[np.logical_and(clt < 8, clt > 1)] = 2  # clt
    echo_id[np.ma.getmaskarray(clt)] = 1  # noise
    echo_id[clt == 1] = 1  # missing Zh

    id_field = pyart.config.get_metadata('radar_echo_id')
    id_field['data'] = echo_id

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field('radar_echo_id', id_field)

    return new_dataset, ind_rad


def process_hydro_mf_to_hydro(procstatus, dscfg, radar_list=None):
    """
    Converts the hydrometeor classification from Météo France to
    that of MeteoSwiss

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
        if datatype == 'hydroMF':
            field = get_fieldname_pyart(datatype)
            break

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if field not in radar.fields:
        warn('hydroMF not present. Unable to obtain hydro')
        return None, None

    hydro = np.zeros((radar.nrays, radar.ngates), dtype=np.uint8)
    hydroMF = radar.fields[field]['data']

    # BRUIT, ZH_MQT, SOL, INSECTES, OISEAUX, MER_CHAFF, PARASITES,
    # ROND_CENTRAL, TYPE_INCONNU, SIMPLE_POLAR are classified as NC
    hydro[hydroMF < 8] = 1
    hydro[hydroMF == 30] = 1
    hydro[hydroMF == 31] = 1
    # PRECIP_INDIFFERENCIEE, PLUIE, PRECIP are classified as RN
    hydro[hydroMF == 8] = 6
    hydro[hydroMF == 9] = 6
    hydro[hydroMF == 32] = 6
    hydro[hydroMF == 10] = 8  # NEIGE_MOUILLEE is WS
    hydro[hydroMF == 11] = 2  # NEIGE_SECHE is AG
    hydro[hydroMF == 12] = 3  # GLACE is CR
    hydro[hydroMF == 13] = 5  # PETITE_GRELE is RP
    # MOYENNE_GRELE, GROSSE_GRELE is IH/HDG
    hydro[hydroMF == 14] = 10
    hydro[hydroMF == 15] = 10
    # Light rain (LR), vertically oriented ice (VI) and melting hail (MH) have
    # no equivalent in the Météo France classification

    hydro_field = pyart.config.get_metadata('radar_echo_classification')
    hydro_field['data'] = hydro

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(
        'radar_echo_classification', hydro_field)

    return new_dataset, ind_rad


def process_echo_filter(procstatus, dscfg, radar_list=None):
    """
    Masks all echo types that are not of the class specified in
    keyword echo_type

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        echo_type : int or list of ints
            The type of echoes to keep: 1 noise, 2 clutter, 3 precipitation.
            Default 3
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

    echoid_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'echoID':
            echoid_field = get_fieldname_pyart(datatype)
            break
    if echoid_field is None:
        warn('echoID field required to filter data')
        return None, None

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if echoid_field not in radar.fields:
        warn('Unable to filter data. Missing echo ID field')
        return None, None

    echo_type = dscfg.get('echo_type', 3)
    mask = np.ma.isin(
        radar.fields[echoid_field]['data'], echo_type, invert=True)

    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'echoID':
            continue

        field_name = get_fieldname_pyart(datatype)
        if field_name not in radar.fields:
            warn('Unable to filter '+field_name+' according to echo ID. ' +
                 'No valid input fields')
            continue
        radar_field = deepcopy(radar.fields[field_name])
        radar_field['data'] = np.ma.masked_where(
            mask, radar_field['data'])

        if field_name.startswith('corrected_'):
            new_field_name = field_name
        elif field_name.startswith('uncorrected_'):
            new_field_name = field_name.replace(
                'uncorrected_', 'corrected_', 1)
        else:
            new_field_name = 'corrected_'+field_name
        new_dataset['radar_out'].add_field(new_field_name, radar_field)

    if not new_dataset['radar_out'].fields:
        return None, None

    return new_dataset, ind_rad


def process_cdf(procstatus, dscfg, radar_list=None):
    """
    Collects the fields necessary to compute the Cumulative Distribution
    Function

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

    echoid_field = None
    hydro_field = None
    vis_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'echoID':
            echoid_field = get_fieldname_pyart(datatype)
        elif datatype == 'hydro':
            hydro_field = get_fieldname_pyart(datatype)
        elif datatype == 'VIS':
            vis_field = get_fieldname_pyart(datatype)
        else:
            field_name = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if field_name not in radar.fields:
        warn('Unable to compute CDF. Missing field')
        return None, None

    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()

    new_dataset['radar_out'].add_field(field_name, radar.fields[field_name])
    if echoid_field is not None:
        if echoid_field not in radar.fields:
            warn('Missing echo ID field. Clutter can not be filtered')
        else:
            new_dataset['radar_out'].add_field(
                echoid_field, radar.fields[echoid_field])
    if hydro_field is not None:
        if hydro_field not in radar.fields:
            warn('Missing hydrometeor type field. ' +
                 'Filtration according to hydrometeor type not possible')
        else:
            new_dataset['radar_out'].add_field(
                hydro_field, radar.fields[hydro_field])
    if vis_field is not None:
        if vis_field not in radar.fields:
            warn('Missing visibility field. Blocked gates can not be filtered')
        else:
            new_dataset['radar_out'].add_field(
                vis_field, radar.fields[vis_field])

    return new_dataset, ind_rad


def process_filter_snr(procstatus, dscfg, radar_list=None):
    """
    filters out low SNR echoes

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        SNRmin : float. Dataset keyword
            The minimum SNR to keep the data.
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
        if datatype in ('SNRh', 'SNRv'):
            snr_field = get_fieldname_pyart(datatype)
            break

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()

    if snr_field not in radar.fields:
        warn('Unable to filter dataset according to SNR. Missing SNR field')
        return None, None

    gatefilter = pyart.filters.snr_based_gate_filter(
        radar, snr_field=snr_field, min_snr=dscfg['SNRmin'])
    is_low_snr = gatefilter.gate_excluded == 1

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)

        if datatype in ('SNRh', 'SNRv'):
            continue

        field_name = get_fieldname_pyart(datatype)
        if field_name not in radar.fields:
            warn('Unable to filter '+field_name +
                 ' according to SNR. '+'No valid input fields')
            continue

        radar_field = deepcopy(radar.fields[field_name])
        radar_field['data'] = np.ma.masked_where(
            is_low_snr, radar_field['data'])

        if field_name.startswith('corrected_'):
            new_field_name = field_name
        elif field_name.startswith('uncorrected_'):
            new_field_name = field_name.replace(
                'uncorrected_', 'corrected_', 1)
        else:
            new_field_name = 'corrected_'+field_name
        new_dataset['radar_out'].add_field(new_field_name, radar_field)

    if not new_dataset['radar_out'].fields:
        return None, None

    return new_dataset, ind_rad


def process_filter_vel_diff(procstatus, dscfg, radar_list=None):
    """
    filters out range gates that could not be used for Doppler velocity
    estimation

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
        if datatype == 'diffV':
            vel_diff_field = get_fieldname_pyart(datatype)
            break

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()

    if vel_diff_field not in radar.fields:
        warn('Unable to filter dataset according to valid velocity. ' +
             'Missing velocity differences field')
        return None, None

    mask = np.ma.getmaskarray(radar.fields[vel_diff_field]['data'])

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)

        if datatype == 'diffV':
            continue

        field_name = get_fieldname_pyart(datatype)
        if field_name not in radar.fields:
            warn('Unable to filter '+field_name +
                 ' according to SNR. '+'No valid input fields')
            continue

        radar_field = deepcopy(radar.fields[field_name])
        radar_field['data'] = np.ma.masked_where(mask, radar_field['data'])

        if field_name.find('corrected_') != -1:
            new_field_name = field_name
        elif field_name.startswith('uncorrected_'):
            new_field_name = field_name.replace(
                'uncorrected_', 'corrected_', 1)
        else:
            new_field_name = 'corrected_'+field_name
        new_dataset['radar_out'].add_field(new_field_name, radar_field)

    if not new_dataset['radar_out'].fields:
        return None, None

    return new_dataset, ind_rad


def process_filter_visibility(procstatus, dscfg, radar_list=None):
    """
    filters out rays gates with low visibility and corrects the reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        VISmin : float. Dataset keyword
            The minimum visibility to keep the data.
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
        if datatype == 'VIS':
            vis_field = get_fieldname_pyart(datatype)
            break

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()

    if vis_field not in radar.fields:
        warn('Unable to filter dataset according to visibility. ' +
             'Missing visibility field')
        return None, None

    gatefilter = pyart.filters.visibility_based_gate_filter(
        radar, vis_field=vis_field, min_vis=dscfg['VISmin'])
    is_lowVIS = gatefilter.gate_excluded == 1

    for datatypedescr in dscfg['datatype']:
        _, _, datatype, _, _ = get_datatype_fields(
            datatypedescr)

        if datatype == 'VIS':
            continue
        field_name = get_fieldname_pyart(datatype)
        if field_name not in radar.fields:
            warn('Unable to filter '+field_name +
                 ' according to visibility. No valid input fields')
            continue

        radar_aux = deepcopy(radar)
        radar_aux.fields[field_name]['data'] = np.ma.masked_where(
            is_lowVIS, radar_aux.fields[field_name]['data'])

        if datatype in ('dBZ', 'dBZc', 'dBuZ', 'dBZv', 'dBZvc', 'dBuZv'):
            radar_field = pyart.correct.correct_visibility(
                radar_aux, vis_field=vis_field, field_name=field_name)
        else:
            radar_field = radar_aux.fields[field_name]

        if field_name.startswith('corrected_'):
            new_field_name = field_name
        elif field_name.startswith('uncorrected_'):
            new_field_name = field_name.replace(
                'uncorrected_', 'corrected_', 1)
        else:
            new_field_name = 'corrected_'+field_name
        new_dataset['radar_out'].add_field(new_field_name, radar_field)

    if not new_dataset['radar_out'].fields:
        return None, None

    return new_dataset, ind_rad


def process_outlier_filter(procstatus, dscfg, radar_list=None):
    """
    filters out gates which are outliers respect to the surrounding

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        threshold : float. Dataset keyword
            The distance between the value of the examined range gate and the
            median of the surrounding gates to consider the gate an outlier
        nb : int. Dataset keyword
            The number of neighbours (to one side) to analyse. i.e. 2 would
            correspond to 24 gates
        nb_min : int. Dataset keyword
            Minimum number of neighbouring gates to consider the examined gate
            valid
        percentile_min, percentile_max : float. Dataset keyword
            gates below (above) these percentiles (computed over the sweep) are
            considered potential outliers and further examined
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

    radarnr, _, datatype, _, _ = get_datatype_fields(
        dscfg['datatype'][0])

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    field_name = get_fieldname_pyart(datatype)
    if field_name not in radar.fields:
        warn('Unable to perform outlier removal. No valid data')
        return None, None

    threshold = dscfg.get('threshold', 10.)
    nb = dscfg.get('nb', 2)
    nb_min = dscfg.get('nb_min', 3)
    percentile_min = dscfg.get('percentile_min', 5.)
    percentile_max = dscfg.get('percentile_max', 95.)

    field = radar.fields[field_name]
    field_out = deepcopy(field)
    for sweep in range(radar.nsweeps):
        # find gates suspected to be outliers
        sweep_start = radar.sweep_start_ray_index['data'][sweep]
        sweep_end = radar.sweep_end_ray_index['data'][sweep]
        nrays_sweep = radar.rays_per_sweep['data'][sweep]
        data_sweep = field['data'][sweep_start:sweep_end+1, :]

        # check if all elements in array are masked
        if np.all(np.ma.getmaskarray(data_sweep)):
            continue

        percent_vals = np.nanpercentile(
            data_sweep.filled(fill_value=np.nan),
            (percentile_min, percentile_max))
        ind_rays, ind_rngs = np.ma.where(
            np.ma.logical_or(
                data_sweep < percent_vals[0], data_sweep > percent_vals[1]))

        for i, ind_ray in enumerate(ind_rays):
            ind_rng = ind_rngs[i]
            # find neighbours of suspected outlier gate
            data_cube = []
            for ray_nb in range(-nb, nb+1):
                for rng_nb in range(-nb, nb+1):
                    if ray_nb == 0 and rng_nb == 0:
                        continue
                    if ((ind_ray+ray_nb >= 0) and
                            (ind_ray+ray_nb < nrays_sweep) and
                            (ind_rng+rng_nb >= 0) and
                            (ind_rng+rng_nb < radar.ngates)):
                        if (data_sweep[ind_ray+ray_nb, ind_rng+rng_nb] is not
                                np.ma.masked):
                            data_cube.append(
                                data_sweep[ind_ray+ray_nb, ind_rng+rng_nb])

            # remove data far from median of neighbours or with not enough
            # valid neighbours
            if len(data_cube) < nb_min:
                field_out['data'][
                    sweep_start+ind_ray, ind_rng] = np.ma.masked
            elif (abs(np.ma.median(data_cube) -
                      data_sweep[ind_ray, ind_rng]) > threshold):
                field_out['data'][sweep_start+ind_ray, ind_rng] = np.ma.masked

    if field_name.startswith('corrected_'):
        new_field_name = field_name
    elif field_name.startswith('uncorrected_'):
        new_field_name = field_name.replace(
            'uncorrected_', 'corrected_', 1)
    else:
        new_field_name = 'corrected_'+field_name

    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(new_field_name, field_out)

    return new_dataset, ind_rad


def process_hydroclass(procstatus, dscfg, radar_list=None):
    """
    Classifies precipitation echoes

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        HYDRO_METHOD : string. Dataset keyword
            The hydrometeor classification method. One of the following:
            SEMISUPERVISED, UKMO
        centroids_file : string or None. Dataset keyword
            Used with HYDRO_METHOD SEMISUPERVISED. The name of the .csv file
            that stores the centroids. The path is given by
            [configpath]/centroids_hydroclass/
            If None is provided default centroids are going to be used
        compute_entropy : bool. Dataset keyword
            Used with HYDRO_METHOD SEMISUPERVISED. If true the entropy is
            computed and the field hydroclass_entropy is output
        output_distances : bool. Dataset keyword
            Used with HYDRO_METHOD SEMISUPERVISED. If true the de-mixing
            algorithm based on the distances to the centroids is computed and
            the field proportions of each hydrometeor in the radar range gate
            is output
        vectorize : bool. Dataset keyword
            Used with HYDRO_METHOD SEMISUPERVISED. If true a vectorized
            version of the algorithm is used
        weights : array of floats. Dataset keyword
            Used with HYDRO_METHOD SEMISUPERVISED. The list of weights given
            to each variable
        hydropath : string. Dataset keyword
            Used with HYDRO_METHOD UKMO. Directory of the UK MetOffice
            hydrometeor classification code
        mf_dir : string. Dataset keyword
            Used with HYDRO_METHOD UKMO. Directory where the UK MetOffice
            hydrometeor classification membership functions are stored
        ml_depth: float. Dataset keyword
            Used with HYDRO_METHOD UKMO. Depth of the melting layer [km].
            Default 500.
        perturb_ml_depth: float. Dataset keyword
            Used with HYDRO_METHOD UKMO. if specified, the depth of the
            melting layer can be varied by +/- this value [km], allowing a
            less-rigidly defined melting layer. Default 0.
        freezing_level: float or None. Dataset keyword
            Used with HYDRO_METHOD UKMO. if desired, a single freezing level
            height can be specified for the entire PPI domain - this will
            over-ride any field found within the input file. Default None
        use_dualpol: Bool. Dataset keyword
            Used with HYDRO_METHOD UKMO. If false no radar data is used and
            the classification is performed using temperature information
            only. Default True
        use_temperature: Bool. Dataset keyword
            Used with HYDRO_METHOD UKMO. If false no temperature information
            is used and the classification is performed using radar data only.
            Default True
        use_interpolation: Bool. Dataset keyword
            Used with HYDRO_METHOD UKMO. If True gaps in the classification
            are filled using a nearest-neighbour interpolation. Default False
        map_to_semisupervised: Bool. Dataset keyword
            Used with HYDRO_METHOD UKMO. If True the output is map to the same
            categories as the semi-supervised classification. Default True
        append_all_fields: Bool. Dataset keyword
            Used with HYDRO_METHOD UKMO. If True auxiliary fields such as
            confidence and probability for each class are going to be added to
            the output

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

    if 'HYDRO_METHOD' not in dscfg:
        raise Exception(
            "ERROR: Undefined parameter 'HYDRO_METHOD' for dataset {}".format(
                dscfg['dsname']))

    temp_field = None
    iso0_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'ZDR':
            zdr_field = 'differential_reflectivity'
        if datatype == 'ZDRc':
            zdr_field = 'corrected_differential_reflectivity'
        if datatype == 'RhoHV':
            rhv_field = 'cross_correlation_ratio'
        if datatype == 'uRhoHV':
            rhv_field = 'uncorrected_cross_correlation_ratio'
        if datatype == 'RhoHVc':
            rhv_field = 'corrected_cross_correlation_ratio'
        if datatype == 'KDP':
            kdp_field = 'specific_differential_phase'
        if datatype == 'KDPc':
            kdp_field = 'corrected_specific_differential_phase'
        if datatype == 'TEMP':
            temp_field = 'temperature'
        if datatype == 'H_ISO0':
            iso0_field = 'height_over_iso0'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if dscfg['HYDRO_METHOD'] == 'SEMISUPERVISED':
        if temp_field is None and iso0_field is None:
            warn('iso0 or temperature fields needed to create hydrometeor ' +
                 'classification field')
            return None, None

        if temp_field is not None and (temp_field not in radar.fields):
            warn('Unable to create hydrometeor classification field. ' +
                 'Missing temperature field')
            return None, None

        if iso0_field is not None and (iso0_field not in radar.fields):
            warn('Unable to create hydrometeor classification field. ' +
                 'Missing height over iso0 field')
            return None, None

        temp_ref = 'temperature'
        if iso0_field is not None:
            temp_ref = 'height_over_iso0'

        if ((refl_field not in radar.fields) or
                (zdr_field not in radar.fields) or
                (rhv_field not in radar.fields) or
                (kdp_field not in radar.fields)):
            warn('Unable to create hydrometeor classification field. ' +
                 'Missing data')
            return None, None

        # user defined parameters
        centroids_file = dscfg.get('centroids_file', None)
        compute_entropy = dscfg.get('compute_entropy', False)
        output_distances = dscfg.get('output_distances', False)
        vectorize = dscfg.get('vectorize', False)
        weights = dscfg.get('weights', np.array((1., 1., 1., 0.75, 0.5)))

        # load centroids
        if centroids_file is not None:
            mass_centers, hydro_names, var_names = read_centroids(
                dscfg['configpath']+'centroids_hydroclass/'+centroids_file)
            if mass_centers is None:
                warn(
                    'Unable to read centroids file. ' +
                    'Default centroids will be used in classification.')
                hydro_names = (
                    'AG', 'CR', 'LR', 'RP', 'RN', 'VI', 'WS', 'MH', 'IH/HDG')
                var_names = ('dBZ', 'ZDR', 'KDP', 'RhoHV', 'H_ISO0')
        else:
            warn(
                'No centroids were specified. ' +
                'Default centroids will be used in classification.')
            mass_centers = None
            hydro_names = (
                'AG', 'CR', 'LR', 'RP', 'RN', 'VI', 'WS', 'MH', 'IH/HDG')
            var_names = ('dBZ', 'ZDR', 'KDP', 'RhoHV', 'H_ISO0')

        fields_dict = pyart.retrieve.hydroclass_semisupervised(
            radar, hydro_names=hydro_names, var_names=var_names,
            mass_centers=mass_centers, weights=weights, refl_field=refl_field,
            zdr_field=zdr_field, rhv_field=rhv_field, kdp_field=kdp_field,
            temp_field=temp_field, iso0_field=iso0_field, hydro_field=None,
            entropy_field=None, temp_ref=temp_ref,
            compute_entropy=compute_entropy,
            output_distances=output_distances, vectorize=vectorize)

        # prepare for exit
        new_dataset = {'radar_out': deepcopy(radar)}
        new_dataset['radar_out'].fields = dict()
        new_dataset['radar_out'].add_field(
            'radar_echo_classification', fields_dict['hydro'])

        if compute_entropy:
            new_dataset['radar_out'].add_field(
                'hydroclass_entropy', fields_dict['entropy'])
            if output_distances:
                for hydro_name in hydro_names:
                    field_name = 'proportion_'+hydro_name
                    new_dataset['radar_out'].add_field(
                        field_name, fields_dict[field_name])

    elif dscfg['HYDRO_METHOD'] == 'UKMO':
        if dscfg['initialized'] == 0:
            # load the UKMO algorithm
            hydropath = dscfg.get(
                'hydropath',
                os.path.expanduser('~')+'/hydrometeor-classification/code')
            if not os.path.isdir(hydropath):
                str1 = "ERROR: Wrong UKMO hydrometeor classification path {}"
                raise Exception(str1.format(hydropath))
            sys.path.append(hydropath)
            try:
                import classify
            except ImportError:
                raise Exception(
                    "ERROR: Unable to load UKMO hydrometeor classification"
                    " module")
            dscfg['initialized'] = 1

        if dscfg['initialized'] == 0:
            return None, None

        if ((refl_field not in radar.fields) or
                (zdr_field not in radar.fields) or
                (rhv_field not in radar.fields) or
                (kdp_field not in radar.fields)):
            warn('Unable to create hydrometeor classification field. '
                 'Missing data')
            return None, None

        mf_dir = dscfg.get(
            'mf_dir',
            os.path.expanduser('~') +
            '/hydrometeor-classification/membership_functions/')
        if not os.path.isdir(mf_dir):
            str1 = (
                'ERROR: Unable to load hydrometeor MF.'
                ' Path {} not a directory')
            raise Exception(str1.format(mf_dir))

        ml_depth = dscfg.get('ml_depth', 0.5)
        perturb_ml_depth = dscfg.get('perturb_ml_depth', 0)
        freezing_level = dscfg.get('freezing_level', None)

        use_dualpol = dscfg.get('use_dualpol', True)
        use_temperature = dscfg.get('use_temperature', True)
        use_interpolation = dscfg.get('use_interpolation', False)
        map_to_semisupervised = dscfg.get('map_to_semisupervised', True)
        append_all_fields = dscfg.get('append_all_fields', True)

        dualpol_vars = (rhv_field, zdr_field, kdp_field)
        dualpol_vars_int = ('RHOHV', 'ZDR', 'KDP')

        # prepare for exit
        hydro_field = 'radar_echo_classification'
        new_dataset = {'radar_out': deepcopy(radar)}
        new_dataset['radar_out'].fields = dict()
        hydro = pyart.config.get_metadata(hydro_field)
        hydro['data'] = np.ma.masked_all(
            (radar.nrays, radar.ngates), dtype=np.uint8)
        new_dataset['radar_out'].add_field(hydro_field, hydro)

        if append_all_fields:
            confidence_field = 'hydroclass_confidence'
            confidence = pyart.config.get_metadata(confidence_field)
            confidence['data'] = np.ma.masked_all((radar.nrays, radar.ngates))
            new_dataset['radar_out'].add_field(confidence_field, confidence)

            prob_fields = [
                'probability_RN', 'probability_WS', 'probability_AG',
                'probability_IH', 'probability_LR', 'probability_IC',
                'probability_RP']
            prob_keys = [
                    'RAIN', 'WET_SNOW', 'DRY_SNOW', 'HAIL', 'DRIZZLE', 'ICE',
                    'GRAUPEL']
            for prob_field in prob_fields:
                prob = pyart.config.get_metadata(prob_field)
                prob['data'] = np.ma.masked_all((radar.nrays, radar.ngates))
                new_dataset['radar_out'].add_field(prob_field, prob)

        for sweep in range(radar.nsweeps):
            scan_object = ScanObject(
                radar.extract_sweeps([sweep]), refl_field, iso0_field,
                dualpol_vars, dualpol_vars_int=dualpol_vars_int)

            hydro_object = classify.classify_hydrometeors(
                scan_object, ml_depth, mf_dir, use_dualpol=use_dualpol,
                use_temperature=use_temperature, freezing_level=freezing_level,
                perturb_ml_depth=perturb_ml_depth,
                dualpol_vars=dualpol_vars_int,
                use_interpolation=use_interpolation,
                append_all_fields=append_all_fields)

            ind_start = radar.sweep_start_ray_index['data'][sweep]
            ind_end = radar.sweep_end_ray_index['data'][sweep]

            if map_to_semisupervised:
                # re-mapping into semi-supervised classes
                hydro_ukmo_aux = hydro_object.hydro_class
                hydro_ukmo = np.ma.masked_all(
                    hydro_ukmo_aux.shape, dtype=np.uint8)
                hydro_ukmo[hydro_ukmo_aux == -1] = 1  # No data to NC
                hydro_ukmo[hydro_ukmo_aux == 1] = 6  # RAIN to RN
                hydro_ukmo[hydro_ukmo_aux == 2] = 8  # WET_SNOW to WS
                hydro_ukmo[hydro_ukmo_aux == 3] = 2  # DRY_SNOW to AG
                hydro_ukmo[hydro_ukmo_aux == 5] = 10  # HAIL to IH/HDG
                hydro_ukmo[hydro_ukmo_aux == 6] = 4  # DRIZZLE to LR
                hydro_ukmo[hydro_ukmo_aux == 7] = 3  # ICE to CR
                hydro_ukmo[hydro_ukmo_aux == 8] = 5  # GRAUPEL to RP
            else:
                hydro_ukmo = hydro_object.hydro_class
            new_dataset['radar_out'].fields[hydro_field]['data'][
                ind_start:ind_end+1] = hydro_ukmo

            if append_all_fields:
                new_dataset['radar_out'].fields[confidence_field]['data'][
                    ind_start:ind_end+1] = hydro_object.confidence

                for prob_field, prob_key in zip(prob_fields, prob_keys):
                    new_dataset['radar_out'].fields[prob_field]['data'][
                        ind_start:ind_end+1] = (
                        hydro_object.probability[prob_key])
    else:
        raise Exception(
            "ERROR: Unknown hydrometeor classification method {}".format(
                dscfg['HYDRO_METHOD']))

    return new_dataset, ind_rad


def process_centroids(procstatus, dscfg, radar_list=None):
    """
    Computes centroids for the semi-supervised hydrometeor classification

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        samples_per_vol : int. Dataset keyword
            Maximum number of samples per volume kept for further analysis.
            Default 20000
        nbins : int.
            Number of bins of the histogram used to make the data platykurtic.
            Default 110
        pdf_zh_max : int
            Multiplicative factor to the Guassian function used to make the
            distribution of the reflectivity platykurtic that determines the
            number of samples for each bin. Default 10000
        pdf_relh_max : int
            Multiplicative factor to the Guassian function used to make the
            distribution of the height relative to the iso-0 platykurtic that
            determines the number of samples for each bin. Default 20000
        sigma_zh, sigma_relh : float
            sigma of the respective Gaussian functions. Defaults 0.75 and 1.5
        randomize : bool
            If True the data is randomized to avoid the effects of the
            quantization. Default True
        platykurtic_dBZ : bool
            If True makes the reflectivity distribution platykurtic. Default
            True
        platykurtic_H_ISO0 : bool
            If True makes the height respect to the iso-0 distribution
            platykurtic. Default True
        relh_slope : float. Dataset keyword
            The slope used to transform the height relative to the iso0 into
            a sigmoid function. Default 0.001
        external_iterations : int. Dataset keywords
            Number of iterations of the external loop. This number will
            determine how many medoids are computed for each hydrometeor
            class. Default 30
        internal_iterations : int. Dataset keyword
            Maximum number of iterations of the internal loop. Default 10
        sample_data : Bool.
            If True the data is going to be sampled prior to each external
            iteration. Default False
        nsamples_iter : int.
            Number of samples per iteration. Default 20000
        alpha : float
            Minimum value to accept the cluster according to p. Default 0.01
        cv_approach : Bool
            If true it is used a critical value approach to reject or accept
            similarity between observations and reference. If false it is used
            a p-value approach. Default True
        n_samples_syn : int
            Number of samples drawn from reference to compare it with
            observations in the KS test. Default 50
        num_samples_arr : array of int
            Number of observation samples used in the KS test to choose from.
            Default (30, 35, 40)
        acceptance_threshold : float. Dataset keyword
            Threshold on the inter-quantile coefficient of dispersion of the
            medoids above which the medoid of the class is not acceptable.
            Default 0.5
        nmedoids_min : int
            Minimum number of intermediate medoids to compute the final
            result. Default 1
        var_names : tupple
            The names of the features. Default ('dBZ', 'ZDR', 'KDP', 'RhoHV',
            'H_ISO0')
        hydro_names: tupple
            The name of the hydrometeor types. Default ('AG', 'CR', 'LR',
            'RP', 'RN', 'VI', 'WS', 'MH', 'IH/HDG')
        weight : tupple
            The weight given to each feature when comparing to the reference.
            It is in the same order as var_names. Default (1., 1., 1., 1.,
            0.75)
        parallelized : bool
            If True the centroids search is going to be parallelized. Default
            False
        kmax_iter : int
            Maximum number of iterations of the k-medoids algorithm. Default
            100
        nsamples_small : int
            Maximum number before using the k-medoids CLARA algorithm. If this
            number is exceeded the CLARA algorithm will be used. Default 40000
        sampling_size_clara : int
            Number of samples used in each iteration of the k-medoids CLARA
            algorithm. Default 10000
        niter_clara : int
            Number of iterations performed by the k-medoids CLARA algorithm.
            Default 5
        keep_labeled_data : bool
            If True the labeled data is going to be kept for storage. Default
            True
        use_median : bool
            If True the intermediate centroids are computed as the median
            of the observation variables and the final centroids are computed
            as the median of the intermediate centroids. If false they are
            computed using the kmedoids algorithm. Default false
        allow_label_duplicates : bool
            If True allow to label multiple clusters with the same label.
            Default True

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

    temp_field = None
    iso0_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'ZDR':
            zdr_field = 'differential_reflectivity'
        if datatype == 'ZDRc':
            zdr_field = 'corrected_differential_reflectivity'
        if datatype == 'RhoHV':
            rhv_field = 'cross_correlation_ratio'
        if datatype == 'uRhoHV':
            rhv_field = 'uncorrected_cross_correlation_ratio'
        if datatype == 'RhoHVc':
            rhv_field = 'corrected_cross_correlation_ratio'
        if datatype == 'KDP':
            kdp_field = 'specific_differential_phase'
        if datatype == 'KDPc':
            kdp_field = 'corrected_specific_differential_phase'
        if datatype == 'TEMP':
            temp_field = 'temperature'
        if datatype == 'H_ISO0':
            iso0_field = 'height_over_iso0'

    ind_rad = int(radarnr[5:8])-1

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if temp_field is None and iso0_field is None:
            warn('iso0 or temperature fields needed to compute centroids')
            return None, None

        if temp_field is not None and (temp_field not in radar.fields):
            warn('temperature field needed to compute centroids')
            return None, None

        if iso0_field is not None and (iso0_field not in radar.fields):
            warn('iso0 field needed to compute centroids')
            return None, None

        temp_ref = 'temperature'
        if iso0_field is not None:
            temp_ref = 'height_over_iso0'

        if ((refl_field not in radar.fields) or
                (zdr_field not in radar.fields) or
                (rhv_field not in radar.fields) or
                (kdp_field not in radar.fields)):
            warn('Unable to compute centroids. Missing data')
            return None, None

        samples_per_vol = dscfg.get('samples_per_vol', 20000)

        (refl, zdr, kdp, rhohv, relh) = pyart.retrieve.data_for_centroids(
            radar, refl_field=refl_field, zdr_field=zdr_field,
            rhv_field=rhv_field, kdp_field=kdp_field, temp_field=temp_field,
            iso0_field=iso0_field, temp_ref=temp_ref,
            nsamples_max=samples_per_vol)

        # first volume: initialize global data
        if dscfg['initialized'] == 0:
            if radar.instrument_parameters is None:
                warn('Unknown radar frequency. C-band assumed')
                band = 'C'
            elif 'frequency' not in radar.instrument_parameters.keys():
                warn('Unknown radar frequency. C-band assumed')
                band = 'C'
            else:
                band = pyart.retrieve.get_freq_band(
                    radar.instrument_parameters['frequency']['data'][0])
                print('Radar frequency band: {}'.format(band))

            data_dict = {
                'dBZ': refl,
                'ZDR': zdr,
                'KDP': kdp,
                'RhoHV': rhohv,
                'H_ISO0': relh,
                'npoints': np.array([refl.size], dtype=int),
                'timeinfo': np.array([dscfg['timeinfo']]),
                'final': False,
                'band': band
            }
            dscfg['global_data'] = data_dict
            dscfg['initialized'] = 1

            new_dataset = {'data_dict': data_dict}
            return new_dataset, ind_rad

        if dscfg['initialized'] == 0:
            return None, None

        dscfg['global_data']['dBZ'] = np.ma.append(
            dscfg['global_data']['dBZ'], refl)
        dscfg['global_data']['ZDR'] = np.ma.append(
            dscfg['global_data']['ZDR'], zdr)
        dscfg['global_data']['KDP'] = np.ma.append(
            dscfg['global_data']['KDP'], kdp)
        dscfg['global_data']['RhoHV'] = np.ma.append(
            dscfg['global_data']['RhoHV'], rhohv)
        dscfg['global_data']['H_ISO0'] = np.ma.append(
            dscfg['global_data']['H_ISO0'], relh)
        dscfg['global_data']['npoints'] = np.ma.append(
            dscfg['global_data']['npoints'], refl.size)
        dscfg['global_data']['timeinfo'] = np.ma.append(
            dscfg['global_data']['timeinfo'], dscfg['timeinfo'])

        new_dataset = {'data_dict': dscfg['global_data']}
        return new_dataset, ind_rad

    if dscfg['initialized'] == 0:
        warn('Unable to compute centroids. No valid data found')
        return None, None

    dscfg['global_data']['final'] = True
    new_dataset = {'data_dict': dscfg['global_data']}

    if not _SKLEARN_AVAILABLE:
        warn(
            'Unable to compute centroids. scikit-learn package not available')
        return new_dataset, ind_rad

    # select data to be used to determine centroids
    nbins = dscfg.get('nbins', 110)
    pdf_zh_max = dscfg.get('pdf_zh_max', 10000)
    pdf_relh_max = dscfg.get('pdf_relh_max', 20000)
    sigma_zh = dscfg.get('sigma_zh', 0.75)
    sigma_relh = dscfg.get('sigma_relh', 1.5)
    randomize = dscfg.get('randomize', True)
    platykurtic_dBZ = dscfg.get('platykurtic_dBZ', True)
    platykurtic_H_ISO0 = dscfg.get('platykurtic_H_ISO0', True)

    fm = np.transpose(np.array([
        dscfg['global_data']['dBZ'],
        dscfg['global_data']['ZDR'],
        dscfg['global_data']['KDP'],
        dscfg['global_data']['RhoHV'],
        dscfg['global_data']['H_ISO0']], dtype=np.float32))

    fm_sample = pyart.retrieve.select_samples(
        fm, np.random.default_rng(seed=0), nbins=nbins,
        pdf_zh_max=pdf_zh_max, pdf_relh_max=pdf_relh_max,
        sigma_zh=sigma_zh, sigma_relh=sigma_relh, randomize=randomize,
        platykurtic_dBZ=platykurtic_dBZ,
        platykurtic_H_ISO0=platykurtic_H_ISO0)

    dscfg['global_data']['dBZ'] = fm_sample[:, 0]
    dscfg['global_data']['ZDR'] = fm_sample[:, 1]
    dscfg['global_data']['KDP'] = fm_sample[:, 2]
    dscfg['global_data']['RhoHV'] = fm_sample[:, 3]
    dscfg['global_data']['H_ISO0'] = fm_sample[:, 4]

    # user selected parameters
    relh_slope = dscfg.get('relh_slope', 0.001)
    external_iterations = dscfg.get('external_iterations', 30)
    internal_iterations = dscfg.get('internal_iterations', 10)
    sample_data = dscfg.get('sample_data', False)
    nsamples_iter = dscfg.get('nsamples_iter', 20000)
    alpha = dscfg.get('alpha', 0.01)
    cv_approach = dscfg.get('cv_approach', True)
    n_samples_syn = dscfg.get('nsamples_syn', 50)
    num_samples_arr = dscfg.get('num_samples_arr', (30, 35, 40))
    acceptance_threshold = dscfg.get('acceptance_threshold', 0.5)
    nmedoids_min = dscfg.get('nmedoids_min', 1)
    var_names = dscfg.get(
        'var_names', ('dBZ', 'ZDR', 'KDP', 'RhoHV', 'H_ISO0'))
    hydro_names = dscfg.get(
        'hydro_names',
        ('AG', 'CR', 'LR', 'RP', 'RN', 'VI', 'WS', 'MH', 'IH/HDG'))
    weight = dscfg.get('weight', (1., 1., 1., 1., 1.))
    parallelized = dscfg.get('parallelized', False)
    kmax_iter = dscfg.get('kmax_iter', 100)
    nsamples_small = dscfg.get('nsamples_small', 40000)
    sampling_size_clara = dscfg.get('sampling_size_clara', 10000)
    niter_clara = dscfg.get('niter_clara', 5)
    keep_labeled_data = dscfg.get('keep_labeled_data', True)
    use_median = dscfg.get('use_median', True)
    allow_label_duplicates = dscfg.get('allow_label_duplicates', True)

    (labeled_data, labels, medoids_dict,
     final_medoids_dict) = pyart.retrieve.compute_centroids(
        fm_sample, weight=weight, var_names=var_names,
        hydro_names=hydro_names, nsamples_iter=nsamples_iter,
        external_iterations=external_iterations,
        internal_iterations=internal_iterations, alpha=alpha,
        cv_approach=cv_approach, num_samples_arr=num_samples_arr,
        n_samples_syn=n_samples_syn, nmedoids_min=nmedoids_min,
        acceptance_threshold=acceptance_threshold,
        band=dscfg['global_data']['band'], relh_slope=relh_slope,
        parallelized=parallelized, sample_data=sample_data,
        kmax_iter=kmax_iter, nsamples_small=nsamples_small,
        sampling_size_clara=sampling_size_clara,
        niter_clara=niter_clara, keep_labeled_data=keep_labeled_data,
        use_median=use_median, allow_label_duplicates=allow_label_duplicates)

    if not medoids_dict:
        return new_dataset, ind_rad

    labeled_data_dict = {
        'hydro_names': hydro_names,
        'var_names': var_names,
        'band': dscfg['global_data']['band'],
        'medoids_dict': medoids_dict,
        'final_medoids_dict': final_medoids_dict,
        'timeinfo': dscfg['global_data']['timeinfo']
    }
    if keep_labeled_data:
        labeled_data_dict.update({
            'dBZ': labeled_data[:, 0],
            'ZDR': labeled_data[:, 1],
            'KDP': labeled_data[:, 2],
            'RhoHV': labeled_data[:, 3],
            'H_ISO0': labeled_data[:, 4],
            'labels': labels})
    new_dataset.update({'labeled_data_dict': labeled_data_dict})

    return new_dataset, ind_rad


def process_melting_layer(procstatus, dscfg, radar_list=None):
    """
    Detects the melting layer

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

    if 'ML_METHOD' not in dscfg:
        str1 = "ERROR: Undefined parameter 'ML_METHOD' for dataset {}"
        raise Exception(str1.format(dscfg['dsname']))

    if dscfg['ML_METHOD'] == 'GIANGRANDE':

        temp_ref = None
        temp_field = None
        iso0_field = None
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            if datatype == 'dBZ':
                refl_field = 'reflectivity'
            if datatype == 'dBZc':
                refl_field = 'corrected_reflectivity'
            if datatype == 'ZDR':
                zdr_field = 'differential_reflectivity'
            if datatype == 'ZDRc':
                zdr_field = 'corrected_differential_reflectivity'
            if datatype == 'RhoHV':
                rhv_field = 'cross_correlation_ratio'
            if datatype == 'RhoHVc':
                rhv_field = 'corrected_cross_correlation_ratio'
            if datatype == 'TEMP':
                temp_field = 'temperature'
            if datatype == 'H_ISO0':
                iso0_field = 'height_over_iso0'

        ind_rad = int(radarnr[5:8])-1
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        # Check which should be the reference field for temperature
        if iso0_field is not None:
            if iso0_field not in radar.fields:
                warn('Unable to detect melting layer. '
                     'Missing height over iso0 field')
                return None, None
            temp_ref = 'height_over_iso0'

        if temp_field is not None:
            if temp_field not in radar.fields:
                warn('Unable to detect melting layer. '
                     'Missing temperature field')
                return None, None
            temp_ref = 'temperature'
            iso0_field = 'height_over_iso0'

        if temp_ref is None:
            iso0_field = 'height_over_iso0'

        if ((refl_field not in radar.fields) or
                (zdr_field not in radar.fields) or
                (rhv_field not in radar.fields)):
            warn('Unable to detect melting layer. Missing data')
            return None, None

        # User defined variables
        nVol = dscfg.get('nVol', 3)
        maxh = dscfg.get('maxh', 6000.)
        hres = dscfg.get('hres', 50.)

        rmin = dscfg.get('rmin', 1000.)
        elmin = dscfg.get('elmin', 4.)
        elmax = dscfg.get('elmax', 10.)
        rhomin = dscfg.get('rhomin', 0.75)
        rhomax = dscfg.get('rhomax', 0.94)
        zhmin = dscfg.get('zhmin', 20.)
        hwindow = dscfg.get('hwindow', 500.)
        mlzhmin = dscfg.get('mlzhmin', 30.)
        mlzhmax = dscfg.get('mlzhmax', 50.)
        mlzdrmin = dscfg.get('mlzdrmin', 1.)
        mlzdrmax = dscfg.get('mlzdrmax', 5.)
        htol = dscfg.get('htol', 500.)
        ml_bottom_diff_max = dscfg.get('ml_bottom_diff_max', 1000.)

        time_accu_max = dscfg.get('time_accu_max', 1800.)
        nml_points_min = dscfg.get('nml_points_min', None)
        wlength = dscfg.get('wlength', 20.)
        percentile_bottom = dscfg.get('percentile_bottom', 0.3)
        percentile_top = dscfg.get('percentile_top', 0.9)
        interpol = dscfg.get('interpol', True)
        time_nodata_allowed = dscfg.get('time_nodata_allowed', 3600.)

        get_iso0 = dscfg.get('get_iso0', True)

        if not dscfg['initialized']:
            # initialize dataset
            ml_obj, ml_dict, iso0_dict, ml_global = (
                pyart.retrieve.melting_layer_giangrande(
                    radar, nVol=nVol, maxh=maxh, hres=hres, rmin=rmin,
                    elmin=elmin, elmax=elmax, rhomin=rhomin, rhomax=rhomax,
                    zhmin=zhmin, hwindow=hwindow, mlzhmin=mlzhmin,
                    mlzhmax=mlzhmax, mlzdrmin=mlzdrmin, mlzdrmax=mlzdrmax,
                    htol=htol, ml_bottom_diff_max=ml_bottom_diff_max,
                    time_accu_max=time_accu_max, nml_points_min=nml_points_min,
                    wlength=wlength, percentile_bottom=percentile_bottom,
                    percentile_top=percentile_top, interpol=interpol,
                    time_nodata_allowed=time_nodata_allowed,
                    refl_field=refl_field, zdr_field=zdr_field,
                    rhv_field=rhv_field, temp_field=temp_field,
                    iso0_field=iso0_field, ml_field='melting_layer',
                    ml_pos_field='melting_layer_height',
                    temp_ref=temp_ref, get_iso0=get_iso0, ml_global=None))
            dscfg['initialized'] = True
        else:
            # use previous detection
            ml_obj, ml_dict, iso0_dict, ml_global = (
                pyart.retrieve.melting_layer_giangrande(
                    radar, nVol=nVol, maxh=maxh, hres=hres, rmin=rmin,
                    elmin=elmin, elmax=elmax, rhomin=rhomin, rhomax=rhomax,
                    zhmin=zhmin, hwindow=hwindow, mlzhmin=mlzhmin,
                    mlzhmax=mlzhmax, mlzdrmin=mlzdrmin, mlzdrmax=mlzdrmax,
                    htol=htol, ml_bottom_diff_max=ml_bottom_diff_max,
                    time_accu_max=time_accu_max, nml_points_min=nml_points_min,
                    wlength=wlength, percentile_bottom=percentile_bottom,
                    percentile_top=percentile_top, interpol=interpol,
                    time_nodata_allowed=time_nodata_allowed,
                    refl_field=refl_field, zdr_field=zdr_field,
                    rhv_field=rhv_field, temp_field=temp_field,
                    iso0_field=iso0_field, ml_field='melting_layer',
                    ml_pos_field='melting_layer_height',
                    temp_ref=temp_ref, get_iso0=get_iso0,
                    ml_global=dscfg['global_data']))

        # update global stack
        dscfg['global_data'] = ml_global

    elif dscfg['ML_METHOD'] == 'WOLFENSBERGER':
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            if datatype == 'dBZ':
                refl_field = 'reflectivity'
            if datatype == 'dBZc':
                refl_field = 'corrected_reflectivity'
            if datatype == 'RhoHV':
                rhohv_field = 'cross_correlation_ratio'
            if datatype == 'RhoHVc':
                rhohv_field = 'corrected_cross_correlation_ratio'

        ind_rad = int(radarnr[5:8])-1
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if ((refl_field not in radar.fields) or
                (rhohv_field not in radar.fields)):
            warn('Unable to detect melting layer. Missing data')
            return None, None

        # User defined parameters
        max_range = dscfg.get('max_range', 20000.)
        detect_threshold = dscfg.get('detect_threshold', 0.02)
        interp_holes = dscfg.get('interp_holes', False)
        max_length_holes = dscfg.get('max_length_holes', 250)
        check_min_length = dscfg.get('check_min_length', True)
        get_iso0 = dscfg.get('get_iso0', True)

        ml_obj, ml_dict, iso0_dict, _ = pyart.retrieve.detect_ml(
            radar, refl_field=refl_field, rhohv_field=rhohv_field,
            ml_field='melting_layer', ml_pos_field='melting_layer_height',
            iso0_field='height_over_iso0', max_range=max_range,
            detect_threshold=detect_threshold, interp_holes=interp_holes,
            max_length_holes=max_length_holes,
            check_min_length=check_min_length, get_iso0=get_iso0)

    elif dscfg['ML_METHOD'] == 'MF':
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            if datatype == 'dBZ':
                refl_field = 'reflectivity'
            if datatype == 'dBZc':
                refl_field = 'corrected_reflectivity'
            if datatype == 'RhoHV':
                rhohv_field = 'cross_correlation_ratio'
            if datatype == 'RhoHVc':
                rhohv_field = 'corrected_cross_correlation_ratio'

        ind_rad = int(radarnr[5:8])-1
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if ((refl_field not in radar.fields) or
                (rhohv_field not in radar.fields)):
            warn('Unable to detect melting layer. Missing data')
            return None, None

        # User defined parameters
        max_range = dscfg.get('max_range', 20000.)
        detect_threshold = dscfg.get('detect_threshold', 0.02)
        interp_holes = dscfg.get('interp_holes', False)
        max_length_holes = dscfg.get('max_length_holes', 250)
        check_min_length = dscfg.get('check_min_length', True)
        get_iso0 = dscfg.get('get_iso0', True)

        ml_obj, ml_dict, iso0_dict, _ = pyart.retrieve.melting_layer_mf(
            radar, refl_field=refl_field, rhohv_field=rhohv_field,
            ml_field='melting_layer', ml_pos_field='melting_layer_height',
            iso0_field='height_over_iso0', max_range=max_range,
            detect_threshold=detect_threshold, interp_holes=interp_holes,
            max_length_holes=max_length_holes,
            check_min_length=check_min_length, get_iso0=get_iso0)

        ml_memory_max = dscfg.get('ml_memory_max', 0.)
        datatypedescr = dscfg.get('ml_datatype', None)

        # read the retrieved ml from the past X hours
        ml_thickness_arr = np.ma.array([])
        ml_bottom_arr = np.ma.array([])
        age_arr = np.ma.array([])
        ang_arr = np.ma.array([])
        ml_memory = None
        if ml_memory_max > 0:
            if (datatypedescr is None or dscfg['loadbasepath'] is None
                    or dscfg['loadname'] is None):
                warn('unable to find files containing'
                     ' melting layer information')
            else:
                flist = get_file_list(
                    datatypedescr,
                    [dscfg['timeinfo']
                     - datetime.timedelta(hours=ml_memory_max)],
                    [dscfg['timeinfo']], dscfg)
                if not flist:
                    warn('No files with melting information found')
                else:
                    for fname in flist:
                        radar_ml = pyart.io.read_cfradial(fname)
                        if radar_ml is None:
                            warn('Unable to use retrieved melting layer data')
                            continue

                        ml_top = (
                            radar_ml.fields['melting_layer_height']['data'][:, 1])
                        ml_bottom = (
                            radar_ml.fields['melting_layer_height']['data'][:, 0])
                        ml_thickness = ml_top-ml_bottom
                        ml_bottom_arr = np.ma.append(
                            ml_bottom_arr,
                            np.ma.zeros(radar_ml.nsweeps)+ml_bottom)
                        ml_thickness_arr = np.ma.append(
                            ml_thickness_arr,
                            np.ma.zeros(radar_ml.nsweeps)+ml_thickness)
                        ang_arr = np.ma.append(
                            ang_arr, radar_ml.elevation['data'])
                        age_arr = np.ma.append(
                            age_arr,
                            np.ma.zeros(radar_ml.nrays)+(
                                dscfg['timeinfo']
                                - get_datetime(fname,
                                            datatypedescr)).seconds/3600.)
                    ml_memory = {
                        'ml_bottom': ml_bottom_arr,
                        'ml_thickness': ml_thickness_arr,
                        'ang': ang_arr,
                        'age': age_arr}

        (ml_obj, ml_dict, iso0_dict,
         ml_retrieved) = pyart.retrieve.melting_layer_mf(
            radar, nvalid_min=nvalid_min, ml_thickness_min=ml_thickness_min,
            ml_thickness_max=ml_thickness_max,
            ml_thickness_step=ml_thickness_step, iso0_max=iso0_max,
            ml_top_diff_max=ml_top_diff_max, ml_top_step=ml_top_step,
            rhohv_snow=rhohv_snow, rhohv_rain=rhohv_rain, rhohv_ml=rhohv_ml,
            zh_snow=zh_snow, zh_rain=zh_rain, zh_ml=zh_ml, zv_snow=zv_snow,
            zv_rain=zv_rain, zv_ml=zv_ml, h_max=h_max, h_res=h_res,
            beam_factor=beam_factor, npts_diagram=npts_diagram,
            rng_bottom_max=rng_bottom_max, ns_factor=ns_factor,
            rhohv_corr_min=rhohv_corr_min, rhohv_nash_min=rhohv_nash_min,
            ang_iso0=ang_iso0, age_iso0=age_iso0,
            ml_thickness_iso0=ml_thickness_iso0, ml_memory=ml_memory,
            rhohv_field_obs=rhohv_field_obs, temp_field=temp_field,
            iso0_field=iso0_field,
            rhohv_field_theo='theoretical_cross_correlation_ratio',
            temp_ref=temp_ref, get_iso0=get_iso0)

    elif dscfg['ML_METHOD'] == 'FROM_HYDROCLASS':
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            if datatype == 'hydro':
                hydro_field = get_fieldname_pyart(datatype)

        ind_rad = int(radarnr[5:8])-1
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if hydro_field not in radar.fields:
            warn('Unable to detect melting layer. Missing data')
            return None, None

        # User defined parameters
        force_continuity = dscfg.get('force_continuity', True)
        dist_max = dscfg.get('dist_max', 350.)
        get_iso0 = dscfg.get('get_iso0', False)

        ml_obj, ml_dict, iso0_dict = pyart.retrieve.melting_layer_hydroclass(
            radar, hydro_field=hydro_field, ml_field='melting_layer',
            ml_pos_field='melting_layer_height',
            iso0_field='height_over_iso0', force_continuity=force_continuity,
            dist_max=dist_max, get_iso0=get_iso0)

    else:
        raise Exception(
            "ERROR: Unknown melting layer retrieval method {}".format(
                dscfg['ML_METHOD']))

    # prepare for exit
    if ml_dict is None:
        return None, None

    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field('melting_layer', ml_dict)
    if iso0_dict is not None:
        new_dataset['radar_out'].add_field('height_over_iso0', iso0_dict)
    new_dataset.update({'ml_obj': ml_obj})

    return new_dataset, ind_rad


def process_zdr_column(procstatus, dscfg, radar_list=None):
    """
    Detects ZDR columns

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

    temp_field = None
    iso0_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'ZDR':
            zdr_field = 'differential_reflectivity'
        if datatype == 'ZDRc':
            zdr_field = 'corrected_differential_reflectivity'
        if datatype == 'RhoHV':
            rhv_field = 'cross_correlation_ratio'
        if datatype == 'RhoHVc':
            rhv_field = 'corrected_cross_correlation_ratio'
        if datatype == 'TEMP':
            temp_field = 'temperature'
        if datatype == 'H_ISO0':
            iso0_field = 'height_over_iso0'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    # Check which should be the reference field for temperature
    if iso0_field is not None and (iso0_field not in radar.fields):
        warn('Unable to detect melting layer. '
             'Missing height over iso0 field')
        return None, None
    temp_ref = 'height_over_iso0'

    if temp_field is not None and (temp_field not in radar.fields):
        warn('Unable to detect melting layer. Missing temperature field')
        return None, None
    temp_ref = 'temperature'
    iso0_field = 'height_over_iso0'

    if ((zdr_field not in radar.fields) or
            (rhv_field not in radar.fields)):
        warn('Unable to detect melting layer. Missing data')
        return None, None

    rhohv_min = dscfg.get('rhohv_min', 0.8)
    zdr_min = dscfg.get('zdr_min', 1.)
    smooth_window = dscfg.get('smooth_window', 0.)
    latlon_tol = dscfg.get('latlon_tol', 0.025)  # approx 3x2 km
    if smooth_window == 0:
        smooth_window_len = 0
    else:
        smooth_window_len = int(
            smooth_window/(radar.range['data'][1]-radar.range['data'][0]))

    zdr_dict = deepcopy(radar.fields[zdr_field])

    if smooth_window_len > 0:
        zdr_dict['data'] = pyart.correct.smooth_masked(
            zdr_dict['data'], wind_len=smooth_window_len, min_valid=1,
            wind_type='mean')

    zdr_dict['data'][
        radar.fields[rhv_field]['data'] < rhohv_min] = np.ma.masked
    zdr_dict['data'][zdr_dict['data'] < zdr_min] = np.ma.masked
    zdr_dict['data'][radar.fields[temp_field]['data'] > 0.] = np.ma.masked
    zdr_valid = np.logical_not(np.ma.getmaskarray(zdr_dict['data']))

    hlowerleft, hupperright = pyart.retrieve._get_res_vol_sides(radar)
    ind_ang_sorted = np.argsort(radar.fixed_angle['data'])

    # get number of suspected ZDR columns
    lat_cols = np.array([], dtype=int)
    lon_cols = np.array([], dtype=int)
    zdr_cols = np.array([], dtype=int)

    g_lat = radar.gate_latitude['data']
    g_lon = radar.gate_longitude['data']
    for ind_ray in range(radar.nrays):
        # Get bins with negative temperatures
        ind_rngs = np.where(
            radar.fields[temp_field]['data'][ind_ray, :] < 0.)[0]
        if ind_rngs.size == 0:
            continue

        # Segment negative temperatures and get start of each segment
        cons_list = np.split(ind_rngs, np.where(np.diff(ind_rngs) != 1)[0]+1)
        for ind_rngs_cell in cons_list:
            if not zdr_valid[ind_ray, ind_rngs_cell[0]]:
                continue

            ind_ray_col = ind_ray
            ind_rng_col = ind_rngs_cell[0]

            # extract data around point:
            ind_rays, ind_rngs = np.where(np.logical_and.reduce((
                np.logical_and(
                    g_lat >= g_lat[ind_ray_col, ind_rng_col]-latlon_tol,
                    g_lat <= g_lat[ind_ray_col, ind_rng_col]+latlon_tol),
                np.logical_and(
                    g_lon >= g_lon[ind_ray_col, ind_rng_col]-latlon_tol,
                    g_lon <= g_lon[ind_ray_col, ind_rng_col]+latlon_tol),
                zdr_valid)))

            # get ZDR column height for each radar sweep
            h_low = np.ma.masked_all(radar.nsweeps)
            h_high = np.ma.masked_all(radar.nsweeps)
            for sweep in range(radar.nsweeps):
                ind = np.where(np.logical_and(
                    ind_rays >= radar.sweep_start_ray_index['data'][sweep],
                    ind_rays <= radar.sweep_end_ray_index['data'][sweep]))[0]
                if ind.size == 0:
                    continue

                h_low[sweep] = np.min(
                    hlowerleft[ind_rays[ind], ind_rngs[ind]])
                h_high[sweep] = np.max(
                    hupperright[ind_rays[ind], ind_rngs[ind]])

            # order data by elevation angle
            h_low = h_low[ind_ang_sorted]
            h_high = h_high[ind_ang_sorted]

            # get the first segment of continuous ZDR valid values
            ind_valid = np.where(np.ma.getmaskarray(h_low) == 0)[0]
            ind_valid = np.split(
                ind_valid, np.where(np.diff(ind_valid) != 1)[0]+1)[0]

            # compute ZDR column
            zdr_col = h_high[ind_valid[-1]]-h_low[ind_valid[0]]

            # put data in output array
            lat_cols = np.append(
                lat_cols,
                radar.gate_latitude['data'][ind_ray_col, ind_rng_col])
            lon_cols = np.append(
                lon_cols,
                radar.gate_longitude['data'][ind_ray_col, ind_rng_col])
            zdr_cols = np.append(zdr_cols, zdr_col)

    zdr_col_dict = pyart.config.get_metadata(
        'differential_reflectivity_column_height')
    zdr_col_dict['data'] = zdr_cols/1000.
    new_dataset = {
        'field_limits': [
            np.min(radar.gate_longitude['data']),
            np.max(radar.gate_longitude['data']),
            np.min(radar.gate_latitude['data']),
            np.max(radar.gate_latitude['data'])],
        'lat': lat_cols,
        'lon': lon_cols,
        'fields': {'differential_reflectivity_column_height': zdr_col_dict}}

    return new_dataset, ind_rad


class ScanObject:
    """generates scan object containing all required radar parameters"""

    def __init__(self, radar, refl_field, iso0_field, dp_fields,
                 dualpol_vars_int=('RHOHV', 'ZDR', 'KDP')):
        """initialises required radar fields / parameters"""

        # (currently pixels are only considered where their flag value
        # matches the integer specified by classify.MET_FLAG_VAL)
        if radar.instrument_parameters is None:
            warn('Radar beamwidth not specified. Assumed 1. deg')
            self.beam_width_degrees = 1.
        elif 'radar_beam_width_h' not in radar.instrument_parameters:
            warn('Radar beamwidth not specified. Assumed 1. deg')
            self.beam_width_degrees = 1.
        else:
            self.beam_width_degrees = (
                radar.instrument_parameters['radar_beam_width_h']['data'][0])
        self.site_altitude_metres = radar.altitude['data'][0]
        self.n_rays = radar.nrays
        self.n_bins = radar.ngates
        self.scan_elevation_degrees = radar.fixed_angle['data'][0]
        self.bin_length_km = (
            (radar.range['data'][1]-radar.range['data'][0])/1000.)

        # from height over iso0 to iso0 height
        self.freezing_level = (
            radar.gate_altitude['data']-radar.fields[iso0_field]['data'])

        self.zh_lin = np.power(10., 0.1*radar.fields[refl_field]['data'])

        # assumes only precipitation has valid reflectivity
        self.flags = np.ma.getmaskarray(radar.fields[refl_field]['data'])

        # 'secondary' (dualpol) fields paired with Zh in membership functions:
        dualpol_fields = {}
        for field_name, field_name_int in zip(dp_fields, dualpol_vars_int):
            dualpol_fields[field_name_int] = radar.fields[field_name]['data']
        self.dualpol_fields = dualpol_fields
