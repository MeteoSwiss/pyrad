"""
pyrad.proc.process_echoclass
===============================

Functions for echo classification and filtering

.. autosummary::
    :toctree: generated/

    process_echo_id
    process_birds_id
    process_clt_to_echo_id
    process_hydro_mf_to_hydro
    process_echo_filter
    process_cdf
    process_filter_snr
    process_filter_vel_diff
    process_filter_visibility
    process_outlier_filter
    process_hydroclass
    process_melting_layer
    process_zdr_column

"""

from copy import deepcopy
from warnings import warn

import numpy as np

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart


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
    hydro[hydroMF<8] = 1
    hydro[hydroMF==30] = 1
    hydro[hydroMF==31] = 1
    # PRECIP_INDIFFERENCIEE, PLUIE, PRECIP are classified as RN
    hydro[hydroMF==8] = 6
    hydro[hydroMF==9] = 6
    hydro[hydroMF==32] = 6
    hydro[hydroMF==10] = 8  # NEIGE_MOUILLEE is WS
    hydro[hydroMF==11] = 2  # NEIGE_SECHE is AG
    hydro[hydroMF==12] = 3  # GLACE is CR
    hydro[hydroMF==13] = 5  # PETITE_GRELE is RP
    # MOYENNE_GRELE, GROSSE_GRELE is IH/HDG
    hydro[hydroMF==14] = 10
    hydro[hydroMF==15] = 10
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
            SEMISUPERVISED
        RADARCENTROIDS : string. Dataset keyword
            Used with HYDRO_METHOD SEMISUPERVISED. The name of the radar of
            which the derived centroids will be used. One of the following: A
            Albis, L Lema, P Plaine Morte, DX50
        compute_entropy : bool. Dataset keyword
            If true the entropy is computed and the field hydroclass_entropy
            is output
        output_distances : bool. Dataset keyword
            If true the de-mixing algorithm based on the distances to the
            centroids is computed and the field proportions of each
            hydrometeor in the radar range gate is output
        vectorize : bool. Dataset keyword
            If true a vectorized version of the algorithm is used
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
            "ERROR: Undefined parameter 'HYDRO_METHOD' for dataset '%s'"
            % dscfg['dsname'])

    if dscfg['HYDRO_METHOD'] == 'SEMISUPERVISED':
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

        mass_centers = np.zeros((9, 5))
        if dscfg['RADARCENTROIDS'] == 'A':
            #      Zh      ZDR     kdp   RhoHV   delta_Z
            mass_centers[0, :] = [
                13.5829, 0.4063, 0.0497, 0.9868, 1330.3]  # AG
            mass_centers[1, :] = [
                02.8453, 0.2457, 0.0000, 0.9798, 0653.8]  # CR
            mass_centers[2, :] = [
                07.6597, 0.2180, 0.0019, 0.9799, -1426.5]  # LR
            mass_centers[3, :] = [
                31.6815, 0.3926, 0.0828, 0.9978, 0535.3]  # RP
            mass_centers[4, :] = [
                39.4703, 1.0734, 0.4919, 0.9876, -1036.3]  # RN
            mass_centers[5, :] = [
                04.8267, -0.5690, 0.0000, 0.9691, 0869.8]  # VI
            mass_centers[6, :] = [
                30.8613, 0.9819, 0.1998, 0.9845, -0066.1]  # WS
            mass_centers[7, :] = [
                52.3969, 2.1094, 2.4675, 0.9730, -1550.2]  # MH
            mass_centers[8, :] = [
                50.6186, -0.0649, 0.0946, 0.9904, 1179.9]  # IH/HDG
        elif dscfg['RADARCENTROIDS'] == 'L':
            #       Zh      ZDR     kdp   RhoHV   delta_Z
            mass_centers[0, :] = [
                13.8231, 0.2514, 0.0644, 0.9861, 1380.6]  # AG
            mass_centers[1, :] = [
                03.0239, 0.1971, 0.0000, 0.9661, 1464.1]  # CR
            mass_centers[2, :] = [
                04.9447, 0.1142, 0.0000, 0.9787, -0974.7]  # LR
            mass_centers[3, :] = [
                34.2450, 0.5540, 0.1459, 0.9937, 0945.3]  # RP
            mass_centers[4, :] = [
                40.9432, 1.0110, 0.5141, 0.9928, -0993.5]  # RN
            mass_centers[5, :] = [
                03.5202, -0.3498, 0.0000, 0.9746, 0843.2]  # VI
            mass_centers[6, :] = [
                32.5287, 0.9751, 0.2640, 0.9804, -0055.5]  # WS
            mass_centers[7, :] = [
                52.6547, 2.7054, 2.5101, 0.9765, -1114.6]  # MH
            mass_centers[8, :] = [
                46.4998, 0.1978, 0.6431, 0.9845, 1010.1]  # IH/HDG
        elif dscfg['RADARCENTROIDS'] == 'D':
            #       Zh      ZDR     kdp   RhoHV   delta_Z
            mass_centers[0, :] = [
                12.567, 0.18934, 0.041193, 0.97693, 1328.1]  # AG
            mass_centers[1, :] = [
                3.2115, 0.13379, 0.0000, 0.96918, 1406.3]  # CR
            mass_centers[2, :] = [
                10.669, 0.18119, 0.0000, 0.97337, -1171.9]  # LR
            mass_centers[3, :] = [
                34.941, 0.13301, 0.090056, 0.9979, 898.44]  # RP
            mass_centers[4, :] = [
                39.653, 1.1432, 0.35013, 0.98501, -859.38]  # RN
            mass_centers[5, :] = [
                2.8874, -0.46363, 0.0000, 0.95653, 1015.6]  # VI
            mass_centers[6, :] = [
                34.122, 0.87987, 0.2281, 0.98003, -234.37]  # WS
            mass_centers[7, :] = [
                53.134, 2.0888, 2.0055, 0.96927, -1054.7]  # MH
            mass_centers[8, :] = [
                46.715, 0.030477, 0.16994, 0.9969, 976.56]  # IH/HDG
        elif dscfg['RADARCENTROIDS'] == 'P':
            #       Zh      ZDR     kdp   RhoHV   delta_Z
            mass_centers[0, :] = [
                13.9882, 0.2470, 0.0690, 0.9939, 1418.1]  # AG
            mass_centers[1, :] = [
                00.9834, 0.4830, 0.0043, 0.9834, 0950.6]  # CR
            mass_centers[2, :] = [
                05.3962, 0.2689, 0.0000, 0.9831, -0479.5]  # LR
            mass_centers[3, :] = [
                35.3411, 0.1502, 0.0940, 0.9974, 0920.9]  # RP
            mass_centers[4, :] = [
                35.0114, 0.9681, 0.1106, 0.9785, -0374.0]  # RN
            mass_centers[5, :] = [
                02.5897, -0.3879, 0.0282, 0.9876, 0985.5]  # VI
            mass_centers[6, :] = [
                32.2914, 0.7789, 0.1443, 0.9075, -0153.5]  # WS
            mass_centers[7, :] = [
                53.2413, 1.8723, 0.3857, 0.9454, -0470.8]  # MH
            mass_centers[8, :] = [
                44.7896, 0.0015, 0.1349, 0.9968, 1116.7]  # IH/HDG
        elif dscfg['RADARCENTROIDS'] == 'W':
            #       Zh      ZDR     kdp   RhoHV   delta_Z
            mass_centers[0, :] = [
                16.7650, 0.3754, 0.0442, 0.9866, 1409.0]  # AG
            mass_centers[1, :] = [
                01.4418, 0.3786, 0.0000, 0.9490, 1415.8]  # CR
            mass_centers[2, :] = [
                16.0987, 0.3238, 0.0000, 0.9871, -0818.7]  # LR
            mass_centers[3, :] = [
                36.5465, 0.2041, 0.0731, 0.9952, 0745.4]  # RP
            mass_centers[4, :] = [
                43.4011, 0.6658, 0.3241, 0.9894, -0778.5]  # RN
            mass_centers[5, :] = [
                00.9077, -0.4793, 0.0000, 0.9502, 1488.6]  # VI
            mass_centers[6, :] = [
                36.8091, 0.7266, 0.1284, 0.9924, -0071.1]  # WS
            mass_centers[7, :] = [
                53.8402, 0.8922, 0.5306, 0.9890, -1017.6]  # MH
            mass_centers[8, :] = [
                45.9686, 0.0845, 0.0963, 0.9940, 0867.4]  # IH/HDG
        elif dscfg['RADARCENTROIDS'] == 'DX50':
            #       Zh      ZDR     kdp   RhoHV   delta_Z
            mass_centers[0, :] = [
                19.0770, 0.4139, 0.0099, 0.9841, 1061.7]  # AG
            mass_centers[1, :] = [
                03.9877, 0.5040, 0.0000, 0.9642, 0856.6]  # CR
            mass_centers[2, :] = [
                20.7982, 0.3177, 0.0004, 0.9858, -1375.1]  # LR
            mass_centers[3, :] = [
                34.7124, -0.3748, 0.0988, 0.9828, 1224.2]  # RP
            mass_centers[4, :] = [
                33.0134, 0.6614, 0.0819, 0.9802, -1169.8]  # RN
            mass_centers[5, :] = [
                08.2610, -0.4681, 0.0000, 0.9722, 1100.7]  # VI
            mass_centers[6, :] = [
                35.1801, 1.2830, 0.1322, 0.9162, -0159.8]  # WS
            mass_centers[7, :] = [
                52.4539, 2.3714, 1.1120, 0.9382, -1618.5]  # MH
            mass_centers[8, :] = [
                44.2216, -0.3419, 0.0687, 0.9683, 1272.7]  # IH/HDG
        else:
            warn(
                ' Unknown radar. ' +
                'Default centroids will be used in classification.')
            mass_centers = None

        compute_entropy = dscfg.get('compute_entropy', False)
        output_distances = dscfg.get('output_distances', False)
        vectorize = dscfg.get('vectorize', False)

        fields_dict = pyart.retrieve.hydroclass_semisupervised(
            radar, mass_centers=mass_centers,
            weights=np.array([1., 1., 1., 0.75, 0.5]), refl_field=refl_field,
            zdr_field=zdr_field, rhv_field=rhv_field, kdp_field=kdp_field,
            temp_field=temp_field, iso0_field=iso0_field, hydro_field=None,
            entropy_field=None, temp_ref=temp_ref,
            compute_entropy=compute_entropy,
            output_distances=output_distances, vectorize=vectorize)
    else:
        raise Exception(
            "ERROR: Unknown hydrometeor classification method " +
            dscfg['HYDRO_METHOD'])

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(
        'radar_echo_classification', fields_dict['hydro'])

    if compute_entropy:
        new_dataset['radar_out'].add_field(
            'hydroclass_entropy', fields_dict['entropy'])

        if output_distances:
            new_dataset['radar_out'].add_field(
                'proportion_AG', fields_dict['prop_AG'])
            new_dataset['radar_out'].add_field(
                'proportion_CR', fields_dict['prop_CR'])
            new_dataset['radar_out'].add_field(
                'proportion_LR', fields_dict['prop_LR'])
            new_dataset['radar_out'].add_field(
                'proportion_RP', fields_dict['prop_RP'])
            new_dataset['radar_out'].add_field(
                'proportion_RN', fields_dict['prop_RN'])
            new_dataset['radar_out'].add_field(
                'proportion_VI', fields_dict['prop_VI'])
            new_dataset['radar_out'].add_field(
                'proportion_WS', fields_dict['prop_WS'])
            new_dataset['radar_out'].add_field(
                'proportion_MH', fields_dict['prop_MH'])
            new_dataset['radar_out'].add_field(
                'proportion_IH', fields_dict['prop_IH'])

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
        raise Exception(
            "ERROR: Undefined parameter 'ML_METHOD' for dataset '%s'"
            % dscfg['dsname'])

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
                warn('Unable to detect melting layer. ' +
                     'Missing height over iso0 field')
                return None, None
            temp_ref = 'height_over_iso0'

        if temp_field is not None:
            if temp_field not in radar.fields:
                warn('Unable to detect melting layer. ' +
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
            "ERROR: Unknown melting layer retrieval method " +
            dscfg['ML_METHOD'])

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
        warn('Unable to detect melting layer. ' +
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
