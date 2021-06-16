"""
pyrad.prod.process_product
==========================

Functions for obtaining Pyrad products from the datasets

.. autosummary::
    :toctree: generated/

    generate_occurrence_products
    generate_cosmo_coord_products
    generate_cosmo_to_radar_products
    generate_sun_hits_products
    generate_qvp_products
    generate_ml_products
    generate_centroids_products

"""

from copy import deepcopy
from warnings import warn
import os

import numpy as np

import pyart

from .process_vol_products import generate_vol_products

from ..io.io_aux import get_fieldname_pyart
from ..io.io_aux import get_save_dir, make_filename

from ..io.read_data_sun import read_sun_retrieval
from ..io.read_data_other import read_ml_ts

from ..io.write_data import write_sun_hits, write_sun_retrieval
from ..io.write_data import write_excess_gates, write_ts_ml, write_histogram
from ..io.write_data import write_timeseries_point, write_centroids

from ..graph.plots import plot_sun_hits, plot_histogram2, plot_scatter
from ..graph.plots import plot_centroids
from ..graph.plots_timeseries import plot_sun_retrieval_ts, plot_ml_ts
from ..graph.plots_vol import plot_fixed_rng, plot_fixed_rng_sun

from ..util.radar_utils import create_sun_hits_field, compute_histogram
from ..util.radar_utils import create_sun_retrieval_field

import matplotlib as mpl
mpl.use('Agg')

# Increase a bit font size
mpl.rcParams.update({'font.size': 16})
mpl.rcParams.update({'font.family':  "sans-serif"})

import matplotlib.pyplot as plt


def generate_occurrence_products(dataset, prdcfg):
    """
    generates occurrence products. Accepted product types:
        'WRITE_EXCESS_GATES': Write the data that identifies radar gates
            with clutter that has a frequency of occurrence above a certain
            threshold.
            User defined parameters:
                quant_min: float
                    Minimum frequency of occurrence in percentage to keep the
                    gate as valid. Default 95.
        All the products of the 'VOL' dataset group

    Parameters
    ----------
    dataset : tuple
        radar object and metadata dictionary

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    filename : str
        the name of the file created. None otherwise

    """
    instant = False
    if 'instant' in prdcfg:
        instant = prdcfg['instant']

    if not instant and not dataset['occu_final']:
        return None

    if prdcfg['type'] == 'WRITE_EXCESS_GATES':
        if not dataset['occu_final']:
            return None

        radar = dataset['radar_out']
        if (('frequency_of_occurrence' not in radar.fields) or
                ('occurrence' not in radar.fields) or
                ('number_of_samples' not in radar.fields)):
            warn('Unable to create quantile excess gates file. '
                 'Missing data')
            return None

        dssavedir = prdcfg['dsname']
        if 'dssavename' in prdcfg:
            dssavedir = prdcfg['dssavename']

        quant_min = 95.
        if 'quant_min' in prdcfg:
            quant_min = prdcfg['quant_min']

        # get index of gates exceeding quantile
        freq_occu = radar.fields['frequency_of_occurrence'][
            'data']
        ind_ray, ind_rng = np.where(freq_occu > quant_min)
        if ind_ray.size == 0:
            warn('No data exceeds the frequency of occurrence ' +
                 str(quant_min)+' %')
            return None

        excess_dict = {
            'starttime': dataset['starttime'],
            'endtime': dataset['endtime'],
            'quant_min': quant_min,
            'ray_ind': ind_ray,
            'rng_ind': ind_rng,
            'ele': radar.elevation['data'][ind_ray],
            'azi': radar.azimuth['data'][ind_ray],
            'rng': radar.range['data'][ind_rng],
            'nsamples': (
                radar.fields['number_of_samples']['data'][ind_ray, ind_rng]),
            'occurrence': (
                radar.fields['occurrence']['data'][ind_ray, ind_rng]),
            'freq_occu': freq_occu[ind_ray, ind_rng]
        }
        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=dataset['endtime'])

        fname = make_filename(
            'excess_gates', prdcfg['dstype'], prdcfg['prdname'], ['csv'],
            prdcfginfo='quant'+'{:.1f}'.format(quant_min),
            timeinfo=dataset['endtime'])

        fname = savedir+fname[0]

        fname = write_excess_gates(excess_dict, fname)

        if fname is not None:
            print('saved excess gates file: {}'.format(fname))

        return fname

    field_name = get_fieldname_pyart(prdcfg['voltype'])
    if ((field_name == 'frequency_of_occurrence') and
            (not dataset['occu_final'])):
        return None
    if dataset['occu_final']:
        prdcfg['timeinfo'] = dataset['endtime']

    return generate_vol_products(dataset, prdcfg)


def generate_cosmo_coord_products(dataset, prdcfg):
    """
    generates COSMO coordinates products. Accepted product types:
        'SAVEVOL': Save an object containing the index of the COSMO model grid
            that corresponds to each radar gate in a C/F radial file.
            User defined parameters:
                file_type: str
                    The type of file used to save the data. Can be 'nc' or
                    'h5'. Default 'nc'
                physical: Bool
                    If True the data will be saved in physical units (floats).
                    Otherwise it will be quantized and saved as binary
                compression: str
                    For ODIM file formats, the type of compression. Can be any
                    of the allowed compression types for hdf5 files. Default
                    gzip
                compression_opts: any
                    The compression options allowed by the hdf5. Depends on
                    the type of compression. Default 6 (The gzip compression
                    level).

    Parameters
    ----------
    dataset : tuple
        radar object containing the COSMO coordinates

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    filename : str
        the name of the file created. None otherwise

    """
    if prdcfg['type'] == 'SAVEVOL':
        radar_obj = dataset['radar_out']
        ind_rad = dataset['ind_rad']

        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in radar_obj.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        file_type = prdcfg.get('file_type', 'nc')
        physical = prdcfg.get('physical', True)
        compression = prdcfg.get('compression', 'gzip')
        compression_opts = prdcfg.get('compression_opts', 6)

        new_dataset = deepcopy(radar_obj)
        new_dataset.fields = dict()
        new_dataset.add_field(field_name, radar_obj.fields[field_name])

        savedir = prdcfg['cosmopath'][ind_rad]+'rad2cosmo/'
        fname = 'rad2cosmo_'+prdcfg['voltype']+'_'+prdcfg['procname']+'.nc'

        if file_type == 'nc':
            pyart.io.cfradial.write_cfradial(
                savedir+fname, new_dataset, physical=physical)
        elif file_type == 'h5':
            pyart.aux_io.write_odim_h5(
                savedir+fname, new_dataset, physical=physical,
                compression=compression, compression_opts=compression_opts)
        else:
            warn('Data could not be saved. ' +
                 'Unknown saving file type '+file_type)
            return None

        print('saved file: {}'.format(savedir+fname))

        return fname

    warn(' Unsupported product type: ' + prdcfg['type'])
    return None


def generate_cosmo_to_radar_products(dataset, prdcfg):
    """
    generates COSMO data in radar coordinates products. Accepted product
    types:
        'SAVEVOL': Save an object containing the COSMO data in radar
            coordinatesin a C/F radial or ODIM file.
                User defined parameters:
                file_type: str
                    The type of file used to save the data. Can be 'nc' or
                    'h5'. Default 'nc'
                physical: Bool
                    If True the data will be saved in physical units (floats).
                    Otherwise it will be quantized and saved as binary
                compression: str
                    For ODIM file formats, the type of compression. Can be any
                    of the allowed compression types for hdf5 files. Default
                    gzip
                compression_opts: any
                    The compression options allowed by the hdf5. Depends on
                    the type of compression. Default 6 (The gzip compression
                    level).
        All the products of the 'VOL' dataset group

    Parameters
    ----------
    dataset : tuple
        radar object containing the COSMO coordinates

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    filename : str
        the name of the file created. None otherwise

    """
    time_index = prdcfg.get('cosmo_time_index', 0)
    if time_index > len(dataset)-1:
        warn(
            'COSMO time index larger than available. Skipping product ' +
            prdcfg['type'])
        return None

    radar_dataset = dataset[time_index]
    if prdcfg['type'] == 'SAVEVOL':
        radar_obj = radar_dataset['radar_out']
        ind_rad = radar_dataset['ind_rad']

        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in radar_obj.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        file_type = prdcfg.get('file_type', 'nc')
        physical = prdcfg.get('physical', True)
        compression = prdcfg.get('compression', 'gzip')
        compression_opts = prdcfg.get('compression_opts', 6)

        new_dataset = deepcopy(radar_obj)
        new_dataset.fields = dict()
        new_dataset.add_field(field_name, radar_obj.fields[field_name])

        savedir = (
            prdcfg['cosmopath'][ind_rad]+prdcfg['voltype']+'/radar/' +
            prdcfg['timeinfo'].strftime('%Y-%m-%d')+'/'+prdcfg['procname']+'/')
        fname = (
            prdcfg['voltype']+'_RUN' +
            prdcfg['timeinfo'].strftime('%Y%m%d%H%M%S')+'_' +
            radar_dataset['dtcosmo'].strftime('%Y%m%d%H%M%S')+'.nc')

        if not os.path.isdir(savedir):
            os.makedirs(savedir)

        if file_type == 'nc':
            pyart.io.cfradial.write_cfradial(
                savedir+fname, new_dataset, physical=physical)
        elif file_type == 'h5':
            pyart.aux_io.write_odim_h5(
                savedir+fname, new_dataset, physical=physical,
                compression=compression, compression_opts=compression_opts)
        else:
            warn('Data could not be saved. ' +
                 'Unknown saving file type '+file_type)
            return None

        print('saved file: {}'.format(savedir+fname))

        return fname

    return generate_vol_products(radar_dataset, prdcfg)


def generate_sun_hits_products(dataset, prdcfg):
    """
    generates sun hits products. Accepted product types:
        'PLOT_SUN_HITS': Plots in a sun-radar azimuth difference-sun-radar
            elevation difference grid the values of all sun hits obtained
            during the processing period
        'PLOT_SUN_RETRIEVAL': Plots in a sun-radar azimuth difference-sun-
            radar elevation difference grid the retrieved sun pattern
        'PLOT_SUN_RETRIEVAL_TS': Plots time series of the retrieved sun
            pattern parameters
            User defined parameters:
                dpi: int
                    The pixel density of the plot. Default 72
                add_date_in_fname: Bool
                    If true the year is added in the plot file name
        'PLOT_SUNSCAN': Plots a constant range radar azimuth-elevation of the
            sunscan field data
        'WRITE_SUN_HITS': Writes the information concerning possible sun hits
            in a csv file
        'WRITE_SUN_RETRIEVAL': Writes the retrieved sun pattern parameters in
            a csv file.
            User defined parameters:
                add_date_in_fname: Bool
                    If true the year is added in the csv file name
        'WRITE_SUNSCAN': Writes the sunscan parameters in a csv file

        All the products of the 'VOL' dataset group

    Parameters
    ----------
    dataset : tuple
        radar object and sun hits dictionary

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    filename : str
        the name of the file created. None otherwise

    """

    dssavedir = prdcfg['dsname']
    if 'dssavename' in prdcfg:
        dssavedir = prdcfg['dssavename']

    prdcfg['timeinfo'] = dataset['timeinfo']

    if prdcfg['type'] == 'WRITE_SUN_HITS':
        if 'sun_hits' not in dataset:
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=dataset['timeinfo'])

        fname = make_filename(
            'info', prdcfg['dstype'], 'detected', ['csv'],
            timeinfo=dataset['timeinfo'], timeformat='%Y%m%d')[0]

        fname = savedir+fname

        write_sun_hits(dataset['sun_hits'], fname)

        print('saved sun hits file: {}'.format(fname))

        return fname[0]

    if prdcfg['type'] == 'PLOT_SUN_HITS':
        if 'sun_hits_final' not in dataset:
            return None

        field_name = get_fieldname_pyart(prdcfg['voltype'])

        if prdcfg['voltype'] not in dataset['sun_hits_final']:
            warn(
                ' Field type ' + prdcfg['voltype'] +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=dataset['timeinfo'])

        fname_list = make_filename(
            'detected', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], timeinfo=dataset['timeinfo'],
            timeformat='%Y%m%d')

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        field = create_sun_hits_field(
            dataset['sun_hits_final']['rad_el'],
            dataset['sun_hits_final']['rad_az'],
            dataset['sun_hits_final']['sun_el'],
            dataset['sun_hits_final']['sun_az'],
            dataset['sun_hits_final'][prdcfg['voltype']],
            prdcfg['sunhitsImageConfig'])

        if field is None:
            warn(
                'Unable to create field '+prdcfg['voltype'] +
                ' Skipping product ' + prdcfg['type'])
            return None

        plot_sun_hits(field, field_name, fname_list, prdcfg)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'WRITE_SUN_RETRIEVAL':
        if 'sun_retrieval' not in dataset:
            return None

        timeinfo = None
        timeformat = None
        if prdcfg.get('add_date_in_fname', False):
            timeinfo = dataset['timeinfo']
            timeformat = '%Y'

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=None)

        fname = make_filename(
            'info', prdcfg['dstype'], 'retrieval', ['csv'], timeinfo=timeinfo,
            timeformat=timeformat, runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname

        write_sun_retrieval(dataset['sun_retrieval'], fname)

        print('saved sun retrieval file: {}'.format(fname))

        return fname

    if prdcfg['type'] == 'PLOT_SUN_RETRIEVAL':
        if 'sun_retrieval' not in dataset:
            return None

        field_name = get_fieldname_pyart(prdcfg['voltype'])
        par = None
        if field_name == 'sun_est_power_h':
            par = 'par_h'
        elif field_name == 'sun_est_power_v':
            par = 'par_v'
        elif field_name == 'sun_est_differential_reflectivity':
            par = 'par_zdr'

        if par not in dataset['sun_retrieval']:
            warn(
                ' Field type ' + prdcfg['voltype'] +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=dataset['timeinfo'])

        fname_list = make_filename(
            'retrieval', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], timeinfo=dataset['timeinfo'],
            timeformat='%Y%m%d')

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        if dataset['sun_retrieval'][par] is None:
            warn(
                ' Invalid retrieval parameters. Skipping product ' +
                prdcfg['type'])
            return None

        field = create_sun_retrieval_field(
            dataset['sun_retrieval'][par], field_name,
            prdcfg['sunhitsImageConfig'],
            lant=dataset['sun_retrieval']['lant'])

        if field is not None:
            plot_sun_hits(field, field_name, fname_list, prdcfg)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'PLOT_SUN_RETRIEVAL_TS':
        if 'sun_retrieval' not in dataset:
            return None

        dpi = prdcfg.get('dpi', 72)
        timeinfo = None
        timeformat = None
        if prdcfg.get('add_date_in_fname', False):
            timeinfo = dataset['timeinfo']
            timeformat = '%Y'

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdid'], timeinfo=None)

        fname = make_filename(
            'info', prdcfg['dstype'], 'retrieval', ['csv'], timeinfo=timeinfo,
            timeformat=timeformat, runinfo=prdcfg['runinfo'])

        fname = savedir + fname[0]

        sun_retrieval = read_sun_retrieval(fname)

        if sun_retrieval[0] is None:
            warn(
                'Unable to read sun retrieval file '+fname)
            return None

        if len(sun_retrieval[0]) < 2:
            warn(
                'Unable to plot sun retrieval time series. ' +
                'Not enough data points.')
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=None)

        fname_list = make_filename(
            'retrieval_ts', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], timeinfo=timeinfo,
            timeformat=timeformat, runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        titl = (prdcfg['runinfo']+' Sun Retrieval ' +
                sun_retrieval[1][0].strftime('%Y%m%d')+'-' +
                sun_retrieval[1][-1].strftime('%Y%m%d'))
        figfname = plot_sun_retrieval_ts(
            sun_retrieval, prdcfg['voltype'], fname_list, titl=titl, dpi=dpi)

        if figfname is None:
            return None

        print('----- save to '+' '.join(fname_list))
        return fname_list

    if prdcfg['type'] == 'WRITE_SUNSCAN':
        if 'sun_retrieval' not in dataset:
            return None

        text = [
            "SunScan info",
            "sun_az:             [deg] Azimuth sun position ",
            "sun_el:             [deg] Elevation sun position",
            "noise_pwr:          [dBm] Noise power",
            "sun_maxpwr_noise:   [dBm]"
            " sun maximal power sample (including noise)",
            "sun_maxpwr_nonoise: [dBm]"
            " sun maximal power sample without noise",
            "sun_maxpwr_fit:     [dBm]"
            " sun maximal fitted power (without noise)",
            "sun_maxpwr_toa:     [dBm]"
            " sun maximal power at top of atmosphere",
            "az_offset:          [deg]"
            " Azimuth shift of fitted maxima to sun azimuth",
            "el_offset:          [deg]"
            " Elevation shift of fitted maxima to sun elevation",
            "az_phi3db:          [deg]"
            " Half-power beam width in azimuth",
            "el_phi3db:          [deg]"
            " Half-power beam width in elevation",
            "fit_stddev:         [dBm]"
            " Standard deviation (fit to samples)",
            "num_samples:        [#]"
            " Number of samples used for the sun power fitting"
            ]

        sunRdata = dataset['sun_retrieval']

        if dataset['field_name'] == 'noisedBm_hh':
            data = {
                'dstype': prdcfg['dstype'],
                'unit': 'dBm',
                'time': sunRdata['sunscan_time'],
                'label': [
                    "sun_az", "sun_el", "noise_pwr", "sun_maxpwr_noise",
                    "sun_maxpwr_nonoise", "sun_maxpwr_fit", "sun_maxpwr_toa",
                    "az_offset", "el_offset", "az_phi3db", "el_phi3db",
                    "fit_stddev", "num_samples"],
                'value': [
                    sunRdata['sunpos_az'], sunRdata['sunpos_el'],
                    sunRdata['noise_pwr'], sunRdata['sun_maxpwr_noise'],
                    sunRdata['sun_maxpwr_nonoise'], sunRdata['dBm_sun_est'],
                    sunRdata['dBm_sun_est_toa'], sunRdata['az_bias_h'],
                    sunRdata['el_bias_h'], sunRdata['az_width_h'],
                    sunRdata['el_width_h'], sunRdata['std(dBm_sun_est)'],
                    sunRdata['nhits_h']]
                }
        elif dataset['field_name'] == 'noisedBm_vv':
            data = {
                'dstype': prdcfg['dstype'],
                'unit': 'dBm',
                'time': sunRdata['sunscan_time'],
                'label': [
                    "sun_az", "sun_el", "noise_pwr", "sun_maxpwr_noise",
                    "sun_maxpwr_nonoise", "sun_maxpwr_fit", "sun_maxpwr_toa",
                    "az_offset", "el_offset", "az_phi3db", "el_phi3db",
                    "fit_stddev", "num_samples"],
                'value': [
                    sunRdata['sunpos_az'], sunRdata['sunpos_el'],
                    sunRdata['noise_pwr'], sunRdata['sun_maxpwr_noise'],
                    sunRdata['sun_maxpwr_nonoise'], sunRdata['dBmv_sun_est'],
                    sunRdata['dBmv_sun_est_toa'], sunRdata['az_bias_v'],
                    sunRdata['el_bias_v'], sunRdata['az_width_v'],
                    sunRdata['el_width_v'], sunRdata['std(dBmv_sun_est)'],
                    sunRdata['nhits_v']]
                }
        else:
            warn('ERROR: No valid datatype for WRITE_SUNSCAN product.')

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], prdcfg['timeinfo'])

        fname1 = make_filename(
            'ts', prdcfg['dstype'], dataset['field_name'], ['csv'],
            timeinfo=prdcfg['timeinfo'], timeformat='%Y%m%d',
            runinfo=prdcfg['runinfo'])[0]

        fname1 = savedir+fname1
        write_timeseries_point(fname1, data, prdcfg['dstype'], text)

        print('saved sunscan file: {}'.format(fname1))

        return fname1

    if prdcfg['type'] == 'PLOT_SUNSCAN':
        radar = dataset['radar_out']
        sun_hits = dataset['sun_hits']

        field_name = dataset['field_name']
        if field_name not in radar.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined parameters
        azi_res = prdcfg.get('azi_res', None)
        ele_res = prdcfg.get('ele_res', None)
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)
        angtol = prdcfg.get('ang_tol', 0.5)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], prdcfg['timeinfo'])

        fname_list = make_filename(
            'constr', prdcfg['dstype'], prdcfg['dsname'],
            prdcfg['imgformat'],
            prdcfginfo='rng'+'{:.1f}'.format(
                dataset['radar_out'].range['data'][0]),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_fixed_rng_sun(
            radar, field_name, sun_hits, prdcfg, fname_list, azi_res=None,
            ele_res=None, ang_tol=angtol, vmin=vmin, vmax=vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if 'radar_out' in dataset:
        return generate_vol_products(dataset, prdcfg)

    return None


def generate_qvp_products(dataset, prdcfg):
    """
    Generates quasi vertical profile-like products. Quasi vertical profiles
    come from azimuthal averaging of polarimetric radar data. With the
    variable 'qvp_type' the user decides if the product has to be generated
    at the end of the processing period ('final') or instantaneously
    ('instant')
    Accepted product types:
        All the products of the 'VOL' dataset group

    Parameters
    ----------
    dataset : dict
        dictionary containing the radar object and a keyword stating the
        status of the processing

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    filename : str
        the name of the file created. None otherwise

    """
    qvp_type = 'final'
    if 'qvp_type' in prdcfg:
        qvp_type = prdcfg['qvp_type']

    if qvp_type == 'final' and dataset['radar_type'] != 'final':
        return None

    prdcfg['timeinfo'] = dataset['start_time']
    return generate_vol_products(dataset, prdcfg)


def generate_ml_products(dataset, prdcfg):
    """
    Generates melting layer products. Accepted product types:
        'ML_TS': Plots and writes a time series of the melting layer, i.e.
            the evolution of the average and standard deviation of the melting
            layer top and thickness and the the number of rays used in the
            retrieval.
            User defined parameters:
                dpi: int
                    The pixel density of the plot. Default 72
        'SAVE_ML': Saves an object containing the melting layer retrieval
            information in a C/F radial file
        All the products of the 'VOL' dataset group

    Parameters
    ----------
    dataset : dict
        dictionary containing the radar object and a keyword stating the
        status of the processing

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    filename : str
        the name of the file created. None otherwise

    """

    dssavedir = prdcfg['dsname']
    if 'dssavename' in prdcfg:
        dssavedir = prdcfg['dssavename']

    if prdcfg['type'] == 'ML_TS':
        dpi = prdcfg.get('dpi', 72)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        csvfname = make_filename(
            'ts', prdcfg['dstype'], 'ml', ['csv'],
            timeinfo=prdcfg['timeinfo'], timeformat='%Y%m%d')[0]

        csvfname = savedir+csvfname

        ml_bottom = dataset['ml_obj'].fields['melting_layer_height']['data'][
            :, 0]
        ml_top = dataset['ml_obj'].fields['melting_layer_height']['data'][:, 1]

        ml_top_avg = np.ma.asarray(np.ma.mean(ml_top))
        ml_top_std = np.ma.asarray(np.ma.std(ml_top))
        thick = ml_top-ml_bottom
        thick_avg = np.ma.asarray(np.ma.mean(thick))
        thick_std = np.ma.asarray(np.ma.std(thick))
        nrays_valid = thick.compressed().size
        nrays_total = thick.size

        write_ts_ml(
            prdcfg['timeinfo'], ml_top_avg, ml_top_std, thick_avg, thick_std,
            nrays_valid, nrays_total, csvfname)

        print('saved CSV file: {}'.format(csvfname))

        (dt_ml_arr, ml_top_avg_arr, ml_top_std_arr, thick_avg_arr,
         thick_std_arr, nrays_valid_arr, nrays_total_arr) = (
             read_ml_ts(csvfname))

        if dt_ml_arr is None:
            warn(
                'Unable to plot time series. No valid data')
            return None

        figfname_list = make_filename(
            'ts', prdcfg['dstype'], 'ml', prdcfg['imgformat'],
            timeinfo=dt_ml_arr[0], timeformat='%Y%m%d')

        for i, figfname in enumerate(figfname_list):
            figfname_list[i] = savedir+figfname

        titl = dt_ml_arr[0].strftime('%Y-%m-%d')+' melting layer time series'

        plot_ml_ts(
            dt_ml_arr, ml_top_avg_arr, ml_top_std_arr, thick_avg_arr,
            thick_std_arr, nrays_valid_arr, nrays_total_arr, figfname_list,
            labelx='Time UTC', titl=titl, dpi=dpi)
        print('----- save to '+' '.join(figfname_list))

        return figfname_list

    if prdcfg['type'] == 'SAVE_ML':
        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'saveml', prdcfg['dstype'], 'ml_h', ['nc'],
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname
        pyart.io.cfradial.write_cfradial(fname, dataset['ml_obj'])
        print('saved file: {}'.format(fname))

        return fname

    return generate_vol_products(dataset, prdcfg)


def generate_centroids_products(dataset, prdcfg):
    """
    Generates centroids products. Accepted product types:
        'HISTOGRAM': Plots the histogram of one of the variables used for
            centroids computation.
            User defined parameters:
                voltype : str
                    The name of the variable to plot. Can be dBZ, ZDR, KDP,
                    RhoHV, H_ISO0 and its standardized form (e.g. dBZ_std)
                write_data : Bool
                    If true writes the histogram in a .csv file. Default True
                step : float
                    bin size. Default 0.1
        'HISTOGRAM2D': Plots the 2D- histogram of two of the variables used
            for centroids computation.
            User defined parameters:
                voltype_x, voltype_y : str
                    The name of the variables to plot. Can be dBZ, ZDR, KDP,
                    RhoHV, H_ISO0 and its standardized form (e.g. dBZ_std)
                step_x, step_y : float
                    bin size. Default 0.1
        'HISTOGRAM_LABELED': Plots the histogram of one of the variables used
            for centroids computation. Only plots labeled data.
            User defined parameters:
                voltype : str
                    The name of the variable to plot. Can be dBZ, ZDR, KDP,
                    RhoHV, H_ISO0 and its standardized form (e.g. dBZ_std)
                write_data : Bool
                    If true writes the histogram in a .csv file. Default True
                step : float
                    bin size. Default 0.1
        'HISTOGRAM_CENTROIDS': Plots the histogram of one of the variables
            used for centroids computation corresponding to a particular
            hydrometeor type, the intermediate centroids and the final
            centroid
            User defined parameters:
                voltype : str
                    The name of the variable to plot. Can be dBZ, ZDR, KDP,
                    RhoHV, H_ISO0 and its standardized form (e.g. dBZ_std)
                hydro_type : str
                    The name of the hydrometeor type.
                write_data : Bool
                    If true writes the histogram in a .csv file. Default True
                step : float
                    bin size. Default 0.1
        'HISTOGRAM2D_CENTROIDS': Plots the 2D- histogram of two of the
            variables used for centroids computatio ncorresponding to a
            particular hydrometeor type, the intermediate centroids and the
            final centroid
            User defined parameters:
                voltype_x, voltype_y : str
                    The name of the variables to plot. Can be dBZ, ZDR, KDP,
                    RhoHV, H_ISO0 and its standardized form (e.g. dBZ_std)
                hydro_type : str
                    The name of the hydrometeor type.
                step_x, step_y : float
                    bin size. Default 0.1
        'WRITE_CENTROIDS': Writes the final centroids in a .csv file.
        'SAVE_DATA': Saves the data used to compute the centroids in an .npz
            file
        'SAVE_LABELED_DATA': Saves the labeled data, the intermediate
            centroids and the final centroids in an .npz file


    Parameters
    ----------
    dataset : dict
        dictionary containing the radar object and a keyword stating the
        status of the processing

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    filename : str
        the name of the file created. None otherwise

    """
    dssavedir = prdcfg['dsname']
    if 'dssavename' in prdcfg:
        dssavedir = prdcfg['dssavename']

    if prdcfg['type'] == 'HISTOGRAM':
        data_dict = dataset['data_dict']

        # check if we have to plot the data
        hist_type = prdcfg.get('hist_type', 'cumulative')
        if hist_type == 'cumulative' and not data_dict['final']:
            return None
        if hist_type == 'instant' and data_dict['final']:
            return None

        voltype = prdcfg['voltype']
        if voltype not in data_dict.keys():
            warn(
                ' Field type ' + voltype +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        write_data = prdcfg.get('write_data', 1)
        step = prdcfg.get('step', 0.1)

        timeformat = '%Y%m%d'
        timeinfo = data_dict['timeinfo'][0]
        titl = timeinfo.strftime('%Y-%m-%d')+'\n'+voltype
        data_vals = data_dict[voltype]
        if hist_type == 'instant':
            timeformat = '%Y%m%d%H%M%S'
            timeinfo = data_dict['timeinfo'][-1]
            titl = timeinfo.strftime('%Y-%m-%d %H:%M:%S')+'\n'+voltype
            nvols = data_dict['npoints'].size
            ind_start = np.sum(data_dict['npoints'][0:nvols-2])
            data_vals = data_vals[ind_start:]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=timeinfo)

        fname_list = make_filename(
            'histogram', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=timeinfo, timeformat=timeformat)

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        if '_std' in voltype:
            bin_edges = np.arange(-1.-step/2., 1.+step/2+step, step)
            values = data_vals
        else:
            field_name = get_fieldname_pyart(voltype)
            bin_edges, values = compute_histogram(
                data_vals, field_name, step=step)
        bin_centers = bin_edges[1:]-step/2
        hist, bin_edges = np.histogram(values, bins=bin_edges)
        plot_histogram2(
            bin_centers, hist, fname_list, labelx=voltype,
            labely='Number of Samples', titl=titl)

        print('----- save to '+' '.join(fname_list))

        if write_data:
            fname = savedir+make_filename(
                'histogram', prdcfg['dstype'], prdcfg['voltype'],
                ['csv'], timeinfo=timeinfo, timeformat=timeformat)[0]

            write_histogram(bin_edges, hist, fname, step=step)
            print('----- save to {}'.format(fname))

            return fname

        return fname_list

    if prdcfg['type'] == 'HISTOGRAM2D':
        data_dict = dataset['data_dict']

        # check if we have to plot the data
        hist_type = prdcfg.get('hist_type', 'cumulative')
        if hist_type == 'cumulative' and not data_dict['final']:
            return None
        if hist_type == 'instant' and data_dict['final']:
            return None

        voltype_x = prdcfg['voltype_x']
        voltype_y = prdcfg['voltype_y']
        if voltype_x not in data_dict.keys():
            warn(
                ' Field type ' + voltype_x +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None
        if voltype_y not in data_dict.keys():
            warn(
                ' Field type ' + voltype_y +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        step_x = prdcfg.get('step_x', 0.1)
        step_y = prdcfg.get('step_y', 0.1)

        timeformat = '%Y%m%d'
        timeinfo = data_dict['timeinfo'][0]
        titl = (
            timeinfo.strftime('%Y-%m-%d')+'\n'+voltype_x+'-'+voltype_y)
        data_vals_x = data_dict[voltype_x]
        data_vals_y = data_dict[voltype_y]
        if hist_type == 'instant':
            timeformat = '%Y%m%d%H%M%S'
            timeinfo = data_dict['timeinfo'][-1]
            titl = (
                timeinfo.strftime('%Y-%m-%d %H:%M:%S')+'\n'+voltype_x+'-' +
                voltype_y)
            nvols = data_dict['npoints'].size
            ind_start = np.sum(data_dict['npoints'][0:nvols-2])
            data_vals_x = data_vals_x[ind_start:]
            data_vals_y = data_vals_y[ind_start:]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=timeinfo)

        fname_list = make_filename(
            '2Dhistogram', prdcfg['dstype'], voltype_x+'-'+voltype_y,
            prdcfg['imgformat'],
            timeinfo=timeinfo, timeformat=timeformat)

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        if '_std' in voltype_x:
            bin_edges_x = np.arange(-1.-step_x/2., 1.+step_x/2+step_x, step_x)
            values_x = data_vals_x
        else:
            field_name_x = get_fieldname_pyart(voltype_x)
            bin_edges_x, values_x = compute_histogram(
                data_vals_x, field_name_x, step=step_x)

        if '_std' in voltype_y:
            bin_edges_y = np.arange(-1.-step_y/2., 1.+step_y/2+step_y, step_y)
            values_y = data_vals_y
        else:
            field_name_y = get_fieldname_pyart(voltype_y)
            bin_edges_y, values_y = compute_histogram(
                data_vals_y, field_name_y, step=step_y)

        hist_2d, bin_edges_x, bin_edges_y = np.histogram2d(
            values_x, values_y, bins=[bin_edges_x, bin_edges_y])

        plot_scatter(
            bin_edges_x, bin_edges_y, np.ma.asarray(hist_2d), voltype_x,
            voltype_y, fname_list, prdcfg, rad1_name='', rad2_name='',
            titl=titl, cmap='viridis')

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'HISTOGRAM_LABELED':
        if 'labeled_data_dict' not in dataset:
            return None

        data_dict = dataset['labeled_data_dict']

        voltype = prdcfg['voltype']
        if voltype not in data_dict.keys():
            warn(
                ' Field type ' + voltype +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        write_data = prdcfg.get('write_data', 1)
        step = prdcfg.get('step', 0.1)

        timeformat = '%Y%m%d'
        timeinfo = data_dict['timeinfo'][0]
        titl = timeinfo.strftime('%Y-%m-%d')+'\n'+voltype
        data_vals = data_dict[voltype]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=timeinfo)

        fname_list = make_filename(
            'histogram', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=timeinfo, timeformat=timeformat)

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        if '_std' in voltype:
            bin_edges = np.arange(-1.-step/2., 1.+step/2+step, step)
            values = data_vals
        else:
            field_name = get_fieldname_pyart(voltype)
            bin_edges, values = compute_histogram(
                data_vals, field_name, step=step)
        bin_centers = bin_edges[1:]-step/2
        hist, bin_edges = np.histogram(values, bins=bin_edges)
        plot_histogram2(
            bin_centers, hist, fname_list, labelx=voltype,
            labely='Number of Samples', titl=titl)

        print('----- save to '+' '.join(fname_list))

        if write_data:
            fname = savedir+make_filename(
                'histogram', prdcfg['dstype'], prdcfg['voltype'],
                ['csv'], timeinfo=timeinfo, timeformat=timeformat)[0]

            write_histogram(bin_edges, hist, fname, step=step)
            print('----- save to {}'.format(fname))

            return fname

        return fname_list

    if prdcfg['type'] == 'HISTOGRAM_CENTROIDS':
        if 'labeled_data_dict' not in dataset:
            return None

        data_dict = dataset['labeled_data_dict']

        voltype = prdcfg['voltype']
        hydro_type = prdcfg['hydro_type']
        if voltype not in data_dict.keys():
            warn(
                ' Field type ' + voltype +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None
        if hydro_type not in data_dict['medoids_dict']:
            warn('No medoids where found for hydrometeor class '+hydro_type)
            return None
        ind_hydro = np.where(
            np.array(data_dict['hydro_names']) == hydro_type)[0]
        ind_medoid = np.where(
            np.array(data_dict['var_names']) == voltype)[0]

        write_data = prdcfg.get('write_data', 1)
        step = prdcfg.get('step', 0.1)
        dpi = 72
        if 'dpi' in prdcfg['ppiImageConfig']:
            dpi = prdcfg['ppiImageConfig']['dpi']

        timeformat = '%Y%m%d'
        timeinfo = data_dict['timeinfo'][0]
        titl = timeinfo.strftime('%Y-%m-%d')+'\n'+voltype+' '+hydro_type
        data_vals = data_dict[voltype]
        data_vals = data_vals[data_dict['labels'] == ind_hydro]

        medoids = np.array(data_dict['medoids_dict'][hydro_type])
        medoids = medoids[:, ind_medoid]

        if hydro_type not in data_dict['final_medoids_dict']:
            warn('No medoid for hydrometeor class '+hydro_type)
            fmedoid = None
        else:
            fmedoid = data_dict['final_medoids_dict'][hydro_type][ind_medoid]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=timeinfo)

        fname_list = make_filename(
            'histogram', prdcfg['dstype'], hydro_type+_+prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=timeinfo, timeformat=timeformat)

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        if '_std' in voltype:
            bin_edges = np.arange(-1.-step/2., 1.+step/2+step, step)
            values = data_vals
        else:
            field_name = get_fieldname_pyart(voltype)
            bin_edges, values = compute_histogram(
                data_vals, field_name, step=step)
        bin_centers = bin_edges[1:]-step/2
        hist, bin_edges = np.histogram(values, bins=bin_edges)
        pos_medoids = []
        for medoid in medoids:
            ind = np.where(bin_edges <= medoid)[0][-1]
            pos_medoids.append(hist[ind])
        pos_medoids = np.array(pos_medoids)
        if fmedoid is not None:
            ind = np.where(bin_edges <= fmedoid)[0][-1]
            pos_fmedoid = hist[ind]

        fig, ax, = plot_histogram2(
            bin_centers, hist, fname_list, labelx=voltype,
            labely='Number of Samples', titl=titl, save_fig=False)

        ax.scatter(medoids, pos_medoids, c='g', marker='o')
        if fmedoid is not None:
            ax.plot(fmedoid, pos_fmedoid, c='r', marker='D')

        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close(fig)

        print('----- save to '+' '.join(fname_list))

        if write_data:
            fname = savedir+make_filename(
                'histogram', prdcfg['dstype'], prdcfg['voltype'],
                ['csv'], timeinfo=timeinfo, timeformat=timeformat)[0]

            write_histogram(bin_edges, hist, fname, step=step)
            print('----- save to {}'.format(fname))

            return fname

        return fname_list

    if prdcfg['type'] == 'HISTOGRAM2D_CENTROIDS':
        if 'labeled_data_dict' not in dataset:
            return None

        labeled_data_dict = dataset['labeled_data_dict']

        voltype_x = prdcfg['voltype_x']
        voltype_y = prdcfg['voltype_y']
        hydro_type = prdcfg['hydro_type']
        if voltype_x not in labeled_data_dict.keys():
            warn(
                ' Field type ' + voltype_x +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None
        if voltype_y not in labeled_data_dict.keys():
            warn(
                ' Field type ' + voltype_y +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None
        if hydro_type not in labeled_data_dict['medoids_dict']:
            warn('No medoids where found for hydrometeor class '+hydro_type)
            return None

        ind_hydro = np.where(
            np.array(labeled_data_dict['hydro_names']) == hydro_type)[0]
        ind_medoid_x = np.where(
            np.array(labeled_data_dict['var_names']) == voltype_x)[0]
        ind_medoid_y = np.where(
            np.array(labeled_data_dict['var_names']) == voltype_y)[0]

        step_x = prdcfg.get('step_x', 0.1)
        step_y = prdcfg.get('step_y', 0.1)

        timeformat = '%Y%m%d'
        timeinfo = labeled_data_dict['timeinfo'][0]
        titl = (
            timeinfo.strftime('%Y-%m-%d')+'\n'+voltype_x+'-'+voltype_y +
            ' '+hydro_type)

        data_vals_x = labeled_data_dict[voltype_x]
        data_vals_y = labeled_data_dict[voltype_y]

        data_vals_x = data_vals_x[labeled_data_dict['labels'] == ind_hydro]
        data_vals_y = data_vals_y[labeled_data_dict['labels'] == ind_hydro]

        medoids = np.array(labeled_data_dict['medoids_dict'][hydro_type])
        medoids_x = medoids[:, ind_medoid_x]
        medoids_y = medoids[:, ind_medoid_y]

        if hydro_type not in labeled_data_dict['final_medoids_dict']:
            warn('No medoid for hydrometeor class '+hydro_type)
            fmedoid_x = None
            fmedoid_y = None
        else:
            fmedoid_x = labeled_data_dict['final_medoids_dict'][hydro_type][
                ind_medoid_x]
            fmedoid_y = labeled_data_dict['final_medoids_dict'][hydro_type][
                ind_medoid_y]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=timeinfo)

        fname_list = make_filename(
            '2Dhistogram', prdcfg['dstype'],
            hydro_type+'_'+voltype_x+'-'+voltype_y, prdcfg['imgformat'],
            timeinfo=timeinfo, timeformat=timeformat)

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        if '_std' in voltype_x:
            bin_edges_x = np.arange(-1.-step_x/2., 1.+step_x/2+step_x, step_x)
            values_x = data_vals_x
        else:
            field_name_x = get_fieldname_pyart(voltype_x)
            bin_edges_x, values_x = compute_histogram(
                data_vals_x, field_name_x, step=step_x)

        if '_std' in voltype_y:
            bin_edges_y = np.arange(-1.-step_y/2., 1.+step_y/2+step_y, step_y)
            values_y = data_vals_y
        else:
            field_name_y = get_fieldname_pyart(voltype_y)
            bin_edges_y, values_y = compute_histogram(
                data_vals_y, field_name_y, step=step_y)

        hist_2d, bin_edges_x, bin_edges_y = np.histogram2d(
            values_x, values_y, bins=[bin_edges_x, bin_edges_y])

        plot_centroids(
            bin_edges_x, bin_edges_y, np.ma.asarray(hist_2d), voltype_x,
            voltype_y, fname_list, prdcfg, titl=titl, medoids_x=medoids_x,
            medoids_y=medoids_y, fmedoid_x=fmedoid_x, fmedoid_y=fmedoid_y,
            cmap='viridis')

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'WRITE_CENTROIDS':
        if 'labeled_data_dict' not in dataset:
            return None

        labeled_data_dict = dataset['labeled_data_dict']
        timeformat = '%Y%m%d'
        timeinfo = labeled_data_dict['timeinfo'][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=timeinfo)

        fname = make_filename(
            'centroids', prdcfg['dstype'], 'centroids',
            ['csv'], timeinfo=timeinfo, timeformat=timeformat)[0]

        fname = savedir+fname

        write_centroids(
            fname, labeled_data_dict['final_medoids_dict'],
            labeled_data_dict['var_names'])
        print('----- save to {}'.format(fname))

        return fname

    if prdcfg['type'] == 'SAVE_DATA':
        data_dict = dataset['data_dict']

        # check if we have to save the data
        if not data_dict['final']:
            return None

        timeinfo = data_dict['timeinfo'][0]
        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=timeinfo)

        fname = make_filename(
            'data', prdcfg['dstype'], 'centroids_data', ['npz'],
            timeinfo=timeinfo)[0]

        fname = savedir+fname

        np.savez_compressed(fname, centroids_data=data_dict)

        print('----- save to {}'.format(fname))

        return fname

    if prdcfg['type'] == 'SAVE_LABELED_DATA':
        if 'labeled_data_dict' not in dataset:
            return None

        data_dict = dataset['labeled_data_dict']

        timeinfo = data_dict['timeinfo'][0]
        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=timeinfo)

        fname = make_filename(
            'labeled_data', prdcfg['dstype'], 'centroids_data', ['npz'],
            timeinfo=timeinfo)[0]

        fname = savedir+fname

        np.savez_compressed(fname, centroids_data=data_dict)

        print('----- save to {}'.format(fname))

        return fname

    warn(' Unsupported product type: ' + prdcfg['type'])
    return None
