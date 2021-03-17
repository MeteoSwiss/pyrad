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

from ..io.write_data import write_sun_hits, write_sun_retrieval, write_timeseries_point
from ..io.write_data import write_excess_gates, write_ts_ml

from ..graph.plots import plot_sun_hits
from ..graph.plots_timeseries import plot_sun_retrieval_ts, plot_ml_ts
from ..graph.plots_vol import plot_fixed_rng, plot_fixed_rng_sun

from ..util.radar_utils import create_sun_hits_field
from ..util.radar_utils import create_sun_retrieval_field


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
            print('saved excess gates file: '+fname)

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

        print('saved file: '+savedir+fname)

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

        print('saved file: '+savedir+fname)

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

        print('saved sun hits file: '+fname)

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

        print('saved sun retrieval file: '+fname)

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

        text = ["SunScan info",
                "sun_az:             [deg] Azimuth sun position ",
                "sun_el:             [deg] Elevation sun position",
                "noise_pwr:          [dBm] Noise power",
                "sun_maxpwr_noise:   [dBm] sun maximal power sample (including noise)",
                "sun_maxpwr_nonoise: [dBm] sun maximal power sample without noise",
                "sun_maxpwr_fit:     [dBm] sun maximal fitted power (without noise)",
                "sun_maxpwr_toa:     [dBm] sun maximal power at top of atmosphere",
                "az_offset:          [deg] Azimuth shift of fitted maxima to sun azimuth",
                "el_offset:          [deg] Elevation shift of fitted maxima to sun elevation",
                "az_phi3db:          [deg] Half-power beam width in azimuth",
                "el_phi3db:          [deg] Half-power beam width in elevation",
                "fit_stddev:         [dBm] Standard deviation (fit to samples)",
                "num_samples:        [#] Number of samples used for the sun power fitting"
                ]

        sunRdata = dataset['sun_retrieval']

        if dataset['field_name'] == 'noisedBm_hh':
            data = {'dstype': prdcfg['dstype'],
                    'unit': 'dBm',
                    'time': sunRdata['sunscan_time'],
                    'label': ["sun_az", "sun_el", "noise_pwr",
                              "sun_maxpwr_noise", "sun_maxpwr_nonoise", "sun_maxpwr_fit",
                              "sun_maxpwr_toa", "az_offset", "el_offset",
                              "az_phi3db", "el_phi3db", "fit_stddev",
                              "num_samples"],
                    'value': [sunRdata['sunpos_az'], sunRdata['sunpos_el'],
                              sunRdata['noise_pwr'], sunRdata['sun_maxpwr_noise'],
                              sunRdata['sun_maxpwr_nonoise'], sunRdata['dBm_sun_est'],
                              sunRdata['dBm_sun_est_toa'], sunRdata['az_bias_h'],
                              sunRdata['el_bias_h'], sunRdata['az_width_h'],
                              sunRdata['el_width_h'], sunRdata['std(dBm_sun_est)'],
                              sunRdata['nhits_h']]
                    }
        elif dataset['field_name'] == 'noisedBm_vv':
            data = {'dstype': prdcfg['dstype'],
                    'unit': 'dBm',
                    'time': sunRdata['sunscan_time'],
                    'label': ["sun_az", "sun_el", "noise_pwr",
                              "sun_maxpwr_noise", "sun_maxpwr_nonoise", "sun_maxpwr_fit",
                              "sun_maxpwr_toa", "az_offset", "el_offset",
                              "az_phi3db", "el_phi3db", "fit_stddev",
                              "num_samples"],
                    'value': [sunRdata['sunpos_az'], sunRdata['sunpos_el'],
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
            'ts', prdcfg['dstype'], dataset['field_name'], ['csv'], timeinfo=prdcfg['timeinfo'],
            timeformat='%Y%m%d', runinfo=prdcfg['runinfo'])[0]

        fname1 = savedir+fname1
        write_timeseries_point(fname1, data, prdcfg['dstype'], text)

        print('saved sunscan file: ' +fname1)

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

        plot_fixed_rng_sun(radar, field_name, sun_hits, prdcfg, fname_list, azi_res=None,
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

        ml_bottom = dataset['ml_obj'].fields['melting_layer_height']['data'][:, 0]
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

        print('saved CSV file: '+csvfname)

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
        print('saved file: '+fname)

        return fname

    return generate_vol_products(dataset, prdcfg)
