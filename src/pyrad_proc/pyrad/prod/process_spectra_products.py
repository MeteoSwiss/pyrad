"""
pyrad.prod.process_spectra_products
===================================

Functions for obtaining Pyrad products from spectra datasets

.. autosummary::
    :toctree: generated/

    generate_spectra_products

"""

from warnings import warn
from copy import deepcopy

import numpy as np

import pyart
from pyart.util import datetime_from_radar

from ..io.io_aux import get_fieldname_pyart
from ..io.io_aux import get_save_dir, make_filename

from ..graph.plots_spectra import plot_range_Doppler, plot_Doppler
from ..graph.plots_spectra import plot_complex_range_Doppler
from ..graph.plots_spectra import plot_amp_phase_range_Doppler
from ..graph.plots_spectra import plot_complex_Doppler, plot_amp_phase_Doppler
from ..graph.plots_spectra import plot_time_Doppler, plot_complex_time_Doppler
from ..graph.plots_spectra import plot_amp_phase_time_Doppler
from ..graph.plots_spectra import plot_angle_Doppler
from ..graph.plots_spectra import plot_complex_angle_Doppler
from ..graph.plots_spectra import plot_amp_phase_angle_Doppler

from ..util.radar_utils import find_ray_index, find_rng_index


def generate_spectra_products(dataset, prdcfg):
    """
    generates spectra products. Accepted product types:
        'AMPLITUDE_PHASE_ANGLE_DOPPLER': Makes an angle Doppler plot of
            complex spectra or IQ data. The plot can be along azimuth or along
            range. It is plotted separately the module and the phase of the
            signal.
            User defined parameters:
                along_azi : bool
                    If true the plot is performed along azimuth, otherwise
                    along elevation. Default true
                ang : float
                    The fixed angle (deg). Default 0.
                rng : float
                    The fixed range (m). Default 0.
                ang_tol : float
                    The fixed angle tolerance (deg). Default 1.
                rng_tol : float
                    The fixed range tolerance (m). Default 50.
                xaxis_info : str
                    The xaxis type. Can be 'Doppler_velocity',
                    'Doppler_frequency' or 'pulse_number'
                ampli_vmin, ampli_vmax, phase_vmin, phase_vmax : float or None
                    Minimum and maximum of the color scale for the module and
                    phase
        'AMPLITUDE_PHASE_DOPPLER': Plots a complex Doppler spectrum or IQ data
            making two separate plots for the module and phase of the signal
            User defined parameters:
                azi, ele, rng : float
                    azimuth and elevation (deg) and range (m) of the ray to
                    plot
                azi_to, ele_tol, rng_tol : float
                    azimuth and elevation (deg) and range (m) tolerance
                    respect to nominal position to plot. Default 1, 1, 50.
                ind_ray, ind_rng : int
                    index of the ray and range to plot. Alternative to
                    defining its antenna coordinates
                xaxis_info : str
                    The xaxis type. Can be 'Doppler_velocity',
                    'Doppler_frequency' or 'pulse_number'
                ampli_vmin, ampli_vmax, phase_vmin, phase_vmax : float or None
                    Minimum and maximum of the color scale for the module and
                    phase
        'AMPLITUDE_PHASE_RANGE_DOPPLER': Plots a complex spectra or IQ data
            range-Doppler making two separate plots for the module and phase
            of the signal User defined parameters:
                azi, ele : float
                    azimuth and elevation (deg) of the ray to plot
                azi_to, ele_tol : float
                    azimuth and elevation (deg) tolerance respect to nominal
                    position to plot. Default 1, 1.
                ind_ray : int
                    index of the ray to plot. Alternative to
                    defining its antenna coordinates
                xaxis_info : str
                    The xaxis type. Can be 'Doppler_velocity',
                    'Doppler_frequency' or 'pulse_number'
                ampli_vmin, ampli_vmax, phase_vmin, phase_vmax : float or None
                    Minimum and maximum of the color scale for the module and
                    phase
        'AMPLITUDE_PHASE_TIME_DOPPLER': Plots a complex spectra or IQ data
            time-Doppler making two separate plots for the module and phase of
            the signal
            User defined parameters:
                xaxis_info : str
                    The xaxis type. Can be 'Doppler_velocity' or
                    'Doppler frequency'
                ampli_vmin, ampli_vmax, phase_vmin, phase_vmax : float or None
                    Minimum and maximum of the color scale for the module and
                    phase
                plot_type : str
                    Can be 'final' or 'temporal'. If final the data is only
                    plotted at the end of the processing
        'ANGLE_DOPPLER': Makes an angle Doppler plot. The plot can be along
            azimuth or along range
            User defined parameters:
                along_azi : bool
                    If true the plot is performed along azimuth, otherwise
                    along elevation. Default true
                ang : float
                    The fixed angle (deg). Default 0.
                rng : float
                    The fixed range (m). Default 0.
                ang_tol : float
                    The fixed angle tolerance (deg). Default 1.
                rng_tol : float
                    The fixed range tolerance (m). Default 50.
                xaxis_info : str
                    The xaxis type. Can be 'Doppler_velocity',
                    'Doppler_frequency' or 'pulse_number'
                vmin, vmax : float or None
                    Minimum and maximum of the color scale
        'COMPLEX_ANGLE_DOPPLER': Makes an angle Doppler plot of complex
            spectra or IQ data. The plot can be along azimuth or along range.
            The real and imaginary parts are plotted separately
            User defined parameters:
                along_azi : bool
                    If true the plot is performed along azimuth, otherwise
                    along elevation. Default true
                ang : float
                    The fixed angle (deg). Default 0.
                rng : float
                    The fixed range (m). Default 0.
                ang_tol : float
                    The fixed angle tolerance (deg). Default 1.
                rng_tol : float
                    The fixed range tolerance (m). Default 50.
                xaxis_info : str
                    The xaxis type. Can be 'Doppler_velocity',
                    'Doppler_frequency' or 'pulse_number'
                vmin, vmax : float or None
                    Minimum and maximum of the color scale
        'COMPLEX_DOPPLER': Plots a complex Doppler spectrum or IQ data making
            two separate plots for the real and imaginary parts
            User defined parameters:
                azi, ele, rng : float
                    azimuth and elevation (deg) and range (m) of the ray to
                    plot
                azi_to, ele_tol, rng_tol : float
                    azimuth and elevation (deg) and range (m) tolerance
                    respect to nominal position to plot. Default 1, 1, 50.
                ind_ray, ind_rng : int
                    index of the ray and range to plot. Alternative to
                    defining its antenna coordinates
                xaxis_info : str
                    The xaxis type. Can be 'Doppler_velocity',
                    'Doppler_frequency' or 'pulse_number'
                vmin, vmax : float or None
                    Minimum and maximum of the color scale
        'COMPLEX_RANGE_DOPPLER': Plots the complex spectra or IQ data
            range-Doppler making two separate plots for the real and imaginary
            parts
            User defined parameters:
                azi, ele : float
                    azimuth and elevation (deg) of the ray to plot
                azi_to, ele_tol : float
                    azimuth and elevation (deg) tolerance respect to nominal
                    position to plot. Default 1, 1.
                ind_ray : int
                    index of the ray to plot. Alternative to
                    defining its antenna coordinates
                xaxis_info : str
                    The xaxis type. Can be 'Doppler_velocity',
                    'Doppler_frequency' or 'pulse_number'
                vmin, vmax : float or None
                    Minimum and maximum of the color scale
        'COMPLEX_TIME_DOPPLER': Plots the complex spectra or IQ data
            time-Doppler making two separate plots for the real and imaginary
            parts
            User defined parameters:
                xaxis_info : str
                    The xaxis type. Can be 'Doppler_velocity' or
                    'Doppler frequency'
                vmin, vmax : float or None
                    Minimum and maximum of the color scale
                plot_type : str
                    Can be 'final' or 'temporal'. If final the data is only
                    plotted at the end of the processing
        'DOPPLER': Plots a Doppler spectrum variable or IQ data variable
            User defined parameters:
                azi, ele, rng : float
                    azimuth and elevation (deg) and range (m) of the ray to
                    plot
                azi_to, ele_tol, rng_tol : float
                    azimuth and elevation (deg) and range (m) tolerance
                    respect to nominal position to plot. Default 1, 1, 50.
                ind_ray, ind_rng : int
                    index of the ray and range to plot. Alternative to
                    defining its antenna coordinates
                xaxis_info : str
                    The xaxis type. Can be 'Doppler_velocity',
                    'Doppler_frequency' or 'pulse_number'
                vmin, vmax : float or None
                    Minimum and maximum of the color scale
        'RANGE_DOPPLER': Makes a range-Doppler plot of spectral or IQ data
            User defined parameters:
                azi, ele : float
                    azimuth and elevation (deg) of the ray to plot
                azi_to, ele_tol : float
                    azimuth and elevation (deg) tolerance respect to nominal
                    position to plot. Default 1, 1.
                ind_ray : int
                    index of the ray to plot. Alternative to
                    defining its antenna coordinates
                xaxis_info : str
                    The xaxis type. Can be 'Doppler_velocity',
                    'Doppler_frequency' or 'pulse_number'
                vmin, vmax : float or None
                    Minimum and maximum of the color scale
        'SAVEALL': Saves radar spectra or IQ volume data including all or a
            list of userdefined fields in a netcdf file
            User defined parameters:
                datatypes: list of str or None
                    The list of data types to save. If it is None, all fields
                    in the radar object will be saved
                physical: Bool
                    If True the data will be saved in physical units (floats).
                    Otherwise it will be quantized and saved as binary
        'SAVEVOL': Saves one field of a radar spectra or IQ volume data in a
            netcdf file
            User defined parameters:
                physical: Bool
                    If True the data will be saved in physical units (floats).
                    Otherwise it will be quantized and saved as binary
        'TIME_DOPPLER': Makes a time-Doppler plot of spectral or IQ data at a
            point of interest.
            User defined parameters:
                xaxis_info : str
                    The xaxis type. Can be 'Doppler_velocity',
                    'Doppler_frequency' or 'pulse_number'
                vmin, vmax : float or None
                    Minimum and maximum of the color scale
                plot_type : str
                    Can be 'final' or 'temporal'. If final the data is only
                    plotted at the end of the processing

    Parameters
    ----------
    dataset : spectra
        spectra object

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    None or name of generated files

    """
    dssavedir = prdcfg['dsname']
    if 'dssavename' in prdcfg:
        dssavedir = prdcfg['dssavename']

    if prdcfg['type'] == 'RANGE_DOPPLER':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        azi = prdcfg.get('azi', None)
        ele = prdcfg.get('ele', None)
        azi_tol = prdcfg.get('azi_tol', 1.)
        ele_tol = prdcfg.get('ele_tol', 1.)

        if azi is None or ele is None:
            ind_ray = prdcfg.get('ind_ray', 0)
            azi = dataset['radar_out'].azimuth['data'][ind_ray]
            ele = dataset['radar_out'].elevation['data'][ind_ray]
        else:
            ind_ray = find_ray_index(
                dataset['radar_out'].elevation['data'],
                dataset['radar_out'].azimuth['data'], ele, azi,
                ele_tol=ele_tol, azi_tol=azi_tol)

        if ind_ray is None:
            warn('Ray azi='+str(azi)+', ele='+str(ele) +
                 ' out of radar coverage')
            return None

        gateinfo = 'az'+'{:.1f}'.format(azi)+'el'+'{:.1f}'.format(ele)

        xaxis_info = prdcfg.get('xaxis_info', 'Doppler_velocity')
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'range_Doppler', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        if dataset['radar_out'].ngates == 1:
            plot_Doppler(
                dataset['radar_out'], field_name, ind_ray, 0, prdcfg,
                fname_list, xaxis_info=xaxis_info, vmin=vmin, vmax=vmax)
        else:
            plot_range_Doppler(
                dataset['radar_out'], field_name, ind_ray, prdcfg, fname_list,
                xaxis_info=xaxis_info, vmin=vmin, vmax=vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'ANGLE_DOPPLER':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        along_azi = prdcfg.get('along_azi', True)
        ang = prdcfg.get('ang', 0)
        rng = prdcfg.get('rng', 0)
        ang_tol = prdcfg.get('ang_tol', 1.)
        rng_tol = prdcfg.get('rng_tol', 50.)

        ind_rng = find_rng_index(
            dataset['radar_out'].range['data'], rng, rng_tol=rng_tol)

        if ind_rng is None:
            warn('No data at rng='+str(rng))
            return None

        if along_azi:
            ind_rays = np.where(np.logical_and(
                dataset['radar_out'].elevation['data'] <= ang+ang_tol,
                dataset['radar_out'].elevation['data'] >= ang-ang_tol))[0]
        else:
            ind_rays = np.where(np.logical_and(
                dataset['radar_out'].azimuth['data'] <= ang+ang_tol,
                dataset['radar_out'].azimuth['data'] >= ang-ang_tol))[0]

        if ind_rays.size == 0:
            warn('No data for angle '+str(ang))
            return None

        # sort angles
        if along_azi:
            ang_selected = dataset['radar_out'].azimuth['data'][ind_rays]

        else:
            ang_selected = dataset['radar_out'].elevation['data'][ind_rays]
        ind_rays = ind_rays[np.argsort(ang_selected)]

        if along_azi:
            gateinfo = 'azi'+'{:.1f}'.format(ang)+'rng'+'{:.1f}'.format(rng)
        else:
            gateinfo = 'ele'+'{:.1f}'.format(ang)+'rng'+'{:.1f}'.format(rng)

        xaxis_info = prdcfg.get('xaxis_info', 'Doppler_velocity')
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'range_Doppler', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        if ind_rays.size == 1:
            plot_Doppler(
                dataset['radar_out'], field_name, ind_rays, ind_rng, prdcfg,
                fname_list, xaxis_info=xaxis_info, vmin=vmin, vmax=vmax)
        else:
            plot_angle_Doppler(
                dataset['radar_out'], field_name, ang, ind_rays, ind_rng,
                prdcfg, fname_list, xaxis_info=xaxis_info,
                along_azi=along_azi, vmin=vmin, vmax=vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'TIME_DOPPLER':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        xaxis_info = prdcfg.get('xaxis_info', 'Doppler_velocity')
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)
        xmin = prdcfg.get('xmin', None)
        xmax = prdcfg.get('xmax', None)
        ymin = prdcfg.get('ymin', None)
        ymax = prdcfg.get('ymax', None)
        plot_type = prdcfg.get('plot_type', 'final')

        if plot_type == 'final' and not dataset['final']:
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

        time_info = datetime_from_radar(dataset['radar_out'])

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=time_info)

        fname_list = make_filename(
            'time_Doppler', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=time_info, runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        if dataset['radar_out'].nrays == 1:
            plot_Doppler(
                dataset['radar_out'], field_name, 0, 0, prdcfg, fname_list,
                xaxis_info=xaxis_info, vmin=vmin, vmax=vmax)
        else:
            plot_time_Doppler(
                dataset['radar_out'], field_name, prdcfg, fname_list,
                xaxis_info=xaxis_info, vmin=vmin, vmax=vmax, xmin=xmin,
                xmax=xmax, ymin=ymin, ymax=ymax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'DOPPLER':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        azi = prdcfg.get('azi', None)
        ele = prdcfg.get('ele', None)
        rng = prdcfg.get('rng', None)
        azi_tol = prdcfg.get('azi_tol', 1.)
        ele_tol = prdcfg.get('ele_tol', 1.)
        rng_tol = prdcfg.get('rng_tol', 50.)

        if azi is None or ele is None or rng is None:
            ind_ray = prdcfg.get('ind_ray', 0)
            ind_rng = prdcfg.get('ind_rng', 0)
            azi = dataset['radar_out'].azimuth['data'][ind_ray]
            ele = dataset['radar_out'].elevation['data'][ind_ray]
            rng = dataset['radar_out'].range['data'][ind_rng]
        else:
            ind_ray = find_ray_index(
                dataset['radar_out'].elevation['data'],
                dataset['radar_out'].azimuth['data'], ele, azi,
                ele_tol=ele_tol, azi_tol=azi_tol)
            ind_rng = find_rng_index(
                dataset['radar_out'].range['data'], rng, rng_tol=rng_tol)

        if ind_rng is None or ind_ray is None:
            warn('Point azi='+str(azi)+', ele='+str(ele)+', rng='+str(rng) +
                 ' out of radar coverage')
            return None

        gateinfo = (
            'az'+'{:.1f}'.format(azi)+'el'+'{:.1f}'.format(ele) +
            'r'+'{:.1f}'.format(rng))

        xaxis_info = prdcfg.get('xaxis_info', 'Doppler_velocity')
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'Doppler', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_Doppler(
            dataset['radar_out'], field_name, ind_ray, ind_rng, prdcfg,
            fname_list, xaxis_info=xaxis_info, vmin=vmin, vmax=vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'COMPLEX_RANGE_DOPPLER':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        azi = prdcfg.get('azi', None)
        ele = prdcfg.get('ele', None)
        azi_tol = prdcfg.get('azi_tol', 1.)
        ele_tol = prdcfg.get('ele_tol', 1.)

        if azi is None or ele is None:
            ind_ray = prdcfg.get('ind_ray', 0)
            azi = dataset['radar_out'].azimuth['data'][ind_ray]
            ele = dataset['radar_out'].elevation['data'][ind_ray]
        else:
            ind_ray = find_ray_index(
                dataset['radar_out'].elevation['data'],
                dataset['radar_out'].azimuth['data'], ele, azi,
                ele_tol=ele_tol, azi_tol=azi_tol)

        if ind_ray is None:
            warn('Ray azi='+str(azi)+', ele='+str(ele) +
                 ' out of radar coverage')
            return None

        gateinfo = 'az'+'{:.1f}'.format(azi)+'el'+'{:.1f}'.format(ele)

        xaxis_info = prdcfg.get('xaxis_info', 'Doppler_velocity')
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'c_range_Doppler', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        if dataset['radar_out'].ngates == 1:
            plot_complex_Doppler(
                dataset['radar_out'], field_name, ind_ray, 0, prdcfg,
                fname_list, xaxis_info=xaxis_info, vmin=vmin, vmax=vmax)
        else:
            plot_complex_range_Doppler(
                dataset['radar_out'], field_name, ind_ray, prdcfg, fname_list,
                xaxis_info=xaxis_info, vmin=vmin, vmax=vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'COMPLEX_ANGLE_DOPPLER':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        along_azi = prdcfg.get('along_azi', True)
        ang = prdcfg.get('ang', 0)
        rng = prdcfg.get('rng', 0)
        ang_tol = prdcfg.get('ang_tol', 1.)
        rng_tol = prdcfg.get('rng_tol', 50.)

        ind_rng = find_rng_index(
            dataset['radar_out'].range['data'], rng, rng_tol=rng_tol)

        if ind_rng is None:
            warn('No data at rng='+str(rng))
            return None

        if along_azi:
            ind_rays = np.where(np.logical_and(
                dataset['radar_out'].elevation['data'] <= ang+ang_tol,
                dataset['radar_out'].elevation['data'] >= ang-ang_tol))[0]
        else:
            ind_rays = np.where(np.logical_and(
                dataset['radar_out'].azimuth['data'] <= ang+ang_tol,
                dataset['radar_out'].azimuth['data'] >= ang-ang_tol))[0]

        if ind_rays.size == 0:
            warn('No data for angle '+str(ang))
            return None

        # sort angles
        if along_azi:
            ang_selected = dataset['radar_out'].azimuth['data'][ind_rays]

        else:
            ang_selected = dataset['radar_out'].elevation['data'][ind_rays]
        ind_rays = ind_rays[np.argsort(ang_selected)]

        if along_azi:
            gateinfo = 'azi'+'{:.1f}'.format(ang)+'rng'+'{:.1f}'.format(rng)
        else:
            gateinfo = 'ele'+'{:.1f}'.format(ang)+'rng'+'{:.1f}'.format(rng)

        xaxis_info = prdcfg.get('xaxis_info', 'Doppler_velocity')
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'range_Doppler', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        if ind_rays.size == 1:
            plot_complex_Doppler(
                dataset['radar_out'], field_name, ind_rays, ind_rng, prdcfg,
                fname_list, xaxis_info=xaxis_info, vmin=vmin, vmax=vmax)
        else:
            plot_complex_angle_Doppler(
                dataset['radar_out'], field_name, ang, ind_rays, ind_rng,
                prdcfg, fname_list, xaxis_info=xaxis_info,
                along_azi=along_azi, vmin=vmin, vmax=vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'COMPLEX_TIME_DOPPLER':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        xaxis_info = prdcfg.get('xaxis_info', 'Doppler_velocity')
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)
        plot_type = prdcfg.get('plot_type', 'final')

        if plot_type == 'final' and not dataset['final']:
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

        time_info = datetime_from_radar(dataset['radar_out'])

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=time_info)

        fname_list = make_filename(
            'c_time_Doppler', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=time_info, runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        if dataset['radar_out'].nrays == 1:
            plot_complex_Doppler(
                dataset['radar_out'], field_name, 0, 0, prdcfg, fname_list,
                xaxis_info=xaxis_info, vmin=vmin, vmax=vmax)
        else:
            plot_complex_time_Doppler(
                dataset['radar_out'], field_name, prdcfg, fname_list,
                xaxis_info=xaxis_info, vmin=vmin, vmax=vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'COMPLEX_DOPPLER':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        azi = prdcfg.get('azi', None)
        ele = prdcfg.get('ele', None)
        rng = prdcfg.get('rng', None)
        azi_tol = prdcfg.get('azi_tol', 1.)
        ele_tol = prdcfg.get('ele_tol', 1.)
        rng_tol = prdcfg.get('rng_tol', 50.)

        if azi is None or ele is None or rng is None:
            ind_ray = prdcfg.get('ind_ray', 0)
            ind_rng = prdcfg.get('ind_rng', 0)
            azi = dataset['radar_out'].azimuth['data'][ind_ray]
            ele = dataset['radar_out'].elevation['data'][ind_ray]
            rng = dataset['radar_out'].range['data'][ind_rng]
        else:
            ind_ray = find_ray_index(
                dataset['radar_out'].elevation['data'],
                dataset['radar_out'].azimuth['data'], ele, azi,
                ele_tol=ele_tol, azi_tol=azi_tol)
            ind_rng = find_rng_index(
                dataset['radar_out'].range['data'], rng, rng_tol=rng_tol)

        if ind_rng is None or ind_ray is None:
            warn('Point azi='+str(azi)+', ele='+str(ele)+', rng='+str(rng) +
                 ' out of radar coverage')
            return None

        gateinfo = (
            'az'+'{:.1f}'.format(azi)+'el'+'{:.1f}'.format(ele) +
            'r'+'{:.1f}'.format(rng))

        xaxis_info = prdcfg.get('xaxis_info', 'Doppler_velocity')
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'c_Doppler', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_complex_Doppler(
            dataset['radar_out'], field_name, ind_ray, ind_rng, prdcfg,
            fname_list, xaxis_info=xaxis_info, vmin=vmin, vmax=vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'AMPLITUDE_PHASE_DOPPLER':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        azi = prdcfg.get('azi', None)
        ele = prdcfg.get('ele', None)
        rng = prdcfg.get('rng', None)
        azi_tol = prdcfg.get('azi_tol', 1.)
        ele_tol = prdcfg.get('ele_tol', 1.)
        rng_tol = prdcfg.get('rng_tol', 50.)

        if azi is None or ele is None or rng is None:
            ind_ray = prdcfg.get('ind_ray', 0)
            ind_rng = prdcfg.get('ind_rng', 0)
            azi = dataset['radar_out'].azimuth['data'][ind_ray]
            ele = dataset['radar_out'].elevation['data'][ind_ray]
            rng = dataset['radar_out'].range['data'][ind_rng]
        else:
            ind_ray = find_ray_index(
                dataset['radar_out'].elevation['data'],
                dataset['radar_out'].azimuth['data'], ele, azi,
                ele_tol=ele_tol, azi_tol=azi_tol)
            ind_rng = find_rng_index(
                dataset['radar_out'].range['data'], rng, rng_tol=rng_tol)

        if ind_rng is None or ind_ray is None:
            warn('Point azi='+str(azi)+', ele='+str(ele)+', rng='+str(rng) +
                 ' out of radar coverage')
            return None

        gateinfo = (
            'az'+'{:.1f}'.format(azi)+'el'+'{:.1f}'.format(ele) +
            'r'+'{:.1f}'.format(rng))

        xaxis_info = prdcfg.get('xaxis_info', 'Doppler_velocity')
        ampli_vmin = prdcfg.get('ampli_vmin', None)
        ampli_vmax = prdcfg.get('ampli_vmax', None)
        phase_vmin = prdcfg.get('phase_vmin', None)
        phase_vmax = prdcfg.get('phase_vmax', None)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'ap_Doppler', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_amp_phase_Doppler(
            dataset['radar_out'], field_name, ind_ray, ind_rng, prdcfg,
            fname_list, xaxis_info=xaxis_info, ampli_vmin=ampli_vmin,
            ampli_vmax=ampli_vmax, phase_vmin=phase_vmin,
            phase_vmax=phase_vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'AMPLITUDE_PHASE_RANGE_DOPPLER':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        azi = prdcfg.get('azi', None)
        ele = prdcfg.get('ele', None)
        azi_tol = prdcfg.get('azi_tol', 1.)
        ele_tol = prdcfg.get('ele_tol', 1.)

        if azi is None or ele is None:
            ind_ray = prdcfg.get('ind_ray', 0)
            azi = dataset['radar_out'].azimuth['data'][ind_ray]
            ele = dataset['radar_out'].elevation['data'][ind_ray]
        else:
            ind_ray = find_ray_index(
                dataset['radar_out'].elevation['data'],
                dataset['radar_out'].azimuth['data'], ele, azi,
                ele_tol=ele_tol, azi_tol=azi_tol)

        if ind_ray is None:
            warn('Ray azi='+str(azi)+', ele='+str(ele) +
                 ' out of radar coverage')
            return None

        gateinfo = 'az'+'{:.1f}'.format(azi)+'el'+'{:.1f}'.format(ele)

        xaxis_info = prdcfg.get('xaxis_info', 'Doppler_velocity')
        ampli_vmin = prdcfg.get('ampli_vmin', None)
        ampli_vmax = prdcfg.get('ampli_vmax', None)
        phase_vmin = prdcfg.get('phase_vmin', None)
        phase_vmax = prdcfg.get('phase_vmax', None)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'ap_range_Doppler', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        if dataset['radar_out'].ngates == 1:
            plot_amp_phase_Doppler(
                dataset['radar_out'], field_name, ind_ray, 0, prdcfg,
                fname_list, xaxis_info=xaxis_info, ampli_vmin=ampli_vmin,
                ampli_vmax=ampli_vmax, phase_vmin=phase_vmin,
                phase_vmax=phase_vmax)
        else:
            plot_amp_phase_range_Doppler(
                dataset['radar_out'], field_name, ind_ray, prdcfg, fname_list,
                xaxis_info=xaxis_info, ampli_vmin=ampli_vmin,
                ampli_vmax=ampli_vmax, phase_vmin=phase_vmin,
                phase_vmax=phase_vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'AMPLITUDE_PHASE_ANGLE_DOPPLER':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        along_azi = prdcfg.get('along_azi', True)
        ang = prdcfg.get('ang', 0)
        rng = prdcfg.get('rng', 0)
        ang_tol = prdcfg.get('ang_tol', 1.)
        rng_tol = prdcfg.get('rng_tol', 50.)

        ind_rng = find_rng_index(
            dataset['radar_out'].range['data'], rng, rng_tol=rng_tol)

        if ind_rng is None:
            warn('No data at rng='+str(rng))
            return None

        if along_azi:
            ind_rays = np.where(np.logical_and(
                dataset['radar_out'].elevation['data'] <= ang+ang_tol,
                dataset['radar_out'].elevation['data'] >= ang-ang_tol))[0]
        else:
            ind_rays = np.where(np.logical_and(
                dataset['radar_out'].azimuth['data'] <= ang+ang_tol,
                dataset['radar_out'].azimuth['data'] >= ang-ang_tol))[0]

        if ind_rays.size == 0:
            warn('No data for angle '+str(ang))
            return None

        # sort angles
        if along_azi:
            ang_selected = dataset['radar_out'].azimuth['data'][ind_rays]

        else:
            ang_selected = dataset['radar_out'].elevation['data'][ind_rays]
        ind_rays = ind_rays[np.argsort(ang_selected)]

        if along_azi:
            gateinfo = 'azi'+'{:.1f}'.format(ang)+'rng'+'{:.1f}'.format(rng)
        else:
            gateinfo = 'ele'+'{:.1f}'.format(ang)+'rng'+'{:.1f}'.format(rng)

        xaxis_info = prdcfg.get('xaxis_info', 'Doppler_velocity')
        ampli_vmin = prdcfg.get('ampli_vmin', None)
        ampli_vmax = prdcfg.get('ampli_vmax', None)
        phase_vmin = prdcfg.get('phase_vmin', None)
        phase_vmax = prdcfg.get('phase_vmax', None)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'range_Doppler', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        if ind_rays.size == 1:
            plot_amp_phase_Doppler(
                dataset['radar_out'], field_name, ind_rays, ind_rng, prdcfg,
                fname_list, xaxis_info=xaxis_info, ampli_vmin=ampli_vmin,
                ampli_vmax=ampli_vmax, phase_vmin=phase_vmin,
                phase_vmax=phase_vmax)
        else:
            plot_amp_phase_angle_Doppler(
                dataset['radar_out'], field_name, ang, ind_rays, ind_rng,
                prdcfg, fname_list, xaxis_info=xaxis_info,
                along_azi=along_azi, ampli_vmin=ampli_vmin,
                ampli_vmax=ampli_vmax, phase_vmin=phase_vmin,
                phase_vmax=phase_vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'AMPLITUDE_PHASE_TIME_DOPPLER':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        xaxis_info = prdcfg.get('xaxis_info', 'Doppler_velocity')
        ampli_vmin = prdcfg.get('ampli_vmin', None)
        ampli_vmax = prdcfg.get('ampli_vmax', None)
        phase_vmin = prdcfg.get('phase_vmin', None)
        phase_vmax = prdcfg.get('phase_vmax', None)
        plot_type = prdcfg.get('plot_type', 'final')

        if plot_type == 'final' and not dataset['final']:
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

        time_info = datetime_from_radar(dataset['radar_out'])

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=time_info)

        fname_list = make_filename(
            'ap_time_Doppler', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=time_info, runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        if dataset['radar_out'].nrays == 1:
            plot_amp_phase_Doppler(
                dataset['radar_out'], field_name, 0, 0, prdcfg, fname_list,
                xaxis_info=xaxis_info, ampli_vmin=ampli_vmin,
                ampli_vmax=ampli_vmax, phase_vmin=phase_vmin,
                phase_vmax=phase_vmax)
        else:
            plot_amp_phase_time_Doppler(
                dataset['radar_out'], field_name, prdcfg, fname_list,
                xaxis_info=xaxis_info, ampli_vmin=ampli_vmin,
                ampli_vmax=ampli_vmax, phase_vmin=phase_vmin,
                phase_vmax=phase_vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'SAVEVOL':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        file_type = prdcfg.get('file_type', 'nc')
        physical = prdcfg.get('physical', True)

        new_dataset = deepcopy(dataset['radar_out'])
        new_dataset.fields = dict()
        new_dataset.add_field(
            field_name, dataset['radar_out'].fields[field_name])

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'savevol', prdcfg['dstype'], prdcfg['voltype'], [file_type],
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname

        pyart.aux_io.write_spectra(fname, new_dataset)

        print('saved file: '+fname)

        return fname

    if prdcfg['type'] == 'SAVEALL':
        file_type = prdcfg.get('file_type', 'nc')
        datatypes = prdcfg.get('datatypes', None)
        physical = prdcfg.get('physical', True)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'savevol', prdcfg['dstype'], 'all_fields', [file_type],
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname

        field_names = None
        if datatypes is not None:
            field_names = []
            for datatype in datatypes:
                field_names.append(get_fieldname_pyart(datatype))

        if field_names is not None:
            radar_aux = deepcopy(dataset['radar_out'])
            radar_aux.fields = dict()
            for field_name in field_names:
                if field_name not in dataset['radar_out'].fields:
                    warn(field_name+' not in radar object')
                else:
                    radar_aux.add_field(
                        field_name,
                        dataset['radar_out'].fields[field_name])
        else:
            radar_aux = dataset['radar_out']
        pyart.aux_io.write_spectra(fname, radar_aux, )

        print('saved file: '+fname)

        return fname

    warn(' Unsupported product type: ' + prdcfg['type'])
    return None
