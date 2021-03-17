"""
pyrad.io.read_data_hzt
========================

Functions for reading HZT data

.. autosummary::
    :toctree: generated/

    hzt2radar_data
    hzt2radar_coord
    get_iso0_field
    read_hzt_data
    _prepare_for_interpolation

"""
from warnings import warn
import datetime
import platform
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from scipy.spatial import cKDTree

import pyart
from pyart.config import get_metadata, get_field_name

# check existence of METRANET library
try:
    METRANET_LIB = pyart.aux_io.get_library(momentms=False)
    if platform.system() == 'Linux':
        METRANET_LIB = pyart.aux_io.get_library(momentms=True)
    _METRANETLIB_AVAILABLE = True
except:
    _METRANETLIB_AVAILABLE = False

from ..io.read_data_cosmo import _put_radar_in_swiss_coord


def hzt2radar_data(radar, hzt_coord, hzt_data, slice_xy=True,
                   field_name='height_over_iso0'):
    """
    get the HZT value corresponding to each radar gate using nearest
    neighbour interpolation

    Parameters
    ----------
    radar : Radar
        the radar object containing the information on the position of the
        radar gates
    hzt_coord : dict
        dictionary containing the HZT coordinates
    hzt_data : dict
        dictionary containing the HZT data
    slice_xy : boolean
        if true the horizontal plane of the COSMO field is cut to the
        dimensions of the radar field
    field_name : str
        name of HZT fields to convert (default height_over_iso0)

    Returns
    -------
    hzt_fields : list of dict
        list of dictionary with the HZT fields and metadata

    """
    x_radar, y_radar, z_radar = _put_radar_in_swiss_coord(radar)

    x_hzt, y_hzt, ind_xmin, ind_ymin, ind_xmax, ind_ymax = (
        _prepare_for_interpolation(
            x_radar, y_radar, hzt_coord, slice_xy=slice_xy))

    values = hzt_data['HZT']['data'][
        ind_ymin:ind_ymax+1, ind_xmin:ind_xmax+1].flatten()
    # find interpolation function
    interp_func = NearestNDInterpolator((y_hzt, x_hzt), values)

    # interpolate
    data_interp = interp_func((y_radar, x_radar))

    # put field
    field_dict = get_metadata(field_name)
    field_dict['data'] = (z_radar-data_interp).astype(float)

    return field_dict


def hzt2radar_coord(radar, hzt_coord, slice_xy=True, field_name=None):
    """
    Given the radar coordinates find the nearest HZT pixel

    Parameters
    ----------
    radar : Radar
        the radar object containing the information on the position of the
        radar gates
    hzt_coord : dict
        dictionary containing the HZT coordinates
    slice_xy : boolean
        if true the horizontal plane of the HZT field is cut to the
        dimensions of the radar field
    field_name : str
        name of the field

    Returns
    -------
    hzt_ind_field : dict
        dictionary containing a field of HZT indices and metadata

    """
    # parse the field parameters
    if field_name is None:
        field_name = get_field_name('hzt_index')

    x_radar, y_radar, _ = _put_radar_in_swiss_coord(radar)

    x_hzt, y_hzt, ind_xmin, ind_ymin, ind_xmax, _ = (
        _prepare_for_interpolation(
            x_radar, y_radar, hzt_coord, slice_xy=slice_xy))

    tree = cKDTree(np.transpose((y_hzt, x_hzt)))
    _, ind_vec = tree.query(
        np.transpose((y_radar.flatten(), x_radar.flatten())), k=1)

    # put the index in the original cosmo coordinates
    nx_hzt = len(hzt_coord['x']['data'])

    nx = ind_xmax-ind_xmin+1

    ind_y = (ind_vec/nx).astype(int)+ind_ymin
    ind_x = (ind_vec % nx).astype(int)+ind_xmin
    ind_hzt = (ind_x+nx_hzt*ind_y).astype(int)

    hzt_ind_field = get_metadata(field_name)
    hzt_ind_field['data'] = ind_hzt.reshape(radar.nrays, radar.ngates)

    return hzt_ind_field


def get_iso0_field(hzt_data, hzt_ind, z_radar, field_name='height_over_iso0'):
    """
    Get the height over iso0 data corresponding to each radar gate
    using a precomputed look up table of the nearest neighbour

    Parameters
    ----------
    hzt_data : dict
        dictionary containing the HZT data and metadata
    hzt_ind : dict
        dictionary containing a field of HZT indices and metadata
    z_radar : ndarray
        gates altitude [m MSL]
    field_name : str
        names of HZT parameters (default height_over_iso0)

    Returns
    -------
    iso0_field : list of dict
        dictionary with the height over iso0 field and metadata

    """
    nrays, ngates = np.shape(hzt_ind['data'])
    values = hzt_data['HZT']['data'][:, :].flatten()
    field_dict = get_metadata(field_name)
    field_dict['data'] = z_radar-values[hzt_ind['data'].flatten()].reshape(
        nrays, ngates).astype(float)

    return field_dict


def read_hzt_data(fname, chy0=255., chx0=-160., read_lib='C'):
    """
    Reads iso-0 degree data from an HZT file

    Parameters
    ----------
    fname : str
        name of the file to read
    chy0, chx0: float
        south west point of grid in Swiss coordinates [km]
    read_lib : str
        Type of METRANET read library used. Can be 'C' or 'python'

    Returns
    -------
    hzt_data : dictionary
        dictionary with the data and metadata

    """
    if read_lib == 'C' and _METRANETLIB_AVAILABLE:
        ret = pyart.aux_io.read_product_c(
            fname, physic_value=True, masked_array=True)
    elif read_lib == 'python':
        ret = pyart.aux_io.read_product_py(
            fname, physic_value=True, masked_array=True)
    else:
        warn('METRANET C-library reader not available or unknown ' +
             'library type. Python library will be used')
        ret = pyart.aux_io.read_product_py(
            fname, physic_value=True, masked_array=True)

    if ret is None:
        warn('Unable to read HZT file '+fname)
        return None

    var_data = {
        'units': 'meters_above_mean_sea_level',
        'data': ret.data,
        'long_name': 'height_over_iso0'
    }
    run_time = datetime.datetime.strptime(ret.header['time'], '%y%j%H%M0')
    time_data = {
        'standard_name': 'time',
        'long_name': 'time',
        'units': 'seconds since '+run_time.strftime('%Y-%m-%d %H:%M:%S'),
        'calendar': 'gregorian',
        'data': [float(ret.header['usr_forecast_hour'])*3600.]
    }

    x_1 = {
        'axis': "X",
        'long_name': "x-coordinate in Swiss coordinate system",
        'standard_name': "projection_x_coordinate",
        'units': "m",
        'data': ((np.arange(int(ret.header['column'])) *
                  float(ret.header['rect_xres'])+chy0 +
                  float(ret.header['rect_xres'])/2.)*1000.)
    }
    y_1 = {
        'axis': "Y",
        'long_name': "y-coordinate in Swiss coordinate system",
        'standard_name': "projection_y_coordinate",
        'units': "m",
        'data': ((np.arange(int(ret.header['row'])) *
                  float(ret.header['rect_yres'])+chx0 +
                  float(ret.header['rect_yres'])/2.)*1000.)
    }

    hzt_data = {
        'HZT': var_data,
        'metadata': ret.header,
        'time': time_data,
        'x': x_1,
        'y': y_1
    }

    return hzt_data


def _prepare_for_interpolation(x_radar, y_radar, hzt_coord, slice_xy=True):
    """
    prepares the HZT 2D volume for interpolation:
        1. if set slices the cosmo data to the area covered by the radar
        2. creates the x, y grid for the interpolation

    Parameters
    ----------
    x_radar, y_radar : arrays
        The Swiss coordinates of the radar
    hzt_coord : dict
        dictionary containing the HZT coordinates
    slice_xy : boolean
        if true the horizontal plane of the HZT field is cut to the
        dimensions of the radar field

    Returns
    -------
    x_hzt, y_hzt : 1D arrays
        arrays containing the flatten swiss coordinates of the HZT data in
        the area of interest [m]
    ind_xmin, ind_ymin, ind_xmax, ind_ymax : ints
        the minimum and maximum indices of each dimension

    """
    nx_hzt = len(hzt_coord['x']['data'])
    ny_hzt = len(hzt_coord['y']['data'])

    if slice_xy:
        # get the COSMO data within the radar range
        xmin = np.min(x_radar)
        xmax = np.max(x_radar)
        ymin = np.min(y_radar)
        ymax = np.max(y_radar)

        ind_xmin = np.where(hzt_coord['x']['data'] < xmin)[0]
        if ind_xmin.size == 0:
            ind_xmin = 0
        else:
            ind_xmin = ind_xmin[-1]

        ind_xmax = np.where(hzt_coord['x']['data'] > xmax)[0]
        if ind_xmax.size == 0:
            ind_xmax = nx_hzt-1
        else:
            ind_xmax = ind_xmax[0]

        ind_ymin = np.where(hzt_coord['y']['data'] < ymin)[0]
        if ind_ymin.size == 0:
            ind_ymin = 0
        else:
            ind_ymin = ind_ymin[-1]

        ind_ymax = np.where(hzt_coord['y']['data'] > ymax)[0]
        if ind_ymax.size == 0:
            ind_ymax = ny_hzt-1
        else:
            ind_ymax = ind_ymax[0]
    else:
        ind_xmin = 0
        ind_xmax = nx_hzt-1
        ind_ymin = 0
        ind_ymax = ny_hzt-1

    nx = ind_xmax-ind_xmin+1
    ny = ind_ymax-ind_ymin+1

    x_hzt = hzt_coord['x']['data'][ind_xmin:ind_xmax+1]
    y_hzt = hzt_coord['y']['data'][ind_ymin:ind_ymax+1]

    x_hzt = (
        np.broadcast_to(x_hzt.reshape(1, nx), (ny, nx))).flatten()
    y_hzt = (
        np.broadcast_to(y_hzt.reshape(ny, 1), (ny, nx))).flatten()

    return x_hzt, y_hzt, ind_xmin, ind_ymin, ind_xmax, ind_ymax
