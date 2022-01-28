"""
pyrad.io.read_data_mxpol
========================
Functions for reading radar mxpol data files
.. autosummary::
    :toctree: generated/
    classes - MXPOL:
        pyrad_MXPOL
    classes - MCH:
        pyrad_MCH
    utilities - read:
        row_stack
        findTimes
        int2float_radar
        readMXPOLRadData
        readCHRadData
    utilities - config:
        load_myconfig
        get_mymetadata
        get_elevation_metadata
        generate_radar_table
        generate_polvar_metadata
        convert_polvar_name
"""

import time
import importlib.util
import re
import os
import datetime
from copy import deepcopy
import warnings

import numpy as np
import netCDF4
try:
    import h5py
    _H5PY_AVAILABLE = True
except ImportError:
    _H5PY_AVAILABLE = False

class MissingOptionalDependency(Exception):
    """ Exception raised when a optional dependency is needed but not found. """
    pass

import pyart

# -------------------------- classes - MXPOL ------------------------------ #


class pyrad_MXPOL(pyart.core.Radar):
    def __init__(self, filename, field_names=None, max_range=np.Inf,
                 min_range=10000, pyrad_names=True):
        # find information based on filename
        all_files = [filename]
        fname_basename = os.path.basename(filename)

        if 'PPI' in fname_basename:
            scan_type = 'ppi'
        elif 'RHI' in fname_basename:
            scan_type = 'rhi'

        strdate = re.findall(r"([0-9]{8}-[0-9]{6})", fname_basename)[0]
        date = datetime.datetime.strptime(strdate, '%Y%m%d-%H%M%S')

        # if field name is None, take all available fields

        if field_names is None:
            field_names = ['Zh', 'Zdr', 'Kdp', 'Phidp', 'Rhohv', 'ZhCorr',
                           'ZdrCorr', 'RVel', 'Rvel', 'Sw', 'SNRh', 'SNRv', 'Psidp']

        # convert fieldname if necessary
        varnames = []
        for fieldname in field_names:
            newname = convert_polvar_name('LTE', fieldname)
            varnames.append(newname)

        # get labels, units etc
        long_names = []
        standard_names = []
        units = []
        vmin = []
        vmax = []

        for varname in varnames:
            metadata = generate_polvar_metadata(varname)
            standard_names.append(metadata['standard_name'])
            long_names.append(metadata['long_name'])
            units.append(metadata['units'])
            vmin.append(metadata['valid_min'])
            vmax.append(metadata['valid_max'])

        # initiate empty vectors
        N_sweeps = len(all_files)
        fields = {}
        fixed_angle = {}
        fixed_angle['data'] = np.zeros(N_sweeps, )

        sweep_start_ray_index = {}
        sweep_start_ray_index['data'] = []
        sweep_stop_ray_index = {}
        sweep_stop_ray_index['data'] = []

        for i, k in enumerate(varnames):
            fields[k] = {}
            fields[k]['data'] = []
            fields[k]['long_name'] = long_names[i]
            fields[k]['standard_name'] = standard_names[i]
            fields[k]['units'] = units[i]
            fields[k]['valid_min'] = vmin[i]
            fields[k]['valid_max'] = vmax[i]

        idx_start = 0
        idx_stop = 0
        elevations = []
        azimuths = []
        ranges = []
        nyquist = []

        # read data and create dictionaries
        for i in range(N_sweeps):
            metadata, data = readMXPOLRadData(
                all_files[i], varnames, max_range)
            if scan_type == 'rhi':
                fixed_angle['data'][i] = np.round(np.mean(data['azimuth']))
            elif scan_type == 'ppi':
                fixed_angle['data'][i] = np.round(np.mean(data['elevation']))

            [N_az, N_ranges] = data[varnames[0]].shape
            idx_stop = idx_start + N_az - 1
            sweep_start_ray_index['data'].append(idx_start)
            sweep_stop_ray_index['data'].append(idx_stop)
            idx_start = idx_stop + 1
            elevations.extend(list(data['elevation']))
            nyquist.extend([data['nyquist_vel']]*N_az)
            azimuths.extend(list(data['azimuth']))
            ranges.extend(list(data['range']))

            for j, v in enumerate(varnames):
                if v in data.keys():
                    if not(len(fields[v]['data'])):
                        fields[v]['data'] = data[v]
                    else:
                        fields[v]['data'] = row_stack(
                            fields[v]['data'], data[v])
                else:
                    print('Variable '+v+' was not found in file!')

        # mask NaNs

        for v in varnames:
            if not len(fields[v]['data']):
                # Remove variable
                fields.pop(v)
            else:
                fields[v]['data'] = np.ma.masked_equal(
                    fields[v]['data'], -99900.0)

        [a, N_ranges] = fields[varnames[0]]['data'].shape

        # create dictionaries according to pyART standard
        latitude = {'data': np.asarray([data['latitude']]),
                    'units': data['lat_units']}
        longitude = {'data': np.asarray([data['longitude']]),
                     'units': data['lon_units']}
        altitude = {'data': np.asarray([data['altitude']]),
                    'units': data['alt_units']}
        sweep_number = {'data': np.arange(0, len(all_files))}
        sweep_mode = {'data': np.asarray([scan_type]*N_sweeps)}
        instrument_parameters = {
            'nyquist_velocity': {'data': np.asarray(nyquist)}}
        azimuth = {'data': np.asarray(azimuths), 'units': data['azim_units']}
        rrange = {'data': np.asarray(ranges),
                  'units': data['range_units']}
        elevation = {'data': np.asarray(elevations),
                     'units': data['elev_units']}
        sweep_start_ray_index['data'] = np.asarray(
            sweep_start_ray_index['data'])
        sweep_stop_ray_index['data'] = np.asarray(
            sweep_stop_ray_index['data'])

        time_units = 'seconds since ' + str(date)
        time_data = {'data': data['time'], 'units': time_units}

        # change keys to match pyART metranet keys
        if pyrad_names:
            fields_copy = deepcopy(fields)
            for keys in fields_copy:
                newkey = fields[keys]['standard_name']
                fields[newkey] = fields.pop(keys)

        # Create PyART instance
        pyart.core.Radar.__init__(
            self, time_data, rrange, fields, metadata, scan_type, latitude,
            longitude, altitude, sweep_number, sweep_mode, fixed_angle,
            sweep_start_ray_index, sweep_stop_ray_index, azimuth, elevation,
            instrument_parameters=instrument_parameters)


# -------------------------- classes - IDL --------------------------- #

class pyrad_IDL(pyart.core.Radar):
    def __init__(self, filename, field_names=None, max_range=np.Inf,
                 min_range=10000):

        # find information based on filename
        all_files = [filename]
        fname_basename = os.path.basename(filename)

        fname = netCDF4.Dataset(filename)

        if 'PPI' in fname_basename:
            scan_type = 'ppi'
        elif 'RHI' in fname_basename:
            scan_type = 'rhi'
        strdate = re.findall(r"([0-9]{8}-[0-9]{6})", fname_basename)[0]
        date = datetime.datetime.strptime(strdate, '%Y%m%d-%H%M%S')

        # if field name is None, take all available fields

        if field_names is None:
            field_names = list(fname.variables.keys())

        # convert fieldname if necessary
        varnames = []
        for fieldname in field_names:
            newname = convert_polvar_name('IDL', fieldname)
            varnames.append(newname)

        # get labels, units etc
        long_names = []
        standard_names = []
        units = []
        vmin = []
        vmax = []

        for varname in varnames:
            metadata = generate_polvar_metadata(varname)
            standard_names.append(metadata['standard_name'])
            long_names.append(metadata['long_name'])
            units.append(metadata['units'])
            vmin.append(metadata['valid_min'])
            vmax.append(metadata['valid_max'])

        # initiate empty vectors
        N_sweeps = len(all_files)
        fields = {}
        fixed_angle = {}
        fixed_angle['data'] = np.zeros(N_sweeps, )

        sweep_start_ray_index = {}
        sweep_start_ray_index['data'] = []
        sweep_stop_ray_index = {}
        sweep_stop_ray_index['data'] = []

        for i, k in enumerate(varnames):
            fields[k] = {}
            fields[k]['data'] = []
            fields[k]['long_name'] = long_names[i]
            fields[k]['standard_name'] = standard_names[i]
            fields[k]['units'] = units[i]
            fields[k]['valid_min'] = vmin[i]
            fields[k]['valid_max'] = vmax[i]

        idx_start = 0
        idx_stop = 0
        elevations = []
        azimuths = []
        nyquist = []

        # read data and create dictionaries
        for i in range(N_sweeps):
            metadata, data = readIDLRadData(
                all_files[i], varnames, max_range)
            if scan_type == 'rhi':
                fixed_angle['data'][i] = np.round(np.mean(data['azimuth']))
            elif scan_type == 'ppi':
                fixed_angle['data'][i] = data['elevation'][0]

            [N_az, N_ranges] = data[varnames[0]].shape
            idx_stop = idx_start + N_az - 1
            sweep_start_ray_index['data'].append(idx_start)
            sweep_stop_ray_index['data'].append(idx_stop)
            idx_start = idx_stop + 1
            elevations.extend([data['elevation'][0]]*N_az)
            nyquist.extend([data['nyquist_vel']]*N_az)
            azimuths.extend(list(data['azimuth']))
            warnings.warn("Warning, sweep rank could not be found, using first rank")

            starttime, endtime = findTimes(1)
            interval = ((endtime-starttime)/N_az)
            #time_lapse = np.arange(starttime+(0.5*interval), endtime, interval)
            # because this is a single sweep
            time_lapse = np.around(
                np.arange(0.+(0.5*interval), endtime-starttime, interval))

            for j, v in enumerate(varnames):
                if v in data.keys():
                    if fields[v]['data'].size == 0:
                        fields[v]['data'] = data[v]
                    else:
                        fields[v]['data'] = row_stack(
                            fields[v]['data'], data[v])
                else:
                    print('Variable '+v+' was not found in file!')

        # mask NaNs

        for v in varnames:
            fields[v]['data'] = np.ma.masked_equal(
                fields[v]['data'], -99900.0)

        [a, N_ranges] = fields[varnames[0]]['data'].shape

        # create dictionaries according to pyART standard
        latitude = {'data': np.asarray([data['latitude']]),
                    'units': data['lat_units']}
        longitude = {'data': np.asarray([data['longitude']]),
                     'units': data['lon_units']}
        altitude = {'data': np.asarray([data['altitude']]),
                    'units': data['alt_units']}
        sweep_number = {'data': np.arange(0, len(all_files))}
        sweep_mode = {'data': np.asarray([scan_type]*N_sweeps)}
        instrument_parameters = {
            'nyquist_velocity': {'data': np.asarray(nyquist)}}
        azimuth = {'data': np.asarray(azimuths), 'units': data['azim_units']}
        rrange = {'data': np.arange(N_ranges)*data['resolution'],
                  'units': data['range_units']}
        elevation = {'data': np.asarray(elevations),
                     'units': data['elev_units']}
        sweep_start_ray_index['data'] = np.asarray(
            sweep_start_ray_index['data'])
        sweep_stop_ray_index['data'] = np.asarray(
            sweep_stop_ray_index['data'])

        time_units = 'seconds since ' + str(date)
        time_lapse = np.asarray(time_lapse)
        time_data = {'data': time_lapse, 'units': time_units}

        # change keys to match pyART metranet keys
        fields_copy = deepcopy(fields)
        for keys in fields_copy:
            newkey = fields[keys]['standard_name']
            fields[newkey] = fields.pop(keys)

        # Create PyART instance
        pyart.core.Radar.__init__(
            self, time_data, rrange, fields, metadata, scan_type, latitude,
            longitude, altitude, sweep_number, sweep_mode, fixed_angle,
            sweep_start_ray_index, sweep_stop_ray_index, azimuth, elevation,
            instrument_parameters=instrument_parameters)

# -------------------------- classes - MCH --------------------------- #

class pyrad_MCH(pyart.core.Radar):
    def __init__(self, filename, field_names=None, max_range=np.Inf):

        # find information based on filename
        all_files = [filename]
        N_sweeps = len(all_files)

        fname_basename = os.path.basename(filename)

        # Get name of radar
        index_letter = fname_basename[2]

        radar_info = generate_radar_table(index_letter)
        radar_name = radar_info['radarID']

        # Get radar resolution
        if fname_basename[1] == 'L':
            rres = 500.
        else:
            rres = 83.3

        scan_type = 'ppi'

        scandate = datetime.datetime.strptime(
            fname_basename[3:12], '%y%j%H%M')
        self.scan_date = scandate.timetuple()

        # if field name is None, take all available fields
        if field_names is None:
            field_names = ['Z', 'ZDR', 'ZV', 'V', 'W', 'RHO', 'CLUT', 'PHIDP']

        # convert fieldname if necessary
        varnames = []
        for fieldname in field_names:
            newname = convert_polvar_name('MCH', fieldname)
            varnames.append(newname)

        # get labels, units etc
        long_names = []
        standard_names = []
        units = []
        vmin = []
        vmax = []

        for varname in varnames:
            metadata = generate_polvar_metadata(varname)
            standard_names.append(metadata['standard_name'])
            long_names.append(metadata['long_name'])
            units.append(metadata['units'])
            vmin.append(metadata['valid_min'])
            vmax.append(metadata['valid_max'])

        # initiate empty vectors
        fields = {}
        fixed_angle = {}
        fixed_angle['data'] = np.zeros(N_sweeps, )

        sweep_start_ray_index = {}
        sweep_start_ray_index['data'] = []
        sweep_stop_ray_index = {}
        sweep_stop_ray_index['data'] = []

        for i, k in enumerate(varnames):
            fields[k] = {}
            fields[k]['data'] = []
            fields[k]['long_name'] = long_names[i]
            fields[k]['standard_name'] = standard_names[i]
            fields[k]['units'] = units[i]
            fields[k]['valid_min'] = vmin[i]
            fields[k]['valid_max'] = vmax[i]

        # Initialize
        idx_start = 0
        idx_stop = 0
        elevations = []
        azimuths = []
        nyquist = []
        time_lapse = []

        # read and organise data
        for i in range(N_sweeps):
            data = readCHRadData(
                all_files[i], radar_name, varnames, rres, max_range)
            fixed_angle['data'][i] = data['elevation']
            [N_ranges, N_az] = data[varnames[0]].shape
            idx_stop = idx_start + N_az - 1
            sweep_start_ray_index['data'].append(idx_start)
            sweep_stop_ray_index['data'].append(idx_stop)
            idx_start = idx_stop + 1
            elevations.extend([data['elevation']]*N_az)
            nyquist.extend([data['nyquist_vel']]*N_az)
            azimuths.extend(list(data['azimuth']))
            # create list of times at the center of each ray
            sweep_rank = 1
            print()
            starttime, endtime = findTimes(sweep_rank)
            interval = ((endtime-starttime)/len(list(data['azimuth'])))
            time_lapse.extend(np.arange(
                starttime+(0.5*interval), endtime, interval))
            for j, v in enumerate(varnames):
                if fields[v]['data'].size == 0:
                    fields[v]['data'] = data[v].T
                else:
                    fields[v]['data'] = row_stack(
                        fields[v]['data'], data[v].T)

        # mask nans
        for v in varnames:
            fields[v]['data'] = np.ma.array(
                fields[v]['data'], mask=np.isnan(fields[v]['data']))

        sweep_start_ray_index['data'] = np.asarray(
            sweep_start_ray_index['data'])
        sweep_stop_ray_index['data'] = np.asarray(
            sweep_stop_ray_index['data'])
        metadata = {}

        [a, N_ranges] = fields[varnames[0]]['data'].shape

        latitude = {'data': np.array([radar_info['coordinates'][0]]),
                    'units': "DegreesNorth"}
        longitude = {'data': np.array([radar_info['coordinates'][1]]),
                     'units': "DegreesEast"}
        altitude = {'data': np.array([radar_info['altitude']]),
                    'units': "MetersAboveSeaLevel"}
        sweep_number = {'data': np.arange(0, len(all_files))}
        sweep_mode = {'data': np.asarray(['ppi']*N_sweeps)}
        instrument_parameters = {
            'nyquist_velocity': {'data': np.array(nyquist)}}

        metadata['Source'] = (
            "Operational radar data processed at MeteoSwiss Locarno-Monti")
        metadata['Institution'] = (
            "MeteoSwiss, MDR, Locarno-Monti, Switzerland")
        metadata['History'] = [
            "created: %s, " % time.ctime(os.path.getctime(filename)) +
            "last modified: %s" % time.ctime(os.path.getmtime(filename))]
        metadata['ContactInformation'] = "marc.schneebeli@meteosvizzera.ch"

        azimuth = {'data': np.array(azimuths), 'units': "Degrees"}
        rrange = {'data': np.arange(N_ranges)*data['resolution'],
                  'units': "Meters"}
        elevation = {'data': np.array(elevations), 'units': "Degrees"}

        time_units = 'seconds since '+str(scandate)
        time_lapse = np.asarray(time_lapse)
        scantime = {'data': time_lapse, 'units': time_units}

        # change keys to match pyART metranet keys
        fields_copy = deepcopy(fields)
        for keys in fields_copy:
            newkey = fields[keys]['standard_name']
            fields[newkey] = fields.pop(keys)

        # Create PyART instance
        pyart.core.Radar.__init__(
            self, scantime, rrange, fields, metadata, scan_type, latitude,
            longitude, altitude, sweep_number, sweep_mode, fixed_angle,
            sweep_start_ray_index, sweep_stop_ray_index, azimuth, elevation,
            instrument_parameters=instrument_parameters)


# ----------------------- utilities - read --------------------- #

def row_stack(a1, a2):
    """
    Stacks data from subsequent sweeps, while padding "empty" columns from
    subsequent sweeps.
    Inputs
    ------
    a1: np.array
        destination array
    a2: np.array
        array which is added onto the first array
    Returns
    -------
    out: np.array
        stacked destination and additional array, with uniform shape
    """
    [N1, M1] = a1.shape
    [N2, M2] = a2.shape

    if M1 > M2:
        a2 = np.pad(a2, ((0, 0), (0, M1-M2)), mode='constant',
                    constant_values=-9999999)
    elif M2 < M1:
        a1 = np.pad(a2, ((0, 0), (0, M2-M1)), mode='constant',
                    constant_values=-9999999)

    out = np.vstack((a1, a2))
    out[out == -9999999] = np.nan

    return out


def findTimes(num_sweep):
    """
    Finds the times at the beginning and at the end of each sweep. Information
    comes from the elapsed time since the beginning of the volume scan, from
    the Rad4Alp: Specifications/ Request for Proposal (RFP) document.
    Inputs
    ------
    num_sweep: int
        rank of the sweep
    Returns
    -------
    elapsed_times[num_sweep][0]: float
        the elapsed time since the beginning of the volume scan at the
        beginning of the sweep
    elapsed_times[num_sweep][1]: float
        the elapsed time since the beginning of the volume scan at the end of
        the sweep
    """

    elapsed_times = {9: [0, 11.4],
                     7: [11.4, 22.8],
                     5: [22.8, 39.2],
                     3: [39.3, 60.5],
                     1: [60.5, 84.7],
                     19: [84.7, 97.2],
                     17: [97.2, 109.6],
                     15: [109.6, 121.6],
                     13: [121.6, 133.1],
                     11: [133.1, 144.4],
                     10: [144.4, 155.8],
                     8: [155.8, 172.2],
                     6: [172.2, 188.6],
                     4: [188.6, 204.9],
                     2: [204.9, 229.4],
                     20: [229.4, 241.9],
                     18: [241.9, 254.4],
                     16: [254.4, 266.6],
                     14: [266.6, 278.3],
                     12: [278.3, 289.9]}

    return elapsed_times[num_sweep][0], elapsed_times[num_sweep][1]


def int2float_radar(data, varname, index_angle):
    """
    Converts radar moments from bit to float
    Inputs
    ------
    data: np.array
        moment data as loaded from h5 file
    varname: str
        name of the moment (i.e. 'ZH')
    index_angle: int
        rank of the sweep-1 (converted to base 0)
    Returns
    -------
    output: np.array
        moment data converted to float
    """
    varname = convert_polvar_name('metranet', varname)
    NYQUIST_VEL = get_mymetadata('nyq_vel')

    output = np.zeros(data.shape)
    if varname in ['ZH', 'ZV', 'Z', 'ZHC']:
        output[data != 0] = (data[data != 0]-64)*0.5
        output[data == 0] = float('nan')
    elif varname == 'VEL':
        output[data != 0] = (data[data != 0]-128)/127*NYQUIST_VEL[index_angle]
        output[data == 0] = float('nan')
    elif varname == 'WID':
        output = data/255*NYQUIST_VEL[index_angle]
    elif varname in ['ZDR', 'ZDRC']:
        output[data != 0] = data[data != 0]*1.0/16.1259842 - 7.9375
        output[data == 0] = float('nan')
    elif varname == 'RHO':
        output[data != 0] = 1.003-10**(-(data[data != 0]-1.0)/100)
        output[data == 0] = float('nan')
    elif varname == 'PHI':
        output[data != 0] = (data[data != 0]-32768)/32767*180
        output[data == 0] = float('nan')
    elif varname == 'CLUT':
        output = data
    else:
        output = data
        warnings.warn(
            ("Warning, %s was not found and could not be converted")
            % (varname))

    return output


def readMXPOLRadData(filename, variableList, max_range=np.Inf, min_range=0):
    """
    Reads a netcdf containing processed radar data in polar coordinates
    Parameters
    ----------
    filename: str
        complete path of the file
    variableList: list
        list of variables to be read
    Returns
    -------
    varPol: dict
        dictionary containing the variables, the azimuth and the range
    metadata: dict
        dictionary containing the metadata of the file
    """

    varPol = {}
    metadata = {}
    ncid = netCDF4.Dataset(filename)

    time_data = ncid.variables['Time']
    time_data -= time_data[0]  # To get time in seconds from beginning of scan

    rrange = ncid.variables['Range'][:]

    # Get indexes between min_range and max_range
    idx2keep = np.where(np.logical_and(
        rrange < max_range, rrange > min_range))[0]
    rrange = rrange[idx2keep]

    # Get variables in polar coordinates
    for varname in variableList:
        try:
            varPol[varname] = ncid.variables[varname][:].T
        except:
            pass

    varPol['resolution'] = ncid.__dict__['RangeResolution-value']
    varPol['range'] = rrange
    varPol['range_units'] = ncid.__dict__['RangeResolution-unit']
    varPol['azimuth'] = ncid.variables['Azimuth'][:]
    try:
        varPol['azim_units'] = ncid.__dict__['Azimuth-unit']
    except KeyError:
        varPol['azim_units'] = ncid.variables['Azimuth'].Units
    varPol['elevation'] = ncid.variables['Elevation'][:]
    try:
        varPol['elev_units'] = ncid.__dict__['Elevation-unit']
    except KeyError:
        varPol['elev_units'] = ncid.variables['Elevation'].Units
    varPol['nyquist_vel'] = ncid.__dict__['NyquistVelocity-value']
    varPol['longitude'] = ncid.__dict__['Longitude-value']
    varPol['lon_units'] = ncid.__dict__['Longitude-unit']
    varPol['latitude'] = ncid.__dict__['Latitude-value']
    varPol['lat_units'] = ncid.__dict__['Latitude-unit']
    varPol['altitude'] = ncid.__dict__['Altitude-value']
    varPol['alt_units'] = ncid.__dict__['Altitude-unit']
    varPol['time'] = time_data

    metadata['Source'] = ncid.__dict__['Source']
    metadata['Institution'] = ncid.__dict__['Institution']
    metadata['History'] = ncid.__dict__['History']
    metadata['ContactInformation'] = ncid.__dict__['ContactInformation']

    # Close netcdf
    ncid.close()

    return metadata, varPol

def readIDLRadData(filename, variableList, max_range=np.Inf, min_range=0):
    """
    Reads a netcdf containing IDL processed radar data in polar coordinates
    Parameters
    ----------
    filename: str
        complete path of the file
    variableList: list
        list of variables to be read
    Returns
    -------
    varPol: dict
        dictionary containing the variables, the azimuth and the range
    metadata: dict
        dictionary containing the metadata of the file
    """

    varPol = {}
    metadata = {}
    ncid = netCDF4.Dataset(filename)

    time_data = ncid.variables['Time']
    time_data -= time_data[0]  # To get time in seconds from beginning of scan

    rrange = ncid.variables['Range'][:]

    # Get indexes between min_range and max_range
    idx2keep = np.where(np.logical_and(
        rrange < max_range, rrange > min_range))[0]
    rrange = rrange[idx2keep]

    # Get variables in polar coordinates
    for varname in variableList:
        try:
            varPol[varname] = ncid.variables[varname][:].T
        except:
            pass

    varPol['resolution'] = ncid.__dict__['RangeResolution-value']
    varPol['range'] = rrange
    varPol['range_units'] = ncid.__dict__['RangeResolution-unit']
    # because this data seems to be on -180 to 180
    #varPol['azimuth'] = (ncid.variables['Azimuth'][:] + 180)%360
    varPol['azimuth'] = ncid.variables['Azimuth'][:]
    try:
        varPol['azim_units'] = ncid.__dict__['Azimuth-unit']
    except KeyError:
        varPol['azim_units'] = ncid.variables['Azimuth'].Units
    varPol['elevation'] = ncid.variables['Elevation'][:]
    try:
        varPol['elev_units'] = ncid.__dict__['Elevation-unit']
    except KeyError:
        varPol['elev_units'] = ncid.variables['Elevation'].Units
    varPol['nyquist_vel'] = ncid.__dict__['NyquistVelocity-value']
    varPol['longitude'] = ncid.__dict__['Longitude-value']
    varPol['lon_units'] = ncid.__dict__['Longitude-unit']
    varPol['latitude'] = ncid.__dict__['Latitude-value']
    varPol['lat_units'] = ncid.__dict__['Latitude-unit']
    varPol['altitude'] = ncid.__dict__['Altitude-value']
    varPol['alt_units'] = ncid.__dict__['Altitude-unit']
    varPol['time'] = time_data

    metadata['Source'] = ncid.__dict__['Source']
    metadata['Institution'] = ncid.__dict__['Institution']
    metadata['History'] = ncid.__dict__['History']
    metadata['ContactInformation'] = ncid.__dict__['ContactInformation']

    # Close netcdf
    ncid.close()

    return metadata, varPol


def readCHRadData(filename, radar_name, variableList, radial_resolution,
                  max_range=np.Inf, min_range=0):
    """
    Reads a HDF5 file containing processed radar data in polar coordinates
    Parameters
    ----------
    filename: str
        complete path of the file
    radar_name: str
        name of MCH radar
    variableList: list
        list of variables to be read
    radial_resolution: float
        resolution of the radar in metres (i.e. high: 83.3, low: 500.)
    max_range: float
        maximum range upto which to read data
    min_range: float
        mimimum range from which to read data
    Returns
    -------
    varPol: dict
        the projected variables, the azimuth and the range
    """
    # check that h5py library is available
    if not _H5PY_AVAILABLE:
        raise MissingOptionalDependency(
            "h5py is required to use readCHRadData but is not installed")

    varPol = {}
    h5id = h5py.File(filename, 'r')

    ELEVATION_ANGLES = get_elevation_metadata(radar_name)
    radar_info = generate_radar_table(radar_name)
    ANG_RES = radar_info['dbbeam']
    NYQUIST_VEL = get_mymetadata('nyq_vel')

    # Get dimensions
    siz = h5id['moments']['Z'].shape
    rng = np.arange(0, siz[1])*radial_resolution
    idx2keep = np.where(np.logical_and(
        rng < max_range, rng > min_range))[0]
    rng = rng[idx2keep]
    azimuth = np.arange(0, siz[0])*ANG_RES
    index_angle = int(re.findall(r"\.([0-9]{3})\.", filename)[0])-1
    elevation = ELEVATION_ANGLES[index_angle]
    # Get variables in polar coordinates
    for varname in variableList:
        varname = convert_polvar_name('MCH', varname)
        data = []
        data = h5id['moments'][varname][:].T
        data = np.asarray(data)
        data = data.astype(float)
        clut = h5id['moments']['CLUT'][:].T
        data[clut >= 100] = float('nan')  # Remove clutter
        data = data[idx2keep, :]
        varPol[varname] = int2float_radar(data, varname, index_angle)

    varPol['resolution'] = rng[3]-rng[2]
    varPol['range'] = rng
    varPol['azimuth'] = azimuth
    varPol['elevation'] = elevation
    varPol['nyquist_vel'] = NYQUIST_VEL[index_angle]
    # Close netcdf
    h5id.close()

    return varPol


# ------------------------ utilities - config ------------------------- #

_dirname = os.path.dirname(__file__)
_DEFAULT_CONFIG_FILE = os.path.join(_dirname, 'mxpol_config.py')


def load_myconfig(filename=None):
    """
    Load configuration from a config file.
    Parameters
    ----------
    filename: str
        Filename of the configuration file. If None the default configuration
        file is loaded from the directory.
    Returns
    -------
    _DEFAULT_METADATA: dict
        Dictionary with metadata
    """

    if filename is None:
        filename = _DEFAULT_CONFIG_FILE

    # private:

    global cfile
    global _DEFAULT_POLARNAMES
    global _DEFAULT_METADATA
    global _DEFAULT_RADAR_INFO

    spec = importlib.util.spec_from_file_location("metadata_config", 
                                                  filename)
    cfile = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cfile)
    
    _DEFAULT_METADATA = cfile.MY_METADATA
    _DEFAULT_POLARNAMES = cfile.MY_POLARNAMES
    _DEFAULT_RADAR_INFO = cfile.RADAR_INFO

    return _DEFAULT_METADATA


def get_mymetadata(p, filename=None):
    """
    Return a dictionary of metadata for a given parameter, p.
    An empty dictionary will be returned if no metadata dictionary exists for
    parameter p.
    Parameters
    ----------
    p: str
        parameter name (i.e. Polvar) for which to return metadata
    filename: str
        Filename of the configuration file. If None the default configuration
        file is loaded from the directory.
    Returns
    -------
    _DEFAULT_METADATA[p].copy(): dict
        a copy of the parameter of interest from the metadata dictionary
    """
    load_myconfig(filename=filename)

    if p in _DEFAULT_METADATA:
        return _DEFAULT_METADATA[p].copy()

    return {}


def get_elevation_metadata(radarname, filename=None):
    """
    Gets the elevation angles for each sweep from the configuration file
    Inputs
    ------
    radarname: str
        name of the radar for which to retrieve elevation angles
    filename: str
        name of the configuration file, if None, the default configuration
        file is used
    Returns
    -------
    _DEFAULT_RADAR_INFO['elevations'][radarname]: list
        list of elevation angles in degrees
    or None if not available
    """
    load_myconfig(filename=filename)

    if radarname in _DEFAULT_RADAR_INFO['elevations']:
        return _DEFAULT_RADAR_INFO['elevations'][radarname]
    else:
        print(("no elevation angles in configfile for radar %s") % (radarname))


def generate_radar_table(radarname, filename=None):
    """
    Generates a table with basic radar info, based on the given (or default)
    configfile
    Parameters
    ----------
    radarname: str
        name of the radar (i.e. 'ALB' or 'A', 'MXPOL' etc)
    filename: str
        path and name of the configfile, if None, the default configfile is
        used
    Returns
    -------
    radar_table: dict
        table containing basic radar info
    """
    load_myconfig(filename=filename)

    if radarname in _DEFAULT_RADAR_INFO['radarID']:
        radarname = _DEFAULT_RADAR_INFO['radarID'][radarname]
        radar_table = get_mymetadata('Radar_info', filename=filename)
        for key in radar_table:
            if key in _DEFAULT_RADAR_INFO:
                radar_table[key] = _DEFAULT_RADAR_INFO[key][radarname]
            else:
                radar_table[key] = None
        return radar_table

    return None


def generate_polvar_metadata(polvar, filename=None):
    """
    Generates a dictionary with metadata for a polarimetric variable
    Parameters
    ----------
    polvar: str
        polatimetric variable of interest
    filename: str
        Filename of the configuration file. If None the default configuration
        file is loaded from the directory.
    Returns
    -------
    polvar_metadata: dict
        dictionary with metatdata for polarimetric variable of interest
    """
    load_myconfig(filename=filename)
    polvar = convert_polvar_name('LTE', polvar)

    if polvar in _DEFAULT_POLARNAMES:
        (standard_name, long_name, units, valid_min, valid_max,
         plot_interval) = _DEFAULT_POLARNAMES[polvar]
    else:
        (standard_name, long_name, units, valid_min, valid_max,
         plot_interval) = None, None, None, None, None, None

    polvar_metadata = get_mymetadata('Polvar', filename)
    polvar_metadata['units'] = units
    polvar_metadata['standard_name'] = standard_name
    polvar_metadata['short_name'] = convert_polvar_name('MCH', polvar)
    polvar_metadata['long_name'] = long_name
    polvar_metadata['valid_min'] = valid_min
    polvar_metadata['valid_max'] = valid_max
    polvar_metadata['plot_interval'] = plot_interval

    return polvar_metadata


def convert_polvar_name(convention, polvar):
    """
    Finds the correct variable name for a given convention (MXPOL, MCH) and
    a given variable name which was spelled with a different case or
    according to a different convention. For example, MXPOL convention uses
    'Z' for the reflectivity variable, but if a user inserted 'Zh' this
    function will convert it to 'Z'.
    Parameters
    ----------
    convention : str, destination convention; either MCH or LTE
    polvar : str, key of polarimetric variable to be converted
    Returns
    -------
    mykey : str, polarimertric variable key as used within the ProfileLab
        toolbox context
    """
    # Generate dictionary for the conversion
    metranet_list = [
        'ZH', 'ZV', 'ZDR', 'PHI', 'VEL', 'VEL', 'WID', 'RHO', 'CLUT', 'MPH',
        'STA1', 'STA2', 'WBN', 'ZHC', 'ZDRC', 'ZDRP', 'Kdpc', 'Rhohvc']
    MCH_list = [
        'Z', 'ZV', 'ZDR', 'PHIDP', 'V', 'V', 'W', 'RHO', 'CLUT', 'MPH', 'STA1',
        'STA2', 'WBN', 'Zhc', 'Zdrc', 'Hydrometeor_type_from_Besic1', 'Kdpc', 'RHOC']
    # ZhCorr and ZdrCorr have been changed to Zhc and Zdrc!
    LTE_list = [
        'Zh', 'Zv', 'Zdr', 'Phidp', 'RVel', 'Rvel', 'Sw', 'Rhohv', 'Clut', 'mph',
        'sta1', 'sta2', 'wbn', 'Zhc', 'Zdrc', 'Hydroclass', 'Kdpc', 'Rhohvc']
    IDL_list = [
        'Zh', 'Zv', 'Zdr', 'Phidp_raw', 'V', 'V', 'W', 'uRhohv', 'CLUT', 'MPH', 'STA1',
        'STA2', 'WBN', 'Zhc', 'Zdrc', 'TYPECLUS2', 'Kdpc', 'Rhohvc']
    pyrad_list = [
        'reflectivity', 'reflectivity_vv', 'differential_reflectivity',
        'differential_phase', 'velocity', 'velocity', 'spectrum_width',
        'uncorrected_cross_correlation_ratio', 'radar_echo_id', 'MPH',
        'STA1', 'STA2', 'WBN', 'corrected_reflectivity',
        'corrected_differential_reflectivity', 'radar_echo_classification',
        'corrected_specific_differential_phase',
        'corrected_cross_correlation_ratio']

    convertkeys = {}
    convertkeys['MCH'] = {}
    convertkeys['LTE'] = {}
    convertkeys['metranet'] = {}
    convertkeys['IDL'] = {}
    convertkeys['pyrad'] = {}

    for i, MCH in enumerate(MCH_list):
        convertkeys['MCH'][MCH] = [
            LTE_list[i], metranet_list[i], IDL_list[i], pyrad_list[i]]

    convertkeys['LTE'] = {}
    for i, LTE in enumerate(LTE_list):
        convertkeys['LTE'][LTE] = [
            MCH_list[i], metranet_list[i], IDL_list[i], pyrad_list[i]]

    for i, metranet in enumerate(metranet_list):
        convertkeys['metranet'][metranet] = [
            MCH_list[i], LTE_list[i], IDL_list[i], pyrad_list[i]]

    for i, IDL in enumerate(IDL_list):
        convertkeys['IDL'][IDL] = [
            metranet_list[i], LTE_list[i], MCH_list[i], pyrad_list[i]]

    for i, pyrad in enumerate(pyrad_list):
        convertkeys['pyrad'][pyrad] = [
            metranet_list[i], LTE_list[i], MCH_list[i], IDL_list[i]]

    # translate between conventions
    mykey = polvar

    for key, value in convertkeys[convention].items():
        if polvar in value:
            mykey = key
            break

    return mykey
