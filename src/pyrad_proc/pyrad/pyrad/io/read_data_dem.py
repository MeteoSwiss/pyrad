"""
pyrad.io.read_data_dem
========================

Functions for reading data derived from Digital Elevation Models (DEM)

.. autosummary::
    :toctree: generated/

    dem2radar_data
    read_idrisi_data
    read_idrisi_metadata
    _prepare_for_interpolation


"""

from warnings import warn
import numpy as np
import pandas as pd
import pathlib
from scipy.interpolate import RegularGridInterpolator

# check existence of gdal
try:
    from osgeo import gdal, osr
    _GDAL_AVAILABLE = True
except ImportError:
    try:
        import gdal
        import osr
        _GDAL_AVAILABLE = True
    except ImportError:
        _GDAL_AVAILABLE = False
from pyproj import CRS
import pyart

from pyart.config import get_metadata
from ..io.read_data_cosmo import _put_radar_in_swiss_coord

# from memory_profiler import profile

# import time


def dem2radar_data(radar, dem_data, slice_xy=True, field_name='visibility'):
    """
    get the DEM value corresponding to each radar gate using nearest
    neighbour interpolation

    Parameters
    ----------
    radar : Radar
        the radar object containing the information on the position of the
        radar gates
    dem_data : dict
        dictionary containing the DEM data
    slice_xy : boolean
        if true the horizontal plane of the DEM field is cut to the
        dimensions of the radar field
    field_names : str
        names of DEM fields to convert

    Returns
    -------
    dem_field : dict
        Dictionary with the DEM fields and metadata

    """
    # debugging
    # start_time = time.time()

    x_radar, y_radar, _ = _put_radar_in_swiss_coord(radar)

    (x_dem, y_dem, ind_xmin, ind_ymin, ind_xmax, ind_ymax) = (
        _prepare_for_interpolation(
            x_radar, y_radar, dem_data, slice_xy=slice_xy))

    if field_name not in dem_data:
        warn('DEM field '+field_name+' data not available')
        return None

    values = dem_data[field_name]['data'][
        ind_xmin:ind_xmax+1, ind_ymin:ind_ymax+1]
    

    # Note RegularGridInterpolator is 10x faster than NDNearestInterpolator
    # and has the advantage of not extrapolating outside of grid domain

    # replace masked values with nans
    values = np.ma.filled(values, np.nan)
    interp_func = RegularGridInterpolator(
        (x_dem, y_dem), values, bounds_error = False)
    

    # interpolate
    data_interp = interp_func((x_radar, y_radar))
    
    del values
    # restore mask
    data_interp = np.ma.masked_equal(data_interp, np.nan)

    # put field
    field_dict = get_metadata(field_name)
    field_dict['data'] = data_interp.astype(float)

    del data_interp

    return field_dict

# @profile
def read_dem(fname, field_name = 'terrain_altitude', fill_value=None,
             projparams = None):
    """
    Generic reader that reads DEM data from any format, will infer the proper
    reader from filename extension

    Parameters
    ----------
    fname : str
        name of the file to read
    field_name : str
        name of the readed variable
    fill_value : float
        The fill value, if not provided will be infered from metadata 
        if possible
    projparams : projection transform as can be used by pyproj, either a 
        OGC WKT or Proj4 string, see epsg.io for a list, if not provided 
        will be infered from the file, or for ASCII, LV1903 will be used
        
    Returns
    -------
    dem_data : dictionary
        dictionary with the data and metadata

    """
    extension = pathlib.Path(fname).suffix
  
    if extension in ['.tif','.tiff','.gtif']:
        return  read_geotiff_data(fname, field_name, fill_value, projparams)
    elif extension in ['.asc','.dem','.txt']:
        return read_ascii_data(fname, field_name, fill_value, projparams)
    elif extension in ['.rst']:
        return read_idrisi_data(fname, field_name, fill_value, projparams)
    else:
        warn('Unable to read file %s, extension must be .tif .tiff .gtif, '+
             '.asc .dem .txt .rst',
             fname)
        return None
    
# @profile
def read_geotiff_data(fname, field_name = 'terrain_altitude', 
                      fill_value = None, projparams = None):
    
    """
    Reads DEM data from a generic geotiff file, unit system is expected to be
    in meters!

    Parameters
    ----------
    fname : str
        name of the file to read
    field_name : str
        name of the readed variable
    fill_value : float
        The fill value, if not provided will be infered from metadata 
        (recommended)
    projparams : projection transform as can be used by pyproj, either a 
        EPSG code integer (21781 for CH1903 for example), or a OGC WKT or 
        Proj4 string, see epsg.io for a list, 
        if not provided will be infered from the idrisi file
        
    Returns
    -------
    dem_data : dictionary
        dictionary with the data and metadata

    """
    if not _GDAL_AVAILABLE:
        warn("gdal is required to use read_geotiff_data but is not installed")
        return None

    if type(projparams) == int: # Retrieve Wkt code from EPSG number
        proj = osr.SpatialReference()
        proj.ImportFromEPSG(projparams)
        projparams = proj.ExportToWkt()
        
    # read the data
    try:
        raster = gdal.Open(fname)

        width = raster.RasterXSize
        height = raster.RasterYSize
        gt = raster.GetGeoTransform()
        minx = gt[0]
        miny = gt[3] + width*gt[4] + height*gt[5] 
        maxx = gt[0] + width*gt[1] + height*gt[2]
        maxy = gt[3] 
        
        metadata = {}
        metadata['resolution'] = np.abs(gt[1])
        metadata['min. X'] = minx
        metadata['min. Y'] = miny
        metadata['rows'] = raster.RasterYSize
        metadata['columns'] = raster.RasterXSize
        metadata['value units'] = 'meters'
        metadata['max. X'] = (metadata['min. X'] + 
                              metadata['resolution'] * metadata['columns'])
        metadata['max. Y'] = (metadata['min. Y'] + 
                              metadata['resolution'] * metadata['rows'])
        metadata['flag value'] = raster.GetRasterBand(1).GetNoDataValue()
        
        if not fill_value: 
            fill_value = metadata['flag value']

        raster_array = raster.ReadAsArray()
        raster_array = np.ma.masked_equal(raster_array, fill_value)
    
    
        field_dict = get_metadata(field_name)
        # Pyart Grid needs 3D array for plotting
        field_dict['data'] = raster_array[::-1,:][None,:,:]
        field_dict['data'] = np.transpose(raster_array)[:, ::-1]
        field_dict['units'] = metadata['value units']
        
        x = get_metadata('x')
        y = get_metadata('y')
        
        x['data'] = (np.linspace(minx,maxx,width))
        y['data'] = (np.linspace(miny,maxy,height))

        dem_data = {
            'metadata': metadata,
            'x': x,
            'y': y,
            field_name: field_dict
        }

        if projparams == None:
            projparams = osr.SpatialReference(wkt=raster.GetProjection())
            projparams = projparams.ExportToWkt()
        if projparams == '':
            warn('No projection info could be found in file, assuming '+\
                 'the coordinate system is LV1903')
            projparams = _get_lv1903_wkt()
            
        time = get_metadata('grid_time')
        time['data'] = np.array([0.0])
        time['units'] = 'seconds since 2000-01-01T00:00:00Z'
    
        # The use of CRS().to_dict() is required to use GridMapDisplay of Pyart
        # which expects a dict for the projection attribute of the grid
        dem_data = pyart.core.Grid(time, {field_name : dem_data[field_name]}, 
             dem_data['metadata'], 0, 0, 0, dem_data['x'], dem_data['y'], 
                     {'data':[0]}, projection = CRS(projparams).to_dict())

        return dem_data
        
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None

# @profile
def read_ascii_data(fname, field_name = 'terrain_altitude', fill_value = None,
                    projparams = None):
    """
    Reads DEM data from an ASCII file in the Swisstopo format, the
    unit system is expected to be in meters!

    Parameters
    ----------
    fname : str
        name of the file to read
    field_name : str
        name of the readed variable
    fill_value : float
        The fill value, if not provided will be infered from metadata 
        (recommended)
    projparams : projection transform as can be used by pyproj, either a 
        EPSG code integer (21781 for CH1903 for example), or a OGC WKT or 
        Proj4 string, see epsg.io for a list, 
        if not provided CH1903 (EPSG:21781) will be used

        
    Returns
    -------
    dem_data : dictionary
        dictionary with the data and metadata

    """
    if type(projparams) == int: # Retrieve Wkt code from EPSG number
        proj = osr.SpatialReference()
        proj.ImportFromEPSG(projparams)
        projparams = proj.ExportToWkt()

    # read the data
    try:
        asciidata = pd.read_csv(fname, header = None)
        
        metadata = {}
        metadata['columns'] = int(asciidata.iloc[0][0].split(' ')[1])
        metadata['rows'] = int(asciidata.iloc[1][0].split(' ')[1])
        metadata['min. X'] = float(asciidata.iloc[2][0].split(' ')[1])
        metadata['min. Y'] = float(asciidata.iloc[3][0].split(' ')[1])
        metadata['resolution'] = float(asciidata.iloc[4][0].split(' ')[1])
        metadata['flag value'] = float(asciidata.iloc[5][0].split(' ')[1])
        metadata['max. X'] = (metadata['min. X'] + 
                              metadata['resolution'] * metadata['columns'])
        metadata['max. Y'] = (metadata['min. Y'] + 
                              metadata['resolution'] * metadata['rows'])
        metadata['value units'] = 'm'
        
        if not fill_value: 
            fill_value = metadata['flag value']
            
        raster_array = pd.read_csv(fname, skiprows = 6, header = None,
                                  sep = ' ')
        raster_array = np.array(raster_array)
        raster_array = raster_array[np.isfinite(raster_array)]
        raster_array = np.reshape(raster_array,
                                 (metadata['rows'],metadata['columns']))
        raster_array = np.ma.masked_equal(raster_array, fill_value)
        
        field_dict = get_metadata(field_name)
        # Pyart Grid needs 3D array for plotting
        field_dict['data'] = raster_array[::-1,:][None,:,:]
        rasterarray = pd.read_csv(fname, skiprows = 6, header = None,
                                  sep = ' ')
        rasterarray = np.array(rasterarray)
        rasterarray = rasterarray[np.isfinite(rasterarray)]
        rasterarray = np.reshape(rasterarray,
                                 (metadata['rows'],metadata['columns']))
        rasterarray = np.ma.masked_equal(rasterarray, fill_value)
        
        field_dict = get_metadata(field_name)
        field_dict['data'] = np.transpose(rasterarray)[:, ::-1]
        field_dict['units'] = metadata['value units']
            
        x = get_metadata('x')
        y = get_metadata('y')
        x['data'] = (
            np.arange(metadata['columns'])*metadata['resolution'] +
            metadata['resolution']/2.+metadata['min. X'])

        y['data'] = (
            np.arange(metadata['rows'])*metadata['resolution'] +
            metadata['resolution']/2.+metadata['min. Y'])
        
        dem_data = {
            'metadata': metadata,
            'x': x,
            'y': y,
            field_name: field_dict
        }
        if projparams == None:
            projparams = _get_lv1903_wkt()
        
        time = get_metadata('grid_time')
        time['data'] = np.array([0.0])
        time['units'] = 'seconds since 2000-01-01T00:00:00Z'
    
        # The use of CRS().to_dict() is required to use GridMapDisplay of Pyart
        # which expects a dict for the projection attribute of the grid
        dem_data = pyart.core.Grid(time, {field_name : dem_data[field_name]}, 
             dem_data['metadata'], 0, 0, 0, dem_data['x'], dem_data['y'], 
                     {'data':[0]}, projection = CRS(projparams).to_dict())


        return dem_data
        
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None
    
# @profile
def read_idrisi_data(fname, field_name = 'terrain_altitude', fill_value = None,
                     projparams = None):
    """
    Reads DEM data from an IDRISI .rst file

    Parameters
    ----------
    fname : str
        name of the file to read
    field_name : str
        name of the readed variable
    fill_value : float
        The fill value
    projparams : projection transform as can be used by pyproj, either a 
        EPSG code integer (21781 for CH1903 for example), or a OGC WKT or 
        Proj4 string, see epsg.io for a list, 
        if not provided will be infered from the idrisi file


    Returns
    -------
    dem_data : dictionary
        dictionary with the data and metadata

    """
    if not _GDAL_AVAILABLE:
        warn("gdal is required to use read_idrisi_data but is not installed")
        return None
    
    if type(projparams) == int: # Retrieve Wkt code from EPSG number
        proj = osr.SpatialReference()
        proj.ImportFromEPSG(projparams)
        projparams = proj.ExportToWkt()
        
    # read the data
    try:
        if fill_value == None:
            fill_value = -99.
            
        raster = gdal.Open(fname)
        raster_array = raster.ReadAsArray()
        raster_array = np.ma.masked_equal(raster_array, fill_value)

        metadata = read_idrisi_metadata(fname)

        if metadata is None:
            return None

        field_dict = get_metadata(field_name)
        # Pyart Grid needs 3D array for plotting
        field_dict['data'] = raster_array[::-1,:][None,:,:]
        field_dict['units'] = metadata['value units']

        x = get_metadata('x')
        y = get_metadata('y')
        x['data'] = (
            np.arange(raster.RasterXSize)*metadata['resolution'] +
            metadata['resolution']/2.+metadata['min. X'])

        y['data'] = (
            np.arange(raster.RasterYSize)*metadata['resolution'] +
            metadata['resolution']/2.+metadata['min. Y'])

        dem_data = {
            'metadata': metadata,
            'x': x,
            'y': y,
            field_name: field_dict
        }
        if projparams == None:
            projparams = osr.SpatialReference(wkt=raster.GetProjection())
            projparams = projparams.ExportToWkt()
        if projparams == '':
            warn('No projection info could be found in file, assuming '+\
                 'the coordinate system is LV1903')
            projparams = _get_lv1903_wkt()
         
        time = get_metadata('grid_time')
        time['data'] = np.array([0.0])
        time['units'] = 'seconds since 2000-01-01T00:00:00Z'
        
        # The use of CRS().to_dict() is required to use GridMapDisplay of Pyart
        # which expects a dict for the projection attribute of the grid
        dem_data = pyart.core.Grid(time, {field_name : dem_data[field_name]}, 
             dem_data['metadata'], 0, 0, 0, dem_data['x'], dem_data['y'], 
                     {'data':[0]}, projection = CRS(projparams).to_dict())

        return dem_data
    
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None


def read_idrisi_metadata(fname):
    """
    Reads DEM metadata from a IDRISI .rdc file

    Parameters
    ----------
    fname : str
        name of the file to read

    Returns
    -------
    metadata : dictionary
        dictionary with the metadata

    """
    # read the data
    fname_rdc = fname.replace('.rst', '.rdc')

    try:
        metadata = dict()
        with open(fname_rdc, 'r', newline='') as txtfile:
            for line in txtfile:
                strs = line.split(':')
                metadata.update({
                    strs[0].strip(): pyart.aux_io.convert_data(strs[1].strip())})

        return metadata
    except EnvironmentError:
        warn('Unable to read file '+fname_rdc)
        return None


def _prepare_for_interpolation(x_radar, y_radar, dem_coord, slice_xy=True):
    """
    prepares the DEM 2D volume for interpolation:
        1. if set slices the DEM data to the area
    covered by the radar
        2. creates the x, y grid for the interpolation

    Parameters
    ----------
    x_radar, y_radar : arrays
        The Swiss coordinates of the radar
    dem_coord : dict
        dictionary containing the DEM coordinates
    slice_xy : boolean
        if true the horizontal plane of the DEM field is cut to the
        dimensions of the radar field

    Returns
    -------
    x_dem, y_dem : 1D arrays
        arrays containing the flatten swiss coordinates of the DEM data in
        the area of interest
    ind_xmin, ind_ymin, ind_xmax, ind_ymax : ints
        the minimum and maximum indices of each dimension

    """
    nx_dem = len(dem_coord['x']['data'])
    ny_dem = len(dem_coord['y']['data'])

    if slice_xy:
        # get the D data within the radar range
        xmin = np.min(x_radar)
        xmax = np.max(x_radar)
        ymin = np.min(y_radar)
        ymax = np.max(y_radar)

        ind_xmin = np.where(dem_coord['x']['data'] < xmin)[0]
        if ind_xmin.size == 0:
            ind_xmin = 0
        else:
            ind_xmin = ind_xmin[-1]

        ind_xmax = np.where(dem_coord['x']['data'] > xmax)[0]
        if ind_xmax.size == 0:
            ind_xmax = nx_dem-1
        else:
            ind_xmax = ind_xmax[0]

        ind_ymin = np.where(dem_coord['y']['data'] < ymin)[0]
        if ind_ymin.size == 0:
            ind_ymin = 0
        else:
            ind_ymin = ind_ymin[-1]

        ind_ymax = np.where(dem_coord['y']['data'] > ymax)[0]
        if ind_ymax.size == 0:
            ind_ymax = ny_dem-1
        else:
            ind_ymax = ind_ymax[0]
    else:
        ind_xmin = 0
        ind_xmax = nx_dem-1
        ind_ymin = 0
        ind_ymax = ny_dem-1


    x_dem = dem_coord['x']['data'][ind_xmin:ind_xmax+1]
    y_dem = dem_coord['y']['data'][ind_ymin:ind_ymax+1]

    # Not used with RegularGridInterpolator
    # nx = ind_xmax-ind_xmin+1
    # ny = ind_ymax-ind_ymin+1
    # x_dem = (
    #     np.broadcast_to(x_dem.reshape(nx, 1), (nx, ny))).flatten()
    # y_dem = (
    #     np.broadcast_to(y_dem.reshape(1, ny), (nx, ny))).flatten()

    return (x_dem, y_dem, ind_xmin, ind_ymin, ind_xmax, ind_ymax)

def _get_lv1903_wkt():
    lv1903 = osr.SpatialReference( )
    lv1903.ImportFromEPSG(21781)
    return lv1903.ExportToWkt()
