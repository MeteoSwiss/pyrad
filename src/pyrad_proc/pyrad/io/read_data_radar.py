"""
pyrad.io.read_data_radar
========================

Functions for reading radar data files

.. autosummary::
    :toctree: generated/

    get_data
    merge_scans_rainbow
    merge_scans_psr
    merge_scans_psr_spectra
    merge_scans_dem
    merge_scans_rad4alp
    merge_scans_odim
    merge_scans_odimgrid
    merge_scans_knmih5_grid
    merge_scans_odimbirds
    merge_scans_gamic
    merge_scans_mfcfradial
    merge_scans_nexrad2
    merge_scans_cfradial
    merge_scans_cfradial2
    merge_scans_skyecho
    merge_scans_cf1
    merge_scans_mxpol
    merge_scans_icon
    merge_scans_icon_rad4alp
    merge_scans_dem_rad4alp
    merge_scans_other_rad4alp
    merge_scans_iq_rad4alp
    merge_fields_rainbow
    merge_fields_psr
    merge_fields_psr_spectra
    merge_fields_rad4alp_grid
    merge_fields_sat_grid
    merge_fields_mf_grid
    merge_fields_pyrad
    merge_fields_pyradicon
    merge_fields_pyradgrid
    merge_fields_pyrad_spectra
    merge_fields_dem
    merge_fields_icon
    get_data_rainbow
    get_data_rad4alp
    get_data_odim
    get_data_odimgrid
    get_data_gamic
    add_field
    interpol_field
    crop_grid
    merge_grids

"""

import glob
import datetime
import platform
import os
from ..util import warn
from copy import deepcopy

import numpy as np

from scipy.interpolate import RegularGridInterpolator

try:
    import wradlib as wrl

    _WRADLIB_AVAILABLE = True
except ImportError:
    _WRADLIB_AVAILABLE = False

import pyart

# check existence of METRANET library
try:
    METRANET_LIB = pyart.aux_io.get_library(momentms=False)
    if platform.system() == "Linux":
        METRANET_LIB = pyart.aux_io.get_library(momentms=True)
    _METRANETLIB_AVAILABLE = True
except BaseException:
    # bare exception needed to capture error
    _METRANETLIB_AVAILABLE = False

try:
    import boto3

    _BOTO3_AVAILABLE = True
except ImportError:
    warn(
        "boto3 is not installed, data cannot be retrieved from S3 buckets!",
        use_debug=False,
    )
    _BOTO3_AVAILABLE = False

from .read_data_other import read_status, read_rad4alp_icon, read_rad4alp_vis
from .read_data_mxpol import pyrad_MXPOL, pyrad_MCH

from .io_aux import get_datatype_metranet, get_fieldname_pyart, get_file_list
from .io_aux import get_datatype_odim, get_datatype_knmi, find_date_in_file_name
from .io_aux import get_datatype_fields, get_datetime, map_hydro, map_Doppler
from .io_aux import find_icon_file, find_rad4alpicon_file
from .io_aux import find_pyradicon_file, get_datatype_skyecho
from .io_aux import get_rad4alp_prod_fname, get_rad4alp_grid_dir
from .io_aux import get_rad4alp_dir, get_scan_files_to_merge
from .io_aux import get_scan_files_to_merge_s3


def _open_s3_client(cfg):
    if len(cfg["s3Certificates"]):
        s3verify = cfg["s3Certificates"]
    else:
        s3verify = cfg["s3Verify"]

    s3_client = boto3.client(
        "s3",
        endpoint_url=cfg["s3EndpointRead"],
        aws_access_key_id=cfg["s3KeyRead"],
        aws_secret_access_key=cfg["s3SecretRead"],
        verify=s3verify,
    )
    return s3_client


def get_data(voltime, datatypesdescr, cfg):
    """
    Reads pyrad input data.

    Parameters
    ----------
    voltime : datetime object
        volume scan time
    datatypesdescr : list
        list of radar field types to read.
        Format : [radarnr]:[datagroup]:[datatype],[dataset],[product]
        'dataset' is only specified for data groups 'ODIM',
        'CFRADIAL', 'CFRADIAL2', 'CF1', 'ODIMPYRAD' 'PYRADGRID' and
        'NETCDFSPECTRA'.
        'product' is only specified for data groups 'CFRADIAL', 'ODIMPYRAD',
        'PYRADGRID' and 'NETCDFSPECTRA'
        The data group specifies the type file from which data is extracted.
        It can be:
            'RAINBOW': Propietary Leonardo format
            'COSMO': COSMO model data saved in Rainbow file format
            'DEM': Visibility data saved in Rainbow file format
            'PSR': Reads PSR data file to extract range gate information
                (Noise and transmitted power)

            'RAD4ALP': METRANET format used for the operational MeteoSwiss
                data. To find out which datatype to use to match a particular
                METRANET field name check the function 'get_datatype_metranet'
                in pyrad/io/io_aux.py
            'RAD4ALPCOSMO': COSMO model data saved in a binary file format.
                Used by operational MeteoSwiss radars
            'RAD4ALPDEM': Visibility data saved in a binary format used by
                operational MeteoSwiss radars
            'RAD4ALPHYDRO': Used to read the MeteoSwiss operational
                hydrometeor classification
            'RAD4ALPDOPPLER': Used to read the MeteoSwiss operational
                dealiased Doppler velocity

            'ODIM': Generic ODIM file format. For such types 'dataset'
                specifies the directory and file name date convention.
                Example: ODIM:dBZ,D{%Y-%m-%d}-F{%Y%m%d%H%M%S}. To find out
                which datatype to use to match a particular ODIM field name
                check the function 'get_datatype_odim' in pyrad/io/io_aux.py

            'ODIMBIRDS': Output of vol2bird algorithm in ODIM convention
                format.For such types 'dataset' specifies the directory and
                file name date convention.
                Example: ODIM:dBZ,D{%Y-%m-%d}-F{%Y%m%d%H%M%S}. To find out
                which datatype to use to match a particular ODIM field name
                check the function 'get_datatype_odim' in pyrad/io/io_aux.py

            'GAMIC': Gamic files.

            'NEXRADII': Nexrad-level II file format.

            'CFRADIAL2': CFRADIAL2 file format. For such types 'dataset'
                specifies the directory and file name date convention.
                Example: ODIM:dBZ,D{%Y-%m-%d}-F{%Y%m%d%H%M%S}. To find out
                which datatype to use to match a particular ODIM field name
                check the function 'get_datatype_odim' in pyrad/io/io_aux.py

            'CFRADIAL': CFRADIAL file format. For such types 'dataset'
                specifies the directory and file name date convention.
                Example: CFRADIAL:dBZ,D{%Y-%m-%d}-F{%Y%m%d%H%M%S}.

            'SKYECHO': SKYECHO netcdf file format. For such types 'dataset'
                specifies the directory and file name date convention.
                Example: SKYECHO:dBZ,D{%Y-%m-%d}-F{%Y%m%d%H%M%S}.

            'CF1': CF1 file format. For such types 'dataset'
                specifies the directory and file name date convention.
                Example: ODIM:dBZ,D{%Y-%m-%d}-F{%Y%m%d%H%M%S}. To find out
                which datatype to use to match a particular ODIM field name
                check the function 'get_datatype_odim' in pyrad/io/io_aux.py

            'MXPOL': MXPOL (EPFL) data written in a netcdf file

            'MFCFRADIAL': radar data from MeteoFrance written in CFRadial

            'CFRADIALPYRAD': CFRadial format with the naming convention and
                directory structure in which Pyrad saves the data. For such
                datatypes 'dataset' specifies the directory where the dataset
                is stored and 'product' specifies the directory where the
                product is stored.
                Example: CFRADIALPYRAD:dBZc,Att_ZPhi,SAVEVOL_dBZc
            'CFRADIALCOSMO': COSMO data in radar coordinates in a CFRadial
                file format.
            'ODIMPYRAD': ODIM file format with the naming convention and
                directory structure in which Pyrad saves the data.  For such
                datatypes 'dataset' specifies the directory where the dataset
                is stored and 'product' specifies the directroy where the
                product is stored.
                Example: ODIMPYRAD:dBZc,Att_ZPhi,SAVEVOL_dBZc

            'RAD4ALPGRID': METRANET format used for the operational MeteoSwiss
                Cartesian products.
            'RAD4ALPGIF': Format used for operational MeteoSwiss Cartesian
                products stored as gif files
            'RAD4ALPBIN': Format used for operational MeteoSwiss Cartesian
                products stored as binary files
            'PYRADGRID': Pyrad generated Cartesian grid products stored in a
                netcdf file. For such datatypes 'dataset' specifies the
                directory where the dataset is stored and 'product' specifies
                the directory where the product is stored.
                Example: PYRADGRID:RR,RZC,SAVEVOL
            'ODIMPYRADGRID': Pyrad generated Cartesian grid products in an
                ODIM HDF5 file. For such datatypes 'dataset' specifies the
                directory where the dataset is stored and 'product' specifies
                the directory where the product is stored.
                Example: ODIMPYRADGRID:RR,RZC,SAVEVOL
            'ODIMGRID': Gridded data in ODIM format. For such types 'dataset'
                specifies the directory and file name date convention.
                Example: ODIMGRID:dBZ,D{%Y-%m-%d}-F{%Y%m%d%H%M%S}.
            'KNMIH5GRID': KNMI gridded data in an H5 file. For such types
                'dataset' specifies the directory and file name date
                convention.
            'SATGRID': CF Netcdf from used for the MeteoSat satellite data
                in the CCS4 (Radar composite) grid.
            'MFBIN': Format used by some MeteoFrance products stored as binary
                files
            'MFPNG': Format used by some MeteoFrance products stored as binary
                files
            'MFGRIB': Format used by some MeteoFrance products stored as GRIB
                files
            'MFDAT': Format used by some MeteoFrance products stored as DAT
                (text) files
            'MFCF': Format used by some MeteoFrance products stored as netcdf
                CF files

            'PSRSPECTRA': Format used to store Rainbow power spectra
                recordings.

            'NETCDFSPECTRA': Format analogous to CFRadial and used to store
                Doppler spectral

            'RAD4ALPIQ': Format used to store rad4alp IQ data

        'RAINBOW', 'RAD4ALP', 'ODIM' 'ODIMBIRDS' CFRADIAL2', 'CF1' 'MFCFRADIAL'
        'GAMIC' and 'MXPOL' are primary data file sources and they cannot be
        mixed for the same radar. It is also the case for their complementary
        data files, i.e. 'COSMO' and 'RAD4ALPCOSMO', etc. 'CFRADIALPYRAD' and
        'ODIMPYRAD' are secondary data file sources and they can be combined
        with any other datagroup type.
        For a list of accepted datatypes and how they map to the Py-ART name
        convention check function 'get_field_name_pyart' in pyrad/io/io_aux.py
    cfg: dictionary of dictionaries
        configuration info to figure out where the data is

    Returns
    -------
    radar : Radar
        radar object

    """
    datatype_rainbow = []
    datatype_rad4alp = []
    datatype_odim = []
    dataset_odim = []
    datatype_gamic = []
    dataset_gamic = []
    datatype_odimbirds = []
    dataset_odimbirds = []
    datatype_mfcfradial = []
    dataset_mfcfradial = []
    datatype_nexrad2 = []
    dataset_nexrad2 = []
    datatype_gecsx = []
    dataset_gecsx = []
    product_gecsx = []
    datatype_cfradialpyrad = []
    dataset_cfradialpyrad = []
    product_cfradialpyrad = []
    datatype_cfradial = []
    dataset_cfradial = []
    datatype_skyecho = []
    dataset_skyecho = []
    datatype_cfradial2 = []
    dataset_cfradial2 = []
    datatype_cf1 = []
    dataset_cf1 = []
    datatype_odimpyrad = []
    dataset_odimpyrad = []
    product_odimpyrad = []
    datatype_icon = []
    datatype_rad4alpicon = []
    datatype_cfradialicon = []
    dataset_cfradialicon = []
    datatype_dem = []
    datatype_rad4alpdem = []
    datatype_rad4alphydro = []
    datatype_rad4alpDoppler = []
    datatype_rad4alpgrid = []
    datatype_rad4alpgif = []
    datatype_rad4alpbin = []
    datatype_mfbin = []
    dataset_mfbin = []
    datatype_mfdat = []
    dataset_mfdat = []
    datatype_mfpng = []
    dataset_mfpng = []
    datatype_mfgrib = []
    dataset_mfgrib = []
    datatype_mfcf = []
    dataset_mfcf = []
    datatype_satgrid = []
    datatype_rad4alpIQ = []
    datatype_mxpol = []
    datatype_pyradgrid = []
    dataset_pyradgrid = []
    product_pyradgrid = []
    datatype_odimpyradgrid = []
    dataset_odimpyradgrid = []
    product_odimpyradgrid = []
    datatype_odimgrid = []
    dataset_odimgrid = []
    datatype_knmih5grid = []
    dataset_knmih5grid = []
    datatype_psr = []
    datatype_psrspectra = []
    datatype_netcdfspectra = []
    dataset_netcdfspectra = []
    product_netcdfspectra = []

    for datatypedescr in datatypesdescr:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr
        )
        if datagroup == "RAINBOW":
            datatype_rainbow.append(datatype)
        elif datagroup == "RAD4ALP":
            datatype_rad4alp.append(datatype)
        elif datagroup == "ODIM":
            datatype_odim.append(datatype)
            dataset_odim.append(dataset)
        elif datagroup == "ODIMBIRDS":
            datatype_odimbirds.append(datatype)
            dataset_odimbirds.append(dataset)
        elif datagroup == "GAMIC":
            datatype_gamic.append(datatype)
            dataset_gamic.append(dataset)
        elif datagroup == "MFCFRADIAL":
            datatype_mfcfradial.append(datatype)
            dataset_mfcfradial.append(dataset)
        elif datagroup == "NEXRADII":
            datatype_nexrad2.append(datatype)
            dataset_nexrad2.append(dataset)
        elif datagroup == "GECSX":
            datatype_gecsx.append(datatype)
            dataset_gecsx.append(dataset)
            product_gecsx.append(product)
        elif datagroup == "CFRADIALPYRAD":
            datatype_cfradialpyrad.append(datatype)
            dataset_cfradialpyrad.append(dataset)
            product_cfradialpyrad.append(product)
        elif datagroup == "CFRADIAL":
            datatype_cfradial.append(datatype)
            dataset_cfradial.append(dataset)
        elif datagroup == "SKYECHO":
            datatype_skyecho.append(datatype)
            dataset_skyecho.append(dataset)
        elif datagroup == "CFRADIAL2":
            datatype_cfradial2.append(datatype)
            dataset_cfradial2.append(dataset)
        elif datagroup == "CF1":
            datatype_cf1.append(datatype)
            dataset_cf1.append(dataset)
        elif datagroup == "ODIMPYRAD":
            datatype_odimpyrad.append(datatype)
            dataset_odimpyrad.append(dataset)
            product_odimpyrad.append(product)
        elif datagroup == "COSMO":
            datatype_icon.append(datatype)
        elif datagroup == "RAD4ALPCOSMO":
            datatype_rad4alpicon.append(datatype)
        elif datagroup == "CFRADIALCOSMO":
            datatype_cfradialicon.append(datatype)
            dataset_cfradialicon.append(dataset)
        elif datagroup == "DEM":
            datatype_dem.append(datatype)
        elif datagroup == "RAD4ALPDEM":
            datatype_rad4alpdem.append(datatype)
        elif datagroup == "RAD4ALPHYDRO":
            datatype_rad4alphydro.append(datatype)
        elif datagroup == "RAD4ALPDOPPLER":
            datatype_rad4alpDoppler.append(datatype)
        elif datagroup == "MXPOL":
            datatype_mxpol.append(datatype)
        elif datagroup == "RAD4ALPGRID":
            datatype_rad4alpgrid.append(datatype)
        elif datagroup == "RAD4ALPGIF":
            datatype_rad4alpgif.append(datatype)
        elif datagroup == "RAD4ALPBIN":
            datatype_rad4alpbin.append(datatype)
        elif datagroup == "MFBIN":
            datatype_mfbin.append(datatype)
            dataset_mfbin.append(dataset)
        elif datagroup == "MFDAT":
            datatype_mfdat.append(datatype)
            dataset_mfdat.append(dataset)
        elif datagroup == "MFPNG":
            datatype_mfpng.append(datatype)
            dataset_mfpng.append(dataset)
        elif datagroup == "MFGRIB":
            datatype_mfgrib.append(datatype)
            dataset_mfgrib.append(dataset)
        elif datagroup == "MFCF":
            datatype_mfcf.append(datatype)
            dataset_mfcf.append(dataset)
        elif datagroup == "RAD4ALPIQ":
            datatype_rad4alpIQ.append(datatype)
        elif datagroup == "PYRADGRID":
            datatype_pyradgrid.append(datatype)
            dataset_pyradgrid.append(dataset)
            product_pyradgrid.append(product)
        elif datagroup == "ODIMPYRADGRID":
            datatype_odimpyradgrid.append(datatype)
            dataset_odimpyradgrid.append(dataset)
            product_odimpyradgrid.append(product)
        elif datagroup == "ODIMGRID":
            datatype_odimgrid.append(datatype)
            dataset_odimgrid.append(dataset)
        elif datagroup == "KNMIH5GRID":
            datatype_knmih5grid.append(datatype)
            dataset_knmih5grid.append(dataset)
        elif datagroup == "SATGRID":
            datatype_satgrid.append(datatype)
        elif datagroup == "PSR":
            datatype_psr.append(datatype)
        elif datagroup == "PSRSPECTRA":
            datatype_psrspectra.append(datatype)
        elif datagroup == "NETCDFSPECTRA":
            datatype_netcdfspectra.append(datatype)
            dataset_netcdfspectra.append(dataset)
            product_netcdfspectra.append(product)

    ind_rad = int(radarnr[5:8]) - 1

    ndatatypes_rainbow = len(datatype_rainbow)
    ndatatypes_rad4alp = len(datatype_rad4alp)
    ndatatypes_odim = len(datatype_odim)
    ndatatypes_odimbirds = len(datatype_odimbirds)
    ndatatypes_gamic = len(datatype_gamic)
    ndatatypes_mfcfradial = len(datatype_mfcfradial)
    ndatatypes_nexrad2 = len(datatype_nexrad2)
    ndatatypes_cfradial = len(datatype_cfradial)
    ndatatypes_skyecho = len(datatype_skyecho)
    ndatatypes_cfradialpyrad = len(datatype_cfradialpyrad)
    ndatatypes_cfradial2 = len(datatype_cfradial2)
    ndatatypes_cf1 = len(datatype_cf1)
    ndatatypes_odimpyrad = len(datatype_odimpyrad)
    ndatatypes_icon = len(datatype_icon)
    ndatatypes_rad4alpicon = len(datatype_rad4alpicon)
    ndatatypes_cfradialicon = len(datatype_cfradialicon)
    ndatatypes_dem = len(datatype_dem)
    ndatatypes_rad4alpdem = len(datatype_rad4alpdem)
    ndatatypes_rad4alphydro = len(datatype_rad4alphydro)
    ndatatypes_rad4alpDoppler = len(datatype_rad4alpDoppler)
    ndatatypes_mxpol = len(datatype_mxpol)
    ndatatypes_rad4alpgrid = len(datatype_rad4alpgrid)
    ndatatypes_rad4alpgif = len(datatype_rad4alpgif)
    ndatatypes_rad4alpbin = len(datatype_rad4alpbin)
    ndatatypes_mfbin = len(datatype_mfbin)
    ndatatypes_mfpng = len(datatype_mfpng)
    ndatatypes_mfgrib = len(datatype_mfgrib)
    ndatatypes_mfcf = len(datatype_mfcf)
    ndatatypes_mfdat = len(datatype_mfdat)
    ndatatypes_satgrid = len(datatype_satgrid)
    ndatatypes_rad4alpIQ = len(datatype_rad4alpIQ)
    ndatatypes_pyradgrid = len(datatype_pyradgrid)
    ndatatypes_odimpyradgrid = len(datatype_odimpyradgrid)
    ndatatypes_odimgrid = len(datatype_odimgrid)
    ndatatypes_knmih5grid = len(datatype_knmih5grid)
    ndatatypes_psr = len(datatype_psr)
    ndatatypes_psrspectra = len(datatype_psrspectra)
    ndatatypes_netcdfspectra = len(datatype_netcdfspectra)
    ndatatypes_gecsx = len(datatype_gecsx)

    rmin = None
    rmax = None
    elmin = None
    elmax = None
    azmin = None
    azmax = None
    if cfg["rmin"] is not None:
        rmin = cfg["rmin"][ind_rad]
    if cfg["rmax"] is not None:
        rmax = cfg["rmax"][ind_rad]
    if cfg["elmin"] is not None:
        elmin = cfg["elmin"][ind_rad]
    if cfg["elmax"] is not None:
        elmax = cfg["elmax"][ind_rad]
    if cfg["azmin"] is not None:
        azmin = cfg["azmin"][ind_rad]
    if cfg["azmax"] is not None:
        azmax = cfg["azmax"][ind_rad]

    radar = None
    if ndatatypes_rainbow > 0 and _WRADLIB_AVAILABLE:
        radar = merge_scans_rainbow(
            cfg["datapath"][ind_rad],
            cfg["ScanList"][ind_rad],
            voltime,
            cfg["ScanPeriod"],
            datatype_rainbow,
            cfg,
            radarnr=radarnr,
        )

    elif ndatatypes_rad4alp > 0:
        radar = merge_scans_rad4alp(
            cfg["datapath"][ind_rad],
            cfg["ScanList"][ind_rad],
            cfg["RadarName"][ind_rad],
            cfg["RadarRes"][ind_rad],
            voltime,
            datatype_rad4alp,
            cfg,
            ind_rad=ind_rad,
        )

    elif ndatatypes_odim > 0:
        try:
            radar_name = cfg["RadarName"][ind_rad]
            radar_res = cfg["RadarRes"][ind_rad]
        except TypeError:
            radar_name = None
            radar_res = None
        radar = merge_scans_odim(
            cfg["datapath"][ind_rad],
            cfg["ScanList"][ind_rad],
            radar_name,
            radar_res,
            voltime,
            datatype_odim,
            dataset_odim,
            cfg,
            ind_rad=ind_rad,
        )

    elif ndatatypes_odimbirds > 0:
        try:
            radar_name = cfg["RadarName"][ind_rad]
            radar_res = cfg["RadarRes"][ind_rad]
        except TypeError:
            radar_name = None
            radar_res = None
        radar = merge_scans_odimbirds(
            cfg["datapath"][ind_rad],
            cfg["ScanList"][ind_rad],
            radar_name,
            radar_res,
            voltime,
            datatype_odimbirds,
            dataset_odimbirds,
            cfg,
            ind_rad=ind_rad,
        )

    elif ndatatypes_gamic > 0:
        try:
            radar_name = cfg["RadarName"][ind_rad]
            radar_res = cfg["RadarRes"][ind_rad]
        except TypeError:
            radar_name = None
            radar_res = None
        radar = merge_scans_gamic(
            cfg["datapath"][ind_rad],
            cfg["ScanList"][ind_rad],
            radar_name,
            radar_res,
            voltime,
            datatype_gamic,
            dataset_gamic,
            cfg,
            ind_rad=ind_rad,
        )

    elif ndatatypes_nexrad2 > 0:
        try:
            radar_name = cfg["RadarName"][ind_rad]
            radar_res = cfg["RadarRes"][ind_rad]
        except TypeError:
            radar_name = None
            radar_res = None
        radar = merge_scans_nexrad2(
            cfg["datapath"][ind_rad],
            cfg["ScanList"][ind_rad],
            voltime,
            datatype_nexrad2,
            dataset_nexrad2,
            cfg,
            ind_rad=ind_rad,
        )

    elif ndatatypes_cfradial > 0:
        try:
            radar_name = cfg["RadarName"][ind_rad]
            radar_res = cfg["RadarRes"][ind_rad]
        except TypeError:
            radar_name = None
            radar_res = None
        radar = merge_scans_cfradial(
            cfg["datapath"][ind_rad],
            cfg["ScanList"][ind_rad],
            radar_name,
            radar_res,
            voltime,
            datatype_cfradial,
            dataset_cfradial,
            cfg,
            ind_rad=ind_rad,
        )

    elif ndatatypes_skyecho > 0:
        try:
            radar_name = cfg["RadarName"][ind_rad]
            radar_res = cfg["RadarRes"][ind_rad]
        except TypeError:
            radar_name = None
            radar_res = None
        radar = merge_scans_skyecho(
            cfg["datapath"][ind_rad],
            cfg["ScanList"][ind_rad],
            voltime,
            datatype_skyecho,
            dataset_skyecho,
            cfg,
            ind_rad=ind_rad,
        )

    elif ndatatypes_cfradial2 > 0:
        try:
            radar_name = cfg["RadarName"][ind_rad]
            radar_res = cfg["RadarRes"][ind_rad]
        except TypeError:
            radar_name = None
            radar_res = None
        radar = merge_scans_cfradial2(
            cfg["datapath"][ind_rad],
            cfg["ScanList"][ind_rad],
            radar_name,
            radar_res,
            voltime,
            datatype_cfradial2,
            dataset_cfradial2,
            cfg,
            ind_rad=ind_rad,
        )

    elif ndatatypes_cf1 > 0:
        try:
            radar_name = cfg["RadarName"][ind_rad]
            radar_res = cfg["RadarRes"][ind_rad]
        except TypeError:
            radar_name = None
            radar_res = None

        radar = merge_scans_cf1(
            cfg["datapath"][ind_rad],
            cfg["ScanList"][ind_rad],
            radar_name,
            radar_res,
            voltime,
            datatype_cf1,
            dataset_cf1,
            cfg,
            ind_rad=ind_rad,
        )

    elif ndatatypes_mxpol > 0:
        radar = merge_scans_mxpol(
            cfg["datapath"][ind_rad],
            cfg["ScanList"][ind_rad],
            voltime,
            datatype_mxpol,
            cfg,
            ind_rad=ind_rad,
        )

    elif ndatatypes_rad4alpgrid > 0:
        radar = merge_fields_rad4alp_grid(
            voltime, datatype_rad4alpgrid, cfg, ind_rad=ind_rad
        )

    elif ndatatypes_psrspectra > 0:
        radar = merge_scans_psr_spectra(
            cfg["datapath"][ind_rad],
            cfg["psrpath"][ind_rad],
            cfg["ScanList"][ind_rad],
            voltime,
            cfg["ScanPeriod"],
            datatype_psrspectra,
            cfg,
            radarnr=radarnr,
        )

    elif ndatatypes_rad4alpIQ > 0:
        radar = merge_scans_iq_rad4alp(
            cfg["datapath"][ind_rad],
            cfg["iqpath"][ind_rad],
            cfg["ScanList"][ind_rad],
            cfg["RadarName"][ind_rad],
            cfg["RadarRes"][ind_rad],
            voltime,
            datatype_rad4alpIQ,
            cfg,
            ind_rad=ind_rad,
        )

    elif ndatatypes_mfcfradial > 0:
        warn("MFCFRADIAL type is kept for legacy purposes but is not needed anymore")
        warn("Please use CFRADIAL type for all CFRadial files,")
        warn("with path_convention ODIM and DataTypeIDInFiles matching your data")
        radar = merge_scans_mfcfradial(
            cfg["datapath"][ind_rad],
            cfg["ScanList"][ind_rad],
            voltime,
            datatype_mfcfradial,
            dataset_mfcfradial,
            cfg,
            ind_rad=ind_rad,
        )

    # add other radar object files
    if ndatatypes_cfradialpyrad > 0:
        radar_aux = merge_fields_pyrad(
            cfg["loadbasepath"][ind_rad],
            cfg["loadname"][ind_rad],
            voltime,
            datatype_cfradialpyrad,
            dataset_cfradialpyrad,
            product_cfradialpyrad,
            rng_min=rmin,
            rng_max=rmax,
            ele_min=elmin,
            ele_max=elmax,
            azi_min=azmin,
            azi_max=azmax,
        )
        radar = add_field(radar, radar_aux)

    if ndatatypes_odimpyrad > 0:
        radar_aux = merge_fields_pyrad(
            cfg["loadbasepath"][ind_rad],
            cfg["loadname"][ind_rad],
            voltime,
            datatype_odimpyrad,
            dataset_odimpyrad,
            product_odimpyrad,
            rng_min=rmin,
            rng_max=rmax,
            ele_min=elmin,
            ele_max=elmax,
            azi_min=azmin,
            azi_max=azmax,
            termination=".h*",
        )
        radar = add_field(radar, radar_aux)

    if ndatatypes_gecsx > 0:
        radar_aux = merge_fields_gecsx(
            cfg["gecsxbasepath"][ind_rad],
            cfg["gecsxname"][ind_rad],
            datatype_gecsx,
            dataset_gecsx,
            product_gecsx,
            rng_min=rmin,
            rng_max=rmax,
            ele_min=elmin,
            ele_max=elmax,
            azi_min=azmin,
            azi_max=azmax,
        )
        radar = add_field(radar, radar_aux)

    if ndatatypes_cfradialicon > 0:
        radar_aux = merge_fields_pyradicon(
            cfg["iconpath"][ind_rad],
            voltime,
            datatype_cfradialicon,
            dataset_cfradialicon,
            cfg,
            rng_min=rmin,
            rng_max=rmax,
            ele_min=elmin,
            ele_max=elmax,
            azi_min=azmin,
            azi_max=azmax,
            termination=".nc",
        )
        radar = add_field(radar, radar_aux)

    # add rainbow ray data from psr files
    if ndatatypes_psr > 0:
        radar_aux = merge_scans_psr(
            cfg["datapath"][ind_rad],
            cfg["psrpath"][ind_rad],
            cfg["ScanList"][ind_rad],
            voltime,
            cfg["ScanPeriod"],
            datatype_psr,
            cfg,
            radarnr=radarnr,
        )
        radar = add_field(radar, radar_aux)

    # add other radar spectra object files
    if ndatatypes_netcdfspectra > 0:
        radar_aux = merge_fields_pyrad_spectra(
            cfg["loadbasepath"][ind_rad],
            cfg["loadname"][ind_rad],
            voltime,
            datatype_netcdfspectra,
            dataset_netcdfspectra,
            product_netcdfspectra,
            rng_min=rmin,
            rng_max=rmax,
            ele_min=elmin,
            ele_max=elmax,
            azi_min=azmin,
            azi_max=azmax,
        )
        radar = add_field(radar, radar_aux)

    # add other grid object files
    if ndatatypes_rad4alpgif > 0:
        radar_aux = merge_fields_rad4alp_grid(
            voltime, datatype_rad4alpgif, cfg, ind_rad=ind_rad, ftype="gif"
        )
        if radar_aux is not None:
            if radar is not None:
                radar = merge_grids(radar, radar_aux)
            else:
                radar = radar_aux

    if ndatatypes_rad4alpbin > 0:
        radar_aux = merge_fields_rad4alp_grid(
            voltime, datatype_rad4alpbin, cfg, ind_rad=ind_rad, ftype="bin"
        )
        if radar_aux is not None:
            if radar is not None:
                radar = merge_grids(radar, radar_aux)
            else:
                radar = radar_aux

    if ndatatypes_pyradgrid > 0:
        radar_aux = merge_fields_pyradgrid(
            cfg["loadbasepath"][ind_rad],
            cfg["loadname"][ind_rad],
            voltime,
            datatype_pyradgrid,
            dataset_pyradgrid,
            product_pyradgrid,
            cfg,
            termination=".nc",
        )
        if radar_aux is not None:
            if radar is not None:
                radar = merge_grids(radar, radar_aux)
            else:
                radar = radar_aux

    if ndatatypes_odimpyradgrid > 0:
        radar_aux = merge_fields_pyradgrid(
            cfg["loadbasepath"][ind_rad],
            cfg["loadname"][ind_rad],
            voltime,
            datatype_odimpyradgrid,
            dataset_odimpyradgrid,
            product_odimpyradgrid,
            cfg,
            termination=".h*",
        )
        if radar_aux is not None:
            if radar is not None:
                radar = merge_grids(radar, radar_aux)
            else:
                radar = radar_aux

    if ndatatypes_odimgrid > 0:
        radar_aux = merge_scans_odimgrid(
            cfg["datapath"][ind_rad],
            cfg["ScanList"][ind_rad],
            voltime,
            datatype_odimgrid,
            dataset_odimgrid,
            cfg,
        )
        if radar_aux is not None:
            if radar is not None:
                radar = merge_grids(radar, radar_aux)
            else:
                radar = radar_aux

    if ndatatypes_knmih5grid > 0:
        radar_aux = merge_scans_knmih5_grid(
            cfg["datapath"][ind_rad],
            cfg["ScanList"][ind_rad],
            voltime,
            datatype_knmih5grid,
            dataset_knmih5grid,
            cfg,
        )
        if radar_aux is not None:
            if radar is not None:
                radar = merge_grids(radar, radar_aux)
            else:
                radar = radar_aux

    if ndatatypes_satgrid > 0:
        radar_aux = merge_fields_sat_grid(
            voltime, datatype_satgrid, cfg, ind_rad=ind_rad
        )
        if radar_aux is not None:
            if radar is not None:
                radar = merge_grids(radar, radar_aux)
            else:
                radar = radar_aux

    if ndatatypes_mfbin > 0:
        radar_aux = merge_fields_mf_grid(
            voltime,
            datatype_mfbin,
            dataset_mfbin,
            cfg["ScanList"][ind_rad],
            cfg,
            ind_rad=ind_rad,
            ftype="bin",
        )
        if radar_aux is not None:
            if radar is not None:
                radar = merge_grids(radar, radar_aux)
            else:
                radar = radar_aux

    if ndatatypes_mfpng > 0:
        radar_aux = merge_fields_mf_grid(
            voltime,
            datatype_mfpng,
            dataset_mfpng,
            cfg["ScanList"][ind_rad],
            cfg,
            ind_rad=ind_rad,
            ftype="png",
        )
        if radar_aux is not None:
            if radar is not None:
                radar = merge_grids(radar, radar_aux)
            else:
                radar = radar_aux

    if ndatatypes_mfgrib > 0:
        radar_aux = merge_fields_mf_grid(
            voltime,
            datatype_mfgrib,
            dataset_mfgrib,
            cfg["ScanList"][ind_rad],
            cfg,
            ind_rad=ind_rad,
            ftype="grib",
        )
        if radar_aux is not None:
            if radar is not None:
                radar = merge_grids(radar, radar_aux)
            else:
                radar = radar_aux

    if ndatatypes_mfdat > 0:
        radar_aux = merge_fields_mf_grid(
            voltime,
            datatype_mfdat,
            dataset_mfdat,
            cfg["ScanList"][ind_rad],
            cfg,
            ind_rad=ind_rad,
            ftype="dat",
        )
        if radar_aux is not None:
            if radar is not None:
                radar = merge_grids(radar, radar_aux)
            else:
                radar = radar_aux

    if ndatatypes_mfcf > 0:
        radar_aux = merge_fields_mf_grid(
            voltime,
            datatype_mfcf,
            dataset_mfcf,
            cfg["ScanList"][ind_rad],
            cfg,
            ind_rad=ind_rad,
            ftype="nc",
        )
        if radar_aux is not None:
            if radar is not None:
                radar = merge_grids(radar, radar_aux)
            else:
                radar = radar_aux

    # add COSMO files to the radar field
    if ndatatypes_icon > 0 and _WRADLIB_AVAILABLE:
        radar_aux = merge_scans_icon(voltime, datatype_icon, cfg, ind_rad=ind_rad)
        radar = add_field(radar, radar_aux)

    elif ndatatypes_rad4alpicon > 0:
        if (cfg["RadarRes"][ind_rad] is None) or (cfg["RadarName"][ind_rad] is None):
            raise ValueError(
                "ERROR: Radar Name and Resolution "
                + "not specified in config file. "
                + "Unable to load rad4alp COSMO data"
            )

        for dt_rad4alpicon in datatype_rad4alpicon:
            radar_aux = merge_scans_icon_rad4alp(
                voltime, dt_rad4alpicon, cfg, ind_rad=ind_rad
            )
            if radar is None:
                radar = radar_aux
                continue
            if radar_aux is None:
                continue
            for field_name in radar_aux.fields.keys():
                try:
                    radar.add_field(field_name, radar_aux.fields[field_name])
                except (ValueError, KeyError) as ee:
                    warn(
                        "Unable to add field '" + field_name + "' to radar object"
                        ": (%s)" % str(ee)
                    )

    # add DEM files to the radar field
    if ndatatypes_dem > 0 and _WRADLIB_AVAILABLE:
        radar_aux = merge_scans_dem(
            cfg["dempath"][ind_rad],
            cfg["ScanList"][ind_rad],
            datatype_dem,
            rng_min=rmin,
            rng_max=rmax,
            ele_min=elmin,
            ele_max=elmax,
            azi_min=azmin,
            azi_max=azmax,
        )
        radar = add_field(radar, radar_aux)

    elif ndatatypes_rad4alpdem > 0:
        if (cfg["RadarRes"][ind_rad] is None) or (cfg["RadarName"][ind_rad] is None):
            raise ValueError(
                "ERROR: Radar Name and Resolution "
                + "not specified in config file. "
                + "Unable to load rad4alp DEM data"
            )

        if cfg["RadarRes"][ind_rad] != "L":
            raise ValueError(
                "ERROR: DEM files only available for rad4alp PL data. "
                + "Current radar "
                + cfg["RadarName"][ind_rad]
                + cfg["RadarRes"][ind_rad]
            )

        for dt_rad4alpdem in datatype_rad4alpdem:
            radar_aux = merge_scans_dem_rad4alp(
                voltime, dt_rad4alpdem, cfg, ind_rad=ind_rad
            )
            if radar is None:
                radar = radar_aux
                continue
            if radar_aux is None:
                continue
            for field_name in radar_aux.fields.keys():
                try:
                    radar.add_field(field_name, radar_aux.fields[field_name])
                except (ValueError, KeyError) as ee:
                    warn(
                        "Unable to add field '" + field_name + "' to radar object"
                        ": (%s)" % str(ee)
                    )

    # add rad4alp radar data
    if ndatatypes_rad4alphydro > 0:
        if (cfg["RadarRes"][ind_rad] is None) or (cfg["RadarName"][ind_rad] is None):
            raise ValueError(
                "ERROR: Radar Name and Resolution "
                + "not specified in config file. "
                + "Unable to load rad4alp hydro data"
            )

        for dt_rad4alphydro in datatype_rad4alphydro:
            radar_aux = merge_scans_other_rad4alp(
                voltime, dt_rad4alphydro, cfg, ind_rad=ind_rad
            )
            if radar is None:
                radar = radar_aux
                continue
            if radar_aux is None:
                continue
            if radar_aux is not None:
                for field_name in radar_aux.fields.keys():
                    try:
                        radar.add_field(field_name, radar_aux.fields[field_name])
                    except (ValueError, KeyError) as ee:
                        warn(
                            "Unable to add field '" + field_name + "' to radar object"
                            ": (%s)" % str(ee)
                        )

    if ndatatypes_rad4alpDoppler > 0:
        if (cfg["RadarRes"][ind_rad] is None) or (cfg["RadarName"][ind_rad] is None):
            raise ValueError(
                "ERROR: Radar Name and Resolution "
                + "not specified in config file. "
                + "Unable to load rad4alp dealiased Doppler data"
            )

        for dt_rad4alpDoppler in datatype_rad4alpDoppler:
            radar_aux = merge_scans_other_rad4alp(
                voltime, dt_rad4alpDoppler, cfg, ind_rad=ind_rad
            )
            if radar is None:
                radar = radar_aux
                continue
            if radar_aux is None:
                continue
            if radar_aux is not None:
                for field_name in radar_aux.fields.keys():
                    try:
                        radar.add_field(field_name, radar_aux.fields[field_name])
                    except (ValueError, KeyError) as ee:
                        warn(
                            "Unable to add field '" + field_name + "' to radar object"
                            ": (%s)" % str(ee)
                        )

    if radar is None:
        return radar

    # if it is specified, get the position from the config file
    if "RadarPosition" in cfg:
        if "latitude" in cfg["RadarPosition"]:
            radar.latitude["data"][0] = cfg["RadarPosition"]["latitude"][ind_rad]
        if "longitude" in cfg["RadarPosition"]:
            radar.longitude["data"][0] = cfg["RadarPosition"]["longitude"][ind_rad]
        if "altitude" in cfg["RadarPosition"]:
            radar.altitude["data"][0] = cfg["RadarPosition"]["altitude"][ind_rad]
        radar.init_gate_longitude_latitude()
        radar.init_gate_altitude()

    # get ray angle resolution
    if "ray_angle_res" in cfg:
        ray_angle_res = pyart.config.get_metadata("ray_angle_res")
        ray_angle_res["data"] = cfg["ray_angle_res"][ind_rad] + np.zeros(
            radar.nsweeps, dtype=np.float32
        )
        radar.ray_angle_res = ray_angle_res

        rays_are_indexed = pyart.config.get_metadata("rays_are_indexed")
        rays_are_indexed["data"] = np.ones(radar.nsweeps, dtype=bool)
        radar.rays_are_indexed = rays_are_indexed

    # get instrument parameters from the config file
    if "frequency" in cfg:
        if radar.instrument_parameters is None:
            frequency = pyart.config.get_metadata("frequency")
            frequency["data"] = np.array([cfg["frequency"][ind_rad]], dtype=np.float32)
            radar.instrument_parameters = {"frequency": frequency}
        elif "frequency" not in radar.instrument_parameters:
            frequency = pyart.config.get_metadata("frequency")
            frequency["data"] = np.array([cfg["frequency"][ind_rad]], dtype=np.float32)
            radar.instrument_parameters.update({"frequency": frequency})
        else:
            radar.instrument_parameters["frequency"]["data"][0] = cfg["frequency"][
                ind_rad
            ]

    if "radar_beam_width_h" in cfg:
        if radar.instrument_parameters is None:
            beamwidth = pyart.config.get_metadata("radar_beam_width_h")
            beamwidth["data"] = np.array(
                [cfg["radar_beam_width_h"][ind_rad]], dtype=np.float32
            )
            radar.instrument_parameters = {"radar_beam_width_h": beamwidth}
        elif "radar_beam_width_h" not in radar.instrument_parameters:
            beamwidth = pyart.config.get_metadata("radar_beam_width_h")
            beamwidth["data"] = np.array(
                [cfg["radar_beam_width_h"][ind_rad]], dtype=np.float32
            )
            radar.instrument_parameters.update({"radar_beam_width_h": beamwidth})
        else:
            radar.instrument_parameters["radar_beam_width_h"]["data"][0] = cfg[
                "radar_beam_width_h"
            ][ind_rad]

    if "radar_beam_width_v" in cfg:
        if radar.instrument_parameters is None:
            beamwidth = pyart.config.get_metadata("radar_beam_width_v")
            beamwidth["data"] = np.array(
                [cfg["radar_beam_width_v"][ind_rad]], dtype=np.float32
            )
            radar.instrument_parameters = {"radar_beam_width_v": beamwidth}
        elif "radar_beam_width_v" not in radar.instrument_parameters:
            beamwidth = pyart.config.get_metadata("radar_beam_width_v")
            beamwidth["data"] = np.array(
                [cfg["radar_beam_width_v"][ind_rad]], dtype=np.float32
            )
            radar.instrument_parameters.update({"radar_beam_width_v": beamwidth})
        else:
            radar.instrument_parameters["radar_beam_width_v"]["data"][0] = cfg[
                "radar_beam_width_v"
            ][ind_rad]

    if "AntennaGainH" in cfg:
        if radar.instrument_parameters is None:
            AntennaGainH = pyart.config.get_metadata("radar_antenna_gain_h")
            AntennaGainH["data"] = np.array(
                [cfg["AntennaGainH"][ind_rad]], dtype=np.float32
            )
            radar.instrument_parameters = {"radar_antenna_gain_h": AntennaGainH}
        elif "radar_antenna_gain_h" not in radar.instrument_parameters:
            AntennaGainH = pyart.config.get_metadata("radar_antenna_gain_h")
            AntennaGainH["data"] = np.array(
                [cfg["AntennaGainH"][ind_rad]], dtype=np.float32
            )
            radar.instrument_parameters.update({"radar_antenna_gain_h": AntennaGainH})
        else:
            radar.instrument_parameters["radar_antenna_gain_h"]["data"][0] = cfg[
                "AntennaGainH"
            ][ind_rad]

    if "AntennaGainV" in cfg:
        if radar.instrument_parameters is None:
            AntennaGainV = pyart.config.get_metadata("radar_antenna_gain_v")
            AntennaGainV["data"] = np.array(
                [cfg["AntennaGainV"][ind_rad]], dtype=np.float32
            )
            radar.instrument_parameters = {"radar_antenna_gain_v": AntennaGainV}
        elif "radar_antenna_gain_v" not in radar.instrument_parameters:
            AntennaGainV = pyart.config.get_metadata("radar_antenna_gain_v")
            AntennaGainV["data"] = np.array(
                [cfg["AntennaGainV"][ind_rad]], dtype=np.float32
            )
            radar.instrument_parameters.update({"radar_antenna_gain_v": AntennaGainV})
        else:
            radar.instrument_parameters["radar_antenna_gain_v"]["data"][0] = cfg[
                "AntennaGainV"
            ][ind_rad]

    # Assumes uniform pulse width in all radar volume
    if "pulse_width" in cfg:
        if radar.instrument_parameters is None:
            pulse_width = pyart.config.get_metadata("pulse_width")
            pulse_width["data"] = cfg["pulse_width"][ind_rad] * np.array(
                radar.nrays, dtype=np.float32
            )
            radar.instrument_parameters = {"pulse_width": pulse_width}
        elif "pulse_width" not in radar.instrument_parameters:
            pulse_width = pyart.config.get_metadata("pulse_width")
            pulse_width["data"] = cfg["pulse_width"][ind_rad] * np.array(
                radar.nrays, dtype=np.float32
            )
            radar.instrument_parameters.update({"pulse_width": pulse_width})
        else:
            radar.instrument_parameters["pulse_width"]["data"] = cfg["pulse_width"][
                ind_rad
            ] * np.array(radar.nrays, dtype=np.float32)

    # Assumes uniform nyquist velocity in all radar volume
    if "nyquist_velocity" in cfg:
        if radar.instrument_parameters is None:
            nyquist_velocity = pyart.config.get_metadata("nyquist_velocity")
            nyquist_velocity["data"] = cfg["nyquist_velocity"][ind_rad] * np.array(
                radar.nrays, dtype=np.float32
            )
            radar.instrument_parameters = {"nyquist_velocity": nyquist_velocity}
        elif "nyquist_velocity" not in radar.instrument_parameters:
            nyquist_velocity = pyart.config.get_metadata("nyquist_velocity")
            nyquist_velocity["data"] = cfg["nyquist_velocity"][ind_rad] * np.array(
                radar.nrays, dtype=np.float32
            )
            radar.instrument_parameters.update({"nyquist_velocity": nyquist_velocity})
        else:
            radar.instrument_parameters["nyquist_velocity"]["data"] = cfg[
                "nyquist_velocity"
            ][ind_rad] * np.array(radar.nrays, dtype=np.float32)

    # Get calibration parameters from config file
    if "radconsth" in cfg:
        if radar.radar_calibration is None:
            radconsth = pyart.config.get_metadata("calibration_constant_hh")
            radconsth["data"] = np.array([cfg["radconsth"][ind_rad]], dtype=np.float32)
            radar.radar_calibration = {"calibration_constant_hh": radconsth}
        elif "calibration_constant_hh" not in radar.radar_calibration:
            radconsth = pyart.config.get_metadata("calibration_constant_hh")
            radconsth["data"] = np.array([cfg["radconsth"][ind_rad]], dtype=np.float32)
            radar.radar_calibration.update({"calibration_constant_hh": radconsth})
        else:
            radar.radar_calibration["calibration_constant_hh"]["data"][0] = cfg[
                "radconsth"
            ][ind_rad]

    if "radconstv" in cfg:
        if radar.radar_calibration is None:
            radconstv = pyart.config.get_metadata("calibration_constant_vv")
            radconstv["data"] = np.array([cfg["radconstv"][ind_rad]], dtype=np.float32)
            radar.radar_calibration = {"calibration_constant_vv": radconstv}
        elif "calibration_constant_vv" not in radar.radar_calibration:
            radconstv = pyart.config.get_metadata("calibration_constant_vv")
            radconstv["data"] = np.array([cfg["radconstv"][ind_rad]], dtype=np.float32)
            radar.radar_calibration.update({"calibration_constant_vv": radconstv})
        else:
            radar.radar_calibration["calibration_constant_vv"]["data"][0] = cfg[
                "radconstv"
            ][ind_rad]

    if "txpwrh" in cfg:
        if radar.radar_calibration is None:
            txpwrh = pyart.config.get_metadata("transmit_power_h")
            txpwrh["data"] = np.array([cfg["txpwrh"][ind_rad]], dtype=np.float32)
            radar.radar_calibration = {"transmit_power_h": txpwrh}
        elif "transmit_power_h" not in radar.radar_calibration:
            txpwrh = pyart.config.get_metadata("transmit_power_h")
            txpwrh["data"] = np.array([cfg["txpwrh"][ind_rad]], dtype=np.float32)
            radar.radar_calibration.update({"transmit_power_h": txpwrh})
        else:
            radar.radar_calibration["transmit_power_h"]["data"][0] = cfg["txpwrh"][
                ind_rad
            ]

    if "txpwrv" in cfg:
        if radar.radar_calibration is None:
            txpwrv = pyart.config.get_metadata("transmit_power_v")
            txpwrv["data"] = np.array([cfg["txpwrv"][ind_rad]], dtype=np.float32)
            radar.radar_calibration = {"transmit_power_v": txpwrv}
        elif "transmit_power_v" not in radar.radar_calibration:
            txpwrv = pyart.config.get_metadata("transmit_power_v")
            txpwrv["data"] = np.array([cfg["txpwrv"][ind_rad]], dtype=np.float32)
            radar.radar_calibration.update({"transmit_power_v": txpwrv})
        else:
            radar.radar_calibration["transmit_power_v"]["data"][0] = cfg["txpwrv"][
                ind_rad
            ]

    if "attg" in cfg:
        if radar.radar_calibration is None:
            attg = pyart.config.get_metadata("path_attenuation")
            attg["data"] = np.array([cfg["attg"][ind_rad]], dtype=np.float32)
            radar.radar_calibration = {"path_attenuation": attg}
        elif "path_attenuation" not in radar.radar_calibration:
            attg = pyart.config.get_metadata("path_attenuation")
            attg["data"] = np.array([cfg["attg"][ind_rad]], dtype=np.float32)
            radar.radar_calibration.update({"path_attenuation": attg})
        else:
            radar.radar_calibration["path_attenuation"]["data"][0] = cfg["attg"][
                ind_rad
            ]

    if "mflossh" in cfg:
        if radar.radar_calibration is None:
            mflossh = pyart.config.get_metadata("matched_filter_loss_h")
            mflossh["data"] = np.array([cfg["mflossh"][ind_rad]], dtype=np.float32)
            radar.radar_calibration = {"matched_filter_loss_h": mflossh}
        elif "matched_filter_loss_h" not in radar.radar_calibration:
            mflossh = pyart.config.get_metadata("matched_filter_loss_h")
            mflossh["data"] = np.array([cfg["mflossh"][ind_rad]], dtype=np.float32)
            radar.radar_calibration.update({"matched_filter_loss_h": mflossh})
        else:
            radar.radar_calibration["matched_filter_loss_h"]["data"][0] = cfg[
                "mflossh"
            ][ind_rad]

    if "mflossv" in cfg:
        if radar.radar_calibration is None:
            mflossv = pyart.config.get_metadata("matched_filter_loss_v")
            mflossv["data"] = np.array([cfg["mflossv"][ind_rad]], dtype=np.float32)
            radar.radar_calibration = {"matched_filter_loss_v": mflossv}
        elif "matched_filter_loss_v" not in radar.radar_calibration:
            mflossv = pyart.config.get_metadata("matched_filter_loss_v")
            mflossv["data"] = np.array([cfg["mflossv"][ind_rad]], dtype=np.float32)
            radar.radar_calibration.update({"matched_filter_loss_v": mflossv})
        else:
            radar.radar_calibration["matched_filter_loss_v"]["data"][0] = cfg[
                "mflossv"
            ][ind_rad]

    if "dBADUtodBmh" in cfg:
        if radar.radar_calibration is None:
            dBADUtodBmh = pyart.config.get_metadata("dBADU_to_dBm_hh")
            dBADUtodBmh["data"] = np.array(
                [cfg["dBADUtodBmh"][ind_rad]], dtype=np.float32
            )
            radar.radar_calibration = {"dBADU_to_dBm_hh": dBADUtodBmh}
        elif "dBADU_to_dBm_hh" not in radar.radar_calibration:
            dBADUtodBmh = pyart.config.get_metadata("dBADU_to_dBm_hh")
            dBADUtodBmh["data"] = np.array(
                [cfg["dBADUtodBmh"][ind_rad]], dtype=np.float32
            )
            radar.radar_calibration.update({"dBADU_to_dBm_hh": dBADUtodBmh})
        else:
            radar.radar_calibration["dBADU_to_dBm_hh"]["data"][0] = cfg["dBADUtodBmh"][
                ind_rad
            ]

    if "dBADUtodBmv" in cfg:
        if radar.radar_calibration is None:
            dBADUtodBmv = pyart.config.get_metadata("dBADU_to_dBm_vv")
            dBADUtodBmv["data"] = np.array(
                [cfg["dBADUtodBmv"][ind_rad]], dtype=np.float32
            )
            radar.radar_calibration = {"dBADU_to_dBm_vv": dBADUtodBmv}
        elif "dBADU_to_dBm_vv" not in radar.radar_calibration:
            dBADUtodBmv = pyart.config.get_metadata("dBADU_to_dBm_vv")
            dBADUtodBmv["data"] = np.array(
                [cfg["dBADUtodBmv"][ind_rad]], dtype=np.float32
            )
            radar.radar_calibration.update({"dBADU_to_dBm_vv": dBADUtodBmv})
        else:
            radar.radar_calibration["dBADU_to_dBm_vv"]["data"][0] = cfg["dBADUtodBmv"][
                ind_rad
            ]

    # If radar name is empty substitude the one from the config file
    if not radar.metadata.get("instrument_name"):
        if "RadarName" in cfg:
            try:
                radar.metadata["instrument_name"] = cfg["RadarName"][ind_rad]
            except (IndexError, TypeError):
                pass

    return radar


def merge_scans_rainbow(
    basepath, scan_list, voltime, scan_period, datatype_list, cfg, radarnr="RADAR001"
):
    """
    merge rainbow scans

    Parameters
    ----------
    basepath : str
        base path of rad4alp radar data
    scan_list : list
        list of scans
    voltime: datetime object
        reference time of the scan
    scan_period : float
        time from reference time where to look for other scans data
    datatype_list : list
        lists of data types to get
    cfg : dict
        configuration dictionary
    radarnr : str
        radar identifier number

    Returns
    -------
    radar : Radar
        radar object

    """

    radar = merge_fields_rainbow(basepath, scan_list[0], voltime, datatype_list)

    # merge scans into a single radar instance
    nscans = len(scan_list)
    if nscans > 1:
        if (datatype_list[0] == "Nh") or (datatype_list[0] == "Nv"):
            datadescriptor = radarnr + ":RAINBOW:dBZ"
        else:
            datadescriptor = radarnr + ":RAINBOW:" + datatype_list[0]
        endtime = voltime + datetime.timedelta(minutes=scan_period)
        for scan in scan_list[1:]:
            filelist = get_file_list(
                datadescriptor, [voltime], [endtime], cfg, scan=scan
            )
            if not filelist:
                warn(
                    "ERROR: No data file found for scan '%s' "
                    "between %s and %s" % (scan, voltime, endtime)
                )
                continue
            scantime = get_datetime(filelist[0], datadescriptor)

            radar_aux = merge_fields_rainbow(basepath, scan, scantime, datatype_list)

            if radar_aux is None:
                continue

            if radar is None:
                radar = radar_aux
            else:
                radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    if radar is None:
        return radar

    ind_rad = int(radarnr[5:8]) - 1
    rmin = None
    rmax = None
    elmin = None
    elmax = None
    azmin = None
    azmax = None
    if cfg["rmin"] is not None:
        rmin = cfg["rmin"][ind_rad]
    if cfg["rmax"] is not None:
        rmax = cfg["rmax"][ind_rad]
    if cfg["elmin"] is not None:
        elmin = cfg["elmin"][ind_rad]
    if cfg["elmax"] is not None:
        elmax = cfg["elmax"][ind_rad]
    if cfg["azmin"] is not None:
        azmin = cfg["azmin"][ind_rad]
    if cfg["azmax"] is not None:
        azmax = cfg["azmax"][ind_rad]

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rmin,
        rng_max=rmax,
        ele_min=elmin,
        ele_max=elmax,
        azi_min=azmin,
        azi_max=azmax,
    )


def merge_scans_psr(
    basepath,
    basepath_psr,
    scan_list,
    voltime,
    scan_period,
    datatype_list,
    cfg,
    radarnr="RADAR001",
):
    """
    merge rainbow scans

    Parameters
    ----------
    basepath : str
        base path of rainbow radar data
    basepath_psr : str
        name of the base path where to find the PSR data
    scan_list : list
        list of scans
    voltime: datetime object
        reference time of the scan
    scan_period : float
        time from reference time where to look for other scans data
    datatype_list : list
        lists of data types to get
    cfg : dict
        configuration dictionary
    radarnr : str
        radar identifier number

    Returns
    -------
    radar : Radar
        radar object

    """
    datadescriptor = radarnr + ":RAINBOW:dBZ"
    endtime = voltime + datetime.timedelta(minutes=scan_period)

    radar = None
    for scan in scan_list:
        filelist = get_file_list(datadescriptor, [voltime], [endtime], cfg, scan=scan)

        if not filelist:
            warn(
                "ERROR: No data file found for scan '%s' "
                "between %s and %s" % (scan, voltime, endtime)
            )
            continue
        scantime = get_datetime(filelist[0], datadescriptor)

        radar_aux = merge_fields_psr(
            basepath,
            basepath_psr,
            scan,
            scantime,
            datatype_list,
            undo_txcorr=cfg["undo_txcorr"],
            cpi=cfg["cpi"],
            ang_tol=cfg["ang_tol"],
            azi_min=cfg["azmin"],
            azi_max=cfg["azmax"],
            ele_min=cfg["elmin"],
            ele_max=cfg["elmax"],
            rng_min=cfg["rmin"],
            rng_max=cfg["rmax"],
        )

        if radar_aux is None:
            continue

        if radar is None:
            radar = radar_aux
        else:
            radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    return radar


def merge_scans_psr_spectra(
    basepath,
    basepath_psr,
    scan_list,
    voltime,
    scan_period,
    datatype_list,
    cfg,
    radarnr="RADAR001",
):
    """
    merge rainbow scans

    Parameters
    ----------
    basepath : str
        base path of rad4alp radar data
    basepath_psr : str
        name of the base path where to find the PSR data
    scan_list : list
        list of scans
    voltime: datetime object
        reference time of the scan
    scan_period : float
        time from reference time where to look for other scans data
    datatype_list : list
        lists of data types to get
    cfg : dict
        configuration dictionary
    radarnr : str
        radar identifier number

    Returns
    -------
    radar : Radar
        radar object

    """
    datadescriptor = radarnr + ":RAINBOW:dBZ"
    endtime = voltime + datetime.timedelta(minutes=scan_period)

    psr = None
    for scan in scan_list:
        filelist = get_file_list(datadescriptor, [voltime], [endtime], cfg, scan=scan)

        if not filelist:
            warn(
                "ERROR: No data file found for scan '%s' "
                "between %s and %s" % (scan, voltime, endtime)
            )
            continue
        scantime = get_datetime(filelist[0], datadescriptor)

        psr_aux = merge_fields_psr_spectra(
            basepath,
            basepath_psr,
            scan,
            scantime,
            datatype_list,
            undo_txcorr=cfg["undo_txcorr"],
            fold=cfg["fold"],
            positive_away=cfg["positive_away"],
            cpi=cfg["cpi"],
            ang_tol=cfg["ang_tol"],
            rng_min=cfg["rmin"],
            rng_max=cfg["rmax"],
            ele_min=cfg["elmin"],
            ele_max=cfg["elmax"],
            azi_min=cfg["azmin"],
            azi_max=cfg["azmax"],
        )

        if psr_aux is None:
            continue

        if psr is None:
            psr = psr_aux
        else:
            psr = pyart.util.radar_utils.join_spectra(psr, psr_aux)

    return psr


def merge_scans_dem(
    basepath,
    scan_list,
    datatype_list,
    rng_min=None,
    rng_max=None,
    azi_min=None,
    azi_max=None,
    ele_min=None,
    ele_max=None,
):
    """
    merge rainbow scans

    Parameters
    ----------
    basepath : str
        base path of rad4alp radar data
    scan_list : list
        list of scans
    datatype_list : list
        lists of data types to get
    radarnr : str
        radar identifier number
    rng_min, rng_max : float
        The range limits [m]. If None the entire coverage of the radar is
        going to be used
    ele_min, ele_max, azi_min, azi_max : float or None
        The limits of the grid [deg]. If None the limits will be the limits
        of the radar volume

    Returns
    -------
    radar : Radar
        radar object

    """
    radar = None
    for scan in scan_list:
        radar_aux = merge_fields_dem(basepath, scan, datatype_list)
        if radar_aux is None:
            continue

        if radar is None:
            radar = radar_aux
            continue

        radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    if radar is None:
        return radar

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rng_min,
        rng_max=rng_max,
        ele_min=ele_min,
        ele_max=ele_max,
        azi_min=azi_min,
        azi_max=azi_max,
    )


def merge_scans_rad4alp(
    basepath, scan_list, radar_name, radar_res, voltime, datatype_list, cfg, ind_rad=0
):
    """
    merge rad4alp data.

    Parameters
    ----------
    basepath : str
        base path of rad4alp radar data
    scan_list : list
        list of scans (001 to 020)
    radar_name : str
        radar_name (A, D, L, ...)
    radar_res : str
        radar resolution (H or L)
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    if (radar_name is None) or (radar_res is None):
        raise ValueError(
            "ERROR: Radar Name and Resolution not specified in config file."
            + " Unable to load rad4alp data"
        )

    timeinfo = voltime.strftime("%H%M")

    radar = None
    for scan in scan_list:
        datapath, basename = get_rad4alp_dir(
            basepath,
            voltime,
            radar_name=radar_name,
            radar_res=radar_res,
            scan=scan,
            path_convention=cfg["path_convention"][ind_rad],
        )

        filename = glob.glob(datapath + basename + timeinfo + "*." + scan + "*")
        if not filename:
            warn("No file found in " + datapath + basename + timeinfo + "*." + scan)
            continue

        radar_aux = get_data_rad4alp(
            filename[0], datatype_list, scan, cfg, ind_rad=ind_rad
        )

        if radar_aux is None:
            continue

        if radar is None:
            radar = radar_aux
        else:
            radar = pyart.util.radar_utils.join_radar(radar, radar_aux)
    if radar is None:
        return radar

    rmin = None
    rmax = None
    elmin = None
    elmax = None
    azmin = None
    azmax = None
    if cfg["rmin"] is not None:
        rmin = cfg["rmin"][ind_rad]
    if cfg["rmax"] is not None:
        rmax = cfg["rmax"][ind_rad]
    if cfg["elmin"] is not None:
        elmin = cfg["elmin"][ind_rad]
    if cfg["elmax"] is not None:
        elmax = cfg["elmax"][ind_rad]
    if cfg["azmin"] is not None:
        azmin = cfg["azmin"][ind_rad]
    if cfg["azmax"] is not None:
        azmax = cfg["azmax"][ind_rad]

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rmin,
        rng_max=rmax,
        ele_min=elmin,
        ele_max=elmax,
        azi_min=azmin,
        azi_max=azmax,
    )


def merge_scans_odim(
    basepath,
    scan_list,
    radar_name,
    radar_res,
    voltime,
    datatype_list,
    dataset_list,
    cfg,
    ind_rad=0,
):
    """
    merge odim data.

    Parameters
    ----------
    basepath : str
        base path of odim radar data
    scan_list : list
        list of scans
    radar_name : str
        radar_name (A, D, L, ...)
    radar_res : str
        radar resolution (H or L)
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    dataset_list : list
        list of datasets. Used to get path
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    # find files to merge
    if "s3BucketRead" in cfg:
        if not _BOTO3_AVAILABLE:
            warn("boto3 not installed")
            return None
        fname_list, scan_list_aux = get_scan_files_to_merge_s3(
            cfg["s3PathRead"],
            scan_list,
            radar_name,
            radar_res,
            voltime,
            dataset_list,
            cfg,
            path_convention=cfg["path_convention"][ind_rad],
            master_scan_time_tol=cfg["MasterScanTimeTol"][ind_rad],
            scan_period=cfg["ScanPeriod"],
        )
    else:
        fname_list, scan_list_aux = get_scan_files_to_merge(
            basepath,
            scan_list,
            radar_name,
            radar_res,
            voltime,
            dataset_list,
            path_convention=cfg["path_convention"][ind_rad],
            master_scan_time_tol=cfg["MasterScanTimeTol"][ind_rad],
            scan_period=cfg["ScanPeriod"],
        )

    if not fname_list:
        warn("No files found")
        return None

    radar = None
    if "s3BucketRead" in cfg:
        s3_client = _open_s3_client(cfg)
        daydir = voltime.strftime("%Y-%m-%d")
        datapath = f"{basepath}{daydir}/"
        if not os.path.isdir(datapath):
            os.makedirs(datapath)

    for fname, scan in zip(fname_list, scan_list_aux):
        if "s3BucketRead" in cfg:
            fname_aux = f"{datapath}{os.path.basename(fname)}"
            s3_client.download_file(cfg["s3BucketRead"], fname, fname_aux)
        else:
            fname_aux = fname

        radar_aux = get_data_odim(fname_aux, datatype_list, scan, cfg, ind_rad=ind_rad)

        if radar_aux is None:
            continue
        if radar is None:
            radar = radar_aux
        else:
            radar = merge_radars(radar, radar_aux)

        if "s3BucketRead" in cfg and "rm_s3_file" in cfg and cfg["rm_s3_file"]:
            os.remove(fname_aux)

    if radar is None:
        return radar

    rmin = None
    rmax = None
    elmin = None
    elmax = None
    azmin = None
    azmax = None
    if cfg["rmin"] is not None:
        rmin = cfg["rmin"][ind_rad]
    if cfg["rmax"] is not None:
        rmax = cfg["rmax"][ind_rad]
    if cfg["elmin"] is not None:
        elmin = cfg["elmin"][ind_rad]
    if cfg["elmax"] is not None:
        elmax = cfg["elmax"][ind_rad]
    if cfg["azmin"] is not None:
        azmin = cfg["azmin"][ind_rad]
    if cfg["azmax"] is not None:
        azmax = cfg["azmax"][ind_rad]
    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rmin,
        rng_max=rmax,
        ele_min=elmin,
        ele_max=elmax,
        azi_min=azmin,
        azi_max=azmax,
    )


def merge_scans_odimgrid(
    basepath, scan_list, voltime, datatype_list, dataset_list, cfg, ind_rad=0
):
    """
    merge odim grid data.

    Parameters
    ----------
    basepath : str
        base path of odim radar data
    scan_list : list
        list of scans
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    dataset_list : list
        list of datasets. Used to get path
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    if cfg["path_convention"][ind_rad] != "ODIM":
        raise ValueError(
            "ERROR: required path convention ODIM for files of type ODIMGRID"
        )

    # find files to merge
    if "s3BucketRead" in cfg:
        if not _BOTO3_AVAILABLE:
            warn("boto3 not installed")
            return None
        fname_list, scan_list_aux = get_scan_files_to_merge_s3(
            cfg["s3PathRead"],
            scan_list,
            None,
            None,
            voltime,
            dataset_list,
            cfg,
            path_convention=cfg["path_convention"][ind_rad],
        )
    else:
        fname_list, scan_list_aux = get_scan_files_to_merge(
            basepath,
            scan_list,
            None,
            None,
            voltime,
            dataset_list,
            path_convention=cfg["path_convention"][ind_rad],
        )

    if not fname_list:
        warn("No files found")
        return None

    grid = None
    if "s3BucketRead" in cfg:
        s3_client = _open_s3_client(cfg)
        daydir = voltime.strftime("%Y-%m-%d")
        datapath = f"{basepath}{daydir}/"
        if not os.path.isdir(datapath):
            os.makedirs(datapath)

    for fname in fname_list:
        if "s3BucketRead" in cfg:
            fname_aux = f"{datapath}{os.path.basename(fname)}"
            s3_client.download_file(cfg["s3BucketRead"], fname, fname_aux)
        else:
            fname_aux = fname

        grid_aux = get_data_odimgrid(
            fname_aux, datatype_list, cfg, mf_scale=cfg["MFScale"], ind_rad=ind_rad
        )

        if grid_aux is None:
            continue
        if grid is None:
            grid = grid_aux
        else:
            for prod_field in grid_aux.fields.keys():
                grid.add_field(prod_field, grid_aux.fields[prod_field])

        if "s3BucketRead" in cfg and "rm_s3_file" in cfg and cfg["rm_s3_file"]:
            os.remove(fname_aux)

    if grid is None:
        return grid

    # Crop the data
    lat_min = cfg.get("latmin", None)
    lat_max = cfg.get("latmax", None)
    lon_min = cfg.get("lonmin", None)
    lon_max = cfg.get("lonmax", None)
    alt_min = cfg.get("altmin", None)
    alt_max = cfg.get("altmax", None)
    ix_min = cfg.get("ixmin", None)
    iy_min = cfg.get("iymin", None)
    iz_min = cfg.get("izmin", None)
    nx = cfg.get("nx", None)
    ny = cfg.get("ny", None)
    nz = cfg.get("nz", None)

    return crop_grid(
        grid,
        lat_min=lat_min,
        lat_max=lat_max,
        lon_min=lon_min,
        lon_max=lon_max,
        alt_min=alt_min,
        alt_max=alt_max,
        nx=nx,
        ny=ny,
        nz=nz,
        ix_min=ix_min,
        iy_min=iy_min,
        iz_min=iz_min,
    )


def merge_scans_knmih5_grid(
    basepath, scan_list, voltime, datatype_list, dataset_list, cfg, ind_rad=0
):
    """
    merge KNMI H5 grid data.

    Parameters
    ----------
    basepath : str
        base path of odim radar data
    scan_list : list
        list of scans
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    dataset_list : list
        list of datasets. Used to get path
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    if cfg["path_convention"][ind_rad] != "ODIM":
        raise ValueError(
            "ERROR: required path convention ODIM for files of type KNMIH5GRID"
        )

    odim_field_names = {}
    for datatype in datatype_list:
        odim_field_names.update(get_datatype_knmi(datatype))

    # find files to merge
    if "s3BucketRead" in cfg:
        if not _BOTO3_AVAILABLE:
            warn("boto3 not installed")
            return None
        fname_list, scan_list_aux = get_scan_files_to_merge_s3(
            cfg["s3PathRead"],
            scan_list,
            None,
            None,
            voltime,
            dataset_list,
            cfg,
            path_convention=cfg["path_convention"][ind_rad],
        )
    else:
        fname_list, scan_list_aux = get_scan_files_to_merge(
            basepath,
            scan_list,
            None,
            None,
            voltime,
            dataset_list,
            path_convention=cfg["path_convention"][ind_rad],
        )

    if not fname_list:
        warn("No files found")
        return None

    grid = None
    if "s3BucketRead" in cfg:
        s3_client = _open_s3_client(cfg)
        daydir = voltime.strftime("%Y-%m-%d")
        datapath = f"{basepath}{daydir}/"
        if not os.path.isdir(datapath):
            os.makedirs(datapath)

    for fname in fname_list:
        if "s3BucketRead" in cfg:
            fname_aux = f"{datapath}{os.path.basename(fname)}"
            s3_client.download_file(cfg["s3BucketRead"], fname, fname_aux)
        else:
            fname_aux = fname

        grid_aux = pyart.aux_io.read_knmi_grid_h5(
            fname_aux, field_names=odim_field_names
        )

        if grid_aux is None:
            continue
        if grid is None:
            grid = grid_aux
        else:
            for prod_field in grid_aux.fields.keys():
                grid.add_field(prod_field, grid_aux.fields[prod_field])

        if "s3BucketRead" in cfg and "rm_s3_file" in cfg and cfg["rm_s3_file"]:
            os.remove(fname_aux)

    if grid is None:
        return grid

    # Crop the data
    lat_min = cfg.get("latmin", None)
    lat_max = cfg.get("latmax", None)
    lon_min = cfg.get("lonmin", None)
    lon_max = cfg.get("lonmax", None)
    alt_min = cfg.get("altmin", None)
    alt_max = cfg.get("altmax", None)
    ix_min = cfg.get("ixmin", None)
    iy_min = cfg.get("iymin", None)
    iz_min = cfg.get("izmin", None)
    nx = cfg.get("nx", None)
    ny = cfg.get("ny", None)
    nz = cfg.get("nz", None)

    return crop_grid(
        grid,
        lat_min=lat_min,
        lat_max=lat_max,
        lon_min=lon_min,
        lon_max=lon_max,
        alt_min=alt_min,
        alt_max=alt_max,
        nx=nx,
        ny=ny,
        nz=nz,
        ix_min=ix_min,
        iy_min=iy_min,
        iz_min=iz_min,
    )


def merge_scans_odimbirds(
    basepath,
    scan_list,
    radar_name,
    radar_res,
    voltime,
    datatype_list,
    dataset_list,
    cfg,
    ind_rad=0,
):
    """
    merge vol2bird odim data.

    Parameters
    ----------
    basepath : str
        base path of odim radar data
    scan_list : list
        list of scans
    radar_name : str
        radar_name (A, D, L, ...)
    radar_res : str
        radar resolution (H or L)
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    dataset_list : list
        list of datasets. Used to get path
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    field_names_dict = {}
    for datatype in datatype_list:
        field_names_dict.update({datatype: get_fieldname_pyart(datatype)})

    # find files to merge
    if "s3BucketRead" in cfg:
        if not _BOTO3_AVAILABLE:
            warn("boto3 not installed")
            return None
        fname_list, scan_list_aux = get_scan_files_to_merge_s3(
            cfg["s3PathRead"],
            scan_list,
            radar_name,
            radar_res,
            voltime,
            dataset_list,
            cfg,
            path_convention=cfg["path_convention"][ind_rad],
            master_scan_time_tol=cfg["MasterScanTimeTol"][ind_rad],
            scan_period=cfg["ScanPeriod"],
        )
    else:
        fname_list, scan_list_aux = get_scan_files_to_merge(
            basepath,
            scan_list,
            radar_name,
            radar_res,
            voltime,
            dataset_list,
            path_convention=cfg["path_convention"][ind_rad],
            master_scan_time_tol=cfg["MasterScanTimeTol"][ind_rad],
            scan_period=cfg["ScanPeriod"],
        )

    if not fname_list:
        warn("No files found")
        return None

    radar = None
    if "s3BucketRead" in cfg:
        s3_client = _open_s3_client(cfg)
        daydir = voltime.strftime("%Y-%m-%d")
        datapath = f"{basepath}{daydir}/"
        if not os.path.isdir(datapath):
            os.makedirs(datapath)

    for fname in fname_list:
        if "s3BucketRead" in cfg:
            fname_aux = f"{datapath}{os.path.basename(fname)}"
            s3_client.download_file(cfg["s3BucketRead"], fname, fname_aux)
        else:
            fname_aux = fname

        radar_aux = pyart.aux_io.read_odim_vp_h5(
            fname_aux, field_names=field_names_dict
        )

        if radar_aux is None:
            continue
        if radar is None:
            radar = radar_aux
        else:
            radar = merge_radars(radar, radar_aux)

        if "s3BucketRead" in cfg and "rm_s3_file" in cfg and cfg["rm_s3_file"]:
            os.remove(fname_aux)

    if radar is None:
        return radar

    rmin = None
    rmax = None
    elmin = None
    elmax = None
    azmin = None
    azmax = None
    if cfg["rmin"] is not None:
        rmin = cfg["rmin"][ind_rad]
    if cfg["rmax"] is not None:
        rmax = cfg["rmax"][ind_rad]
    if cfg["elmin"] is not None:
        elmin = cfg["elmin"][ind_rad]
    if cfg["elmax"] is not None:
        elmax = cfg["elmax"][ind_rad]
    if cfg["azmin"] is not None:
        azmin = cfg["azmin"][ind_rad]
    if cfg["azmax"] is not None:
        azmax = cfg["azmax"][ind_rad]

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rmin,
        rng_max=rmax,
        ele_min=elmin,
        ele_max=elmax,
        azi_min=azmin,
        azi_max=azmax,
    )


def merge_scans_gamic(
    basepath,
    scan_list,
    radar_name,
    radar_res,
    voltime,
    datatype_list,
    dataset_list,
    cfg,
    ind_rad=0,
):
    """
    merge GAMIC data.

    Parameters
    ----------
    basepath : str
        base path of gamic radar data
    scan_list : list
        list of scans
    radar_name : str
        radar_name (A, D, L, ...)
    radar_res : str
        radar resolution (H or L)
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    dataset_list : list
        list of datasets. Used to get path
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    # find files to merge
    if "s3BucketRead" in cfg:
        if not _BOTO3_AVAILABLE:
            warn("boto3 not installed")
            return None
        fname_list, _ = get_scan_files_to_merge_s3(
            cfg["s3PathRead"],
            scan_list,
            radar_name,
            radar_res,
            voltime,
            dataset_list,
            cfg,
            path_convention=cfg["path_convention"][ind_rad],
            master_scan_time_tol=cfg["MasterScanTimeTol"][ind_rad],
            scan_period=cfg["ScanPeriod"],
        )
    else:
        fname_list, _ = get_scan_files_to_merge(
            basepath,
            scan_list,
            radar_name,
            radar_res,
            voltime,
            dataset_list,
            path_convention=cfg["path_convention"][ind_rad],
            master_scan_time_tol=cfg["MasterScanTimeTol"][ind_rad],
            scan_period=cfg["ScanPeriod"],
        )

    if not fname_list:
        warn("No files found")
        return None

    radar = None
    if "s3BucketRead" in cfg:
        s3_client = _open_s3_client(cfg)
        daydir = voltime.strftime("%Y-%m-%d")
        datapath = f"{basepath}{daydir}/"
        if not os.path.isdir(datapath):
            os.makedirs(datapath)

    for fname in fname_list:
        if "s3BucketRead" in cfg:
            fname_aux = f"{datapath}{os.path.basename(fname)}"
            s3_client.download_file(cfg["s3BucketRead"], fname, fname_aux)
        else:
            fname_aux = fname

        radar_aux = get_data_gamic(fname_aux, datatype_list, cfg["pulse_width_gamic"])

        if radar_aux is None:
            continue
        if radar is None:
            radar = radar_aux
        else:
            radar = merge_radars(radar, radar_aux)

        if "s3BucketRead" in cfg and "rm_s3_file" in cfg and cfg["rm_s3_file"]:
            os.remove(fname_aux)

    if radar is None:
        return radar


def merge_scans_mfcfradial(
    basepath, scan_list, voltime, datatype_list, dataset_list, cfg, ind_rad=0
):
    """
    merge CF radial data where one moment is stored per file (mf-moment files)

    Parameters
    ----------
    basepath : str
        base path of odim radar data
    scan_list : list
        list of scans
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    dataset_list : list
        list of datasets. Used to get path
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    # Mapping of MeteoFrance and JMA field names to Py-ART field names
    DataTypeIDInFiles_defaults = {
        "dBZ": "DBZH",  # JMA
        "dBuZ": "TH",  # MF
        "ZDR": "ZDR",  # MF & JMA
        "RhoHV": "RHOHV",  # MF & JMA
        "PhiDP": "PHIDP",  # MF
        "uPhiDP": "PSIDP",  # JMA
        "KDP": "KDP",  # JMA
        "V": "VRAD",  # MF
        "W": "WIDTH",  # JMA
        "sigma_zh": "SIGMA",  # MF
    }

    field_names = {}
    # Use custom name mapping
    for datatype in datatype_list:
        if datatype in cfg["DataTypeIDInFiles"][ind_rad]:
            field_names[
                cfg["DataTypeIDInFiles"][ind_rad][datatype]
            ] = get_fieldname_pyart(datatype)
        else:
            field_names[datatype] = DataTypeIDInFiles_defaults[datatype]

    radar = None

    # get list of files to combine
    flist = []
    scan_list_aux = []

    fpath_strf = dataset_list[0][
        dataset_list[0].find("D") + 2 : dataset_list[0].find("F") - 2
    ]
    fdate_strf = dataset_list[0][dataset_list[0].find("F") + 2 : -1]
    datapath = basepath + voltime.strftime(fpath_strf) + "/"
    filenames = glob.glob(datapath + "*" + scan_list[0] + "*")
    filename = []
    for filename_aux in filenames:
        fdatetime = find_date_in_file_name(filename_aux, date_format=fdate_strf)
        if fdatetime == voltime:
            filename = [filename_aux]
            break

    if not filename:
        warn(
            f"No file found in {datapath} for scan {scan_list[0]}"
            f" and time {voltime}"
        )
    else:
        flist.append(filename[0])
        scan_list_aux.append(scan_list[0])

    if len(scan_list) > 1:
        # merge the elevations into a single radar instance
        for scan in scan_list[1:]:
            filenames = glob.glob(f"{datapath}*{scan}*")
            filename = []
            for filename_aux in filenames:
                fdatetime = find_date_in_file_name(filename_aux, date_format=fdate_strf)
                if cfg["MasterScanTimeTol"][ind_rad] == 0:
                    if fdatetime == voltime:
                        filename = [filename_aux]
                        break
                elif cfg["MasterScanTimeTol"][ind_rad] == 1:
                    if (
                        voltime
                        <= fdatetime
                        < voltime + datetime.timedelta(minutes=cfg["ScanPeriod"])
                    ):
                        filename = [filename_aux]
                        print(os.path.basename(filename[0]))
                        break
                else:
                    if (
                        voltime - datetime.timedelta(minutes=cfg["ScanPeriod"])
                        < fdatetime
                        <= voltime
                    ):
                        filename = [filename_aux]
                        print(os.path.basename(filename[0]))
                        break
            if not filename:
                warn(
                    f"No file found in {datapath} for scan {scan}"
                    f" and time {voltime}"
                )
                continue
            flist.append(filename[0])
            scan_list_aux.append(scan)

    if not flist:
        return radar

    if cfg["DataTypeIDInFiles"] is None:
        for fname, scan in zip(flist, scan_list_aux):
            radar_aux = pyart.io.read_cfradial(fname, field_names=field_names)
            if radar_aux is None:
                continue
            if radar is None:
                radar = radar_aux
                continue
            radar = pyart.util.radar_utils.join_radar(radar, radar_aux)
    else:
        for datatype in datatype_list:
            if datatype not in cfg["DataTypeIDInFiles"].keys():
                warn(f"No file contains data type {datatype}")
                continue
            nscans = 0
            radar_aux = None
            for fname, scan in zip(flist, scan_list_aux):
                if cfg["DataTypeIDInFiles"][datatype] not in os.path.basename(fname):
                    continue
                radar_aux2 = pyart.io.read_cfradial(fname, field_names=field_names)
                if radar_aux2 is None:
                    continue
                if radar_aux is None:
                    radar_aux = radar_aux2
                    nscans += 1
                    continue
                radar_aux = pyart.util.radar_utils.join_radar(radar_aux, radar_aux2)
                nscans += 1
            if radar is None:
                radar = radar_aux
                nscans_expected = nscans
                continue
            if nscans != nscans_expected:
                warn(
                    f"Unable to merge datatype {datatype} into radar object."
                    f" Number of scans containing the datatype {nscans}"
                    f" different from number of scans expected"
                    f" {nscans_expected}"
                )
                continue
            field_name = get_fieldname_pyart(datatype)
            try:
                radar.add_field(field_name, radar_aux.fields[field_name])
            except ValueError:
                if radar.nrays * radar.ngates > radar_aux.nrays * radar_aux.ngates:
                    warn(f"Field {field_name} will be interpolated")
                    radar = add_field(radar, radar_aux)
                else:
                    warn(f"Fields will be adapted to {field_name} field size")
                    radar = add_field(radar_aux, radar)
                print(f"nrays: {radar.nrays} ngates: {radar.ngates}")

    if radar is None:
        return radar

    rmin = None
    rmax = None
    elmin = None
    elmax = None
    azmin = None
    azmax = None
    if cfg["rmin"] is not None:
        rmin = cfg["rmin"][ind_rad]
    if cfg["rmax"] is not None:
        rmax = cfg["rmax"][ind_rad]
    if cfg["elmin"] is not None:
        elmin = cfg["elmin"][ind_rad]
    if cfg["elmax"] is not None:
        elmax = cfg["elmax"][ind_rad]
    if cfg["azmin"] is not None:
        azmin = cfg["azmin"][ind_rad]
    if cfg["azmax"] is not None:
        azmax = cfg["azmax"][ind_rad]

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rmin,
        rng_max=rmax,
        ele_min=elmin,
        ele_max=elmax,
        azi_min=azmin,
        azi_max=azmax,
    )


def merge_scans_nexrad2(
    basepath,
    scan_list,
    voltime,
    datatype_list,
    dataset_list,
    cfg,
    ind_rad=0,
):
    """
    merge NEXRAD level 2 data.

    Parameters
    ----------
    basepath : str
        base path of nexrad radar data
    scan_list : list
        list of scans
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    dataset_list : list
        list of datasets. Used to get path
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    # find files to merge
    if "s3BucketRead" in cfg:
        if not _BOTO3_AVAILABLE:
            warn("boto3 not installed")
            return None
        fname_list, _ = get_scan_files_to_merge_s3(
            cfg["s3PathRead"],
            scan_list,
            None,
            None,
            voltime,
            dataset_list,
            cfg,
            path_convention=cfg["path_convention"][ind_rad],
            master_scan_time_tol=cfg["MasterScanTimeTol"][ind_rad],
            scan_period=cfg["ScanPeriod"],
        )
    else:
        fname_list, _ = get_scan_files_to_merge(
            basepath,
            scan_list,
            None,
            None,
            voltime,
            dataset_list,
            path_convention=cfg["path_convention"][ind_rad],
            master_scan_time_tol=cfg["MasterScanTimeTol"][ind_rad],
            scan_period=cfg["ScanPeriod"],
        )

    if not fname_list:
        warn("No files found")
        return None

    radar = None
    if "s3BucketRead" in cfg:
        s3_client = _open_s3_client(cfg)
        daydir = voltime.strftime("%Y-%m-%d")
        datapath = f"{basepath}{daydir}/"
        if not os.path.isdir(datapath):
            os.makedirs(datapath)

    for fname in fname_list:
        if "s3BucketRead" in cfg:
            fname_aux = f"{datapath}{os.path.basename(fname)}"
            s3_client.download_file(cfg["s3BucketRead"], fname, fname_aux)
        else:
            fname_aux = fname

        radar_aux = pyart.io.read(fname)
        if radar_aux is None:
            continue
        if radar is None:
            radar = radar_aux
        else:
            radar = merge_radars(radar, radar_aux)

    if radar is None:
        return radar

    rmin = None
    rmax = None
    elmin = None
    elmax = None
    azmin = None
    azmax = None
    if cfg["rmin"] is not None:
        rmin = cfg["rmin"][ind_rad]
    if cfg["rmax"] is not None:
        rmax = cfg["rmax"][ind_rad]
    if cfg["elmin"] is not None:
        elmin = cfg["elmin"][ind_rad]
    if cfg["elmax"] is not None:
        elmax = cfg["elmax"][ind_rad]
    if cfg["azmin"] is not None:
        azmin = cfg["azmin"][ind_rad]
    if cfg["azmax"] is not None:
        azmax = cfg["azmax"][ind_rad]

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rmin,
        rng_max=rmax,
        ele_min=elmin,
        ele_max=elmax,
        azi_min=azmin,
        azi_max=azmax,
    )


def merge_scans_cfradial(
    basepath,
    scan_list,
    radar_name,
    radar_res,
    voltime,
    datatype_list,
    dataset_list,
    cfg,
    ind_rad=0,
):
    """
    merge CFRADIAL data.

    Parameters
    ----------
    basepath : str
        base path of CFRADIAL radar data
    scan_list : list
        list of scans
    radar_name : str
        radar name
    radar_res : str
        radar resolution type
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    dataset_list : list
        list of datasets. Used to get path
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    field_names = {}
    # Use custom name mapping
    for datatype in datatype_list:
        if datatype in cfg["DataTypeIDInFiles"][ind_rad]:
            field_names[
                cfg["DataTypeIDInFiles"][ind_rad][datatype]
            ] = get_fieldname_pyart(datatype)
        else:
            field_names[datatype] = get_fieldname_pyart(datatype)

    # find files to merge
    if "s3BucketRead" in cfg:
        if not _BOTO3_AVAILABLE:
            warn("boto3 not installed")
            return None
        fname_list, _ = get_scan_files_to_merge_s3(
            cfg["s3PathRead"],
            scan_list,
            radar_name,
            radar_res,
            voltime,
            dataset_list,
            cfg,
            path_convention=cfg["path_convention"][ind_rad],
            master_scan_time_tol=cfg["MasterScanTimeTol"][ind_rad],
            scan_period=cfg["ScanPeriod"],
        )
    else:
        fname_list, _ = get_scan_files_to_merge(
            basepath,
            scan_list,
            radar_name,
            radar_res,
            voltime,
            dataset_list,
            path_convention=cfg["path_convention"][ind_rad],
            master_scan_time_tol=cfg["MasterScanTimeTol"][ind_rad],
            scan_period=cfg["ScanPeriod"],
        )

    if not fname_list:
        warn("No files found")
        return None

    radar = None
    if "s3BucketRead" in cfg:
        s3_client = _open_s3_client(cfg)
        daydir = voltime.strftime("%Y-%m-%d")
        datapath = f"{basepath}{daydir}/"
        if not os.path.isdir(datapath):
            os.makedirs(datapath)

    for fname in fname_list:
        if "s3BucketRead" in cfg:
            fname_aux = f"{datapath}{os.path.basename(fname)}"
            s3_client.download_file(cfg["s3BucketRead"], fname, fname_aux)
        else:
            fname_aux = fname

        radar_aux = pyart.io.read_cfradial(
            fname_aux, field_names=field_names, include_fields=field_names.values()
        )

        if radar_aux is None:
            continue
        if radar is None:
            radar = radar_aux
        else:
            radar = merge_radars(radar, radar_aux)

        if "s3BucketRead" in cfg and "rm_s3_file" in cfg and cfg["rm_s3_file"]:
            os.remove(fname_aux)

    if radar is None:
        return radar

    rmin = None
    rmax = None
    elmin = None
    elmax = None
    azmin = None
    azmax = None
    if cfg["rmin"] is not None:
        rmin = cfg["rmin"][ind_rad]
    if cfg["rmax"] is not None:
        rmax = cfg["rmax"][ind_rad]
    if cfg["elmin"] is not None:
        elmin = cfg["elmin"][ind_rad]
    if cfg["elmax"] is not None:
        elmax = cfg["elmax"][ind_rad]
    if cfg["azmin"] is not None:
        azmin = cfg["azmin"][ind_rad]
    if cfg["azmax"] is not None:
        azmax = cfg["azmax"][ind_rad]

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rmin,
        rng_max=rmax,
        ele_min=elmin,
        ele_max=elmax,
        azi_min=azmin,
        azi_max=azmax,
    )


def merge_scans_cfradial2(
    basepath,
    scan_list,
    radar_name,
    radar_res,
    voltime,
    datatype_list,
    dataset_list,
    cfg,
    ind_rad=0,
):
    """
    merge CFRADIAL2 data.

    Parameters
    ----------
    basepath : str
        base path of CFRADIAL2 radar data
    scan_list : list
        list of scans
    radar_name : str
        radar name
    radar_res : str
        radar resolution type
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    dataset_list : list
        list of datasets. Used to get path
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    field_names = {}
    # Use custom name mapping
    for datatype in datatype_list:
        if datatype in cfg["DataTypeIDInFiles"][ind_rad]:
            field_names[
                cfg["DataTypeIDInFiles"][ind_rad][datatype]
            ] = get_fieldname_pyart(datatype)
        else:
            field_names[datatype] = get_fieldname_pyart(datatype)

    # find files to merge
    if "s3BucketRead" in cfg:
        if not _BOTO3_AVAILABLE:
            warn("boto3 not installed")
            return None
        fname_list, _ = get_scan_files_to_merge_s3(
            cfg["s3PathRead"],
            scan_list,
            radar_name,
            radar_res,
            voltime,
            dataset_list,
            cfg,
            path_convention=cfg["path_convention"][ind_rad],
            master_scan_time_tol=cfg["MasterScanTimeTol"][ind_rad],
            scan_period=cfg["ScanPeriod"],
        )
    else:
        fname_list, _ = get_scan_files_to_merge(
            basepath,
            scan_list,
            radar_name,
            radar_res,
            voltime,
            dataset_list,
            path_convention=cfg["path_convention"][ind_rad],
            master_scan_time_tol=cfg["MasterScanTimeTol"][ind_rad],
            scan_period=cfg["ScanPeriod"],
        )

    if not fname_list:
        warn("No files found")
        return None

    radar = None
    if "s3BucketRead" in cfg:
        s3_client = _open_s3_client(cfg)
        daydir = voltime.strftime("%Y-%m-%d")
        datapath = f"{basepath}{daydir}/"
        if not os.path.isdir(datapath):
            os.makedirs(datapath)

    for fname in fname_list:
        if "s3BucketRead" in cfg:
            fname_aux = f"{datapath}{os.path.basename(fname)}"
            s3_client.download_file(cfg["s3BucketRead"], fname, fname_aux)
        else:
            fname_aux = fname

        radar_aux = pyart.io.read_cfradial2(
            fname, field_names=field_names, include_fields=field_names.values()
        )
        if radar_aux is None:
            continue
        if radar is None:
            radar = radar_aux
        else:
            radar = merge_radars(radar, radar_aux)

    if radar is None:
        return radar

    rmin = None
    rmax = None
    elmin = None
    elmax = None
    azmin = None
    azmax = None
    if cfg["rmin"] is not None:
        rmin = cfg["rmin"][ind_rad]
    if cfg["rmax"] is not None:
        rmax = cfg["rmax"][ind_rad]
    if cfg["elmin"] is not None:
        elmin = cfg["elmin"][ind_rad]
    if cfg["elmax"] is not None:
        elmax = cfg["elmax"][ind_rad]
    if cfg["azmin"] is not None:
        azmin = cfg["azmin"][ind_rad]
    if cfg["azmax"] is not None:
        azmax = cfg["azmax"][ind_rad]

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rmin,
        rng_max=rmax,
        ele_min=elmin,
        ele_max=elmax,
        azi_min=azmin,
        azi_max=azmax,
    )


def merge_scans_skyecho(
    basepath, scan_list, voltime, datatype_list, dataset_list, cfg, ind_rad=0
):
    """
    merge SKYECHO data.

    Parameters
    ----------
    basepath : str
        base path of CFRADIAL radar data
    scan_list : list
        list of scans
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    dataset_list : list
        list of datasets. Used to get path
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    skyecho_field_names = {}
    for datatype in datatype_list:
        skyecho_field_names.update(get_datatype_skyecho(datatype))

    # find files to merge
    if "s3BucketRead" in cfg:
        if not _BOTO3_AVAILABLE:
            warn("boto3 not installed")
            return None
        fname_list, _ = get_scan_files_to_merge_s3(
            cfg["s3PathRead"],
            scan_list,
            None,
            None,
            voltime,
            dataset_list,
            cfg,
            path_convention="SkyEcho",
        )
    else:
        fname_list, _ = get_scan_files_to_merge(
            basepath,
            scan_list,
            None,
            None,
            voltime,
            dataset_list,
            path_convention="SkyEcho",
        )

    if not fname_list:
        warn("No files found")
        return None

    radar = None
    if "s3BucketRead" in cfg:
        datapath = f"{basepath}"
        if not os.path.isdir(datapath):
            os.makedirs(datapath)

    for fname in fname_list:
        if "s3BucketRead" in cfg:
            fname_aux = f"{datapath}{os.path.basename(fname)}"
        else:
            fname_aux = fname

        radar_aux = pyart.aux_io.read_skyecho(
            fname_aux, sweep_end_time=voltime, field_names=skyecho_field_names
        )

        if radar_aux is None:
            continue
        if radar is None:
            radar = radar_aux
        else:
            radar = merge_radars(radar, radar_aux)

    if radar is None:
        return radar

    rmin = None
    rmax = None
    elmin = None
    elmax = None
    azmin = None
    azmax = None
    if cfg["rmin"] is not None:
        rmin = cfg["rmin"][ind_rad]
    if cfg["rmax"] is not None:
        rmax = cfg["rmax"][ind_rad]
    if cfg["elmin"] is not None:
        elmin = cfg["elmin"][ind_rad]
    if cfg["elmax"] is not None:
        elmax = cfg["elmax"][ind_rad]
    if cfg["azmin"] is not None:
        azmin = cfg["azmin"][ind_rad]
    if cfg["azmax"] is not None:
        azmax = cfg["azmax"][ind_rad]

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rmin,
        rng_max=rmax,
        ele_min=elmin,
        ele_max=elmax,
        azi_min=azmin,
        azi_max=azmax,
    )


def merge_scans_cf1(
    basepath,
    scan_list,
    radar_name,
    radar_res,
    voltime,
    datatype_list,
    dataset_list,
    cfg,
    ind_rad=0,
):
    """
    merge CF1 data.

    Parameters
    ----------
    basepath : str
        base path of CF1 radar data
    scan_list : list
        list of scans
    radar_name : str
        radar name
    radar_res : str
        radar resolution type
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    dataset_list : list
        list of datasets. Used to get path
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    field_names = {}
    # Use custom name mapping
    for datatype in datatype_list:
        if datatype in cfg["DataTypeIDInFiles"][ind_rad]:
            field_names[
                cfg["DataTypeIDInFiles"][ind_rad][datatype]
            ] = get_fieldname_pyart(datatype)
        else:
            field_names[datatype] = get_fieldname_pyart(datatype)

    # find files to merge
    if "s3BucketRead" in cfg:
        if not _BOTO3_AVAILABLE:
            warn("boto3 not installed")
            return None
        fname_list, _ = get_scan_files_to_merge_s3(
            cfg["s3PathRead"],
            scan_list,
            radar_name,
            radar_res,
            voltime,
            dataset_list,
            cfg,
            path_convention=cfg["path_convention"][ind_rad],
            master_scan_time_tol=cfg["MasterScanTimeTol"][ind_rad],
            scan_period=cfg["ScanPeriod"],
        )
    else:
        fname_list, _ = get_scan_files_to_merge(
            basepath,
            scan_list,
            radar_name,
            radar_res,
            voltime,
            dataset_list,
            path_convention=cfg["path_convention"][ind_rad],
            master_scan_time_tol=cfg["MasterScanTimeTol"][ind_rad],
            scan_period=cfg["ScanPeriod"],
        )

    if not fname_list:
        warn("No files found")
        return None

    radar = None
    if "s3BucketRead" in cfg:
        s3_client = _open_s3_client(cfg)
        daydir = voltime.strftime("%Y-%m-%d")
        datapath = f"{basepath}{daydir}/"
        if not os.path.isdir(datapath):
            os.makedirs(datapath)

    for fname in fname_list:
        if "s3BucketRead" in cfg:
            fname_aux = f"{datapath}{os.path.basename(fname)}"
            s3_client.download_file(cfg["s3BucketRead"], fname, fname_aux)
        else:
            fname_aux = fname

        radar_aux = pyart.aux_io.read_cf1(
            fname_aux, field_names=field_names, include_fields=field_names.values()
        )

        if radar_aux is None:
            continue
        if radar is None:
            radar = radar_aux
        else:
            radar = merge_radars(radar, radar_aux)

        if "s3BucketRead" in cfg and "rm_s3_file" in cfg and cfg["rm_s3_file"]:
            os.remove(fname_aux)

    if radar is None:
        return radar

    rmin = None
    rmax = None
    elmin = None
    elmax = None
    azmin = None
    azmax = None
    if cfg["rmin"] is not None:
        rmin = cfg["rmin"][ind_rad]
    if cfg["rmax"] is not None:
        rmax = cfg["rmax"][ind_rad]
    if cfg["elmin"] is not None:
        elmin = cfg["elmin"][ind_rad]
    if cfg["elmax"] is not None:
        elmax = cfg["elmax"][ind_rad]
    if cfg["azmin"] is not None:
        azmin = cfg["azmin"][ind_rad]
    if cfg["azmax"] is not None:
        azmax = cfg["azmax"][ind_rad]

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rmin,
        rng_max=rmax,
        ele_min=elmin,
        ele_max=elmax,
        azi_min=azmin,
        azi_max=azmax,
    )


def merge_scans_mxpol(basepath, scan_list, voltime, datatype_list, cfg, ind_rad=0):
    """
    merge MXPOL data.

    Parameters
    ----------
    basepath : str
        base path of mxpol radar data
    scan_list : list
        list of scans, in the case of mxpol, the elevation or azimuth denoted
        as 005 or 090 (for 5 or 90 degrees elevation) or 330 (for 330 degrees
        azimuth respectively)
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    radar = None
    for scan in scan_list:
        if cfg["path_convention"][ind_rad] == "LTE":
            sub1 = str(voltime.year)
            sub2 = voltime.strftime("%m")
            sub3 = voltime.strftime("%d")
            dayinfo = voltime.strftime("%Y%m%d")
            timeinfo = voltime.strftime("%H%M")
            datapath = basepath + "/" + sub1 + "/" + sub2 + "/" + sub3 + "/"
            scanname = "MXPol-polar-" + dayinfo + "-" + timeinfo + "*-"
            filename = glob.glob(datapath + scanname + scan + "*")
        else:
            daydir = voltime.strftime("%Y-%m-%d")
            dayinfo = voltime.strftime("%Y%m%d")
            timeinfo = voltime.strftime("%H%M")
            datapath = basepath + scan + "/" + daydir + "/"
            if not os.path.isdir(datapath):
                warn("WARNING: Unknown datapath '%s'" % datapath)
                return None
            filename = glob.glob(
                datapath
                + "MXPol-polar-"
                + dayinfo
                + "-"
                + timeinfo
                + "*-"
                + scan
                + ".nc"
            )
        if not filename:
            warn("No file found in " + datapath + scanname + scan)
            continue

        radar_aux = get_data_mxpol(filename[0], datatype_list)

        if radar is None:
            radar = radar_aux
            continue

        radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    if radar is None:
        return radar

    rmin = None
    rmax = None
    elmin = None
    elmax = None
    azmin = None
    azmax = None
    if cfg["rmin"] is not None:
        rmin = cfg["rmin"][ind_rad]
    if cfg["rmax"] is not None:
        rmax = cfg["rmax"][ind_rad]
    if cfg["elmin"] is not None:
        elmin = cfg["elmin"][ind_rad]
    if cfg["elmax"] is not None:
        elmax = cfg["elmax"][ind_rad]
    if cfg["azmin"] is not None:
        azmin = cfg["azmin"][ind_rad]
    if cfg["azmax"] is not None:
        azmax = cfg["azmax"][ind_rad]

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rmin,
        rng_max=rmax,
        ele_min=elmin,
        ele_max=elmax,
        azi_min=azmin,
        azi_max=azmax,
    )


def merge_scans_icon(voltime, datatype_list, cfg, ind_rad=0):
    """
    merge rainbow scans

    Parameters
    ----------
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    radar = None
    for scan in cfg["ScanList"][ind_rad]:
        filename_list = []
        for datatype in datatype_list:
            filename = find_icon_file(voltime, datatype, cfg, scan, ind_rad=ind_rad)
            if filename is not None:
                filename_list.append(filename)

        nfiles_valid = len(filename_list)
        if nfiles_valid > 0:
            radar_aux = merge_fields_icon(filename_list)

            if radar is None:
                radar = radar_aux
            else:
                radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    if radar is None:
        return radar

    rmin = None
    rmax = None
    elmin = None
    elmax = None
    azmin = None
    azmax = None
    if cfg["rmin"] is not None:
        rmin = cfg["rmin"][ind_rad]
    if cfg["rmax"] is not None:
        rmax = cfg["rmax"][ind_rad]
    if cfg["elmin"] is not None:
        elmin = cfg["elmin"][ind_rad]
    if cfg["elmax"] is not None:
        elmax = cfg["elmax"][ind_rad]
    if cfg["azmin"] is not None:
        azmin = cfg["azmin"][ind_rad]
    if cfg["azmax"] is not None:
        azmax = cfg["azmax"][ind_rad]

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rmin,
        rng_max=rmax,
        ele_min=elmin,
        ele_max=elmax,
        azi_min=azmin,
        azi_max=azmax,
    )


def merge_scans_icon_rad4alp(voltime, datatype, cfg, ind_rad=0):
    """
    merge icon rad4alp scans. If data for all the scans cannot be retrieved
    returns None

    Parameters
    ----------
    voltime: datetime object
        reference time of the scan
    datatype : str
        name of the data type to read
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    # look for rad4alp ICON data. Data must be present in all scans
    # to consider the volume valid
    radar = None
    for scan in cfg["ScanList"][ind_rad]:
        # create the radar object where to store the data
        # taking as reference the metranet polar file
        # The radar cutting is going to be done at the end
        cfg_aux = deepcopy(cfg)
        cfg_aux["rmin"] = None
        cfg_aux["rmax"] = None
        cfg_aux["elmin"] = None
        cfg_aux["elmax"] = None
        cfg_aux["azmin"] = None
        cfg_aux["azmax"] = None
        radar_aux = merge_scans_rad4alp(
            cfg["datapath"][ind_rad],
            [scan],
            cfg["RadarName"][ind_rad],
            cfg["RadarRes"][ind_rad],
            voltime,
            ["dBZ"],
            cfg_aux,
            ind_rad=ind_rad,
        )

        if radar_aux is None:
            return None

        # read the icon file
        filename = find_rad4alpicon_file(voltime, datatype, cfg, scan)
        if filename is None:
            return None
        icon_dict = read_rad4alp_icon(filename, datatype)
        radar_aux.add_field(get_fieldname_pyart(datatype), icon_dict)

        if radar is None:
            radar = radar_aux
        else:
            radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    if radar is None:
        return radar

    rmin = None
    rmax = None
    elmin = None
    elmax = None
    azmin = None
    azmax = None
    if cfg["rmin"] is not None:
        rmin = cfg["rmin"][ind_rad]
    if cfg["rmax"] is not None:
        rmax = cfg["rmax"][ind_rad]
    if cfg["elmin"] is not None:
        elmin = cfg["elmin"][ind_rad]
    if cfg["elmax"] is not None:
        elmax = cfg["elmax"][ind_rad]
    if cfg["azmin"] is not None:
        azmin = cfg["azmin"][ind_rad]
    if cfg["azmax"] is not None:
        azmax = cfg["azmax"][ind_rad]

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rmin,
        rng_max=rmax,
        ele_min=elmin,
        ele_max=elmax,
        azi_min=azmin,
        azi_max=azmax,
    )


def merge_scans_dem_rad4alp(voltime, datatype, cfg, ind_rad=0):
    """
    merge DEM rad4alp scans. If data for all the scans cannot be retrieved
    returns None

    Parameters
    ----------
    voltime: datetime object
        reference time of the scan
    datatype : str
        name of the data type to read
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    # read visibility data file
    vis_list = read_rad4alp_vis(
        cfg["dempath"][ind_rad] + cfg["RadarName"][ind_rad] + "_visib_volume_40",
        datatype,
    )
    if vis_list is None:
        return None

    radar = None
    for scan in cfg["ScanList"][ind_rad]:
        # create the radar object where to store the data
        # taking as reference the metranet polar file
        # The radar cutting is going to be done at the end
        cfg_aux = deepcopy(cfg)
        cfg_aux["rmin"] = None
        cfg_aux["rmax"] = None
        cfg_aux["elmin"] = None
        cfg_aux["elmax"] = None
        cfg_aux["azmin"] = None
        cfg_aux["azmax"] = None
        radar_aux = merge_scans_rad4alp(
            cfg["datapath"][ind_rad],
            [scan],
            cfg["RadarName"][ind_rad],
            cfg["RadarRes"][ind_rad],
            voltime,
            ["dBZ"],
            cfg_aux,
            ind_rad=ind_rad,
        )

        if radar_aux is None:
            return None

        # add visibility data
        radar_aux.fields = {}
        radar_aux.add_field(get_fieldname_pyart(datatype), vis_list[int(scan) - 1])

        if radar is None:
            radar = radar_aux
        else:
            radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    if radar is None:
        return radar

    rmin = None
    rmax = None
    elmin = None
    elmax = None
    azmin = None
    azmax = None
    if cfg["rmin"] is not None:
        rmin = cfg["rmin"][ind_rad]
    if cfg["rmax"] is not None:
        rmax = cfg["rmax"][ind_rad]
    if cfg["elmin"] is not None:
        elmin = cfg["elmin"][ind_rad]
    if cfg["elmax"] is not None:
        elmax = cfg["elmax"][ind_rad]
    if cfg["azmin"] is not None:
        azmin = cfg["azmin"][ind_rad]
    if cfg["azmax"] is not None:
        azmax = cfg["azmax"][ind_rad]

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rmin,
        rng_max=rmax,
        ele_min=elmin,
        ele_max=elmax,
        azi_min=azmin,
        azi_max=azmax,
    )


def merge_scans_other_rad4alp(voltime, datatype, cfg, ind_rad=0):
    """
    merge other rad4alp polar products not contained in the basic M or P
    files, i.e. hydro, dealiased velocity or precip. If data for all the
    scans cannot be retrieved returns None

    Parameters
    ----------
    voltime: datetime object
        reference time of the scan
    datatype : str
        name of the data type to read
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    radar_name = cfg["RadarName"][ind_rad]
    radar_res = cfg["RadarRes"][ind_rad]
    basepath = cfg["datapath"][ind_rad]
    scan_list = cfg["ScanList"][ind_rad]
    dayinfo = voltime.strftime("%y%j")
    timeinfo = voltime.strftime("%H%M")

    acronym, _ = get_rad4alp_prod_fname(datatype)
    prod_field = get_fieldname_pyart(datatype)
    prod_dict = pyart.config.get_metadata(prod_field)
    basename_prod = acronym + radar_name + dayinfo

    radar = None
    for scan in scan_list:
        # read product data file
        if cfg["path_convention"][ind_rad] == "LTE":
            yy = dayinfo[0:2]
            dy = dayinfo[2:]
            subf = acronym + radar_name + yy + "hdf" + dy
            datapath_prod = basepath + subf + "/"
        elif cfg["path_convention"][ind_rad] == "MCH":
            datapath_prod = basepath + dayinfo + "/" + basename_prod + "/"
        else:
            datapath_prod = basepath + acronym + radar_name + "/"

        filename_prod = glob.glob(
            datapath_prod + basename_prod + timeinfo + "*." + str(800 + int(scan)) + "*"
        )
        if not filename_prod:
            warn(
                "No file found in "
                + datapath_prod
                + basename_prod
                + timeinfo
                + "*."
                + str(800 + int(scan))
            )
            return None

        filename_prod = filename_prod[0]
        if cfg["metranet_read_lib"] == "C" and _METRANETLIB_AVAILABLE:
            prod_obj = pyart.aux_io.read_product_c(
                filename_prod, physic_value=False, masked_array=True
            )
        elif cfg["metranet_read_lib"] == "python":
            prod_obj = pyart.aux_io.read_product_py(
                filename_prod, physic_value=False, masked_array=True
            )
        else:
            warn(
                "METRANET C-library reader not available or unknown "
                + "library type. Python library will be used"
            )
            prod_obj = pyart.aux_io.read_product_py(
                filename_prod, physic_value=False, masked_array=True
            )
        if prod_obj is None:
            warn("Unable to read file " + filename_prod)
            return None

        if datatype == "hydro":
            prod_dict["data"] = map_hydro(prod_obj.data)
        elif datatype == "dealV":
            prod_dict["data"] = map_Doppler(
                prod_obj.data, float(prod_obj.header["nyquist"])
            )

        # create the radar object where to store the data
        # taking as reference the metranet polar file
        # The radar cutting is going to be done at the end
        cfg_aux = deepcopy(cfg)
        cfg_aux["rmin"] = None
        cfg_aux["rmax"] = None
        cfg_aux["elmin"] = None
        cfg_aux["elmax"] = None
        cfg_aux["azmin"] = None
        cfg_aux["azmax"] = None
        radar_aux = merge_scans_rad4alp(
            basepath,
            [scan],
            radar_name,
            radar_res,
            voltime,
            ["dBZ"],
            cfg_aux,
            ind_rad=ind_rad,
        )

        if radar_aux is None:
            return None

        # add product data
        radar_aux.fields = {}
        radar_aux.add_field(prod_field, prod_dict)

        if radar is None:
            radar = radar_aux
        else:
            radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    if radar is None:
        return radar

    rmin = None
    rmax = None
    elmin = None
    elmax = None
    azmin = None
    azmax = None
    if cfg["rmin"] is not None:
        rmin = cfg["rmin"][ind_rad]
    if cfg["rmax"] is not None:
        rmax = cfg["rmax"][ind_rad]
    if cfg["elmin"] is not None:
        elmin = cfg["elmin"][ind_rad]
    if cfg["elmax"] is not None:
        elmax = cfg["elmax"][ind_rad]
    if cfg["azmin"] is not None:
        azmin = cfg["azmin"][ind_rad]
    if cfg["azmax"] is not None:
        azmax = cfg["azmax"][ind_rad]

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rmin,
        rng_max=rmax,
        ele_min=elmin,
        ele_max=elmax,
        azi_min=azmin,
        azi_max=azmax,
    )


def merge_scans_iq_rad4alp(
    basepath,
    basepath_iq,
    scan_list,
    radar_name,
    radar_res,
    voltime,
    datatype_list,
    cfg,
    ang_tol=0.1,
    ang_step=0.01,
    ind_rad=0,
):
    """
    merge rad4alp IQ scans

    Parameters
    ----------
    basepath : str
        base path of rad4alp radar data
    basepath_iq : str
        base path of rad4alp IQ data
    scan_list : list
        list of scans (001 to 020)
    radar_name : str
        radar_name (A, D, L, ...)
    radar_res : str
        radar resolution (H or L)
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    cfg : dict
        configuration dictionary
    ang_tol : float
        Tolerance between nominal elevation and actual elevation
    ang_step : float
        The elevation angular step used when checking valid ray files
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    if (radar_name is None) or (radar_res is None):
        raise ValueError(
            "ERROR: Radar Name and Resolution not specified in config file."
            + " Unable to load rad4alp data"
        )

    timeinfo = voltime.strftime("%H%M")
    dayinfo = voltime.strftime("%y%j")

    ele_vec = [
        -0.2,
        0.4,
        1.0,
        1.6,
        2.5,
        3.5,
        4.5,
        5.5,
        6.5,
        7.5,
        8.5,
        9.5,
        11.0,
        13.0,
        16.0,
        20.0,
        25.0,
        30.0,
        35.0,
        40.0,
    ]
    prfs = [
        600.0,
        700.0,
        600.0,
        900.0,
        800.0,
        900.0,
        1000.0,
        900.0,
        1000.0,
        1200.0,
        1200.0,
        1200.0,
        1500.0,
        1500.0,
        1500.0,
        1500.0,
        1500.0,
        1500.0,
        1500.0,
        1500.0,
    ]

    field_names = []
    for datatype in datatype_list:
        field_names.append(get_fieldname_pyart(datatype))

    # read status
    root = read_status(voltime, cfg, ind_rad=ind_rad)

    radconst_h = None
    radconst_v = None
    if "radconsth" in cfg:
        radconst_h = cfg["radconsth"][ind_rad]
    if "radconstv" in cfg:
        radconst_v = cfg["radconstv"][ind_rad]

    mfloss_h = None
    mfloss_v = None
    if "mflossh" in cfg:
        mfloss_h = cfg["mflossh"][ind_rad]
    if "mflossv" in cfg:
        mfloss_v = cfg["mflossv"][ind_rad]

    radar = None
    for scan in scan_list:
        datapath, basename = get_rad4alp_dir(
            basepath,
            voltime,
            radar_name=radar_name,
            radar_res=radar_res,
            scan=scan,
            path_convention=cfg["path_convention"][ind_rad],
        )

        filename = glob.glob(datapath + basename + timeinfo + "*." + scan + "*")
        if not filename:
            warn("No file found in " + datapath + basename + timeinfo + "*." + scan)
            continue

        ele = ele_vec[int(scan) - 1]
        datapath_iq = basepath_iq + dayinfo + "/IQ" + radar_name + dayinfo + "/"

        filenames_iq = []
        for i in np.arange(-ang_tol, ang_tol + ang_step, ang_step):
            ele_str = "{:04d}".format(int(100.0 * (ele + i)))
            filenames_iq.extend(
                glob.glob(
                    datapath_iq
                    + "IQ20"
                    + dayinfo
                    + timeinfo
                    + "*-E"
                    + ele_str
                    + "*.dat"
                )
            )
        if not filenames_iq:
            warn(
                "No files found in "
                + datapath_iq
                + "IQ20"
                + dayinfo
                + timeinfo
                + "_*.dat"
            )
            continue

        # get metadata from status file
        sweep_number = int(scan) - 1
        noise_h = None
        noise_v = None
        rconst_h = None
        rconst_v = None
        for sweep in root.findall("sweep"):
            sweep_number_file = int(sweep.attrib["name"].split(".")[1]) - 1
            if sweep_number_file == sweep_number:
                noise_h = float(
                    (
                        sweep.find("./RADAR/STAT/CALIB/noisepower_frontend_h_inuse")
                    ).attrib["value"]
                )
                rconst_h = float(
                    (sweep.find("./RADAR/STAT/CALIB/rconst_h")).attrib["value"]
                )
                noise_v = float(
                    (
                        sweep.find("./RADAR/STAT/CALIB/noisepower_frontend_v_inuse")
                    ).attrib["value"]
                )
                rconst_v = float(
                    (sweep.find("./RADAR/STAT/CALIB/rconst_v")).attrib["value"]
                )

        radar_aux = pyart.aux_io.read_iq(
            filename[0],
            filenames_iq,
            field_names=field_names,
            prf=prfs[int(scan) - 1],
            noise_h=noise_h,
            noise_v=noise_v,
            rconst_h=rconst_h,
            rconst_v=rconst_v,
            radconst_h=radconst_h,
            radconst_v=radconst_v,
            mfloss_h=mfloss_h,
            mfloss_v=mfloss_v,
            ang_tol=cfg["ang_tol"],
            rng_min=cfg["rmin"],
            rng_max=cfg["rmax"],
            ele_min=cfg["elmin"],
            ele_max=cfg["elmax"],
            azi_min=cfg["azmin"],
            azi_max=cfg["azmax"],
        )

        if radar_aux is None:
            continue

        if radar is None:
            radar = radar_aux
        else:
            radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    if radar is None:
        return radar

    rmin = None
    rmax = None
    elmin = None
    elmax = None
    azmin = None
    azmax = None
    if cfg["rmin"] is not None:
        rmin = cfg["rmin"][ind_rad]
    if cfg["rmax"] is not None:
        rmax = cfg["rmax"][ind_rad]
    if cfg["elmin"] is not None:
        elmin = cfg["elmin"][ind_rad]
    if cfg["elmax"] is not None:
        elmax = cfg["elmax"][ind_rad]
    if cfg["azmin"] is not None:
        azmin = cfg["azmin"][ind_rad]
    if cfg["azmax"] is not None:
        azmax = cfg["azmax"][ind_rad]

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rmin,
        rng_max=rmax,
        ele_min=elmin,
        ele_max=elmax,
        azi_min=azmin,
        azi_max=azmax,
    )


def merge_fields_rainbow(basepath, scan_name, voltime, datatype_list):
    """
    Merge Rainbow fields into a single radar object.

    Parameters
    ----------
    basepath : str
        Name of the base path where to find the data.
    scan_name: str
        Name of the scan.
    voltime : datetime object
        Reference time of the scan.
    datatype_list : list
        Lists of data types to get.

    Returns
    -------
    radar : Radar
        Radar object.
    """
    datapath = os.path.join(basepath, scan_name, voltime.strftime("%Y-%m-%d") + "/")
    fdatetime = voltime.strftime("%Y%m%d%H%M%S") + "00"

    if (datatype_list[0] != "Nh") and (datatype_list[0] != "Nv"):
        filename = glob.glob(os.path.join(datapath, f"{fdatetime}{datatype_list[0]}.*"))
    elif datatype_list[0] == "Nh":
        filename = glob.glob(os.path.join(datapath, f"{fdatetime}dBZ.*"))
    else:
        filename = glob.glob(os.path.join(datapath, f"{fdatetime}dBZv.*"))

    # create radar object
    radar = None
    if not filename:
        warn(
            "No file found in "
            + os.path.join(datapath, f"{fdatetime}{datatype_list[0]}.*")
        )
    else:
        radar = get_data_rainbow(filename[0], datatype_list[0])

    if len(datatype_list) == 1:
        return radar

    # add other fields in the same scan
    for datatype in datatype_list[1:]:
        if datatype not in ("Nh", "Nv"):
            filename = glob.glob(os.path.join(datapath, f"{fdatetime}{datatype}.*"))
        elif datatype == "Nh":
            filename = glob.glob(os.path.join(datapath, f"{fdatetime}dBZ.*"))
        else:
            filename = glob.glob(os.path.join(datapath, f"{fdatetime}dBZv.*"))

        if not filename:
            warn(
                "No file found in " + os.path.join(datapath, f"{fdatetime}{datatype}.*")
            )
        else:
            radar_aux = get_data_rainbow(filename[0], datatype)
            if radar_aux is None:
                continue

            if radar is None:
                radar = radar_aux
            else:
                for field_name in radar_aux.fields.keys():
                    break
                try:
                    radar.add_field(field_name, radar_aux.fields[field_name])
                except (ValueError, KeyError) as ee:
                    warn(
                        "Unable to add field '" + field_name + "' to radar object"
                        ": (%s)" % str(ee)
                    )

    return radar


def merge_fields_psr_spectra(
    basepath,
    basepath_psr,
    scan_name,
    voltime,
    datatype_list,
    undo_txcorr=True,
    fold=True,
    positive_away=True,
    cpi="low_prf",
    ang_tol=0.5,
    azi_min=None,
    azi_max=None,
    ele_min=None,
    ele_max=None,
    rng_min=None,
    rng_max=None,
):
    """
    merge Rainbow fields into a single radar object.

    Parameters
    ----------
    basepath : str
        name of the base path where to find the data
    basepath_psr : str
        name of the base path where to find the PSR data
    scan_name: str
        name of the scan
    voltime : datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    undo_txcorr: Bool
        If True the correction of the transmitted power is removed from the
        noise signal
    fold: Bool
        If True the spectra is folded so that 0-Doppler is in the middle
    positive_away: Bool
        If True the spectra is reversed so that positive velocities are
        away from the radar
    cpi : str
        The CPI to use. Can be 'low_prf', 'intermediate_prf', 'high_prf' or
        'all'
    ang_tol : float
        Tolerated angle distance between nominal radar angle and angle in
        PSR files
    azi_min, azi_max, ele_min, ele_max : float or None
        The minimum and maximum angles to keep (deg)
    rng_min, rng_max : float or None
        The minimum and maximum ranges to keep (m)

    Returns
    -------
    psr : radar spectra object
        radar spectra object

    """
    psr = None
    # Find reference file
    datapath = basepath + scan_name + voltime.strftime("%Y-%m-%d") + "/"
    fdatetime = voltime.strftime("%Y%m%d%H%M%S") + "00"
    filename = glob.glob(datapath + fdatetime + "dBZ.*")
    if not filename:
        warn("No reference file found in " + datapath + fdatetime + "dBZ.*")
        return psr
    filename = filename[0]

    datapath_psr = basepath_psr + scan_name + voltime.strftime("%Y-%m-%d") + "/"
    for datatype in datatype_list:
        if datatype in ("ShhADUu", "sNADUh"):
            filestr = datapath_psr + fdatetime + "_*.ufh.psr.rd"
        elif datatype in ("SvvADUu", "sNADUv"):
            filestr = datapath_psr + fdatetime + "_*.ufv.psr.rd"
        else:
            warn("Unknown data type " + datatype)
            continue

        filenames_psr = glob.glob(filestr)
        if not filenames_psr:
            warn("No file found in " + filestr)
            continue

        psr_aux = pyart.aux_io.read_rainbow_psr_spectra(
            filename,
            filenames_psr,
            field_names=[get_fieldname_pyart(datatype)],
            undo_txcorr=undo_txcorr,
            fold=fold,
            positive_away=positive_away,
            cpi=cpi,
            ang_tol=ang_tol,
            azi_min=azi_min,
            azi_max=azi_max,
            ele_min=ele_min,
            ele_max=ele_max,
            rng_min=rng_min,
            rng_max=rng_max,
        )

        if psr_aux is None:
            continue

        if psr is None:
            psr = psr_aux
            continue

        for field_name in psr_aux.fields.keys():
            try:
                psr.add_field(field_name, psr_aux.fields[field_name])
            except (ValueError, KeyError) as ee:
                warn(
                    "Unable to add field '" + field_name + "' to radar spectra "
                    "object: (%s)" % str(ee)
                )
    return psr


def merge_fields_psr(
    basepath,
    basepath_psr,
    scan_name,
    voltime,
    datatype_list,
    undo_txcorr=True,
    cpi="low_prf",
    ang_tol=0.5,
    azi_min=None,
    azi_max=None,
    ele_min=None,
    ele_max=None,
    rng_min=None,
    rng_max=None,
):
    """
    merge Rainbow fields into a single radar object.

    Parameters
    ----------
    basepath : str
        name of the base path where to find the data
    basepath_psr : str
        name of the base path where to find the PSR data
    scan_name: str
        name of the scan
    voltime : datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    undo_txcorr : Bool
        If true the correction for transmitted power is undone when
        getting the noise
    cpi : str
        The CPI to use. Can be 'low_prf', 'intermediate_prf', 'high_prf',
        'mean', 'all'. If 'mean' the mean within the angle step is taken
    ang_tol : float
        Tolerated angle distance between nominal radar angle and angle in
        PSR files
    azi_min, azi_max, ele_min, ele_max : float or None
        The minimum and maximum angles to keep (deg)
    rng_min, rng_max : float or None
        The minimum and maximum ranges to keep (m)

    Returns
    -------
    radar : Radar
        radar object

    """
    radar = None

    # Find reference file
    datapath = basepath + scan_name + voltime.strftime("%Y-%m-%d") + "/"
    fdatetime = voltime.strftime("%Y%m%d%H%M%S") + "00"
    filename = glob.glob(datapath + fdatetime + "dBZ.*")
    if not filename:
        warn("No reference file found in " + datapath + fdatetime + "dBZ.*")
        return radar
    filename = filename[0]

    datapath_psr = basepath_psr + scan_name + voltime.strftime("%Y-%m-%d") + "/"
    for datatype in datatype_list:
        if datatype in ("Nh", "NdBADUh", "NdBmh", "TXh"):
            filestr = datapath_psr + fdatetime + "_*.ufh.psr.rd"
        elif datatype in ("Nv", "NdBADUv", "NdBmv", "TXv"):
            filestr = datapath_psr + fdatetime + "_*.ufv.psr.rd"
        else:
            warn("Unknown data type " + datatype)
            continue

        filenames_psr = glob.glob(filestr)
        if not filenames_psr:
            warn("No file found in " + filestr)
            continue

        radar_aux = pyart.aux_io.read_rainbow_psr(
            filename,
            filenames_psr,
            field_names=[get_fieldname_pyart(datatype)],
            undo_txcorr=undo_txcorr,
            cpi=cpi,
            ang_tol=ang_tol,
            azi_min=azi_min,
            azi_max=azi_max,
            ele_min=ele_min,
            ele_max=ele_max,
            rng_min=rng_min,
            rng_max=rng_max,
        )

        if radar_aux is None:
            continue

        if radar is None:
            radar = radar_aux
            continue

        for field_name in radar_aux.fields.keys():
            try:
                radar.add_field(field_name, radar_aux.fields[field_name])
            except (ValueError, KeyError) as ee:
                warn(
                    "Unable to add field '" + field_name + "' to radar object"
                    ": (%s)" % str(ee)
                )

    return radar


def merge_fields_rad4alp_grid(voltime, datatype_list, cfg, ind_rad=0, ftype="METRANET"):
    """
    merge rad4alp Cartesian products

    Parameters
    ----------
    voltime: datetime object
        reference time of the scan
    datatype : str
        name of the data type to read
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index
    ftype : str
        File type. Can be 'METRANET', 'gif' or 'bin'

    Returns
    -------
    radar : Radar
        radar object

    """
    grid = None
    for datatype in datatype_list:
        # read product data file
        acronym, termination = get_rad4alp_prod_fname(datatype)
        if datatype.startswith("d") and datatype not in (
            "dGZC",
            "dACC",
            "dACCH",
            "dARC",
        ):
            dir_day = voltime - datetime.timedelta(days=1)
            timeinfo = "2400"
            dayinfo = dir_day.strftime("%y%j")
        else:
            dir_day = voltime
            timeinfo = voltime.strftime("%H%M")
            dayinfo = voltime.strftime("%y%j")

        basename_prod = acronym + dayinfo
        prod_field = get_fieldname_pyart(datatype)

        datapath_prod = get_rad4alp_grid_dir(
            cfg["datapath"][ind_rad],
            dir_day,
            datatype,
            acronym,
            path_convention=cfg["path_convention"][ind_rad],
        )

        filename_prod = glob.glob(
            datapath_prod + basename_prod + timeinfo + "*" + termination
        )
        if not filename_prod:
            warn(
                "No file found in "
                + datapath_prod
                + basename_prod
                + timeinfo
                + "*"
                + termination
            )
            continue
        filename_prod = filename_prod[0]

        if ftype == "METRANET":
            grid_aux = pyart.aux_io.read_cartesian_metranet(
                filename_prod, reader=cfg["metranet_read_lib"]
            )
        elif ftype == "gif":
            grid_aux = pyart.aux_io.read_gif(filename_prod)
        else:
            grid_aux = pyart.aux_io.read_bin(filename_prod)

        if grid_aux is None:
            continue

        if grid is None:
            grid = grid_aux
        else:
            if not datatype.startswith("OZC"):
                grid.add_field(prod_field, grid_aux.fields[prod_field])
            else:
                # Zh CAPPI product. Merge grids
                grid = merge_grids(grid, grid_aux)

    if grid is None:
        return grid

    # Crop the data
    lat_min = cfg.get("latmin", None)
    lat_max = cfg.get("latmax", None)
    lon_min = cfg.get("lonmin", None)
    lon_max = cfg.get("lonmax", None)
    alt_min = cfg.get("altmin", None)
    alt_max = cfg.get("altmax", None)
    ix_min = cfg.get("ixmin", None)
    iy_min = cfg.get("iymin", None)
    iz_min = cfg.get("izmin", None)
    nx = cfg.get("nx", None)
    ny = cfg.get("ny", None)
    nz = cfg.get("nz", None)

    return crop_grid(
        grid,
        lat_min=lat_min,
        lat_max=lat_max,
        lon_min=lon_min,
        lon_max=lon_max,
        alt_min=alt_min,
        alt_max=alt_max,
        nx=nx,
        ny=ny,
        nz=nz,
        ix_min=ix_min,
        iy_min=iy_min,
        iz_min=iz_min,
    )


def merge_fields_sat_grid(voltime, datatype_list, cfg, ind_rad=0):
    """
    merge satellite Cartesian products

    Parameters
    ----------
    voltime: datetime object
        reference time of the scan
    datatype : str
        name of the data type to read
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    grid = None
    daydir = voltime.strftime("%Y/%m/%d/")
    dayinfo = voltime.strftime("%Y%m%d%H%M")
    datapath = cfg["satpath"][ind_rad] + daydir
    if not os.path.isdir(datapath):
        # warn("WARNING: Unknown datapath '%s'" % datapath)
        return grid

    filename = glob.glob(datapath + "MSG?_ccs4_" + dayinfo + "*_rad_PLAX.nc")
    if not filename:
        warn("No file found in " + datapath + "MSG?_ccs4_" + dayinfo + "*_rad_PLAX.nc")
        return grid

    field_names = []
    for datatype in datatype_list:
        field_names.append(get_fieldname_pyart(datatype))

    grid = pyart.aux_io.read_cf1_cartesian(filename[0], field_names)

    if grid is None:
        return grid

    # Crop the data
    lat_min = cfg.get("latmin", None)
    lat_max = cfg.get("latmax", None)
    lon_min = cfg.get("lonmin", None)
    lon_max = cfg.get("lonmax", None)
    alt_min = cfg.get("altmin", None)
    alt_max = cfg.get("altmax", None)
    ix_min = cfg.get("ixmin", None)
    iy_min = cfg.get("iymin", None)
    iz_min = cfg.get("izmin", None)
    nx = cfg.get("nx", None)
    ny = cfg.get("ny", None)
    nz = cfg.get("nz", None)

    return crop_grid(
        grid,
        lat_min=lat_min,
        lat_max=lat_max,
        lon_min=lon_min,
        lon_max=lon_max,
        alt_min=alt_min,
        alt_max=alt_max,
        nx=nx,
        ny=ny,
        nz=nz,
        ix_min=ix_min,
        iy_min=iy_min,
        iz_min=iz_min,
    )


def merge_fields_mf_grid(
    voltime, datatype_list, dataset_list, scan_list, cfg, ind_rad=0, ftype="bin"
):
    """
    merge MF Cartesian products

    Parameters
    ----------
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        name of the data types to read
    dataset_list : list
        list containing the date format to find the path to each data type
    scan_list : list
        list of scans
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index
    ftype : str
        File type. Can be 'png', 'bin', 'grib', 'dat' or 'nc'

    Returns
    -------
    radar : Radar
        radar object

    """
    grid = None

    fpath_strf = dataset_list[0][
        dataset_list[0].find("D") + 2 : dataset_list[0].find("F") - 2
    ]
    fdate_strf = dataset_list[0][dataset_list[0].find("F") + 2 : -1]
    datapath = cfg["datapath"][ind_rad] + voltime.strftime(fpath_strf) + "/"

    for scan, datatype in zip(scan_list, cfg["BinFileParams"]["datatype"]):
        filename = None
        filenames = glob.glob(datapath + "*" + scan + "*")
        for filename_aux in filenames:
            fdatetime = find_date_in_file_name(filename_aux, date_format=fdate_strf)
            if fdatetime == voltime:
                filename = filename_aux
        if filename is None:
            warn(f"No file found for scan {scan}")
            continue

        prod_field = get_fieldname_pyart(datatype)
        grid_aux = None
        if ftype == "bin":
            grid_aux = pyart.aux_io.read_bin_mf(
                filename,
                field_name=prod_field,
                xres=cfg["BinFileParams"]["xres"],
                yres=cfg["BinFileParams"]["yres"],
                nx=cfg["BinFileParams"]["nx"],
                ny=cfg["BinFileParams"]["ny"],
                nz=cfg["BinFileParams"]["nz"],
                dtype=cfg["BinFileParams"]["dtype"],
                date_format=cfg["BinFileParams"]["date_format"],
                added_time=cfg["BinFileParams"]["added_time"],
                x_offset=cfg["BinFileParams"]["x_offset"],
                y_offset=cfg["BinFileParams"]["y_offset"],
                lat_0=cfg["BinFileParams"]["lat_0"],
                lon_0=cfg["BinFileParams"]["lon_0"],
                proj=cfg["BinFileParams"]["proj"],
            )
        elif ftype == "png":
            grid_aux = pyart.aux_io.read_png(
                filename,
                field_name=prod_field,
                xres=cfg["BinFileParams"]["xres"],
                yres=cfg["BinFileParams"]["yres"],
                nz=cfg["BinFileParams"]["nz"],
                date_format=cfg["BinFileParams"]["date_format"],
                added_time=cfg["BinFileParams"]["added_time"],
                x_offset=cfg["BinFileParams"]["x_offset"],
                y_offset=cfg["BinFileParams"]["y_offset"],
                lat_0=cfg["BinFileParams"]["lat_0"],
                lon_0=cfg["BinFileParams"]["lon_0"],
                proj=cfg["BinFileParams"]["proj"],
            )
        elif ftype == "grib":
            grid_aux = pyart.aux_io.read_grib(filename)
        elif ftype == "dat":
            grid_aux = pyart.aux_io.read_dat_mf(filename)
        elif ftype == "nc":
            grid_aux = pyart.aux_io.read_cf1_cartesian_mf(filename)
        if grid_aux is None:
            continue

        if grid is None:
            grid = grid_aux
        else:
            grid.add_field(prod_field, grid_aux.fields[prod_field])

    if grid is None:
        return grid

    # Crop the data
    lat_min = cfg.get("latmin", None)
    lat_max = cfg.get("latmax", None)
    lon_min = cfg.get("lonmin", None)
    lon_max = cfg.get("lonmax", None)
    alt_min = cfg.get("altmin", None)
    alt_max = cfg.get("altmax", None)
    ix_min = cfg.get("ixmin", None)
    iy_min = cfg.get("iymin", None)
    iz_min = cfg.get("izmin", None)
    nx = cfg.get("nx", None)
    ny = cfg.get("ny", None)
    nz = cfg.get("nz", None)

    return crop_grid(
        grid,
        lat_min=lat_min,
        lat_max=lat_max,
        lon_min=lon_min,
        lon_max=lon_max,
        alt_min=alt_min,
        alt_max=alt_max,
        nx=nx,
        ny=ny,
        nz=nz,
        ix_min=ix_min,
        iy_min=iy_min,
        iz_min=iz_min,
    )


def merge_fields_pyrad(
    basepath,
    loadname,
    voltime,
    datatype_list,
    dataset_list,
    product_list,
    rng_min=None,
    rng_max=None,
    azi_min=None,
    azi_max=None,
    ele_min=None,
    ele_max=None,
    termination=".nc",
):
    """
    merge fields from Pyrad-generated files into a single radar object.
    Accepted file types are CFRadial and ODIM.

    Parameters
    ----------
    basepath : str
        name of the base path where to find the data
    loadname: str
        name of the saving directory
    voltime : datetime object
        reference time of the scan
    datatype_list : list
        list of data types to get
    dataset_list : list
        list of datasets that produced the data type to get.
        Used to get path.
    product_list : list
        list of products. Used to get path
    rng_min, rng_max : float
        The range limits [m]. If None the entire coverage of the radar is
        going to be used
    ele_min, ele_max, azi_min, azi_max : float or None
        The limits of the grid [deg]. If None the limits will be the limits
        of the radar volume
    termination : str
        file termination type. Can be '.nc' or '.h*'

    Returns
    -------
    radar : Radar
        radar object

    """
    fdatetime = voltime.strftime("%Y%m%d%H%M%S")

    radar = None
    for i, dataset in enumerate(dataset_list):
        datapath = (
            basepath
            + loadname
            + "/"
            + voltime.strftime("%Y-%m-%d")
            + "/"
            + dataset
            + "/"
            + product_list[i]
            + "/"
        )
        filename = glob.glob(
            datapath + fdatetime + "*" + datatype_list[i] + termination
        )
        if not filename:
            warn(
                "No file found in "
                + datapath
                + fdatetime
                + "*"
                + datatype_list[i]
                + ".nc"
            )
            continue

        if termination == ".nc":
            try:
                radar_aux = pyart.io.read_cfradial(filename[0])
            except (OSError, KeyError) as ee:
                warn(str(ee))
                warn("Unable to read file " + filename[0])
                radar_aux = None
        else:
            try:
                radar_aux = pyart.aux_io.read_odim_h5(filename[0])
            except OSError as ee:
                warn(str(ee))
                warn("Unable to read file " + filename[0])
                radar_aux = None

        if radar_aux is None:
            continue

        if radar is None:
            radar = radar_aux
            continue

        radar = add_field(radar, radar_aux)

    if radar is None:
        return radar

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rng_min,
        rng_max=rng_max,
        ele_min=ele_min,
        ele_max=ele_max,
        azi_min=azi_min,
        azi_max=azi_max,
    )


def merge_fields_gecsx(
    basepath,
    loadname,
    datatype_list,
    dataset_list,
    product_list,
    rng_min=None,
    rng_max=None,
    azi_min=None,
    azi_max=None,
    ele_min=None,
    ele_max=None,
):
    """
    merge fields from GECSX Pyrad-generated files into a single radar object.
    Accepted file types are CFRadial.

    Parameters
    ----------
    basepath : str
        name of the base path where to find the data
    loadname: str
        name of the saving directory
    datatype_list : list
        list of data types to get
    dataset_list : list
        list of datasets that produced the data type to get.
        Used to get path.
    product_list : list
        list of products. Used to get path
    rng_min, rng_max : float
        The range limits [m]. If None the entire coverage of the radar is
        going to be used
    ele_min, ele_max, azi_min, azi_max : float or None
        The limits of the grid [deg]. If None the limits will be the limits
        of the radar volume
    Returns
    -------
    radar : Radar
        radar object

    """
    radar = None
    for i, dataset in enumerate(dataset_list):
        datapath = basepath + loadname + "/" + dataset + "/" + product_list[i] + "/"
        filename = glob.glob(datapath + "*" + datatype_list[i] + "*.nc")
        if not filename:
            warn("No file found in " + datapath + "*" + datatype_list[i] + "*.nc")
            continue
        try:
            radar_aux = pyart.io.read_cfradial(filename[0])
        except (OSError, KeyError) as ee:
            warn(str(ee))
            warn("Unable to read file " + filename[0])
            radar_aux = None

        if radar_aux is None:
            continue

        if radar is None:
            radar = radar_aux
            continue

        radar = add_field(radar, radar_aux)

    if radar is None:
        return radar

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rng_min,
        rng_max=rng_max,
        ele_min=ele_min,
        ele_max=ele_max,
        azi_min=azi_min,
        azi_max=azi_max,
    )


def merge_fields_pyradicon(
    basepath,
    voltime,
    datatype_list,
    dataset_list,
    cfg,
    rng_min=None,
    rng_max=None,
    azi_min=None,
    azi_max=None,
    ele_min=None,
    ele_max=None,
    termination=".nc",
):
    """
    merge fields from Pyrad-generated files into a single radar object.
    Accepted file types are CFRadial and ODIM.

    Parameters
    ----------
    basepath : str
        name of the base path where to find the data
    voltime : datetime object
        reference time of the scan
    datatype_list : list
        list of data types to get
    dataset_list : list
        list of datasets that produced the data type to get.
        Used to get path.
    cfg : dictionary of dictionaries
        configuration info
    rng_min, rng_max : float
        The range limits [m]. If None the entire coverage of the radar is
        going to be used
    ele_min, ele_max, azi_min, azi_max : float or None
        The limits of the grid [deg]. If None the limits will be the limits
        of the radar volume
    termination : str
        file termination type. Can be '.nc' or '.h*'

    Returns
    -------
    radar : Radar
        radar object

    """
    voltime.strftime("%Y%m%d%H%M%S")

    radar = None
    for datatype, dataset in zip(datatype_list, dataset_list):
        filename = find_pyradicon_file(basepath, voltime, datatype, cfg, dataset)
        if filename is None:
            continue

        if termination == ".nc":
            try:
                radar_aux = pyart.io.read_cfradial(filename)
            except (OSError, KeyError) as ee:
                warn(str(ee))
                warn("Unable to read file " + filename)
                radar_aux = None
        else:
            try:
                radar_aux = pyart.aux_io.read_odim_h5(filename)
            except OSError as ee:
                warn(str(ee))
                warn("Unable to read file " + filename)
                radar_aux = None

        if radar_aux is None:
            continue

        if radar is None:
            radar = radar_aux
            continue

        radar = add_field(radar, radar_aux)

    if radar is None:
        return radar

    return pyart.util.subset_radar(
        radar,
        radar.fields.keys(),
        rng_min=rng_min,
        rng_max=rng_max,
        ele_min=ele_min,
        ele_max=ele_max,
        azi_min=azi_min,
        azi_max=azi_max,
    )


def merge_fields_pyrad_spectra(
    basepath,
    loadname,
    voltime,
    datatype_list,
    dataset_list,
    product_list,
    rng_min=None,
    rng_max=None,
    azi_min=None,
    azi_max=None,
    ele_min=None,
    ele_max=None,
    termination=".nc",
):
    """
    merge fields from Pyrad-generated files into a single radar spectra
    object. Accepted file types are netcdf

    Parameters
    ----------
    basepath : str
        name of the base path where to find the data
    loadname: str
        name of the saving directory
    voltime : datetime object
        reference time of the scan
    datatype_list : list
        list of data types to get
    dataset_list : list
        list of datasets that produced the data type to get.
        Used to get path.
    product_list : list
        list of products. Used to get path
    rng_min, rng_max : float
        The range limits [m]. If None the entire coverage of the radar is
        going to be used
    ele_min, ele_max, azi_min, azi_max : float or None
        The limits of the grid [deg]. If None the limits will be the limits
        of the radar volume
    termination : str
        file termination type. Can be '.nc' or '.h*'

    Returns
    -------
    radar : Radar
        radar object

    """
    fdatetime = voltime.strftime("%Y%m%d%H%M%S")

    radar = None
    for i, dataset in enumerate(dataset_list):
        datapath = (
            basepath
            + loadname
            + "/"
            + voltime.strftime("%Y-%m-%d")
            + "/"
            + dataset
            + "/"
            + product_list[i]
            + "/"
        )
        filename = glob.glob(
            datapath + fdatetime + "*" + datatype_list[i] + termination
        )
        if not filename:
            warn(
                "No file found in "
                + datapath
                + fdatetime
                + "*"
                + datatype_list[i]
                + ".nc"
            )
            continue

        radar_aux = None
        if termination == ".nc":
            try:
                radar_aux = pyart.aux_io.read_spectra(filename[0])
            except OSError as ee:
                warn(str(ee))
                warn("Unable to read file " + filename[0])
                radar_aux = None
        # else:
        #     try:
        #         radar_aux = pyart.aux_io.read_odim_h5(filename[0])
        #     except OSError as ee:
        #         warn(str(ee))
        #         warn('Unable to read file '+filename[0])
        if radar_aux is None:
            continue

        if radar is None:
            radar = radar_aux
            continue

        radar = add_field(radar, radar_aux)

    if radar is None:
        return radar

    return pyart.util.subset_radar_spectra(
        radar,
        radar.fields.keys(),
        rng_min=rng_min,
        rng_max=rng_max,
        ele_min=ele_min,
        ele_max=ele_max,
        azi_min=azi_min,
        azi_max=azi_max,
    )


def merge_fields_pyradgrid(
    basepath,
    loadname,
    voltime,
    datatype_list,
    dataset_list,
    product_list,
    cfg,
    termination=".nc",
):
    """
    merge fields from Pyrad-generated files into a single radar object.
    Accepted file types are CFRadial and ODIM.

    Parameters
    ----------
    basepath : str
        name of the base path where to find the data
    loadname: str
        name of the saving directory
    voltime : datetime object
        reference time of the scan
    datatype_list : list
        list of data types to get
    dataset_list : list
        list of datasets that produced the data type to get.
        Used to get path.
    product_list : list
        list of products. Used to get path
    cfg : dict
        dictionary containing configuration parameters
    termination : str
        file termination type. Can be '.nc' or '.h*'

    Returns
    -------
    grid : Grid
        grid object

    """
    grid = None
    fdatetime = voltime.strftime("%Y%m%d%H%M%S")

    for i, dataset in enumerate(dataset_list):
        datapath = (
            basepath
            + loadname
            + "/"
            + voltime.strftime("%Y-%m-%d")
            + "/"
            + dataset
            + "/"
            + product_list[i]
            + "/"
        )
        filename = glob.glob(
            datapath + fdatetime + "*" + datatype_list[i] + termination
        )
        if not filename:
            warn(
                "No file found in "
                + datapath
                + fdatetime
                + "*"
                + datatype_list[i]
                + termination
            )
            continue

        try:
            if termination == ".nc":
                grid_aux = pyart.io.read_grid(filename[0])
            else:
                grid_aux = pyart.aux_io.read_odim_grid_h5(filename[0])
        except OSError as ee:
            warn(str(ee))
            warn("Unable to read file " + filename[0])
            grid_aux = None

        if grid_aux is None:
            continue

        if grid is None:
            grid = grid_aux
        else:
            for field_name in grid_aux.fields:
                grid.add_field(field_name, grid_aux.fields[field_name])

    if grid is None:
        return grid

    # Crop the data
    lat_min = cfg.get("latmin", None)
    lat_max = cfg.get("latmax", None)
    lon_min = cfg.get("lonmin", None)
    lon_max = cfg.get("lonmax", None)
    alt_min = cfg.get("altmin", None)
    alt_max = cfg.get("altmax", None)
    ix_min = cfg.get("ixmin", None)
    iy_min = cfg.get("iymin", None)
    iz_min = cfg.get("izmin", None)
    nx = cfg.get("nx", None)
    ny = cfg.get("ny", None)
    nz = cfg.get("nz", None)

    return crop_grid(
        grid,
        lat_min=lat_min,
        lat_max=lat_max,
        lon_min=lon_min,
        lon_max=lon_max,
        alt_min=alt_min,
        alt_max=alt_max,
        nx=nx,
        ny=ny,
        nz=nz,
        ix_min=ix_min,
        iy_min=iy_min,
        iz_min=iz_min,
    )


def merge_fields_dem(basepath, scan_name, datatype_list):
    """
    merge DEM fields into a single radar object.

    Parameters
    ----------
    basepath : str
        name of the base path where to find the data
    scan_name: str
        name of the scan
    datatype_list : list
        lists of data types to get

    Returns
    -------
    radar : Radar
        radar object

    """
    scan_name_aux = scan_name.partition("/")[0]
    radar = None
    # add other fields in the same scan
    for datatype in datatype_list:
        datapath = basepath + datatype + "/" + scan_name + "/"
        filename = glob.glob(datapath + datatype + "_" + scan_name_aux)
        if not filename:
            warn("No file found in " + datapath + datatype + "_" + scan_name_aux)
            continue

        radar_aux = get_data_rainbow(filename[0], datatype)
        if radar is None:
            radar = radar_aux
            continue

        for field_name in radar_aux.fields.keys():
            break
        try:
            radar = radar.add_field(field_name, radar_aux.fields[field_name])
        except (ValueError, KeyError):
            warn("Unable to add field " + field_name + " to radar object")

    return radar


def merge_fields_icon(filename_list):
    """
    merge COSMO fields in Rainbow file format

    Parameters
    ----------
    filename_list : str
        list of file paths where to find the data

    Returns
    -------
    radar : Radar
        radar object

    """
    # add other COSMO fields in the same scan
    radar = None
    for filename in filename_list:
        try:
            radar_aux = pyart.aux_io.read_rainbow_wrl(filename)
        except OSError as ee:
            warn(str(ee))
            warn("Unable to read file " + filename)
            continue
        if radar_aux is None:
            continue

        if radar is None:
            radar = radar_aux
            continue

        for field_name in radar_aux.fields.keys():
            break
        radar.add_field(field_name, radar_aux.fields[field_name])

    return radar


def get_data_rainbow(filename, datatype):
    """
    gets rainbow radar data

    Parameters
    ----------
    filename : str
        name of file containing rainbow data
    datatype : str
        field name

    Returns
    -------
    radar : Radar or None
        radar object if the reading of the data has been successful.
        None otherwise

    """
    try:
        radar = pyart.aux_io.read_rainbow_wrl(filename)
    except OSError as ee:
        warn(str(ee))
        warn("Unable to read file " + filename)
        return None
    if radar is None:
        return None

    if datatype in ("Nh", "Nv"):
        try:
            with open(filename, "rb") as fid:
                rbf = wrl.io.read_rainbow(fid, loaddata=True)
                fid.close()
        except OSError as ee:
            warn(str(ee))
            warn("Unable to read file " + filename)
            return None

        # check the number of slices
        nslices = int(rbf["volume"]["scan"]["pargroup"]["numele"])
        if nslices > 1:
            common_slice_info = rbf["volume"]["scan"]["slice"][0]
        else:
            common_slice_info = rbf["volume"]["scan"]["slice"]

        if datatype == "Nh":
            noisedBZ1km_h = float(common_slice_info["noise_power_dbz"])
            noisedBZ_h = pyart.retrieve.compute_noisedBZ(
                radar.nrays,
                noisedBZ1km_h,
                radar.range["data"],
                1.0,
                noise_field="noisedBZ_hh",
            )
            radar.fields = {}
            radar.add_field("noisedBZ_hh", noisedBZ_h)
        else:
            noisedBZ1km_v = float(common_slice_info["noise_power_dbz_dpv"])
            noisedBZ_v = pyart.retrieve.compute_noisedBZ(
                radar.nrays,
                noisedBZ1km_v,
                radar.range["data"],
                1.0,
                noise_field="noisedBZ_vv",
            )
            radar.fields = {}
            radar.add_field("noisedBZ_vv", noisedBZ_v)

    return radar


def get_data_rad4alp(filename, datatype_list, scan_name, cfg, ind_rad=0):
    """
    gets rad4alp radar data

    Parameters
    ----------
    filename : str
        name of file containing rainbow data
    datatype_list : list of strings
        list of data fields to get
    scan_name : str
        name of the elevation (001 to 020)
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object. None if the reading has not been successful

    """
    metranet_field_names = {}
    for datatype in datatype_list:
        if datatype not in ("Nh", "Nv"):
            metranet_field_names.update(get_datatype_metranet(datatype))

    if cfg["path_convention"][ind_rad] == "LTE":
        radar = pyrad_MCH(filename, field_names=metranet_field_names)
    else:
        try:
            radar = pyart.aux_io.read_metranet(
                filename,
                field_names=metranet_field_names,
                reader=cfg["metranet_read_lib"],
            )
        except (ValueError, TypeError) as ee:
            warn("Unable to read file '" + filename + ": (%s)" % str(ee))
            return None

    if ("Nh" not in datatype_list) and ("Nv" not in datatype_list):
        return radar

    # create noise moments
    # read radar information in status file
    voltime = get_datetime(filename, "RAD4ALP:dBZ")
    root = read_status(voltime, cfg, ind_rad=ind_rad)
    if root is None:
        return radar

    sweep_number = int(scan_name) - 1
    if "Nh" in datatype_list:
        found = False
        for sweep in root.findall("sweep"):
            sweep_number_file = int(sweep.attrib["name"].split(".")[1]) - 1
            if sweep_number_file == sweep_number:
                noise_h = sweep.find("./RADAR/STAT/CALIB/noisepower_frontend_h_inuse")
                rconst_h = sweep.find("./RADAR/STAT/CALIB/rconst_h")
                if noise_h is None or rconst_h is None:
                    warn(
                        "Horizontal channel noise power not "
                        + "available for sweep "
                        + scan_name
                    )
                    break

                noisedBADU_h = 10.0 * np.log10(float(noise_h.attrib["value"]))
                rconst_h = float(rconst_h.attrib["value"])

                noisedBZ_h = pyart.retrieve.compute_noisedBZ(
                    radar.nrays,
                    noisedBADU_h + rconst_h,
                    radar.range["data"],
                    100.0,
                    noise_field="noisedBZ_hh",
                )

                radar.add_field("noisedBZ_hh", noisedBZ_h)

                found = True
        if not found:
            warn(
                "Horizontal channel noise power not "
                + "available for sweep "
                + scan_name
            )

    if "Nv" in datatype_list:
        found = False
        for sweep in root.findall("sweep"):
            sweep_number_file = int(sweep.attrib["name"].split(".")[1]) - 1
            if sweep_number_file == sweep_number:
                noise_v = sweep.find("./RADAR/STAT/CALIB/noisepower_frontend_v_inuse")
                rconst_v = sweep.find("./RADAR/STAT/CALIB/rconst_v")
                if noise_v is None or rconst_v is None:
                    warn(
                        "Vertical channel noise power not "
                        + "available for sweep "
                        + scan_name
                    )
                    break

                noisedBADU_v = 10.0 * np.log10(float(noise_v.attrib["value"]))
                rconst_v = float(rconst_v.attrib["value"])

                noisedBZ_v = pyart.retrieve.compute_noisedBZ(
                    radar.nrays,
                    noisedBADU_v + rconst_v,
                    radar.range["data"],
                    100.0,
                    noise_field="noisedBZ_vv",
                )

                radar.add_field("noisedBZ_vv", noisedBZ_v)

                found = True
        if not found:
            warn(
                "Horizontal channel noise power not "
                + "available for sweep "
                + scan_name
            )

    return radar


def get_data_odim(filename, datatype_list, scan_name, cfg, ind_rad=0):
    """
    gets ODIM radar data

    Parameters
    ----------
    filename : str
        name of file containing odim data
    datatype_list : list of strings
        list of data fields to get
    scan_name : str
        name of the elevation (001 to 020)
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object. None if the reading has not been successful

    """

    odim_field_names = {}
    # Use custom name mapping
    for datatype in datatype_list:
        if datatype in cfg["DataTypeIDInFiles"][ind_rad]:
            odim_field_names[
                cfg["DataTypeIDInFiles"][ind_rad][datatype]
            ] = get_fieldname_pyart(datatype)
        else:
            odim_field_names.update(get_datatype_odim(datatype))

    try:
        if cfg["MFScale"]:
            # assumes only a data type per file
            offset = 0.0
            gain = 1.0
            nodata = np.nan
            undetect = np.nan
            use_file_conversion = True
            if datatype_list[0] == "dBuZ":
                # SEMAFOR campaign:
                # offset = -41.25
                # gain = 0.5

                # regular
                offset = -10.5
                gain = 1

                nodata = 0
                undetect = 255
                use_file_conversion = False
            elif datatype_list[0] == "ZDR":
                offset = -9.95
                gain = 0.1
                nodata = 255
                undetect = 255
                use_file_conversion = False
            elif datatype_list[0] == "RhoHV":
                offset = 0.305
                gain = 0.01
                nodata = 255
                undetect = 255
                use_file_conversion = False
            elif datatype_list[0] == "PhiDP":
                offset = 0
                gain = 1
                nodata = 65535
                undetect = 65535
                use_file_conversion = False
            elif datatype_list[0] == "sigma_zh":
                offset = 0.125
                gain = 0.25
                nodata = 255
                undetect = 0
                use_file_conversion = False

            radar = pyart.aux_io.read_odim_h5(
                filename,
                field_names=odim_field_names,
                offset=offset,
                gain=gain,
                nodata=nodata,
                undetect=undetect,
                use_file_conversion=use_file_conversion,
                include_fields=odim_field_names.values(),
            )

            if datatype_list[0] == "PhiDP":
                radar.fields["differential_phase"]["data"][
                    radar.fields["differential_phase"]["data"] > 180.0
                ] -= 360.0
        else:
            radar = pyart.aux_io.read_odim_h5(filename, field_names=odim_field_names)

            if "differential_phase" in radar.fields:
                # make sure that data is within [-180, 180] deg
                radar.fields["differential_phase"]["data"][
                    radar.fields["differential_phase"]["data"] > 180.0
                ] -= 360.0

    except (ValueError, OSError) as ee:
        warn("Unable to read file '" + filename + ": (%s)" % str(ee))
        return None

    if ("Nh" not in datatype_list) and ("Nv" not in datatype_list):
        return radar

    # create noise moments
    # read radar information in status file
    voltime = get_datetime(filename, "ODIM:dBZ")
    root = read_status(voltime, cfg, ind_rad=ind_rad)
    if root is None:
        return radar

    sweep_number = int(scan_name) - 1
    if "Nh" in datatype_list:
        found = False
        for sweep in root.findall("sweep"):
            sweep_number_file = int(sweep.attrib["name"].split(".")[1]) - 1
            if sweep_number_file == sweep_number:
                noise_h = sweep.find("./RADAR/STAT/CALIB/noisepower_frontend_h_inuse")
                rconst_h = sweep.find("./RADAR/STAT/CALIB/rconst_h")
                if noise_h is None or rconst_h is None:
                    warn(
                        "Horizontal channel noise power not "
                        + "available for sweep "
                        + scan_name
                    )
                    break

                noisedBADU_h = 10.0 * np.log10(float(noise_h.attrib["value"]))
                rconst_h = float(rconst_h.attrib["value"])

                noisedBZ_h = pyart.retrieve.compute_noisedBZ(
                    radar.nrays,
                    noisedBADU_h + rconst_h,
                    radar.range["data"],
                    100.0,
                    noise_field="noisedBZ_hh",
                )

                radar.add_field("noisedBZ_hh", noisedBZ_h)

                found = True
        if not found:
            warn(
                "Horizontal channel noise power not "
                + "available for sweep "
                + scan_name
            )

    if "Nv" in datatype_list:
        found = False
        for sweep in root.findall("sweep"):
            sweep_number_file = int(sweep.attrib["name"].split(".")[1]) - 1
            if sweep_number_file == sweep_number:
                noise_v = sweep.find("./RADAR/STAT/CALIB/noisepower_frontend_v_inuse")
                rconst_v = sweep.find("./RADAR/STAT/CALIB/rconst_v")
                if noise_v is None or rconst_v is None:
                    warn(
                        "Vertical channel noise power not "
                        + "available for sweep "
                        + scan_name
                    )
                    break

                noisedBADU_v = 10.0 * np.log10(float(noise_v.attrib["value"]))
                rconst_v = float(rconst_v.attrib["value"])

                noisedBZ_v = pyart.retrieve.compute_noisedBZ(
                    radar.nrays,
                    noisedBADU_v + rconst_v,
                    radar.range["data"],
                    100.0,
                    noise_field="noisedBZ_vv",
                )

                radar.add_field("noisedBZ_vv", noisedBZ_v)

                found = True
        if not found:
            warn(
                "Horizontal channel noise power not "
                + "available for sweep "
                + scan_name
            )

    return radar


def get_data_odimgrid(filename, datatype_list, cfg, mf_scale=False, ind_rad=0):
    """
    gets ODIM grid data

    Parameters
    ----------
    filename : str
        name of file containing odim data
    datatype_list : list of strings
        list of data fields to get
    cfg : dict
        configuration dictionary
    mf_scale : bool
        Whether to impose a MF scale

    Returns
    -------
    grid : Grid
        Grid object. None if the reading has not been successful

    """
    grid = None
    if mf_scale:
        for datatype in datatype_list:
            if datatype in cfg["DataTypeIDInFiles"][ind_rad]:
                odim_field_name = {
                    cfg["DataTypeIDInFiles"][ind_rad][datatype]: get_fieldname_pyart(
                        datatype
                    )
                }
            else:
                odim_field_name = get_datatype_odim(datatype)
            field_name = get_fieldname_pyart(datatype)

            offset = 0.0
            gain = 1.0
            nodata = np.nan
            undetect = np.nan
            use_file_conversion = True
            if field_name == "reflectivity":
                offset = -10.5
                gain = 1
                nodata = 2047
                undetect = 0
                use_file_conversion = False
            elif field_name == "sigma_zh":
                offset = 0.125
                gain = 0.25
                nodata = 255
                undetect = 0
                use_file_conversion = False
            if grid is None:
                try:
                    grid = pyart.aux_io.read_odim_grid_h5(
                        filename,
                        field_names=odim_field_name,
                        offset=offset,
                        gain=gain,
                        nodata=nodata,
                        undetect=undetect,
                        use_file_conversion=use_file_conversion,
                        include_fields=odim_field_name.keys(),
                    )
                except (ValueError, OSError) as ee:
                    warn("Unable to read file '" + filename + ": (%s)" % str(ee))
            else:
                try:
                    grid_aux = pyart.aux_io.read_odim_grid_h5(
                        filename,
                        field_names=odim_field_name,
                        offset=offset,
                        gain=gain,
                        nodata=nodata,
                        undetect=undetect,
                        use_file_conversion=use_file_conversion,
                        include_fields=odim_field_name.keys(),
                    )
                    if grid_aux is not None:
                        grid.add_field(field_name, grid_aux.fields[field_name])
                except (ValueError, OSError) as ee:
                    warn("Unable to read file '" + filename + ": (%s)" % str(ee))
    else:
        odim_field_names = {}
        for datatype in datatype_list:
            if datatype in cfg["DataTypeIDInFiles"][ind_rad]:
                odim_field_name = {
                    cfg["DataTypeIDInFiles"][ind_rad][datatype]: get_fieldname_pyart(
                        datatype
                    )
                }
            else:
                odim_field_name = get_datatype_odim(datatype)
        try:
            grid = pyart.aux_io.read_odim_grid_h5(
                filename,
                field_names=odim_field_names,
                include_fields=odim_field_name.keys(),
            )
        except (ValueError, OSError) as ee:
            warn("Unable to read file '" + filename + ": (%s)" % str(ee))

    return grid


def get_data_gamic(filename, datatype_list, pulse_width):
    """
    gets GAMIC radar data

    Parameters
    ----------
    filename : str
        name of file containing GAMIC data
    datatype_list : list of strings
        list of data fields to get
    pulse_width : list or None
        List of pulse width of the system

    Returns
    -------
    radar : Radar
        radar object. None if the reading has not been successful

    """
    # gamic_field_names = {}
    # for datatype in datatype_list:
    #     # gamic_field_names.update(get_datatype_gamic(datatype))
    try:
        radar = pyart.aux_io.read_gamic(filename, pulse_width=pulse_width)
    except (ValueError, OSError) as ee:
        warn("Unable to read file '" + filename + ": (%s)" % str(ee))
        return None

    return radar


def get_data_mxpol(filename, datatype_list):
    """
    gets MXPol radar data

    Parameters
    ----------
    filename : str
        name of file containing MXPol data
    datatype_list : list of strings
        list of data fields to get

    Returns
    -------
    radar : Radar
        radar object

    """
    field_names = {}
    for datatype in datatype_list:
        if datatype not in ("Nh", "Nv"):
            field_names.update(get_datatype_metranet(datatype))
    radar = pyrad_MXPOL(filename, field_names=field_names)

    # create secondary moments (TODO)
    if ("Nh" in datatype_list) or ("Nv" in datatype_list):
        pass

    return radar


def add_field(radar_dest, radar_orig):
    """
    adds the fields from orig radar into dest radar. If they are not in the
    same grid, interpolates them to dest grid

    Parameters
    ----------
    radar_dest : radar object
        the destination radar
    radar_orig : radar object
        the radar object containing the original field

    Returns
    -------
    field_dest : dict
        interpolated field and metadata

    """
    if radar_dest is None:
        radar_dest = radar_orig
    else:
        if radar_orig is not None:
            if radar_dest.nrays == radar_orig.nrays:
                if (
                    np.allclose(
                        radar_dest.azimuth["data"],
                        radar_orig.azimuth["data"],
                        atol=0.5,
                        equal_nan=True,
                    )
                ) and (
                    np.allclose(
                        radar_dest.elevation["data"],
                        radar_orig.elevation["data"],
                        atol=0.5,
                        equal_nan=True,
                    )
                ):
                    for field_name in radar_orig.fields.keys():
                        radar_dest.add_field(field_name, radar_orig.fields[field_name])
                else:
                    for field_name in radar_orig.fields.keys():
                        field_interp = interpol_field(
                            radar_dest, radar_orig, field_name
                        )
                        radar_dest.add_field(field_name, field_interp)
            else:
                for field_name in radar_orig.fields.keys():
                    field_interp = interpol_field(radar_dest, radar_orig, field_name)
                    radar_dest.add_field(field_name, field_interp)

    return radar_dest


def interpol_field(radar_dest, radar_orig, field_name, fill_value=None, ang_tol=0.5):
    """
    interpolates field field_name contained in radar_orig to the grid in
    radar_dest

    Parameters
    ----------
    radar_dest : radar object
        the destination radar
    radar_orig : radar object
        the radar object containing the original field
    field_name: str
        name of the field to interpolate
    fill_value: float
        The fill value
    ang_tol : float
        angle tolerance to determine whether the radar origin sweep is the
        radar destination sweep

    Returns
    -------
    field_dest : dict
        interpolated field and metadata

    """
    if radar_dest.nsweeps != radar_orig.nsweeps:
        warn(
            "Number of sweeps in destination radar object different from "
            + "origin radar object. Orig: "
            + str(radar_orig.nsweeps)
            + " Dest : "
            + str(radar_dest.nsweeps)
        )

    if fill_value is None:
        fill_value = pyart.config.get_fillvalue()

    field_orig_data = radar_orig.fields[field_name]["data"].filled(
        fill_value=fill_value
    )
    field_dest = deepcopy(radar_orig.fields[field_name])
    field_dest["data"] = np.ma.masked_all(
        (radar_dest.nrays, radar_dest.ngates), dtype=field_orig_data.dtype
    )

    for sweep in range(radar_dest.nsweeps):
        sweep_start_dest = radar_dest.sweep_start_ray_index["data"][sweep]
        sweep_end_dest = radar_dest.sweep_end_ray_index["data"][sweep]
        fixed_angle = radar_dest.fixed_angle["data"][sweep]
        nrays_sweep = radar_dest.rays_per_sweep["data"][sweep]

        # look for nearest angle
        delta_ang = np.absolute(radar_orig.fixed_angle["data"] - fixed_angle)
        ind_sweep_orig = np.argmin(delta_ang)

        if delta_ang[ind_sweep_orig] > ang_tol:
            warn(
                "No fixed angle of origin radar object matches the fixed "
                + "angle of destination radar object for sweep nr "
                + str(sweep)
                + " with fixed angle "
                + str(fixed_angle)
                + "+/-"
                + str(ang_tol)
            )
            field_dest_sweep = np.ma.masked_all(
                (nrays_sweep, radar_dest.ngates), dtype=field_orig_data.dtype
            )
        else:
            sweep_start_orig = radar_orig.sweep_start_ray_index["data"][ind_sweep_orig]
            sweep_end_orig = radar_orig.sweep_end_ray_index["data"][ind_sweep_orig]

            if radar_dest.scan_type == "ppi":
                angle_old = np.sort(
                    radar_orig.azimuth["data"][sweep_start_orig : sweep_end_orig + 1]
                )
                ind_ang = np.argsort(
                    radar_orig.azimuth["data"][sweep_start_orig : sweep_end_orig + 1]
                )
                angle_new = radar_dest.azimuth["data"][
                    sweep_start_dest : sweep_end_dest + 1
                ]
            elif radar_dest.scan_type == "rhi":
                angle_old = np.sort(
                    radar_orig.elevation["data"][sweep_start_orig : sweep_end_orig + 1]
                )
                ind_ang = np.argsort(
                    radar_orig.elevation["data"][sweep_start_orig : sweep_end_orig + 1]
                )
                angle_new = radar_dest.elevation["data"][
                    sweep_start_dest : sweep_end_dest + 1
                ]
            else:
                warn(
                    f"Scan type of destination radar is {radar_dest.scan_type}."
                    f" Only ppi and rhi are supported"
                )

            # Reduce precision to avoid issues with floating numbers
            angle_new = np.around(angle_new, 6)
            angle_old = np.around(angle_old, 6)

            field_orig_sweep_data = field_orig_data[
                sweep_start_orig : sweep_end_orig + 1, :
            ]
            interpol_func = RegularGridInterpolator(
                (angle_old, radar_orig.range["data"]),
                field_orig_sweep_data[ind_ang],
                method="nearest",
                bounds_error=False,
                fill_value=fill_value,
            )

            # interpolate data to radar_dest grid
            angv, rngv = np.meshgrid(angle_new, radar_dest.range["data"], indexing="ij")

            field_dest_sweep = interpol_func((angv, rngv))
            field_dest_sweep = np.ma.masked_where(
                field_dest_sweep == fill_value, field_dest_sweep
            )

        field_dest["data"][sweep_start_dest : sweep_end_dest + 1, :] = field_dest_sweep

    return field_dest


def crop_grid(
    grid,
    lat_min=None,
    lat_max=None,
    lon_min=None,
    lon_max=None,
    alt_min=None,
    alt_max=None,
    ix_min=None,
    iy_min=None,
    iz_min=None,
    nx=None,
    ny=None,
    nz=None,
):
    """
    crops a grid object. The cropping can be done either specifying min and
    max lat, lon and altitude or by specifying the min lat, lon and altitude
    and the length in pixels of each side

    Parameters
    ----------
    grid : grid object
        the grid object to crop
    lat_min, lat_max, lon_min, lon_max : float
        the lat/lon limits of the object (deg)
    alt_min, alt_max : float
        the altitude limits of the object (m MSL)
    ix_min, iy_min, iz_min : int
        the pixel index where to start cropping
    nx, ny ,nz : int
        The number of pixels in each direction

    Returns
    -------
    grid_crop : grid object
        The cropped grid

    """
    grid_crop = deepcopy(grid)
    if (
        lat_min is None
        and lat_max is None
        and lon_min is None
        and lon_max is None
        and alt_min is None
        and alt_max is None
        and nx is None
        and ny is None
        and nz is None
        and ix_min is None
        and iy_min is None
        and iz_min is None
    ):
        return grid_crop

    if iy_min is not None:
        if grid.ny < iy_min:
            warn(
                "pixel index for latitude larger than number of pixels. "
                "Data will not be cropped"
            )
            iy_min = 0
    elif lat_min is not None:
        iz, iy, ix = np.where(grid.point_latitude["data"] >= lat_min)
        if iy.size == 0:
            warn(
                "Min latitude "
                + str(lat_min)
                + " outside of grid. "
                + "The data will not be cropped"
            )
            iy_min = 0
        else:
            iy_min = np.min(iy)
    else:
        iy_min = 0

    if ny is not None:
        iy_max = iy_min + ny
    elif lat_max is not None:
        iz, iy, ix = np.where(grid.point_latitude["data"] <= lat_max)
        if iy.size == 0:
            warn(
                "Max latitude "
                + str(lat_max)
                + " outside of grid. "
                + "The data will not be cropped"
            )
            iy_max = grid.ny
        else:
            iy_max = np.max(iy) + 1
    else:
        iy_max = grid.ny

    if ix_min is not None:
        if grid.nx < ix_min:
            warn(
                "pixel index for longitude larger than number of pixels. "
                "Data will not be cropped"
            )
            ix_min = 0
    elif lon_min is not None:
        iz, iy, ix = np.where(grid.point_longitude["data"] >= lon_min)
        if ix.size == 0:
            warn(
                "Min longitude "
                + str(lon_min)
                + " outside of grid. "
                + "The data will not be cropped"
            )
            ix_min = 0
        else:
            ix_min = np.min(ix)
    else:
        ix_min = 0

    if nx is not None:
        ix_max = ix_min + nx
    elif lon_max is not None:
        iz, iy, ix = np.where(grid.point_longitude["data"] <= lon_max)
        if ix.size == 0:
            warn(
                "Max longitude "
                + str(lon_max)
                + " outside of grid. "
                + "The data will not be cropped"
            )
            ix_max = grid.nx
        else:
            ix_max = np.max(ix) + 1
    else:
        ix_max = grid.nx

    if iz_min is not None:
        if grid.nz < iz_min:
            warn(
                "pixel index for altitude larger than number of pixels. "
                "Data will not be cropped"
            )
            iz_min = 0
    elif alt_min is not None:
        iz, iy, ix = np.where(grid.point_altitude["data"] >= alt_min)
        if iz.size == 0:
            warn(
                "Min altitude "
                + str(alt_min)
                + " outside of grid. "
                + "The data will not be cropped"
            )
            iz_min = 0
        else:
            iz_min = np.min(iz)
    else:
        iz_min = 0

    if nz is not None:
        iz_max = iz_min + nz
    elif alt_max is not None:
        iz, iy, ix = np.where(grid.point_altitude["data"] <= alt_max)
        if iz.size == 0:
            warn(
                "Max longitude "
                + str(lon_max)
                + " outside of grid. "
                + "The data will not be cropped"
            )
            iz_max = grid.nz
        else:
            iz_max = np.max(iz) + 1
    else:
        iz_max = grid.nz

    grid_crop.x["data"] = grid_crop.x["data"][ix_min:ix_max]
    grid_crop.y["data"] = grid_crop.y["data"][iy_min:iy_max]
    grid_crop.z["data"] = grid_crop.z["data"][iz_min:iz_max]
    grid_crop.nx = grid_crop.x["data"].size
    grid_crop.ny = grid_crop.y["data"].size
    grid_crop.nz = grid_crop.z["data"].size

    for field in grid_crop.fields:
        grid_crop.fields[field]["data"] = grid_crop.fields[field]["data"][
            iz_min:iz_max, iy_min:iy_max, ix_min:ix_max
        ]

    grid_crop.init_point_x_y_z()
    grid_crop.init_point_longitude_latitude()
    grid_crop.init_point_altitude()

    return grid_crop


def merge_radars(radar1, radar2):
    """
    Versatile function that merges two radar objects even
    if they do not have the same fields

    Parameters
    ----------
    radar1, radar2 : grid object
        the radar objects to merge

    Returns
    -------
    radar : radar object
        The merged radar

    """

    if set(radar1.fields.keys()) == set(radar2.fields.keys()):
        # Same fields
        radar = pyart.util.join_radar(radar1, radar2)
    elif any([f in radar1.fields.keys() for f in radar2.fields.keys()]):
        # At least one field of radar_aux already in radar
        # create object with common fields
        common_keys = [
            value for value in radar1.fields.keys() if value in radar2.fields.keys()
        ]

        radar1_aux = deepcopy(radar1)
        radar1_aux.fields = {}
        radar2_aux = deepcopy(radar2)
        radar2_aux.fields = {}

        for k in common_keys:
            radar1_aux.fields[k] = radar1.fields[k]
            radar2_aux.fields[k] = radar2.fields[k]

        radar_aux = pyart.util.join_radar(radar1_aux, radar2_aux)
        for key in radar1.fields:
            if key not in common_keys:
                radar1_aux = deepcopy(radar1)
                radar1_aux.fields = {key: radar1.fields[key]}
                radar_aux = add_field(radar_aux, radar1_aux)
        for key in radar2.fields:
            if key not in common_keys:
                radar2_aux = deepcopy(radar2)
                radar2_aux.fields = radar2.fields[key]
                radar_aux = add_field(radar_aux, radar2_aux)
        radar = radar_aux
    else:
        radar = add_field(radar1, radar2)

    return radar


def merge_grids(grid1, grid2):
    """
    Merges two grids

    Parameters
    ----------
    grid1, grid2 : grid object
        the grid objects to merge

    Returns
    -------
    grid : grid object
        The merged grid

    """
    # Check if the projections are the same. Change otherwise
    if grid1.projection != grid2.projection:
        grid2.projection = grid1.projection
        grid2.init_point_longitude_latitude()
        grid2.init_point_altitude()

    # Create new vectors of x, y and z
    x_equal = True
    y_equal = True
    z_equal = True

    x = pyart.config.get_metadata("x")
    if np.allclose(grid1.x["data"], grid2.x["data"]):
        x["data"] = grid1.x["data"]
    else:
        x["data"] = np.sort(np.unique(np.append(grid1.x["data"], grid2.x["data"])))
        x_equal = False

    y = pyart.config.get_metadata("y")
    if np.allclose(grid1.y["data"], grid2.y["data"]):
        y["data"] = grid1.y["data"]
    else:
        y["data"] = np.sort(np.unique(np.append(grid1.y["data"], grid2.y["data"])))
        y_equal = False

    z = pyart.config.get_metadata("z")
    if np.allclose(grid1.z["data"], grid2.z["data"]):
        z["data"] = grid1.z["data"]
    else:
        z["data"] = np.sort(np.unique(np.append(grid1.z["data"], grid2.z["data"])))
        z_equal = False

    nx = x["data"].size
    ny = y["data"].size
    nz = z["data"].size

    # if the grids are identical add the new fields directly
    if x_equal and y_equal and z_equal:
        grid = deepcopy(grid1)
        for field in grid2.fields.keys():
            if field in grid1.fields:
                warn("Field " + field + " already exists")
                continue
            grid.add_field(field, grid2.fields[field])

        return grid

    warn("non-identical grids, interpolating...")

    # create new grid object
    grid = pyart.core.grid.Grid(
        grid1.time,
        {},
        grid1.metadata,
        grid1.origin_latitude,
        grid1.origin_longitude,
        grid1.origin_altitude,
        x,
        y,
        z,
        projection=grid1.projection,
    )

    fields = np.unique(np.append(list(grid1.fields.keys()), list(grid2.fields.keys())))

    for field in fields:
        field1_data = None
        field2_data = None
        field_dict = pyart.config.get_metadata(field)
        field_dict["data"] = np.ma.masked_all((nz, ny, nx))

        if field in grid1.fields:
            field1_data = grid1.fields[field]["data"]
        if field in grid2.fields:
            field2_data = grid2.fields[field]["data"]

        # grids identical in at least two dimensions
        if x_equal and y_equal:
            np.shape(field1_data)
            np.shape(field2_data)
            for i, z_el in enumerate(z["data"]):
                if field1_data is not None:
                    ind_z = np.where(grid1.z["data"] == z_el)[0]
                    if ind_z.size > 0:
                        field_dict["data"][i, :, :] = field1_data[ind_z, :, :]
                if field2_data is not None:
                    ind_z = np.where(grid2.z["data"] == z_el)[0]
                    if ind_z.size > 0:
                        field_dict["data"][i, :, :] = field2_data[ind_z, :, :]
        elif x_equal and z_equal:
            for i, y_el in enumerate(y["data"]):
                if field1_data is not None:
                    ind_y = np.where(grid1.y["data"] == y_el)[0]
                    if ind_y.size > 0:
                        field_dict["data"][:, i, :] = field1_data[:, ind_y, :]
                if field2_data is not None:
                    ind_y = np.where(grid2.y["data"] == y_el)[0]
                    if ind_y.size > 0:
                        field_dict["data"][:, i, :] = field2_data[:, ind_y, :]
        elif y_equal and z_equal:
            for i, x_el in enumerate(x["data"]):
                if field1_data is not None:
                    ind_x = np.where(grid1.x["data"] == x_el)[0]
                    if ind_x.size > 0:
                        field_dict["data"][:, :, i] = field1_data[:, :, ind_x]
                if field2_data is not None:
                    ind_x = np.where(grid2.x["data"] == x_el)[0]
                    if ind_x.size > 0:
                        field_dict["data"][:, i, :] = field2_data[:, :, ind_x]
        else:
            # grids completely different
            for i, z_el in enumerate(z["data"]):
                for j, y_el in enumerate(y["data"]):
                    for k, x_el in enumerate(x["data"]):
                        if field1_data is not None:
                            ind_z = np.where(grid1.z["data"] == z_el)[0]
                            ind_y = np.where(grid1.y["data"] == y_el)[0]
                            ind_x = np.where(grid1.x["data"] == x_el)[0]
                            if ind_z.size > 0 and ind_y.size > 0 and ind_x.size > 0:
                                field_dict["data"][i, j, k] = field1_data[
                                    ind_z, ind_y, ind_x
                                ]
                        if field2_data is not None:
                            ind_z = np.where(grid2.z["data"] == z_el)[0]
                            ind_y = np.where(grid2.y["data"] == y_el)[0]
                            ind_x = np.where(grid2.x["data"] == x_el)[0]
                            if ind_z.size > 0 and ind_y.size > 0 and ind_x.size > 0:
                                field_dict["data"][i, j, k] = field2_data[
                                    ind_z, ind_y, ind_x
                                ]

        grid.add_field(field, field_dict)

    return grid
