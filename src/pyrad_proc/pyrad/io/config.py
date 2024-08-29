"""
pyrad.io.config
===============

Functions for reading pyrad config files

.. autosummary::
    :toctree: generated/

    read_config
    get_num_elements
    string_to_datatype
    get_array
    get_struct
    get_array_type
    init_array

"""

import os
import numpy as np
from warnings import warn

_defaults_loc = {
    "ppiImageConfig": {
        "xsize": 10,
        "ysize": 10,
        "dpi": 72,
        "xmin": -50,
        "xmax": 50,
        "ymin": -50,
        "ymax": 50,
    },
    "ppiMapImageConfig": {"xsize": 10, "ysize": 10, "dpi": 72},
    "rhiImageConfig": {
        "xsize": 10,
        "ysize": 5,
        "dpi": 72,
        "xmin": 0,
        "xmax": 50,
        "ymin": 0,
        "ymax": 12,
    },
    "xsecImageConfig": {
        "xsize": 10,
        "ysize": 5,
        "dpi": 72,
        "xmin": 0,
        "xmax": 50,
        "ymin": 0,
        "ymax": 12,
    },
    "gridMapImageConfig": {"xsize": 10, "ysize": 10, "dpi": 72},
    "sunhitsImageConfig": {
        "xsize": 10,
        "ysize": 5,
        "dpi": 72,
        "azmin": -2.0,
        "azmax": 2.0,
        "elmin": -2.0,
        "elmax": 2.0,
        "azres": 0.1,
        "elres": 0.1,
    },
    "spectraImageConfig": {
        "xsize": 10,
        "ysize": 5,
        "dpi": 72,
        "ymin": None,
        "ymax": None,
        "velmin": None,
        "velmax": None,
    },
    "ScanPeriod": 5,
    "IconRunFreq": 3,
    "IconForecasted": 7,
    "RadarName": None,
    "RadarRes": None,
    "NumRadars": 1,
    "TimeTol": 3600.0,
    "ScanList": None,
}

_defaults_main = {
    "lastStateFile": None,
    "datapath": None,
    "satpath": None,
    "iconpath": None,
    "psrpath": None,
    "iqpath": None,
    "colocgatespath": None,
    "excessgatespath": None,
    "dempath": None,
    "smnpath": None,
    "disdropath": None,
    "solarfluxpath": None,
    "selfconsistencypath": None,
    "loadbasepath": None,
    "loadname": None,
    "gecsxbasepath": None,
    "gecsxname": None,
    "metranet_read_lib": "C",
}

DEFAULT_CONFIG = {"loc": _defaults_loc, "main": _defaults_main}


def read_config(fname, cfg=None, defaults=None):
    """
    Read a pyrad config file.

    Parameters
    ----------
    fname : str
        Name of the configuration file to read.

    cfg : dict of dicts, optional
        dictionary of dictionaries containing configuration parameters where
        the new parameters will be placed

    defaults: dict of dicts, optional
        dictionary of dictionaries containing default values. If a key is
        contained in the defaults dict but not in the configuration file,
        the value from the defaults dict will be assigned

    Returns
    -------
    cfg : dict of dicts
        dictionary of dictionaries containing the configuration parameters

    """

    # check if the file can be read
    try:
        cfgfile = open(fname, "r", encoding="utf-8", errors="ignore")
    except Exception:
        raise Exception("ERROR: Could not find|open config file '" + fname + "'")

    # if config dictionary does not exist yet create it
    if cfg is None:
        cfg = dict()

    # if default does not exist, create it
    if defaults is None:
        defaults = dict()

    # read file contents
    fileend = 0
    while fileend == 0:
        # remove leading and trailing whitespace
        line = os.path.expandvars(cfgfile.readline())
        if line:
            line = line.strip()

            # ignore white lines
            if not line:
                continue

            # ignore comments
            if line.startswith("#"):
                continue

            line = line.partition("#")[0]  # Remove comments
            line = line.strip()

            try:
                vals = line.split()

                fieldname = vals[0]
                nvals = len(vals)

                if nvals < 3:
                    raise Exception(
                        "FILE FORMAT ERROR: file: "
                        + fname
                        + ", variable: "
                        + fieldname
                        + ": Wrong number of elements!"
                    )

                valtype = vals[1]
                valuestr = vals[2:nvals]
                nel, isstruct = get_num_elements(valtype, valuestr)

                if nel > 0:
                    if isstruct:
                        pos = cfgfile.tell()
                        fieldvalue, newpos = get_struct(cfgfile, pos, nel, fname)
                        cfgfile.seek(newpos)
                    else:
                        pos = cfgfile.tell()
                        fieldvalue, newpos = get_array(cfgfile, pos, nel, valtype)
                        cfgfile.seek(newpos)
                else:
                    fieldvalue = string_to_datatype(valtype, valuestr)

                cfg.update({fieldname: fieldvalue})
            except BaseException:
                if line != "":
                    warn("Parsing failed after line {:s}".format(line))
                raise
        else:
            fileend = 1

    cfgfile.close()
    # Verify that all keys in default are in newly created config
    cfg = merge_dicts(cfg, defaults)
    return cfg


def get_num_elements(dtype, nelstr):
    """
    Checks if data type is an array or a structure.

    Parameters
    ----------
    dtype : str
        data type specifier

    nelstr : str
        number of elements

    Returns
    -------
    nel : int
        number of elements if type is *ARR or STRUCT. 0 otherwise

    isstruct : bool
        true if the type is STRUCT

    """

    uptype = dtype.upper()
    isstruct = False
    nel = 0
    narr = uptype.count("ARR")
    if uptype == "STRUCT":
        isstruct = True
    if (narr > 0) or (isstruct is True):
        nel = int(nelstr[0])

    return nel, isstruct


def string_to_datatype(dtype, strval):
    """
    Converts a string containing a value into its Python value

    Parameters
    ----------
    dtype : str
        data type specifier

    strval : str
        string value

    Returns
    -------
    val : scalar
        value contained in the string

    """

    uptype = dtype.upper()

    if uptype == "BYTE":
        return int(strval[0])
    elif uptype == "BOOL":
        return bool(strval[0])
    elif uptype == "INT":
        return int(strval[0])
    elif uptype == "LONG":
        return int(strval[0])
    elif uptype == "HEX":
        return int(strval[0])
    elif uptype == "EXP":
        return float(strval[0])
    elif uptype == "FLOAT":
        return float(strval[0])
    elif uptype == "DOUBLE":
        return float(strval[0])
    elif uptype == "STRING":
        # Replace all environement variables
        return os.path.expandvars(str(strval[0]))
    else:
        raise Exception("ERROR: Unexpected data type " + uptype)


def get_array(cfgfile, pos, nel, valtype):
    """
    reads an array in a config file

    Parameters
    ----------
    cfgfile : file object
        config file

    pos : int
        position in file object

    nel : int
        number of elements of the ray

    valtype : str
        type of array

    Returns
    -------
    arr : array
        array values

    newpos : int
        new position in file object

    """

    arr_type = get_array_type(valtype)
    if arr_type == "STRING":
        arr = []
    else:
        arr = init_array(nel, arr_type)

    newpos = pos
    for i in range(nel):
        pos = cfgfile.seek(newpos)
        line = cfgfile.readline()

        # remove leading and trailing whitespace
        line = line.strip()

        line = line.partition("#")[0]  # Remove comments
        line = line.strip()

        vals = line.split()

        value = string_to_datatype(arr_type, vals)
        if arr_type == "STRING":
            arr.append(value)
        else:
            arr[i] = value
        newpos = cfgfile.tell()

    return arr, newpos


def get_struct(cfgfile, pos, nels, fname):
    """
    reads an struct in a config file

    Parameters
    ----------
    cfgfile : file object
        config file
    pos : int
        position in file object
    nel : int
        number of elements of the ray
    fname : str
        config file name

    Returns
    -------
    struct : dict
        dictionary of struct values

    newpos : int
        new position in file object

    """

    struct = dict()
    newpos = pos
    for i in range(nels):
        pos = cfgfile.seek(newpos)
        line = os.path.expandvars(cfgfile.readline())

        # remove leading and trailing whitespace
        line = line.strip()

        line = line.partition("#")[0]  # Remove comments
        line = line.strip()

        try:
            vals = line.split()

            sfieldname = vals[0]
            nvals = len(vals)

            if nvals < 3:
                raise Exception(
                    "FILE FORMAT ERROR: file: "
                    + fname
                    + ", struct variable: "
                    + sfieldname
                    + ": Wrong number of elements!"
                )

            svaltype = vals[1]
            svaluestr = vals[2:nvals]
            nel, isstruct = get_num_elements(svaltype, svaluestr)

            if nel > 0:
                if isstruct:
                    pos = cfgfile.tell()
                    sfieldvalue, newpos = get_struct(cfgfile, pos, nel, fname)
                    cfgfile.seek(newpos)
                else:
                    pos = cfgfile.tell()
                    sfieldvalue, newpos = get_array(cfgfile, pos, nel, svaltype)
                    cfgfile.seek(newpos)
            else:
                newpos = cfgfile.tell()
                sfieldvalue = string_to_datatype(svaltype, svaluestr)

            struct.update({sfieldname: sfieldvalue})
        except BaseException:
            if line != "":
                warn("Parsing failed after line {:s}".format(line))
            raise

    return struct, newpos


def get_array_type(dtype):
    """
    Determines Python array type from the config file array type

    Parameters
    ----------
    dtype : str
        config file data type

    Returns
    -------
    pytype : str
        Python array type

    """

    uptype = dtype.upper()

    if uptype == "BYTARR":
        return "BYTE"
    elif uptype == "INTARR":
        return "INT"
    elif uptype == "LONARR":
        return "LON"
    elif uptype == "HEXARR":
        return "HEX"
    elif uptype == "EXPARR":
        return "EXP"
    elif uptype == "FLTARR":
        return "FLOAT"
    elif uptype == "DBLARR":
        return "DOUBLE"
    elif uptype == "STRARR":
        return "STRING"
    else:
        raise Exception("ERROR: Unexpected data type " + uptype)


def init_array(nel, dtype):
    """
    Initializes a Python array

    Parameters
    ----------
    nel : int
        number of elements in the array
    dtype : str
        config file data type

    Returns
    -------
    pyarr : array
        Python array

    """

    uptype = dtype.upper()

    if uptype == "BYTE":
        return np.empty(nel, dtype="byte")
    elif uptype == "INT":
        return np.empty(nel, dtype="int32")
    elif uptype == "LONG":
        return np.empty(nel, dtype="int64")
    elif uptype == "HEX":
        return np.empty(nel, dtype="H")
    elif uptype == "EXP":
        return np.empty(nel, dtype="float64")
    elif uptype == "FLOAT":
        return np.empty(nel, dtype="float64")
    elif uptype == "DOUBLE":
        return np.empty(nel, dtype="double")
    elif uptype == "STRING":
        return []
    else:
        raise Exception("ERROR: Unexpected array data type " + uptype)


def merge_dicts(ref, defaults):
    """
    Merge two nested dictionaries recursively.

    Parameters:
    -----------
    ref : dict
        The base dictionary to be merged into.
    defaults : dict
        The dictionary whose keys and values are merged into dictionary ref.

    Returns:
    --------
    ref: dict
        The base dictionary after addition of missing keys from defaults
    """
    for key, value in defaults.items():
        if key in ref:
            if isinstance(ref[key], dict) and isinstance(value, dict):
                merge_dicts(ref[key], value)
        else:
            ref[key] = value
    return ref
