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
import re
import numpy as np


def read_config(fname, cfg=None):
    """
    Read a pyrad config file.

    Parameters
    ----------
    fname : str
        Name of the configuration file to read.

    cfg : dict of dicts, optional
        dictionary of dictionaries containing configuration parameters where
        the new parameters will be placed

    Returns
    -------
    cfg : dict of dicts
        dictionary of dictionaries containing the configuration parameters

    """

    # check if the file can be read
    try:
        cfgfile = open(fname, "r", encoding='utf-8', errors='ignore')
    except:
        raise Exception("ERROR: Could not find|open config file '"+fname+"'")

    # if config dictionary does not exist yet create it
    if cfg is None:
        cfg = dict()

    # read file contents
    fileend = 0
    while fileend == 0:
        # remove leading and trailing whitespace
        line = cfgfile.readline()
        if line:
            line = line.strip()

            # ignore white lines
            if not line:
                continue

            # ignore comments
            if line.startswith('#'):
                continue

            line = line.partition('#')[0]  # Remove comments
            line = line.strip()

            vals = line.split()

            fieldname = vals[0]
            nvals = len(vals)

            if nvals < 3:
                raise Exception(
                    "FILE FORMAT ERROR: file: " + fname +
                    ", variable: " + fieldname +
                    ": Wrong number of elements!")

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
        else:
            fileend = 1

    cfgfile.close()
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
    narr = uptype.count('ARR')
    if uptype == 'STRUCT':
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

    if uptype == 'BYTE':
        return int(strval[0])
    elif uptype == 'BOOL':
        return  bool(strval[0])
    elif uptype == 'INT':
        return int(strval[0])
    elif uptype == 'LONG':
        return int(strval[0])
    elif uptype == 'HEX':
        return int(strval[0])
    elif uptype == 'EXP':
        return float(strval[0])
    elif uptype == 'FLOAT':
        return float(strval[0])
    elif uptype == 'DOUBLE':
        return float(strval[0])
    elif uptype == 'STRING':
        # Replace all environement variables
        return os.path.expandvars(str(strval[0]))
    else:
        raise Exception("ERROR: Unexpected data type "+uptype)


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
    if arr_type == 'STRING':
        arr = []
    else:
        arr = init_array(nel, arr_type)

    newpos = pos
    for i in range(nel):
        pos = cfgfile.seek(newpos)
        line = cfgfile.readline()

        # remove leading and trailing whitespace
        line = line.strip()

        line = line.partition('#')[0]  # Remove comments
        line = line.strip()

        vals = line.split()

        value = string_to_datatype(arr_type, vals)
        if arr_type == 'STRING':
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
        line = cfgfile.readline()

        # remove leading and trailing whitespace
        line = line.strip()

        line = line.partition('#')[0]  # Remove comments
        line = line.strip()

        vals = line.split()

        sfieldname = vals[0]
        nvals = len(vals)

        if nvals < 3:
            raise Exception(
                "FILE FORMAT ERROR: file: " + fname +
                ", struct variable: " + sfieldname +
                ": Wrong number of elements!")

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

    if uptype == 'BYTARR':
        return 'BYTE'
    elif uptype == 'INTARR':
        return 'INT'
    elif uptype == 'LONARR':
        return 'LON'
    elif uptype == 'HEXARR':
        return 'HEX'
    elif uptype == 'EXPARR':
        return 'EXP'
    elif uptype == 'FLTARR':
        return 'FLOAT'
    elif uptype == 'DBLARR':
        return 'DOUBLE'
    elif uptype == 'STRARR':
        return 'STRING'
    else:
        raise Exception("ERROR: Unexpected data type "+uptype)


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

    if uptype == 'BYTE':
        return np.empty(nel, dtype='byte')
    elif uptype == 'INT':
        return np.empty(nel, dtype='int32')
    elif uptype == 'LONG':
        return np.empty(nel, dtype='int64')
    elif uptype == 'HEX':
        return np.empty(nel, dtype='H')
    elif uptype == 'EXP':
        return np.empty(nel, dtype='float64')
    elif uptype == 'FLOAT':
        return np.empty(nel, dtype='float64')
    elif uptype == 'DOUBLE':
        return np.empty(nel, dtype='double')
    elif uptype == 'STRING':
        return []
    else:
        raise Exception("ERROR: Unexpected array data type "+uptype)
