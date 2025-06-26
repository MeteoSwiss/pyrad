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
    "configpath": None,
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
