# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 10:48:31 2016

@author: fvanden

Configuration file for mxpol pyart.core.Radar class. Some information may be
redundant because this file is a copy from the ProfileLab toolkit.

Functions to retrieve data from this file may be found in
 pyrad.io.read_data_mxpol under the utilities section

"""

# radar information

MCH_elev = [
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
NYQUIST_VEL = [
    8.3,
    9.6,
    8.3,
    12.4,
    11.0,
    12.4,
    13.8,
    12.4,
    13.8,
    16.5,
    16.5,
    16.5,
    20.6,
    20.6,
    20.6,
    20.6,
    20.6,
    20.6,
    20.6,
    20.6,
]

RADAR_INFO = {
    "coordinates": {
        "ALB": [47.284, 8.512],
        "DOL": [46.425, 6.099],
        "PPM": [46.371, 7.487],
        "MLE": [46.041, 8.833],
        "DX50": [46.8425, 6.9184],
        "MXPOL": [46.8133, 6.9428],
    },
    "altitude": {
        "ALB": 938,
        "DOL": 1682,
        "PPM": 2937,
        "MLE": 1626,
        "DX50": 451,
        "MXPOL": 489,
    },
    "searchkey": {
        "ALB": "PHA*hdf*",
        "DOL": "PHD*hdf*",
        "PPM": "PHP*hdf*",
        "MLE": "PHL*hdf*",
        "DX50": None,
        "MXPOL": None,
    },
    "radarID": {
        "ALB": "ALB",
        "A": "ALB",
        "DOL": "DOL",
        "D": "DOL",
        "PPM": "PPM",
        "P": "PPM",
        "MLE": "MLE",
        "M": "MLE",
        "DX50": "DX50",
        "MXPOL": "MXPOL",
    },
    "dbbeam": {
        "ALB": 1.0,
        "DOL": 1.0,
        "PPM": 1.0,
        "MLE": 1.0,
        "MXPOL": 1.4,
        "DX50": 1.27,
    },
    "elevations": {
        "ALB": MCH_elev,
        "DOL": MCH_elev,
        "PPM": MCH_elev,
        "MLE": MCH_elev,
        "DX50": None,
        "MXPOL": None,
    },
}

MY_METADATA = {
    "nyq_vel": NYQUIST_VEL,
    # Metadata for instrument tables
    "Radar_info": {
        "searchkey": None,
        "coordinates": None,
        "altitude": None,
        "dbbeam": None,
        "radarID": None,
    },
    "Polvar": {
        "units": None,
        "standard_name": None,
        "short_name": None,
        "long_name": None,
        "valid_min": None,
        "valid_max": None,
        "plot_interval": None,
    },
}

# Metadata for polarimetric short and long names
MY_POLARNAMES = {
    "Zh": ["reflectivity", "reflectivity", "dBZ", 0.0, 55.0, 1.0],
    "Zdr": [
        "differential_reflectivity",
        "Differential reflectivity",
        "dB",
        -1.0,
        5.0,
        0.1,
    ],
    "Kdp": [
        "specific_differential_phase",
        "Specific differential phase",
        "deg/km",
        -2.0,
        7.0,
        0.1,
    ],
    "Phidp": [
        "corrected_differential_phase",
        "Differential phase",
        "deg",
        0.0,
        150.0,
        1.0,
    ],
    "Psidp": [
        "uncorrected_differential_phase",
        "Total differential phase",
        "deg",
        0.0,
        150.0,
        1.0,
    ],
    "Rhohv": [
        "uncorrected_cross_correlation_ratio",
        "Copolar correlation coefficient",
        "-",
        0.57,
        1.0,
        0.05,
    ],
    "ZhCorr": [
        "corrected_unfiltered_reflectivity",
        "Attenuation corrected reflectivity",
        "dBZ",
        0.0,
        55.0,
        1.0,
    ],
    "ZdrCorr": [
        "corrected_differential_reflectivity",
        "Attenuation corrected differential reflectivity",
        "dB",
        0.0,
        3.0,
        0.1,
    ],
    "RVel": ["velocity", "Mean doppler velocity", "m/s", -15.0, 15.0, 0.5],
    "Rvel": ["velocity", "Mean doppler velocity", "m/s", -15.0, 15.0, 0.5],
    "Sw": ["spectrum_width", "Spectral Width", "m2/s2", 0.0, 3.0, 0.1],
    "Zv": ["reflectivity_vv", "Vertical reflectivity", "dBZ", 0.0, 45.0, 1.0],
    "Clut": ["clutter", "Output clutter algorithm", "-", 0.0, 100.0, 10.0],
    "corrected_Z": [
        "corrected_reflectivity",
        "Clutter filtered reflectivity",
        "dBZ",
        0.0,
        55.0,
        1.0,
    ],
    "SNRh": [
        "signal_noise_ratio_h",
        "Signal to noise ratio at hor. pol",
        "-",
        0.0,
        50.0,
        0.5,
    ],
    "SNRv": [
        "signal_noise_ratio_v",
        "Signal to noise ratio at vert. pol",
        "-",
        0.0,
        50.0,
        0.5,
    ],
}
