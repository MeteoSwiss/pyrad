# Note that pull requests should target the dev branch!
**Master branch**

[![Test pyrad](https://github.com/MeteoSwiss/pyrad/actions/workflows/pyrad_tests_base.yml/badge.svg)](https://github.com/MeteoSwiss/pyrad/actions/workflows/pyrad_tests.yml)
[![Test pyrad](https://github.com/MeteoSwiss/pyrad/actions/workflows/pyrad_tests_mch.yml/badge.svg)](https://github.com/MeteoSwiss/pyrad/actions/workflows/pyrad_tests_arm.yml)
[![build-deploy-site](https://github.com/MeteoSwiss/pyrad/actions/workflows/build_docs.yml/badge.svg)](https://github.com/MeteoSwiss/pyrad/actions/workflows/build_docs.yml)

**Dev branch**

[![Test pyrad base dev](https://github.com/MeteoSwiss/pyrad/actions/workflows/pyrad_tests_base_dev.yml/badge.svg?branch=dev)](https://github.com/MeteoSwiss/pyrad/actions/workflows/pyrad_tests_base_dev.yml)
[![Test pyrad mch dev](https://github.com/MeteoSwiss/pyrad/actions/workflows/pyrad_tests_mch_dev.yml/badge.svg?branch=dev)](https://github.com/MeteoSwiss/pyrad/actions/workflows/pyrad_tests_mch_dev.yml)


---

[![Ref doc](https://img.shields.io/badge/docs-users-4088b8.svg)](https://meteoswiss.github.io/pyrad/)


![version](https://img.shields.io/badge/python-3.8-blue.svg)
![version](https://img.shields.io/badge/python-3.9-blue.svg)
![version](https://img.shields.io/badge/python-3.10-blue.svg)
![version](https://img.shields.io/badge/python-3.11-blue.svg)

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pyrad_mch/badges/version.svg)](https://anaconda.org/conda-forge/pyrad_mch)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pyrad_mch/badges/latest_release_date.svg)](https://anaconda.org/conda-forge/pyrad_mch)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pyrad_mch/badges/downloads.svg)](https://anaconda.org/conda-forge/pyrad_mch)

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pyrad_mch/badges/platforms.svg)](https://anaconda.org/conda-forge/pyrad_mch)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pyrad_mch/badges/license.svg)](https://anaconda.org/conda-forge/pyrad_mch)

# pyrad
Python Radar Data Processing

# What is Pyrad?
Pyrad is a real-time data processing framework developed by MeteoSwiss and MeteoFrance. The framework is
aimed at processing and visualizing polar data from individual weather radars as well as composite
Cartesian products both off-line and in real time. It is written in the Python language. The
framework is version controlled and automatic documentation is generated based on doc-strings.
It is capable of ingesting data from all the weather radars in Switzerland, namely the operational
MeteoSwiss C-band rad4alp radar network, the MeteoSwiss X-band DX50 radar and the EPFL MXPol radar
and radar data in the OPERA file format. Additionally, it can ingest C/FRadial and NEXRAD level 2 files.

The processing flow is controlled by 3 simple configuration files. Multiple levels of processing can
be performed. At each level new datasets (e.g. attenuation corrected reflectivity) are created which
can be stored in a file and/or used in the next processing level (for example, creating a rainfall rate
dataset from the corrected reflectivity). Multiple products can be generated from each dataset (e.g.
PPI, RHI images, histograms, etc.). In the off-line mode, data from multiple radars can be ingested
in order to obtain products such as the inter-comparison of reflectivity values at co-located range
gates.

The framework is able to ingest polarimetric and Doppler radar moments as well as auxiliary data
such as numerical weather prediction parameters (e.g. temperature, wind speed, etc.), DEM-based
visibility and data used in the generation of the products such as rain gauge measurements,
disdrometer measurements, solar flux, etc. It can as well work with I/Q data, spectral data and Cartesian
data.

The signal processing and part of the data visualization is performed by a [MeteoSwiss developed version of the Py-ART radar toolkit](https://github.com/meteoswiss-mdr/pyart) which contains enhanced features. MeteoSwiss regularly contributes back to the [main Py-ART branch](https://github.com/ARM-DOE/pyart) once a new functionality has been thoroughly tested and it is considered of interest for the broad weather radar community.

The processing framework has multiple and expanding capabilities, include various forms of echo classification and
filtering, differential phase and specific differential phase estimation, attenuation correction, data
quality monitoring, multiple rainfall rate algorithms, etc. In addition time series of data in points,
regions or trajectories of interest can be extracted and comparisons can be performed with other
sensors. This is particularly useful when performing measurement campaigns where remote
sensing retrievals are validated with in-situ airplane or ground-based measurements.

# Cloning from github
Make sure to also get the submodules by running

```
    git clone --recursive https://github.com/MeteoSwiss/pyrad.git
```

# Installation with conda
To install from the conda repositories simply run

```
    conda install -c conda-forge pyart_mch
    conda install -c conda-forge pyrad_mch
```
Note that you can also install *arm_pyart* instead of *pyart_mch*, which will use the [official Py-ART](https://github.com/ARM-DOE/pyart), but some functionalities of pyrad will be missing!

# Installation
To install Pyrad and its submodules please have a look at the [Pyrad user manual (pdf)](./additional_doc/pyrad_user_manual.pdf).

# Use
Before using it have a look at the [cookbook (pdf)](./additional_doc/pyrad-framework-cookbook/DataProcessing.pdf).

For details on the implemented functions and a list of all pyrad features please check the [pyrad library reference for users](https://meteoswiss.github.io/pyrad/). 

Example configuration files can be found in the repository directory [pyrad/config/processing/](./config/processing) and in the dedicated [examples repository](https://github.com/MeteoSwiss/pyrad-examples). 

To use Pyrad for data quality monitoring check the report [pyrad_monitoring_fvj.pdf](./additional_doc/pyrad_monitoring_fvj.pdf).

# Learning pyrad

Please also check the course that the pyrad team gave at the ASEAN-WMO 2024 workshop: https://github.com/openradar/asean2024-pyrad-course.
It includes a list of 10 configuration files and their descriptions, that encompasses a large number of possible use cases.

# Newsletter

If you would like to be informed about the addition of major features in Pyrad as well as the release of new Pyrad versions, you can subscribe to our [mailing list](https://github.us1.list-manage.com/subscribe?u=2e8561343b026dacadbc70c92&id=9a43ddb927), where we will periodically publish a newsletter. Note that we do not sell, communicate or divulgate your email address to anyone, and you can unsubscribe at any time via a link provided in every newsletter.

# Development
We welcome contributions, suggestions of developments and bug reports.

Suggestions of developments and bug reports should use the [Issues page of the github repository](https://github.com/meteoswiss-mdr/pyrad/issues).

The process to contribute by partners external to MeteoSwiss is described in the [user manual](./additional_doc/pyrad_user_manual.pdf).

# Citation
The core of Pyrad is based on Py-ART. Py-ART was originally developed in the context of the [ARM Research Facility](https://www.arm.gov/). If you use Pyrad for your work, please cite BOTH Py-ART and Pyrad papers in your paper:

J.J. Helmus, S.M. Collis, (2016). The Python ARM Radar Toolkit (Py-ART), a Library for Working with Weather Radar Data in the Python Programming Language. Journal of Open Research Software. 4(1), p.e25. DOI: http://doi.org/10.5334/jors.119

J. Figueras i Ventura, M. Lainer, Z. Schauwecker, J. Grazioli, U. Germann, (2020). Pyrad: A Real-Time Weather Radar Data Processing Framework Based on Py-ART. Journal of Open Research Software, 8(1), p.28. DOI: http://doi.org/10.5334/jors.330 

# Disclaimer
The software is still in a development stage. Please let us know if you would like to test it.

MeteoSwiss cannot be held responsible for errors in the code or problems that could arise from its use.

