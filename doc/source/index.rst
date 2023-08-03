=====================================
The Python ARM Radar Toolkit - Py-ART
=====================================

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: User Guide

   userguide/index.rst

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Reference Guide

   overview/index.rst

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Product overview

   overview/index.rst


About Pyrad
===============

Pyrad is a real-time data processing framework developed by MeteoSwiss and MeteoFrance. The framework is aimed at processing and visualizing polar data from individual weather radars as well as composite Cartesian products both off-line and in real time. It is written in the Python language. The framework is version controlled and automatic documentation is generated based on doc-strings. It is capable of ingesting data from all the weather radars in Switzerland, namely the operational MeteoSwiss C-band rad4alp radar network, the MeteoSwiss X-band DX50 radar and the EPFL MXPol radar and radar data in the OPERA file format. Additionally, it can ingest C/FRadial and NEXRAD level 2 files.

The processing flow is controlled by 3 simple configuration files. Multiple levels of processing can be performed. At each level new datasets (e.g. attenuation corrected reflectivity) are created which can be stored in a file and/or used in the next processing level (for example, creating a rainfall rate dataset from the corrected reflectivity). Multiple products can be generated from each dataset (e.g. PPI, RHI images, histograms, etc.). In the off-line mode, data from multiple radars can be ingested in order to obtain products such as the inter-comparison of reflectivity values at co-located range gates.

The framework is able to ingest polarimetric and Doppler radar moments as well as auxiliary data such as numerical weather prediction parameters (e.g. temperature, wind speed, etc.), DEM-based visibility and data used in the generation of the products such as rain gauge measurements, disdrometer measurements, solar flux, etc. It can as well work with I/Q data, spectral data and Cartesian data.

The signal processing and part of the data visualization is performed by a MeteoSwiss developed version of the Py-ART radar toolkit which contains enhanced features. MeteoSwiss regularly contributes back to the main Py-ART branch once a new functionality has been thoroughly tested and it is considered of interest for the broad weather radar community.

The processing framework has multiple and expanding capabilities, include various forms of echo classification and filtering, differential phase and specific differential phase estimation, attenuation correction, data quality monitoring, multiple rainfall rate algorithms, etc. In addition time series of data in points, regions or trajectories of interest can be extracted and comparisons can be performed with other sensors. This is particularly useful when performing measurement campaigns where remote sensing retrievals are validated with in-situ airplane or ground-based measurements.


Development
=============

We welcome contributions, suggestions of developments and bug reports.

Suggestions of developments and bug reports should use the Issues page of the `github repository <https://github.com/meteoswiss-mdr/pyrad/issues
>`_.


The process to contribute by partners external to MeteoSwiss is described in the `user manual <https://github.com/MeteoSwiss/pyrad/blob/master/additional_doc/pyrad_user_manual.pdf>`_.


Citing Pyrad
=============

The core of Pyrad is based on Py-ART. Py-ART was originally developed in the context of the ARM Research Facility. If you use Pyrad for your work, please cite BOTH Py-ART and Pyrad papers in your paper:

J.J. Helmus, S.M. Collis, (2016). The Python ARM Radar Toolkit (Py-ART), a Library for Working with Weather Radar Data in the Python Programming Language. Journal of Open Research Software. 4(1), p.e25. DOI: http://doi.org/10.5334/jors.119

J. Figueras i Ventura, M. Lainer, Z. Schauwecker, J. Grazioli, U. Germann, (2020). Pyrad: A Real-Time Weather Radar Data Processing Framework Based on Py-ART. Journal of Open Research Software, 8(1), p.28. DOI: http://doi.org/10.5334/jors.330

Disclaimer
=============

The software is still in a development stage. Please let us know if you would like to test it.

MeteoSwiss cannot be held responsible for errors in the code or problems that could arise from its use.