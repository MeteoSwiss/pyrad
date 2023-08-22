===============
What is Pyrad?
===============


Pyrad is a real-time data processing framework developed by MeteoSwiss. The framework is aimed at 
processing and visualizing data from weather radars both off-line and in real time. There is also limited 
support to visualize cloud radar, lidar and satellite data. It is written in the Python language. The 
framework is version controlled and automatic documentation is generated based on doc-strings. It is 
capable of ingesting data from all the weather radars in Switzerland, namely the operational 
MeteoSwiss C-band radar network, the MeteoSwiss X-band METEOR 50DX radar and the EPFL 
MXPol radar. It can also read ODIM complying files, CFRadial and NEXRAD.

The processing flow is controlled by 3 simple configuration files. Multiple levels of processing can be 
performed. At each level new datasets (e.g. attenuation corrected reflectivity) are created which can 
be stored in a file and/or used in the next processing level (e.g. creating a rainfall rate dataset from the 
corrected reflectivity). Multiple products can be generated from each dataset (i.e PPI, RHI images, 
histograms, etc.). In the off-line mode, data from multiple radars can be ingested in order to obtain 
products such as the inter-comparison of reflectivity values at co-located range gates.
The framework is able to ingest polarimetric and Doppler radar moments as well as auxiliary data such 
as numerical weather prediction parameters (e.g. temperature, wind speed, etc.), DEM-based visibility 
and data used in the generation of the products such as rain gauge measurements, disdrometer 
measurements, solar flux, etc.

The signal processing and part of the data visualization is performed by a MeteoSwiss developed 
version of the Py-ART radar toolkit [1] which contains enhanced features. MeteoSwiss regularly 
contributes back to the main Py-ART branch once a new functionality has been thoroughly tested and 
it is considered of interest for the broad weather radar community.

The capabilities of the processing framework include various forms of echo classification and filtering, 
differential phase and specific differential phase estimation, attenuation correction, data quality 
monitoring, multiple rainfall rate algorithms, etc. In addition, time series of data in points, regions or 
trajectories of interest can be extracted and comparisons can be performed with other sensors. This is 
particularly useful when performing measurement campaigns where remote sensing retrievals are 
validated with in-situ airplane or ground-based measurements. The capabilities of the framework are 
expanded on an almost daily basis.

A certain degree of parallelization has been included. The user may choose to parallelize the 
generation of datasets of the same processing level, the generation of all the products of each dataset 
or both.

Radar volumetric data can be stored in C/F radial format or in ODIM format. Other data is typically 
stored as csv files. Plots can be output in any format accepted by Matplotlib.
There exists several software packages that are Py-ART dependent. At the moment Pyrad wraps the 
software package PyTDA for turbulence computation