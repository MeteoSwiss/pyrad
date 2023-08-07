.. _API:

####################
Features
####################

:Release: |version|
:Date: |today|

This guide provides an overview of all available Pyrad features. 
It gives as list of all possible processes and products, as well as all variable names that are currently supported.

Processes are separated by the nature of their output datasets (e.g. VOL, GRID, SPECTRA,...).
For every type of output dataset, a certain number of products can be generated. Hence similarly, all avaialble products are separated by the nature of the dataset to which they apply.

To take an example, the ATTENUATION process returns a dataset of the VOL type. For the VOL type, many products are available, for example CAPPI_IMAGE or PPI_MAP, but not SURFACE_RAW (for example) which is a product that operates on GRID datasets only.

.. toctree::
   :maxdepth: 1

   list_process
   list_products
   list_variables
