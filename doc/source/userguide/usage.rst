=============
Using Pyrad
=============

Configuration files
-----------------------

Pyrad uses 3 different configuration files which are typically stored in the folder::
*pyrad/config/processing/*

The first file specifies the input data, output data and configuration files packages, the second 
specifies radar related parameters (radar name, scan name and frequency, etc.) and the general 
configuration of the various image output, the last file specifies the datasets and products to be 
produced.

The easiest way to start is to copy one of the available config files and modify it according to your 
needs. 

Examples of pyrad config files are available in the repository https://github.com/MeteoSwiss/pyrad-examples/.


Running the programs
---------------------------

To run the programs first you need to activate the conda pyrad environment::

        conda activate pyrad

A number of script are available as executables. You can use them by simply typing::
    
        [name_of_the_program] [arguments]

At the moment there are three main programs:

main_process_data.py
    will process (and optionally post-process) data from a starting point in time to an ending point in time.

main_process_data_period.py
        will process (and optionally post-process) data over several days starting the processing at a given starting and ending time (default 00:00:00 for start and 23:59:59 for 
        the end).

main_process_gecsx.py
    will run the GECSX algorithm (visibility and clutter estimation from a DEM map) 

