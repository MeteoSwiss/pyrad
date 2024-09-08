List of pyrad variables
==============================

.. note::
   Pyrad uses the following mappings to map the variable names in your files to the short names used by pyrad (dBZ, dBZc, RhoHV and else).
   the "Py-ART to Pyrad" mappings are used for CFRadial files and the "ODIM to Pyrad" mappings are used for ODIM HDF5 files. If your files do not contain the standard 
   names listed below, you can use the keyword DataTypeIDInFiles in the loc files to provide them. Please see :doc:`loc`.


ODIM to Pyrad
------------------------------------

.. csv-table:: ODIM to Pyrad mappings
   :file: mappings/pyrad_to_odim.txt
   :header-rows: 1


Py-ART to Pyrad
------------------------------------

.. csv-table::  Py-ART to Pyrad mappings
   :file: mappings/pyrad_to_pyart.txt
   :header-rows: 1


ICON to Pyrad
------------------------------------

.. csv-table:: ICON to Pyrad mappings
   :file: mappings/pyrad_to_icon.txt
   :header-rows: 1


Metranet (Swiss radars) to Pyrad
------------------------------------

.. csv-table:: Metranet (Swiss radars) to Pyrad mappings
   :file: mappings/pyrad_to_metranet.txt
   :header-rows: 1


