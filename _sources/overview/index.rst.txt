Pyrad configuration files
==============================

.. role:: red

The configuration of the data processing in pyrad is divided into three files. The main configuration file, the location configuration file describing the location of the
weather radar and the used scans. Finally, the product configuration file describes the datasets and
products.

Two formats are accepted for configuration files. The classical pyrad format which uses the following syntax::

   ppiMapImageConfig STRUCT 12
      mapres STRING 10m
      alpha FLOAT 0.4
      latmin FLOAT 46.5
      latmax FLOAT 47.1
      lonmin FLOAT 6.5
      lonmax FLOAT 7.5
      lonstep FLOAT 0.1
      latstep FLOAT 0.1
      xsize FLOAT 18.
      ysize FLOAT  10.
      background_zoom INT 10
      maps STRARR 2        # maps to overplot (cartopy)
         OTM
         rivers

**Supported Data Types**

Scalar types:

* BYTE, INT, LONG, HEX: interpreted as integers (HEX supports 0x notation).
* FLOAT, DOUBLE, EXP: interpreted as floating-point numbers.
* BOOL: accepts true/1 (case-insensitive) → True, others → False.
* STRING: interpreted as a string.

Array types:

* BYTARR, INTARR, LONARR, HEXARR: interpreted as a list of integers
* FLTARR, DBLAR: EXPARR: interpreted as a list of floating-point numbers.
* STRARR: interpreted as a list of strings.


Every entry has three elements:

* *fieldname*: indentifier for the field
* *TYPE*: one of the supported scalar or array types
* *value* (for scalar types) or *number of elements* (for STRUCT or arrays)

**The number of elements must be correct or the processing will fail**



:red:`Since v2.1.0 of pyrad, you can also write your config files in the YAML format.`
A `converter script <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/scripts/config_to_yaml.py>`_
is also available to convert your pyrad-style config files to YAML.

.. toctree::
   :maxdepth: 1

   main
   loc
   prod
