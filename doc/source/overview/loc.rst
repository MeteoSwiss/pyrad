Location configuration file
==============================

The location configuration files describes some parameters that are depending on the specific
location of a radar (type of scans we want to measure, radar name, etc ). The location of the
weather radar (its position) itself, is instead usually read from the radar metadata directly and
it is not necessarily defined in this file. The fields are described in the table below.

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   * - Name
     - Type
     - Description
   * - RadarName
     - STRING
     - Short version name of a C-band radar (i.e., A, D, L, P) or DX50, MXPol for the X-band radars.
   * - RadarRes
     - STRING
     - rad4alp radar resolution (H or L). Only necessary if rad4alp (swiss C-band) data is processed.
   * - RadarBeamwidth
     - FLOAT
     - Radar antenna beam width [Deg].
   * - DataTypeIDInFiles
     - STRUCT
     - Structure of strings defining the mapping between pyrad names and variable names in input files. If not provided it will be expected that the name in the files are the standard py-ART names (defined in the py-ART config file). If you process multiple radars, you can append the index of the radar to DataTypeIDInFiles: DataTypeIDInFiles1, will be used for radar 1, DataTypeIDInFiles2, for radar 2, and so on. If you don't add this index, the same DataTypeIDInFiles struc will be used for all radars. An example of such a struct could be:

       - dBZ: DBZ
       - ZDR: ZDR
       - RhoHV: RHOHV
       - PhiDP: PHIDP
       - KDP: KDP
       - V: VEL
       - W: WIDTH
   * - AntennaGaindB
     - FLOAT
     - Antenna gain [dB].
   * - ScanList
     - STRARR
     - A list with the scans used for this data processing. Note that the first scan in this list is used as master scan. The master scan must be the first (temporal) scan of the corresponding rainbow task. In case of composite volumes the master scan is usually a PPI and the following are RHIs. If the radar processed is MCH C-band the scan list consists of the radar elevation (i.e., from 001 to 020). All scan names must have a trailing '/' except if rad4alp data is processed.
   * - ScanPeriod
     - FLOAT
     - Repetition period of each scan in minutes.
   * - MasterScanTimeTol
     - INT
     - Add a tolerance to the scan times when creating radar volumes. If set to 1 will consider all scans which timestamp falls within T0 to T0 + ScanPeriod as belonging to the same volume. If set to -1 will consider all scans which timestamp falls within T0 - ScanPeriod to T0 as belonging to the same volume. T0 = timestamp of masterscan. Default is 0 (= no tolerance)
   * - Azimtol
     - FLOAT
     - Tolerance in azimuth for irregular data. (0.5 is a good value).
   * - clutterMap
     - STRING
     - Clutter map of the data processing. The clutter map is located at ``<configpath>/clutter/<clutterMap>``.
   * - AntennaGain
     - FLOAT
     - Radar antenna gain. Not used for X-band MCH data.
   * - radarconsth(v)
     - FLOAT
     - Radar constant h (v). Not mandatory.
   * - mflossh(v)
     - FLOAT
     - Matched filter losses h (v). Not mandatory.
   * - attg
     - FLOAT
     - Gas attenuation coefficient (units? (1 way attenuation)).
   * - IconRunFreq
     - INT
     - Frequency of a Icon model run in hours.
   * - IconForecasted
     - INT
     - Hours forecasted by the Icon model.
   * - rmax
     - FLOAT
     - For C-band data, the maximum range in [m] to be considered. Useful for speed considerations.
   * - elmax
     - FLOAT
     - Maximum elevation [°] to consider.
   * - ppiImageConfig
     - STRUCT
     - Structure defining the PPI image generating (PPI_IMAGE or PSEUDOPPI_IMAGE products).
   * - rhiImageConfig
     - STRUCT
     - Structure defining the RHI image generating (RHI_IMAGE or PSEUDORHI_IMAGE products).
   * - ppiMapImageConfig
     - STRUCT
     - Structure defining the PPI image overlaid on a map (PPI_MAP product).
   * - gridMapImageConfig
     - STRUCT
     - Structure defining the display of gridded data overlaid on a map (SURFACE_IMAGE product).
   * - xsecImageConfig
     - STRUCT
     - Structure defining the cross-section images generated from gridded data (CROSS_SECTION, LATITUDE_SLICE and LONGITUDE_SLICE products).
   * - spectraImageConfig
     - STRUCT
     - Structure defining the Doppler spectral plots.
   * - sunhitsImageConfig
     - STRUCT
     - Structure defining the sun hits image.
   * - azPatternFile
     - STRING
     - Name of the azimuth pattern file of the antenna. This file and path must be ``<configpath>/antenna/<azPatternFile>``.
   * - elPatternFile
     - STRING
     - Name of the elevation pattern file of the antenna. This file and path must be ``<configpath>/antenna/<elPatternFile>``.
   * - fixed_angle
     - FLOAT
     - Fixed angle of a PAR antenna in degrees. For the PAR azimuth antenna this is the elevation angle. For the elevation antenna it is the azimuth angle.
