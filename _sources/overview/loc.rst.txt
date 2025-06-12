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
     - Structure defining the PPI image generating (PPI_IMAGE or PSEUDOPPI_IMAGE products). It can contain the following keys:

       - dpi: INT: dpi of the output image
       - xsize: FLOAT: width of the image in inches
       - ysize: FLOAT: height of the image in inches
       - ymin: FLOAT: minimum displayed y-coordinate (northwards, in radar centric coordinates)
       - ymax: FLOAT: minimum displayed y-coordinate (northwards, in radar centric coordinates)
       - xmin: FLOAT: minimum displayed x-coordinate (eastwards, in radar centric coordinates)
       - xmax: FLOAT: minimum displayed x-coordinate (northwards, in radar centric coordinates)
       - rngRing: INT: if defined will plot range rings every x kilometers from the radar (x is the specified value)
   * - rhiImageConfig
     - STRUCT
     - Structure defining the RHI image generating (RHI_IMAGE or PSEUDORHI_IMAGE products). It can contain the following keys:

       - dpi: INT: dpi of the output image
       - xsize: FLOAT: width of the image in inches
       - ysize: FLOAT: height of the image in inches
       - ymin: FLOAT: minimum displayed y-coordinate (altitude)
       - ymax: FLOAT: minimum displayed y-coordinate (altitude)
       - xmin: FLOAT: minimum displayed x-coordinate (distance at ground along RHI)
       - xmax: FLOAT: minimum displayed x-coordinate (distance at ground along RHI)
   * - ppiMapImageConfig
     - STRUCT
     - Structure defining the PPI image overlaid on a map (PPI_MAP product). It can contain the following keys:

       - dpi: INT: dpi of the output image
       - xsize: FLOAT: width of the image in inches
       - ysize: FLOAT: height of the image in inches
       - min_lon: FLOAT: minimum displayed longitude
       - max_lon: FLOAT: maximum displayed longitude
       - min_lat: FLOAT: minimum displayed latitude
       - max_lat: FLOAT: maximum displayed latitude
       - resolution: INT: resolution of additional cartopy geodata, only "10m", "50m" or "110m" is supported
       - alpha: FLOAT: transparency value (between 0 and 1) of the displayed radar data, alpha smaller than 1 is required when plotting over a relief or OTM map
       - background_zoom: INT: Zoom level of the additional cartopy raster geodata, the higher the value, the higher the resolution (typically between 8 and 12). Default is 8.
       - maps: list of STRING: list of additional geodata to plot. The following are supported: relief (hillshade), OTM (opentopomaps), provinces, urban_areas, roads, railroads, coastlines, lakes, lakes_europe, rivers, rivers_europe.
       - rngRing: INT: if defined will plot range rings every x kilometers from the radar (x is the specified value)
   * - gridMapImageConfig
     - STRUCT
     - Structure defining the display of gridded data overlaid on a map (SURFACE_IMAGE product). It can contain the following keys:

       - dpi: INT: dpi of the output image
       - xsize: FLOAT: width of the image in inches
       - ysize: FLOAT: height of the image in inches
       - min_lon: FLOAT: minimum displayed longitude
       - max_lon: FLOAT: maximum displayed longitude
       - min_lat: FLOAT: minimum displayed latitude
       - max_lat: FLOAT: maximum displayed latitude
       - resolution: INT: resolution of additional cartopy geodata, only "10m", "50m" or "110m" is supported
       - alpha: FLOAT: transparency value (between 0 and 1) of the displayed radar data, alpha smaller than 1 is required when plotting over a relief or OTM map
       - background_zoom: INT: Zoom level of the additional cartopy raster geodata, the higher the value, the higher the resolution (typically between 8 and 12). Default is 8.
       - maps: list of STRING: list of additional geodata to plot. The following are supported: relief (hillshade), OTM (opentopomaps), provinces, urban_areas, roads, railroads, coastlines, lakes, lakes_europe, rivers, rivers_europe.
   * - xsecImageConfig
     - STRUCT
     - Structure defining the cross-section images generated from gridded data (CROSS_SECTION, LATITUDE_SLICE and LONGITUDE_SLICE products). It can contain the following keys:

       - dpi: INT: dpi of the output image
       - xsize: FLOAT: width of the image in inches
       - ysize: FLOAT: height of the image in inches
       - ymin: FLOAT: minimum displayed y-coordinate (altitude)
       - ymax: FLOAT: minimum displayed y-coordinate (altitude)
       - xmin: FLOAT: minimum displayed x-coordinate (distance at ground along cross-section)
       - xmax: FLOAT: minimum displayed x-coordinate (distance at ground along cross-section)
   * - spectraImageConfig
     - STRUCT
     - Structure defining the Doppler spectral plots (radar observable as a function of altitude and velocity).  It can contain the following keys:

       - dpi: INT: dpi of the output image
       - xsize: FLOAT: width of the image in inches
       - ysize: FLOAT: height of the image in inches
       - ymin: FLOAT: minimum displayed y-coordinate (altitude)
       - ymax: FLOAT: minimum displayed y-coordinate (altitude)
       - velmin: FLOAT: minimum displayed velocity
       - velmax: FLOAT: minimum displayed velocity
   * - sunhitsImageConfig
     - STRUCT
     - Structure defining the sun hits image. It shows the sun observables as a function of azimuth and elevation.  It can contain the following keys:

       - dpi: INT: dpi of the output image
       - xsize: FLOAT: width of the image in inches
       - ysize: FLOAT: height of the image in inches
       - azmin: FLOAT: minimum displayed azimuth angle
       - azmax: FLOAT: minimum displayed azimuth angle
       - elmin: FLOAT: minimum displayed elevation angle
       - elmax: FLOAT: minimum displayed elevation angle
       - azres: FLOAT: azimuth step to use in the plot
       - elres: FLOAT: elevation step to use in the plot
   * - azPatternFile
     - STRING
     - Name of the azimuth pattern file of the antenna. This file and path must be ``<configpath>/antenna/<azPatternFile>``.
   * - elPatternFile
     - STRING
     - Name of the elevation pattern file of the antenna. This file and path must be ``<configpath>/antenna/<elPatternFile>``.
   * - fixed_angle
     - FLOAT
     - Fixed angle of a PAR antenna in degrees. For the PAR azimuth antenna this is the elevation angle. For the elevation antenna it is the azimuth angle.
