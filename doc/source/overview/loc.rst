Location configuration file
==============================

The location configuration files describes some parameters that are depending on the specific
location of a radar (type of scans we want to measure, radar name, etc ). The location of the
weather radar (its position) itself, is instead usually read from the radar metadata directly and
it is not necessarily defined in this file. The fields are described in the table below.

====================  =======  =======================================================================================
Name                  Type     Description
====================  =======  =======================================================================================
RadarName             STRING   Short version name of a C-band radar (i.e., A, D, L, P) or DX50, MXPol for the X-band radars.
RadarRes              STRING   rad4alp radar resolution (H or L). Only necessary if rad4alp (swiss C-band) data is processed.
RadarBeamwidth        FLOAT    Radar antenna beam width [Deg].
AntennaGaindB         FLOAT    Antenna gain [dB].
ScanList              STRARR   A list with the scans used for this data processing. Note that the first scan in this list is used as master scan. The master scan must be the first (temporal) scan of the corresponding rainbow task. In case of composite volumes the master scan is usually a PPI and the following are RHIs. If the radar processed is MCH C-band the scan list consists of the radar elevation (i.e., from 001 to 020). All scan names must have a trailing '/' except if rad4alp data is processed.
ScanPeriod            FLOAT    Repetition period of each scan in minutes.
Azimtol               FLOAT    Tolerance in azimuth for irregular data. (0.5 is a good value).
clutterMap            STRING   Clutter map of the data processing. The clutter map is located at ``<configpath>/clutter/<clutterMap>``.
AntennaGain           FLOAT    Radar antenna gain. Not used for X-band MCH data.
radarconsth(v)        FLOAT    Radar constant h (v). Not mandatory.
mflossh(v)            FLOAT    Matched filter losses h (v). Not mandatory.
attg                  FLOAT    Gas attenuation coefficient (units? (1 way attenuation)).
CosmoRunFreq          INT      Frequency of a COSMO model run in hours.
CosmoForecasted       INT      Hours forecasted by the COSMO model.
rmax                  FLOAT    For C-band data, the maximum range in [m] to be considered. Useful for speed considerations.
elmax                 FLOAT    Maximum elevation [°] to consider.
ppiImageConfig        STRUCT   Structure defining the PPI image generating (PPI_IMAGE or PSEUDOPPI_IMAGE products). The following 6 fields are described below:
- xsize               INT      Number of horizontal pixels of the picture (without frame).
- ysize               INT      Number of vertical pixels of the picture (without frame).
- xmin                FLOAT    Distance of the left image boundary to the radar in km.
- xmax                FLOAT    Distance of the right image boundary to the radar in km.
- ymin                FLOAT    Distance of the lower image boundary to the radar in km.
- ymax                FLOAT    Distance of the upper image boundary to the radar in km.
- dpi				  INT      Resolution of the image in dots per inch.
rhiImageConfig        STRUCT   Structure defining the RHI image generating (RHI_IMAGE or PSEUDORHI_IMAGE products). The following 6 fields are described below:
- xsize               INT      Number of horizontal pixels of the picture (without frame).
- ysize               INT      Number of vertical pixels of the picture (without frame).
- xmin                FLOAT    Distance of the left image boundary to the radar in km.
- xmax                FLOAT    Distance of the right image boundary to the radar in km.
- ymin                FLOAT    Distance of the lower image boundary (vertical direction) to the radar in km.
- ymax                FLOAT    Distance of the upper image boundary (vertical direction) to the radar in km.
- dpi				  INT      Resolution of the image in dots per inch.
ppiMapImageConfig     STRUCT   Structure defining the PPI image overlaid on a map (PPI_MAP product). The following 9 fields are described below:
- rngRing             FLOAT    Distance between range rings (0 means no range ring) [km].
- xsize               FLOAT    Image size (inches) [inch].
- ysize               FLOAT    Image size (inches) [inch].
- lonmin              FLOAT    Minimum WGS84 longitude [°].
- lonmax              FLOAT    Maximum WGS84 longitude [°].
- latmin              FLOAT    Minimum WGS84 latitude [°].
- latmax              FLOAT    Maximum WGS84 latitude [°].
- latstep             FLOAT    Step in latitude [°] used in the map gridlines.
- lonstep             FLOAT    Step in longitude [°] used in the map gridlines.
- exact_limits	      INT      If set to 1 will use exactly the user-specified latmin/latmax, lonmin/lonmax, if set to 0 will round them to the nearest integer.
- mapres              STRING   Map resolution. Accepted strings are: “10m”, “50m”, “110m”.
- maps                STRARR   String array of possible maps to overplot. Accepted entries include: relief, countries, provinces, 
                                 urban_areas, roads, railroads, coastline, lakes, lakes_europe, rivers, rivers_europe.
- dpi				  INT      Resolution of the image in dots per inch.
gridMapImageConfig    STRUCT   Structure defining the display of gridded data overlaid on a map (SURFACE_IMAGE product).
- xsize               FLOAT    Image size (inches) [inch].
- ysize               FLOAT    Image size (inches) [inch].
- lonmin              FLOAT    Minimum WGS84 longitude [°].
- lonmax              FLOAT    Maximum WGS84 longitude [°].
- latmin              FLOAT    Minimum WGS84 latitude [°].
- latmax              FLOAT    Maximum WGS84 latitude [°].
- latstep             FLOAT    Step in latitude [°] used in the map gridlines.
- lonstep             FLOAT    Step in longitude [°] used in the map gridlines.
- exact_limits	      INT      If set to 1 will use exactly the user-specified latmin/latmax, lonmin/lonmax, if set to 0 will round them to the nearest integer.
- mapres              STRING   Map resolution. Accepted strings are: “10m”, “50m”, “110m”.
- maps                STRARR   String array of possible maps to overplot. Accepted entries include: relief, countries, provinces, 
                                 urban_areas, roads, railroads, coastline, lakes, lakes_europe, rivers, rivers_europe
- dpi				  INT      Resolution of the image in dots per inch.
xsecImageConfig       STRUCT   Structure defining the cross-section images generated from gridded data (CROSS_SECTION, LATITUDE_SLICE and LONGITUDE_SLICE products)
- xsize               INT      Number of horizontal pixels of the picture (without frame).
- ysize               INT      Number of vertical pixels of the picture (without frame).
- xmin                FLOAT    Distance of the left image boundary to the radar in km.
- xmax                FLOAT    Distance of the right image boundary to the radar in km.
- ymin                FLOAT    Distance of the lower image boundary (vertical direction) to the radar in km.
- ymax                FLOAT    Distance of the upper image boundary (vertical direction) to the radar in km.
- dpi				  INT      Resolution of the image in dots per inch.
sunhitsImageConfig    STRUCT   Structure defining the sun hits image. The following 6 fields are described below:
- xsize               INT      Number of horizontal pixels of the picture (without frame).
- ysize               INT      Number of vertical pixels of the picture (without frame).
- xmin                FLOAT    Minimum azimuth angle difference (between sun and radar).
- xmax                FLOAT    Maximum azimuth angle difference (between sun and radar).
- ymin                FLOAT    Minimum elevation angle difference (between sun and radar).
- ymax                FLOAT    Maximum azimuth angle difference (between sun and radar).
- dpi				  INT      Resolution of the image in dots per inch.
azPatternFile         STRING   Name of the azimuth pattern file of the antenna. This file and path must be ``<configpath>/antenna/<azPatternFile>``.
elPatternFile         STRING   Name of the elevation pattern file of the antenna. This file and path must be ``<configpath>/antenna/<elPatternFile>``.
fixed_angle           FLOAT    Fixed angle of a PAR antenna in degrees. For the PAR azimuth antenna this is the elevation angle. For the elevation antenna it is the azimuth angle.
====================  =======  =======================================================================================




