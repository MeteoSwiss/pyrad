List of pyrad products
==============================

VOL
-----------------------------
ANTENNA_POS
""""""""""""""""""""""""""""""
description
   Plots time series of the antenna position with respect to elevation or azimuth
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L724>`_
parameters
   dpi float
        dpi of the plot 
   datatype str
        type of data to plot. Can be AZ or EL 

CDF
""""""""""""""""""""""""""""""
description
   plots and writes the cumulative density function of data
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L2943>`_
parameters
   quantiles list
        of floats The quantiles to compute in percent. Default None 
   sector dict
        dictionary defining the sector where to compute the CDF. Default is None and the CDF is computed over all the data May 
   contain rmin
       , 
   rmax float
        min and max range [m] azmin, 
   azmax float
        min and max azimuth angle [deg] elmin, 
   elmax float
        min and max elevation angle [deg] hmin, 
   hmax float
        min and max altitude [m MSL] 
   vismin float
        The minimum visibility to use the data. Default None 
   absolute Bool
        If true the absolute values of the data will be used. Default False 
   use_nans Bool
        If true NaN values will be used. Default False 
   nan_value Bool
        The value by which the NaNs are substituted if NaN values are to be used in the computation 
   filterclt Bool
        If True the gates containing clutter are filtered 
   filterprec list
        of ints The hydrometeor types that are filtered from the analysis. Default empty list.

BSCOPE_IMAGE
""""""""""""""""""""""""""""""
description
   Creates a B-scope image (azimuth, range)
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L2449>`_
parameters
   anglenr  int
        The elevation angle number to use 
   ray_dim  str
        the ray dimension. Can be 'ang' or 'time'. Default 'ang' 
   xaxis_rng  bool
        if True the range will be in the x-axis. Otherwise it will be in the y-axis. Default True vmin, 
   vmax float
        or None The minimum and maximum values of the color scale. If None the scale is going to be set according to the Py-ART config file

CAPPI_IMAGE
""""""""""""""""""""""""""""""
description
   Creates a CAPPI image
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L2080>`_
parameters
   altitude flt
        CAPPI altitude [m MSL] 
   wfunc str
        The function used to produce the CAPPI as defined in pyart.map.grid_from_radars. Default 'NEAREST' 
   cappi_res float
        The CAPPI resolution [m]. Default 500.

CROSS_SECTION
""""""""""""""""""""""""""""""
description
   Plots a cross-section of polar data through arbitrary coordinates
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L2109>`_
parameters
   coordN dict
        The two lat-lon coordinates marking the limits. They have the keywords 'lat' and 'lon' [degree]. 
   step  int
        Step in meters to use between reference points to calculate the cross-section (i.e horizontal resolution). 
   vert_res  int
        Vertical resolution in meters used to calculate the cross-section 
   alt_max  int
        Maximum altitude of the vertical cross-section 
   beamwidth  float
        3dB beamwidth in degrees to be used in the calculations, if not provided will be read from the loc file 
   demfile  str
        Name of the DEM file to use to plot the topography, it must be in the dempath specified in the main config file

FIELD_COVERAGE
""""""""""""""""""""""""""""""
description
   Gets the field coverage over a certain sector
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L2801>`_
parameters
   threshold float
        or None Minimum value to consider the data valid. Default None 
   nvalid_min float
        Minimum number of valid gates in the ray to consider it valid. Default 5 ele_res, 
   azi_res float
        Elevation and azimuth resolution of the sectors [deg]. Default 1. and 2. ele_min, 
   ele_max float
        Min and max elevation angle defining the sector [deg]. Default 0. and 30. 
   ele_step float
        Elevation step [deg]. Default 5. ele_sect_start, 
   ele_sect_stop float
        or None start and stop angles of the sector coverage. Default None 
   quantiles list
        of floats The quantiles to compute in the sector. Default 10. to 90. by steps of 10. 
   AngTol float
        The tolerance in elevation angle when putting the data in a fixed grid

FIXED_RNG_IMAGE
""""""""""""""""""""""""""""""
description
   Plots a fixed range image
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L2188>`_
parameters
   AngTol  float
        The tolerance between the nominal angles and the actual radar angles. Default 1. ele_res, 
   azi_res float
        or None The resolution of the fixed grid [deg]. If None it will be obtained from the separation between angles vmin, 
   vmax  float
        or None Min and Max values of the color scale. If None the values are taken from the Py-ART config file

FIXED_RNG_SPAN_IMAGE
""""""""""""""""""""""""""""""
description
   Plots a user-defined statistic over a fixed range image
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L2227>`_
parameters
   AngTol  float
        The tolerance between the nominal angles and the actual radar angles. Default 1. ele_res, 
   azi_res float
        or None The resolution of the fixed grid [deg]. If None it will be obtained from the separation between angles 
   stat  str
        The statistic to compute. Can be 'min', 'max', 'mean', 'mode'. Default 'max'

HISTOGRAM
""""""""""""""""""""""""""""""
description
   Computes a histogram of the radar volum data
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L2598>`_
parameters
   step float
        or None the data quantization step. If none it will be obtained from the Py-ART configuration file 
   write_data Bool
        If true the histogram data is written in a csv file

PLOT_ALONG_COORD
""""""""""""""""""""""""""""""
description
   Plots the radar volume data along a particular coordinate
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L2265>`_
parameters
   colors list
        of str or None The colors of each ploted line 
   data_on_y  bool
        If True the x-axis is the coordinates, and the y the data values. False swaps the axis. Default True 
   plot_legend  bool
        If True a legend will be plotted. Default True 
   mode str
        Ploting mode. Can be 'ALONG_RNG', 'ALONG_AZI' or 'ALONG_ELE' value_start, 
   value_stop float
        The starting and ending points of the data to plot. According to the mode it may refer to the range, azimuth or elevation. If not specified the minimum and maximum possible values are used fix_elevations, fix_azimuths, 
   fix_ranges list
        of floats The elevations, azimuths or ranges to plot for each mode. 'ALONG_RNG' would use fix_elevations and fix_azimuths 'ALONG_AZI' fix_ranges and fix_elevations 'ALONG_ELE' fix_ranges and fix_azimuths 
   AngTol float
        The tolerance to match the radar angle to the fixed angles Default 1. 
   RngTol float
        The tolerance to match the radar range to the fixed ranges Default 50. 
   use_altitude  bool
        If true and in ALON_RNG mode the coordinate used is the gate altitude. Otherwise is the range. Default False

PLOT_TXH
""""""""""""""""""""""""""""""
description
   Plots the transmitted signal power (H) for a standard sunscan.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L2410>`_
parameters

PPI_CONTOUR
""""""""""""""""""""""""""""""
description
   Plots a PPI countour plot
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L1095>`_
parameters
   contour_values list
        of floats or None The list of contour values to plot. If None the contour values are going to be obtained from the Py-ART config file either with the dictionary key 'contour_values' or from the minimum and maximum values of the field with an assumed division of 10 levels. 
   anglenr float
        The elevation angle number

PPI_CONTOUR_OVERPLOT
""""""""""""""""""""""""""""""
description
   Plots a PPI of a field with another field overplotted as a contour plot.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L1000>`_
parameters
   contour_values list
        of floats or None The list of contour values to plot. If None the contour values are going to be obtained from the Py-ART config file either with the dictionary key 'contour_values' or from the minimum and maximum values of the field with an assumed division of 10 levels. 
   anglenr float
        The elevation angle number

PPI_IMAGE
""""""""""""""""""""""""""""""
description
   Plots a PPI image. It can also plot the histogram and the quantiles of the data in the PPI.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L2080>`_
parameters
   anglenr float
        The elevation angle number 
   plot_type str
        The type of plot to perform. Can be 'PPI', 'QUANTILES' or 'HISTOGRAM' 
   write_data Bool
        If True the histrogram will be also written in a csv file 
   step float
        or None If the plot type is 'HISTOGRAM', the width of the histogram bin. If None it will be obtained from the Py-ART config file 
   quantiles list
        of float or None If the plot type is 'QUANTILES', the list of quantiles to compute. If None a default list of quantiles will be computed vmin, 
   vmax float
        or None The minimum and maximum values of the color scale. If None the scale is going to be set according to the Py-ART config file

PPI_MAP
""""""""""""""""""""""""""""""
description
   Plots a PPI image over a map. The map resolution and the type of maps used are defined in the variables 'mapres' and 'maps' in 'ppiMapImageConfig' in the loc config file.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L2041>`_
parameters
   anglenr float
        The elevation angle number

PPIMAP_ROI_OVERPLOT
""""""""""""""""""""""""""""""
description
   Over plots a polygon delimiting a region of interest on a PPI map. The map resolution and the type of maps used are defined in the variables 'mapres' and 'maps' in 'ppiMapImageConfig' in the loc config file.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L950>`_
parameters
   anglenr float
        The elevation angle number

PROFILE_STATS
""""""""""""""""""""""""""""""
description
   Computes and plots a vertical profile statistics. The statistics are saved in a csv file
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L1575>`_
parameters
   heightResolution float
        The height resolution of the profile [m]. Default 100. heightMin, 
   heightMax float
        or None The minimum and maximum altitude of the profile [m MSL]. If None the values will be obtained from the minimum and maximum gate altitude. 
   quantity str
        The type of statistics to plot. Can be 'quantiles', 'mode', 'reqgression_mean' or 'mean'. 
   quantiles list
        of floats If quantity type is 'quantiles' the list of quantiles to compute. Default 25., 50., 75. 
   nvalid_min int
        The minimum number of valid points to consider the statistic valid. Default 4 
   make_linear Bool
        If true the data is converted from log to linear before computing the stats 
   include_nans Bool
        If true NaN values are included in the statistics 
   fixed_span Bool
        If true the profile plot has a fix X-axis vmin, 
   vmax float
        or None If fixed_span is set, the minimum and maximum values of the X-axis. If None, they are obtained from the Py-ART config file

PSEUDOPPI_CONTOUR
""""""""""""""""""""""""""""""
description
   Plots a pseudo-PPI countour plot
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L1095>`_
parameters
   contour_values list
        of floats or None The list of contour values to plot. If None the contour values are going to be obtained from the Py-ART config file either with the dictionary key 'contour_values' or from the minimum and maximum values of the field with an assumed division of 10 levels. 
   angle float
        The elevation angle at which compute the PPI 
   EleTol float
        The tolerance between the actual radar elevation angle and the nominal pseudo-PPI elevation angle.

PSEUDOPPI_CONTOUR_OVERPLOT
""""""""""""""""""""""""""""""
description
   Plots a pseudo-PPI of a field with another field over-plotted as a contour plot
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L1000>`_
parameters
   contour_values list
        of floats or None The list of contour values to plot. If None the contour values are going to be obtained from the Py-ART config file either with the dictionary key 'contour_values' or from the minimum and maximum values of the field with an assumed division of 10 levels. 
   angle float
        The elevation angle at which compute the PPI 
   EleTol float
        The tolerance between the actual radar elevation angle and the nominal pseudo-PPI elevation angle.

PSEUDOPPI_IMAGE
""""""""""""""""""""""""""""""
description
   Plots a pseudo-PPI image. It can also plot the histogram and the quantiles of the data in the pseudo-PPI.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L814>`_
parameters
   angle float
        The elevation angle of the pseudo-PPI 
   EleTol float
        The tolerance between the actual radar elevation angle and the nominal pseudo-PPI elevation angle. 
   plot_type str
        The type of plot to perform. Can be 'PPI', 'QUANTILES' or 'HISTOGRAM' 
   step float
        or None If the plot type is 'HISTOGRAM', the width of the histogram bin. If None it will be obtained from the Py-ART config file 
   quantiles list
        of float or None If the plot type is 'QUANTILES', the list of quantiles to compute. If None a default list of quantiles will be computed vmin, 
   vmax  float
        or None Min and Max values of the color scale. If None the values are taken from the Py-ART config file

PSEUDOPPI_MAP
""""""""""""""""""""""""""""""
description
   Plots a pseudo-PPI image over a map. The map resolution and the type of maps used are defined in the variables 'mapres' and 'maps' in 'ppiMapImageConfig' in the loc config file.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L2041>`_
parameters
   angle float
        The elevation angle of the pseudo-PPI 
   EleTol float
        The tolerance between the actual radar elevation angle and the nominal pseudo-PPI elevation angle.

PSEUDORHI_CONTOUR
""""""""""""""""""""""""""""""
description
   Plots a pseudo-RHI countour plot
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L1376>`_
parameters
   contour_values list
        of floats or None The list of contour values to plot. If None the contour values are going to be obtained from the Py-ART config file either with the dictionary key 'contour_values' or from the minimum and maximum values of the field with an assumed division of 10 levels. 
   angle float
        The azimuth angle at which to compute the RPI 
   AziTol float
        The tolerance between the actual radar azimuth angle and the nominal pseudo-RHI azimuth angle.

PSEUDORHI_CONTOUR_OVERPLOT
""""""""""""""""""""""""""""""
description
   Plots a pseudo-RHI of a field with another field over-plotted as a contour plot
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L1280>`_
parameters
   contour_values list
        of floats or None The list of contour values to plot. If None the contour values are going to be obtained from the Py-ART config file either with the dictionary key 'contour_values' or from the minimum and maximum values of the field with an assumed division of 10 levels. 
   angle float
        The azimuth angle at which to compute the RPI 
   AziTol float
        The tolerance between the actual radar azimuth angle and the nominal pseudo-RHI azimuth angle.

PSEUDORHI_IMAGE
""""""""""""""""""""""""""""""
description
   Plots a pseudo-RHI image. It can also plot the histogram and the quantiles of the data in the pseudo-RHI.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L1178>`_
parameters
   angle float
        The azimuth angle at which to compute the RPI 
   AziTol float
        The tolerance between the actual radar azimuth angle and the nominal pseudo-RHI azimuth angle. 
   plot_type str
        The type of plot to perform. Can be 'RHI', 'QUANTILES' or 'HISTOGRAM' 
   step float
        or None If the plot type is 'HISTOGRAM', the width of the histogram bin. If None it will be obtained from the Py-ART config file 
   quantiles list
        of float or None If the plot type is 'QUANTILES', the list of quantiles to compute. If None a default list of quantiles will be computed vmin, 
   vmax  float
        or None Min and Max values of the color scale. If None the values are taken from the Py-ART config file

QUANTILES
""""""""""""""""""""""""""""""
description
   Plots and writes the quantiles of a radar volume
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L2659>`_
parameters
   quantiles list
        of floats or None the list of quantiles to compute. If None a default list of quantiles will be computed. 
   write_data Bool
        If True the computed data will be also written in a csv file 
   fixed_span Bool
        If true the quantile plot has a fix Y-axis vmin, 
   vmax float
        or None If fixed_span is set, the minimum and maximum values of the Y-axis. If None, they are obtained from the Py-ART config file

RHI_CONTOUR
""""""""""""""""""""""""""""""
description
   Plots an RHI countour plot
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L1376>`_
parameters
   contour_values list
        of floats or None The list of contour values to plot. If None the contour values are going to be obtained from the Py-ART config file either with the dictionary key 'contour_values' or from the minimum and maximum values of the field with an assumed division of 10 levels. 
   anglenr int
        The azimuth angle number

RHI_CONTOUR_OVERPLOT
""""""""""""""""""""""""""""""
description
   Plots an RHI of a field with another field over-plotted as a contour plot
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L1280>`_
parameters
   contour_values list
        of floats or None The list of contour values to plot. If None the contour values are going to be obtained from the Py-ART config file either with the dictionary key 'contour_values' or from the minimum and maximum values of the field with an assumed division of 10 levels. 
   anglenr int
        The azimuth angle number

RHI_IMAGE
""""""""""""""""""""""""""""""
description
   Plots an RHI image. It can also plot the histogram and the quantiles of the data in the RHI.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L1178>`_
parameters
   anglenr int
        The azimuth angle number 
   plot_type str
        The type of plot to perform. Can be 'RHI', 'QUANTILES' or 'HISTOGRAM' 
   step float
        or None If the plot type is 'HISTOGRAM', the width of the histogram bin. If None it will be obtained from the Py-ART config file 
   quantiles list
        of float or None If the plot type is 'QUANTILES', the list of quantiles to compute. If None a default list of quantiles will be computed vmin, 
   vmax float
        or None The minimum and maximum values of the color scale. If None the scale is going to be set according to the Py-ART config file

RHI_PROFILE
""""""""""""""""""""""""""""""
description
   Computes and plots a vertical profile statistics out of an RHI. The statistics are saved in a csv file
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L1419>`_
parameters
   rangeStop float
        The range start and stop of the data to extract from the RHI to compute the statistics [m]. Default 0., 25000. 
   heightResolution float
        The height resolution of the profile [m]. Default 100. heightMin, 
   heightMax float
        or None The minimum and maximum altitude of the profile [m MSL]. If None the values will be obtained from the minimum and maximum gate altitude. 
   quantity str
        The type of statistics to plot. Can be 'quantiles', 'mode', 'reqgression_mean' or 'mean'. 
   quantiles list
        of floats If quantity type is 'quantiles' the list of quantiles to compute. Default 25., 50., 75. 
   nvalid_min int
        The minimum number of valid points to consider the statistic valid. Default 4 
   make_linear Bool
        If true the data is converted from log to linear before computing the stats 
   include_nans Bool
        If true NaN values are included in the statistics 
   fixed_span Bool
        If true the profile plot has a fix X-axis vmin, 
   vmax float
        or None If fixed_span is set, the minimum and maximum values of the X-axis. If None, they are obtained from the Py-ART config file

SAVEALL
""""""""""""""""""""""""""""""
description
   Saves radar volume data including all or a list of user- defined fields in a C/F radial or ODIM file
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L3513>`_
parameters
   file_type str
        The type of file used to save the data. Can be 'nc' or 'h5'. Default 'nc' 
   datatypes list
        of str or None The list of data types to save. If it is None, all fields in the radar object will be saved 
   physical Bool
        If True the data will be saved in physical units (floats). Otherwise it will be quantized and saved as binary 
   compression str
        For ODIM file formats, the type of compression. Can be any of the allowed compression types for hdf5 files. Default gzip 
   compression_opts any
        The compression options allowed by the hdf5. Depends on the type of compression. Default 6 (The gzip compression level).

SAVEALL_VOL
""""""""""""""""""""""""""""""
description
   Same as before but can be used in a mixed GRID/VOL dataset, as there is no ambiguity with SAVEALL for VOL datasets
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L3513>`_
parameters

SAVESTATE
""""""""""""""""""""""""""""""
description
   Saves the last processed data in a file. Used for real- time data processing
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L3563>`_
parameters

SAVEPSEUDORHI
""""""""""""""""""""""""""""""
description
   Saves one field of a pseudo-RHI computed from a volume scan in C/F radial or ODIM file User defined paraeters: file_type: str The type of file used to save the data. Can be 'nc' or 'h5'. Default 'nc' physical: Bool If True the data will be saved in physical units (floats). Otherwise it will be quantized and saved as binary compression: str For ODIM file formats, the type of compression. Can be any of the allowed compression types for hdf5 files. Default gzip compression_opts: any The compression options allowed by the hdf5. Depends on the type of compression. Default 6 (The gzip compression level).
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L3692>`_
parameters

SAVEPSEUDOPPI
""""""""""""""""""""""""""""""
description
   Saves one field of a pseudo-PPI computed from a volume scan in C/F radial or ODIM file User defined paraeters: file_type: str The type of file used to save the data. Can be 'nc' or 'h5'. Default 'nc' physical: Bool If True the data will be saved in physical units (floats). Otherwise it will be quantized and saved as binary compression: str For ODIM file formats, the type of compression. Can be any of the allowed compression types for hdf5 files. Default gzip compression_opts: any The compression options allowed by the hdf5. Depends on the type of compression. Default 6 (The gzip compression level).
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L3743>`_
parameters

SAVEVOL
""""""""""""""""""""""""""""""
description
   Saves one field of a radar volume data in a C/F radial or ODIM file
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L3469>`_
parameters
   file_type str
        The type of file used to save the data. Can be 'nc' or 'h5'. Default 'nc' 
   physical Bool
        If True the data will be saved in physical units (floats). Otherwise it will be quantized and saved as binary. Default True 
   compression str
        For ODIM file formats, the type of compression. Can be any of the allowed compression types for hdf5 files. Default gzip 
   compression_opts any
        The compression options allowed by the hdf5. Depends on the type of compression. Default 6 (The gzip compression level).

SAVEVOL_CSV
""""""""""""""""""""""""""""""
description
   Saves one field of a radar volume data in a CSV file
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L3391>`_
parameters
   ignore_masked bool
        If True masked values will not be saved. Default False

SAVEVOL_KML
""""""""""""""""""""""""""""""
description
   Saves one field of a radar volume data in a KML file
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L3420>`_
parameters
   ignore_masked bool
        If True masked values will not be saved. Default False 
   azi_res  float
        or None azimuthal resolution of the range bins. If None the antenna beamwidth is going to be used to determine the resolution

SAVEVOL_VOL
""""""""""""""""""""""""""""""
description
   Same as before but can be used in a mixed GRID/VOL dataset, as there is no ambiguity with SAVEVOL for GRID datasets
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L3469>`_
parameters

SAVE_FIXED_ANGLE
""""""""""""""""""""""""""""""
description
   Saves the position of the first fix angle in a csv file
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L3796>`_
parameters

SELFCONSISTENCY
""""""""""""""""""""""""""""""
description
   Plots a ZDR versus KDP/ZH histogram of data.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L3145>`_
parameters
   retrieve_relation  bool
        If True plots also the retrieved relationship. Default True 
   plot_theoretical  bool
        If True plots also the theoretical relationship. Default True 
   normalize  bool
        If True the occurrence density of ZK/KDP for each ZDR bin is going to be represented. Otherwise it will show the number of gates at each bin. Default True

SELFCONSISTENCY2
""""""""""""""""""""""""""""""
description
   Plots a ZH measured versus ZH inferred from a self-consistency relation histogram of data.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L3145>`_
parameters
   normalize  bool
        If True the occurrence density of ZK/KDP for each ZDR bin is going to be represented. Otherwise it will show the number of gates at each bin. Default True

TIME_RANGE
""""""""""""""""""""""""""""""
description
   Plots a time-range/azimuth/elevation plot
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L2493>`_
parameters
   anglenr float
        The number of the fixed angle to plot vmin, 
   vmax float
        or None The minimum and maximum values of the color scale. If None the scale is going to be set according to the Py-ART config file

VOL_TS
""""""""""""""""""""""""""""""
description
   Writes and plots a value corresponding to a time series. Meant primarily for writing and plotting the results of the SELFCONSISTENCY2 algorithm
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L3178>`_
parameters
   ref_value float
        The reference value. Default 0 
   sort_by_date Bool
        If true when reading the csv file containing the statistics the data is sorted by date. Default False 
   rewrite Bool
        If true the csv file containing the statistics is rewritten 
   add_data_in_fname Bool
        If true and the data used is cumulative the year is written in the csv file name and the plot file name 
   npoints_min int
        Minimum number of points to use the data point in the plotting and to send an alarm. Default 0 vmin, 
   vmax float
        or None Limits of the Y-axis (data value). If None the limits are obtained from the Py-ART config file 
   alarm Bool
        If true an alarm is sent 
   tol_abs float
        Margin of tolerance from the reference value. If the current value is above this margin an alarm is sent. If the margin is not specified it is not possible to send any alarm 
   tol_trend float
        Margin of tolerance from the reference value. If the trend of the last X events is above this margin an alarm is sent. If the margin is not specified it is not possible to send any alarm 
   nevents_min int
        Minimum number of events with sufficient points to send an alarm related to the trend. If not specified it is not possible to send any alarm 
   sender str
        The mail of the alarm sender. If not specified it is not possible to send any alarm 
   receiver_list list
        of str The list of emails of the people that will receive the alarm.. If not specified it is not possible to send any alarm

WIND_PROFILE
""""""""""""""""""""""""""""""
description
   Plots vertical profile of wind data (U, V, W components and wind velocity and direction) out of a radar volume containing the retrieved U,V and W components of the wind, the standard deviation of the retrieval and the velocity difference between the estimated radial velocity (assuming the wind to be uniform) and the actual measured radial velocity.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_vol_products.py#L1710>`_
parameters
   heightResolution float
        The height resolution of the profile [m]. Default 100. heightMin, 
   heightMax float
        or None The minimum and maximum altitude of the profile [m MSL]. If None the values will be obtained from the minimum and maximum gate altitude. 
   min_ele float
        The minimum elevation to be used in the computation of the vertical velocities. Default 5. 
   max_ele float
        The maximum elevation to be used in the computation of the horizontal velocities. Default 85. 
   fixed_span Bool
        If true the profile plot has a fix X-axis vmin, 
   vmax float
        or None If fixed_span is set, the minimum and maximum values of the X-axis. If None, they are obtained from the span of the U component defined in the Py-ART config file 

CENTROIDS
-----------------------------
HISTOGRAM
""""""""""""""""""""""""""""""
description
   Plots the histogram of one of the variables used for centroids computation.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L1419>`_
parameters
   voltype  str
        The name of the variable to plot. Can be dBZ, ZDR, KDP, RhoHV, H_ISO0 and its standardized form (e.g. dBZ_std) 
   write_data  Bool
        If true writes the histogram in a .csv file. Default True 
   step  float
        bin size. Default 0.1

HISTOGRAM2D
""""""""""""""""""""""""""""""
description
   Plots the 2D- histogram of two of the variables used for centroids computation.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L1419>`_
parameters
   voltype_y  str
        The name of the variables to plot. Can be dBZ, ZDR, KDP, RhoHV, H_ISO0 and its standardized form (e.g. dBZ_std) step_x, 
   step_y  float
        bin size. Default 0.1

HISTOGRAM_LABELED
""""""""""""""""""""""""""""""
description
   Plots the histogram of one of the variables used for centroids computation. Only plots labeled data.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L1258>`_
parameters
   voltype  str
        The name of the variable to plot. Can be dBZ, ZDR, KDP, RhoHV, H_ISO0 and its standardized form (e.g. dBZ_std) 
   write_data  Bool
        If true writes the histogram in a .csv file. Default True 
   step  float
        bin size. Default 0.1

HISTOGRAM_CENTROIDS
""""""""""""""""""""""""""""""
description
   Plots the histogram of one of the variables used for centroids computation corresponding to a particular hydrometeor type, the intermediate centroids and the final centroid
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L1319>`_
parameters
   voltype  str
        The name of the variable to plot. Can be dBZ, ZDR, KDP, RhoHV, H_ISO0 and its standardized form (e.g. dBZ_std) 
   hydro_type  str
        The name of the hydrometeor type. 
   write_data  Bool
        If true writes the histogram in a .csv file. Default True 
   step  float
        bin size. Default 0.1

HISTOGRAM2D_CENTROIDS
""""""""""""""""""""""""""""""
description
   Plots the 2D- histogram of two of the variables used for centroids computatio ncorresponding to a particular hydrometeor type, the intermediate centroids and the final centroid
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L1419>`_
parameters
   voltype_y  str
        The name of the variables to plot. Can be dBZ, ZDR, KDP, RhoHV, H_ISO0 and its standardized form (e.g. dBZ_std) 
   hydro_type  str
        The name of the hydrometeor type. step_x, 
   step_y  float
        bin size. Default 0.1

WRITE_CENTROIDS
""""""""""""""""""""""""""""""
description
   Writes the final centroids in a .csv file.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L1535>`_
parameters

SAVE_DATA
""""""""""""""""""""""""""""""
description
   Saves the data used to compute the centroids in an .npz file
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L1560>`_
parameters

COLOCATED_GATES
-----------------------------
COSMO_COORD
-----------------------------
SAVEVOL
""""""""""""""""""""""""""""""
description
   Save an object containing the index of the COSMO model grid that corresponds to each radar gate in a C/F radial file.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L201>`_
parameters
   file_type str
        The type of file used to save the data. Can be 'nc' or 'h5'. Default 'nc' 
   physical Bool
        If True the data will be saved in physical units (floats). Otherwise it will be quantized and saved as binary 
   compression str
        For ODIM file formats, the type of compression. Can be any of the allowed compression types for hdf5 files. Default gzip 
   compression_opts any
        The compression options allowed by the hdf5. Depends on the type of compression. Default 6 (The gzip compression level). 

COSMO2RADAR
-----------------------------
SAVEVOL
""""""""""""""""""""""""""""""
description
   Save an object containing the COSMO data in radar coordinatesin a C/F radial or ODIM file.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L291>`_
parameters
   file_type str
        The type of file used to save the data. Can be 'nc' or 'h5'. Default 'nc' 
   physical Bool
        If True the data will be saved in physical units (floats). Otherwise it will be quantized and saved as binary 
   compression str
        For ODIM file formats, the type of compression. Can be any of the allowed compression types for hdf5 files. Default gzip 
   compression_opts any
        The compression options allowed by the hdf5. Depends on the type of compression. Default 6 (The gzip compression level). All the products of the 'VOL' dataset group 

GRID
-----------------------------
CROSS_SECTION
""""""""""""""""""""""""""""""
description
   Plots a cross-section of gridded data
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_grid_products.py#L645>`_
parameters
   dict
        The two lat-lon coordinates marking the limits. They have the keywords 'lat' and 'lon' [degree]. The altitude limits are defined by the parameters in 'xsecImageConfig' in the 'loc' configuration file

HISTOGRAM
""""""""""""""""""""""""""""""
description
   Computes a histogram of the radar volum data
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_grid_products.py#L697>`_
parameters
   step float
        or None the data quantization step. If none it will be obtained from the Py-ART configuration file vmin, 
   vmax float
        or None The minimum and maximum values. If None they will be obtained from the Py-ART configuration file 
   mask_val float
        or None A value to mask. 
   write_data Bool
        If true the histogram data is written in a csv file

LATITUDE_SLICE
""""""""""""""""""""""""""""""
description
   Plots a cross-section of gridded data over a constant latitude.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_grid_products.py#L962>`_
parameters
   lat floats
        The starting point of the cross-section. The ending point is defined by the parameters in 'xsecImageConfig' in the 'loc' configuration file

LONGITUDE_SLICE
""""""""""""""""""""""""""""""
description
   Plots a cross-ection of gridded data over a constant longitude.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_grid_products.py#L897>`_
parameters
   lat floats
        The starting point of the cross-section. The ending point is defined by the parameters in 'xsecImageConfig' in the 'loc' configuration file

SAVEALL
""""""""""""""""""""""""""""""
description
   Saves a gridded data object including all or a list of user-defined fields in a netcdf file
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_grid_products.py#L813>`_
parameters
   datatypes list
        of str or None The list of data types to save. If it is None, all fields in the radar object will be saved

SAVEALL_GRID
""""""""""""""""""""""""""""""
description
   Same as before but can be used in a mixed GRID/VOL dataset, as there is no ambiguity with SAVEALL for VOL datasets
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_grid_products.py#L813>`_
parameters

SAVEVOL
""""""""""""""""""""""""""""""
description
   Saves on field of a gridded data object in a netcdf file.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_grid_products.py#L762>`_
parameters
   file_type str
        The type of file used to save the data. Can be 'nc' or 'h5'. Default 'nc' 
   physical Bool
        If True the data will be saved in physical units (floats). Otherwise it will be quantized and saved as binary. Default True 
   compression str
        For ODIM file formats, the type of compression. Can be any of the allowed compression types for hdf5 files. Default gzip 
   compression_opts any
        The compression options allowed by the hdf5. Depends on the type of compression. Default 6 (The gzip compression level).

SAVEVOL_GRID
""""""""""""""""""""""""""""""
description
   Same as before but can be used in a mixed GRID/VOL dataset, as there is no ambiguity with SAVEVOL for VOL datasets
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_grid_products.py#L762>`_
parameters

STATS
""""""""""""""""""""""""""""""
description
   Computes statistics over the whole images and stores them in a file.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_grid_products.py#L320>`_
parameters
   stat str
        The statistic used. Can be mean, median, min, max

SURFACE_RAW
""""""""""""""""""""""""""""""
description
   Plots a surface image of gridded data without projecting it into a map
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_grid_products.py#L369>`_
parameters
   level int
        The altitude level to plot. The rest of the parameters are defined by the parameters in 'ppiImageConfig' and 'ppiMapImageConfig' in the 'loc' configuration file

SURFACE_IMAGE
""""""""""""""""""""""""""""""
description
   Plots a surface image of gridded data.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_grid_products.py#L400>`_
parameters
   level int
        The altitude level to plot. The rest of the parameters are defined by the parameters in 'ppiImageConfig' and 'ppiMapImageConfig' in the 'loc' configuration file

SURFACE_CONTOUR
""""""""""""""""""""""""""""""
description
   Plots a surface image of contour gridded data.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_grid_products.py#L467>`_
parameters
   level int
        The altitude level to plot. The rest of the parameters are defined by the parameters in 'ppiImageConfig' and 'ppiMapImageConfig' in the 'loc' configuration file 
   contour_values  float
        array or None The contour values. If None the values are taken from the 'boundaries' keyword in the field description in the Py-ART config file. If 'boundaries' is not set the countours are 10 values linearly distributed from vmin to vmax 
   linewidths  float
        width of the contour lines 
   colors  color
        string or sequence of colors The contour colours

SURFACE_CONTOUR_OVERPLOT
""""""""""""""""""""""""""""""
description
   Plots a surface image of gridded data with a contour overplotted.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_grid_products.py#L467>`_
parameters
   level int
        The altitude level to plot. The rest of the parameters are defined by the parameters in 'ppiImageConfig' and 'ppiMapImageConfig' in the 'loc' configuration file 
   contour_values  float
        array or None The contour values. If None the values are taken from the 'boundaries' keyword in the field description in the Py-ART config file. If 'boundaries' is not set the countours are 10 values linearly distributed from vmin to vmax 
   linewidths  float
        width of the contour lines 
   colors  color
        string or sequence of colors The contour colours

SURFACE_OVERPLOT
""""""""""""""""""""""""""""""
description
   Plots on the same surface two images, one on top of the other.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_grid_products.py#L523>`_
parameters
   level int
        The altitude level to plot. The rest of the parameters are defined by the parameters in 'ppiImageConfig' and 'ppiMapImageConfig' in the 'loc' configuration file 
   contour_values  float
        array or None The contour values. If None the values are taken from the 'boundaries' keyword in the field description in the Py-ART config file. If 'boundaries' is not set the countours are 10 values linearly distributed from vmin to vmax

DDA_MAP
""""""""""""""""""""""""""""""
description
   Plots wind vectors obtained from a DDA analysis. The pyDDA package is required
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_grid_products.py#L850>`_
parameters
   level int
        The altitude level to plot. The rest of the parameters are defined by the parameters in 'ppiImageConfig' and 'ppiMapImageConfig' in the 'loc' configuration file 
   display_type str
        Display method for the wind vectors, can be either 'streamline', 'quiver' or 'barbs' 
   bg_ref_rad int
        Which radar to use as reference to display the background voltype. 
   u_vel_contours list
        of int The contours to use for plotting contours of u. Set to None to not display such contours. 
   v_vel_contours list
        of int The contours to use for plotting contours of v. Set to None to not display such contours. 
   w_vel_contours list
        of int The contours to use for plotting contours of w. Set to None to not display such contours. 
   vector_spacing_km float
        Spacing in km between wind vectors in x and y axis (only used for barbs and quiver plots) 
   quiver_len float
        Length to use for the quiver key in m/s. (only used for quiver plots) 
   streamline_arrowsize float
        Factor scale arrow size for streamlines. (only used for streamline plots) 
   linewidth float
        Linewidths for streamlines. (only used for streamline plots)

SPECTRA
-----------------------------
AMPLITUDE_PHASE_ANGLE_DOPPLER
""""""""""""""""""""""""""""""
description
   Makes an angle Doppler plot of complex spectra or IQ data. The plot can be along azimuth or along range. It is plotted separately the module and the phase of the signal.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_spectra_products.py#L959>`_
parameters
   along_azi  bool
        If true the plot is performed along azimuth, otherwise along elevation. Default true 
   ang  float
        The fixed angle (deg). Default 0. 
   rng  float
        The fixed range (m). Default 0. 
   ang_tol  float
        The fixed angle tolerance (deg). Default 1. 
   rng_tol  float
        The fixed range tolerance (m). Default 50. 
   xaxis_info  str
        The xaxis type. Can be 'Doppler_velocity', 'Doppler_frequency' or 'pulse_number' ampli_vmin, ampli_vmax, phase_vmin, 
   phase_vmax  float
        or None Minimum and maximum of the color scale for the module and phase

AMPLITUDE_PHASE_DOPPLER
""""""""""""""""""""""""""""""
description
   Plots a complex Doppler spectrum or IQ data making two separate plots for the module and phase of the signal
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_spectra_products.py#L818>`_
parameters
   rng  float
        azimuth and elevation (deg) and range (m) of the ray to plot azi_to, ele_tol, 
   rng_tol  float
        azimuth and elevation (deg) and range (m) tolerance respect to nominal position to plot. Default 1, 1, 50. ind_ray, 
   ind_rng  int
        index of the ray and range to plot. Alternative to defining its antenna coordinates 
   xaxis_info  str
        The xaxis type. Can be 'Doppler_velocity', 'Doppler_frequency' or 'pulse_number' ampli_vmin, ampli_vmax, phase_vmin, 
   phase_vmax  float
        or None Minimum and maximum of the color scale for the module and phase

AMPLITUDE_PHASE_RANGE_DOPPLER
""""""""""""""""""""""""""""""
description
   Plots a complex spectra or IQ data range-Doppler making two separate plots for the module and phase
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_spectra_products.py#L892>`_
parameters
   ele  float
        azimuth and elevation (deg) of the ray to plot azi_to, 
   ele_tol  float
        azimuth and elevation (deg) tolerance respect to nominal position to plot. Default 1, 1. 
   ind_ray  int
        index of the ray to plot. Alternative to defining its antenna coordinates 
   xaxis_info  str
        The xaxis type. Can be 'Doppler_velocity', 'Doppler_frequency' or 'pulse_number' ampli_vmin, ampli_vmax, phase_vmin, 
   phase_vmax  float
        or None Minimum and maximum of the color scale for the module and phase

AMPLITUDE_PHASE_TIME_DOPPLER
""""""""""""""""""""""""""""""
description
   Plots a complex spectra or IQ data time-Doppler making two separate plots for the module and phase of the signal
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_spectra_products.py#L1046>`_
parameters
   xaxis_info  str
        The xaxis type. Can be 'Doppler_velocity' or 'Doppler frequency' ampli_vmin, ampli_vmax, phase_vmin, 
   phase_vmax  float
        or None Minimum and maximum of the color scale for the module and phase 
   plot_type  str
        Can be 'final' or 'temporal'. If final the data is only plotted at the end of the processing

ANGLE_DOPPLER
""""""""""""""""""""""""""""""
description
   Makes an angle Doppler plot. The plot can be along azimuth or along range
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_spectra_products.py#L959>`_
parameters
   along_azi  bool
        If true the plot is performed along azimuth, otherwise along elevation. Default true 
   ang  float
        The fixed angle (deg). Default 0. 
   rng  float
        The fixed range (m). Default 0. 
   ang_tol  float
        The fixed angle tolerance (deg). Default 1. 
   rng_tol  float
        The fixed range tolerance (m). Default 50. 
   xaxis_info  str
        The xaxis type. Can be 'Doppler_velocity', 'Doppler_frequency' or 'pulse_number' vmin, 
   vmax  float
        or None Minimum and maximum of the color scale

COMPLEX_ANGLE_DOPPLER
""""""""""""""""""""""""""""""
description
   Makes an angle Doppler plot of complex spectra or IQ data. The plot can be along azimuth or along range. The real and imaginary parts are plotted separately
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_spectra_products.py#L608>`_
parameters
   along_azi  bool
        If true the plot is performed along azimuth, otherwise along elevation. Default true 
   ang  float
        The fixed angle (deg). Default 0. 
   rng  float
        The fixed range (m). Default 0. 
   ang_tol  float
        The fixed angle tolerance (deg). Default 1. 
   rng_tol  float
        The fixed range tolerance (m). Default 50. 
   xaxis_info  str
        The xaxis type. Can be 'Doppler_velocity', 'Doppler_frequency' or 'pulse_number' vmin, 
   vmax  float
        or None Minimum and maximum of the color scale

COMPLEX_DOPPLER
""""""""""""""""""""""""""""""
description
   Plots a complex Doppler spectrum or IQ data making two separate plots for the real and imaginary parts
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_spectra_products.py#L748>`_
parameters
   rng  float
        azimuth and elevation (deg) and range (m) of the ray to plot azi_to, ele_tol, 
   rng_tol  float
        azimuth and elevation (deg) and range (m) tolerance respect to nominal position to plot. Default 1, 1, 50. ind_ray, 
   ind_rng  int
        index of the ray and range to plot. Alternative to defining its antenna coordinates 
   xaxis_info  str
        The xaxis type. Can be 'Doppler_velocity', 'Doppler_frequency' or 'pulse_number' vmin, 
   vmax  float
        or None Minimum and maximum of the color scale

COMPLEX_RANGE_DOPPLER
""""""""""""""""""""""""""""""
description
   Plots the complex spectra or IQ data range-Doppler making two separate plots for the real and imaginary parts
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_spectra_products.py#L547>`_
parameters
   ele  float
        azimuth and elevation (deg) of the ray to plot azi_to, 
   ele_tol  float
        azimuth and elevation (deg) tolerance respect to nominal position to plot. Default 1, 1. 
   ind_ray  int
        index of the ray to plot. Alternative to defining its antenna coordinates 
   xaxis_info  str
        The xaxis type. Can be 'Doppler_velocity', 'Doppler_frequency' or 'pulse_number' vmin, 
   vmax  float
        or None Minimum and maximum of the color scale

COMPLEX_TIME_DOPPLER
""""""""""""""""""""""""""""""
description
   Plots the complex spectra or IQ data time-Doppler making two separate plots for the real and imaginary parts
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_spectra_products.py#L689>`_
parameters
   xaxis_info  str
        The xaxis type. Can be 'Doppler_velocity' or 'Doppler frequency' vmin, 
   vmax  float
        or None Minimum and maximum of the color scale 
   plot_type  str
        Can be 'final' or 'temporal'. If final the data is only plotted at the end of the processing

DOPPLER
""""""""""""""""""""""""""""""
description
   Plots a Doppler spectrum variable or IQ data variable
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_spectra_products.py#L1046>`_
parameters
   rng  float
        azimuth and elevation (deg) and range (m) of the ray to plot azi_to, ele_tol, 
   rng_tol  float
        azimuth and elevation (deg) and range (m) tolerance respect to nominal position to plot. Default 1, 1, 50. ind_ray, 
   ind_rng  int
        index of the ray and range to plot. Alternative to defining its antenna coordinates 
   xaxis_info  str
        The xaxis type. Can be 'Doppler_velocity', 'Doppler_frequency' or 'pulse_number' vmin, 
   vmax  float
        or None Minimum and maximum of the color scale

RANGE_DOPPLER
""""""""""""""""""""""""""""""
description
   Makes a range-Doppler plot of spectral or IQ data
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_spectra_products.py#L892>`_
parameters
   ele  float
        azimuth and elevation (deg) of the ray to plot azi_to, 
   ele_tol  float
        azimuth and elevation (deg) tolerance respect to nominal position to plot. Default 1, 1. 
   ind_ray  int
        index of the ray to plot. Alternative to defining its antenna coordinates 
   xaxis_info  str
        The xaxis type. Can be 'Doppler_velocity', 'Doppler_frequency' or 'pulse_number' vmin, 
   vmax  float
        or None Minimum and maximum of the color scale

SAVEALL
""""""""""""""""""""""""""""""
description
   Saves radar spectra or IQ volume data including all or a list of userdefined fields in a netcdf file
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_spectra_products.py#L1144>`_
parameters
   datatypes list
        of str or None The list of data types to save. If it is None, all fields in the radar object will be saved 
   physical Bool
        If True the data will be saved in physical units (floats). Otherwise it will be quantized and saved as binary

SAVEVOL
""""""""""""""""""""""""""""""
description
   Saves one field of a radar spectra or IQ volume data in a netcdf file
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_spectra_products.py#L1111>`_
parameters
   physical Bool
        If True the data will be saved in physical units (floats). Otherwise it will be quantized and saved as binary

TIME_DOPPLER
""""""""""""""""""""""""""""""
description
   Makes a time-Doppler plot of spectral or IQ data at a point of interest.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_spectra_products.py#L1046>`_
parameters
   xaxis_info  str
        The xaxis type. Can be 'Doppler_velocity', 'Doppler_frequency' or 'pulse_number' vmin, 
   vmax  float
        or None Minimum and maximum of the color scale 
   plot_type  str
        Can be 'final' or 'temporal'. If final the data is only plotted at the end of the processing 

GRID_TIMEAVG
-----------------------------
INTERCOMP
-----------------------------
PLOT_AND_WRITE_INTERCOMP_TS
""""""""""""""""""""""""""""""
description
   Writes statistics of radar intercomparison in a file and plots the time series of the statistics.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_intercomp_products.py#L214>`_
parameters
   Bool
        If true adds the year in the csv file containing the statistics. Default False 'sort_by_date': Bool If true sorts the statistics by date when reading the csv file containing the statistics. Default False 'rewrite': Bool If true rewrites the csv file containing the statistics. Default False 'npoints_min'
   int
        The minimum number of points to consider the statistics valid and therefore use the data point in the plotting. Default 0 'corr_min'
   float
        The minimum correlation to consider the statistics valid and therefore use the data point in the plotting. Default 0.

PLOT_SCATTER_INTERCOMP
""""""""""""""""""""""""""""""
description
   Plots a density plot with the points of radar 1 versus the points of radar 2
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_intercomp_products.py#L136>`_
parameters
   float
        The quantization step of the data. If none it will be computed using the Py-ART config file. Default None 'scatter_type'
   str
        Type of scatter plot. Can be a plot for each radar volume ('instant') or at the end of the processing period ('cumulative'). Default is 'cumulative'

WRITE_INTERCOMP
""""""""""""""""""""""""""""""
description
   Writes the instantaneously intercompared data (gate positions, values, etc.) in a csv file.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_intercomp_products.py#L214>`_
parameters

ML
-----------------------------
ML_TS
""""""""""""""""""""""""""""""
description
   Plots and writes a time series of the melting layer, i.e. the evolution of the average and standard deviation of the melting layer top and thickness and the the number of rays used in the retrieval.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L804>`_
parameters
   dpi int
        The pixel density of the plot. Default 72

SAVE_ML
""""""""""""""""""""""""""""""
description
   Saves an object containing the best estimate of the melting layer retrieved in a C/F radial file
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L861>`_
parameters

MONITORING
-----------------------------
ANGULAR_DENSITY
""""""""""""""""""""""""""""""
description
   For a specified elevation angle, plots a 2D histogram with the azimuth angle in the X-axis and the data values in the Y-axis. The reference values and the user defined quantiles are also plot on the same figure
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_monitoring_products.py#L310>`_
parameters
   anglenr int
        The elevation angle number to plot 
   quantiles list
        of floats The quantiles to plot. Default 25., 50., 75. 
   ref_value float
        The reference value vmin, 
   vmax  floats
        or None The minimum and maximum values of the data points. If not specified they are obtained from the Py-ART config file

CUMUL_VOL_TS
""""""""""""""""""""""""""""""
description
   Plots time series of the average of instantaneous quantiles stored in a csv file.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_monitoring_products.py#L560>`_
parameters
   quantiles list
        of 3 floats the quantiles to compute. Default 25., 50., 75. 
   ref_value float
        The reference value. Default 0 
   sort_by_date Bool
        If true when reading the csv file containing the statistics the data is sorted by date. Default False 
   rewrite Bool
        If true the csv file containing the statistics is rewritten 
   add_data_in_fname Bool
        If true and the data used is cumulative the year is written in the csv file name and the plot file name 
   npoints_min int
        Minimum number of points to use the data point in the plotting and to send an alarm. Default 0 vmin, 
   vmax float
        or None Limits of the Y-axis (data value). If None the limits are obtained from the Py-ART config file 
   alarm Bool
        If true an alarm is sent 
   tol_abs float
        Margin of tolerance from the reference value. If the current value is above this margin an alarm is sent. If the margin is not specified it is not possible to send any alarm 
   tol_trend float
        Margin of tolerance from the reference value. If the trend of the last X events is above this margin an alarm is sent. If the margin is not specified it is not possible to send any alarm 
   plot_until_year_end Bool
        If true will always set the xmax of the timeseries to the end of the current year 
   nevents_min int
        Minimum number of events with sufficient points to send an alarm related to the trend. If not specified it is not possible to send any alarm 
   sender str
        The mail of the alarm sender. If not specified it is not possible to send any alarm 
   receiver_list list
        of str The list of emails of the people that will receive the alarm.. If not specified it is not possible to send any alarm

PPI_HISTOGRAM
""""""""""""""""""""""""""""""
description
   Plots a histogram of data at a particular elevation angle.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_monitoring_products.py#L257>`_
parameters
   anglenr int
        The elevation angle number to plot

SAVEVOL
""""""""""""""""""""""""""""""
description
   Saves the monitoring data in a C/F radar file. The data field contains histograms of data for each pair of azimuth and elevation angles
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_monitoring_products.py#L781>`_
parameters

VOL_HISTOGRAM
""""""""""""""""""""""""""""""
description
   Plots a histogram of data collected from all the radar volume.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_monitoring_products.py#L195>`_
parameters
   write_data bool
        If true the resultant histogram is also saved in a csv file. Default True.

VOL_TS
""""""""""""""""""""""""""""""
description
   Computes statistics of the gathered data and writes them in a csv file and plots a time series of those statistics.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_monitoring_products.py#L560>`_
parameters
   quantiles list
        of 3 floats the quantiles to compute. Default 25., 50., 75. 
   ref_value float
        The reference value. Default 0 
   sort_by_date Bool
        If true when reading the csv file containing the statistics the data is sorted by date. Default False 
   rewrite Bool
        If true the csv file containing the statistics is rewritten 
   add_data_in_fname Bool
        If true and the data used is cumulative the year is written in the csv file name and the plot file name 
   npoints_min int
        Minimum number of points to use the data point in the plotting and to send an alarm. Default 0 vmin, 
   vmax float
        or None Limits of the Y-axis (data value). If None the limits are obtained from the Py-ART config file 
   alarm Bool
        If true an alarm is sent 
   tol_abs float
        Margin of tolerance from the reference value. If the current value is above this margin an alarm is sent. If the margin is not specified it is not possible to send any alarm 
   tol_trend float
        Margin of tolerance from the reference value. If the trend of the last X events is above this margin an alarm is sent. If the margin is not specified it is not possible to send any alarm 
   plot_until_year_end Bool
        If true will always set the xmax of the timeseries to the end of the current year 
   nevents_min int
        Minimum number of events with sufficient points to send an alarm related to the trend. If not specified it is not possible to send any alarm 
   sender str
        The mail of the alarm sender. If not specified it is not possible to send any alarm 
   receiver_list list
        of str The list of emails of the people that will receive the alarm.. If not specified it is not possible to send any alarm 

OCCURRENCE
-----------------------------
WRITE_EXCESS_GATES
""""""""""""""""""""""""""""""
description
   Write the data that identifies radar gates with clutter that has a frequency of occurrence above a certain threshold.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L94>`_
parameters
   quant_min float
        Minimum frequency of occurrence in percentage to keep the gate as valid. Default 95. All the products of the 'VOL' dataset group 

QVP
-----------------------------
SPARSE_GRID
-----------------------------
SURFACE_IMAGE
""""""""""""""""""""""""""""""
description
   Generates a surface image
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_grid_products.py#L90>`_
parameters
   list
        of floats The limits of the surface to plot [deg] lon0, lon1, lat0, lat1 

SUN_HITS
-----------------------------
PLOT_SUN_HITS
""""""""""""""""""""""""""""""
description
   Plots in a sun-radar azimuth difference-sun-radar elevation difference grid the values of all sun hits obtained during the processing period
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L416>`_
parameters

PLOT_SUN_RETRIEVAL
""""""""""""""""""""""""""""""
description
   Plots in a sun-radar azimuth difference-sun- radar elevation difference grid the retrieved sun pattern
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L536>`_
parameters

PLOT_SUN_RETRIEVAL_TS
""""""""""""""""""""""""""""""
description
   Plots time series of the retrieved sun pattern parameters
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L536>`_
parameters
   dpi int
        The pixel density of the plot. Default 72 
   add_date_in_fname Bool
        If true the year is added in the plot file name

PLOT_SUNSCAN
""""""""""""""""""""""""""""""
description
   Plots a constant range radar azimuth-elevation of the sunscan field data
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L684>`_
parameters

WRITE_SUN_HITS
""""""""""""""""""""""""""""""
description
   Writes the information concerning possible sun hits in a csv file
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L396>`_
parameters

WRITE_SUN_RETRIEVAL
""""""""""""""""""""""""""""""
description
   Writes the retrieved sun pattern parameters in a csv file.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L461>`_
parameters
   add_date_in_fname Bool
        If true the year is added in the csv file name

TIMEAVG
-----------------------------
TIMESERIES
-----------------------------
COMPARE_CUMULATIVE_POINT
""""""""""""""""""""""""""""""
description
   Plots in the same graph 2 time series of data accumulation (tipically rainfall rate). One time series is a point measurement of radar data while the other is from a co-located instrument (rain gauge or disdrometer)
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_timeseries_products.py#L472>`_
parameters
   dpi int
        The pixel density of the plot. Default 72 vmin, 
   vmax float
        The limits of the Y-axis. If none they will be obtained from the Py-ART config file. 
   sensor str
        The sensor type. Can be 'rgage' or 'disdro' 
   sensorid str
        The sensor ID. 
   location str
        A string identifying the location of the disdrometer 
   freq float
        The frequency used to retrieve the polarimetric variables of a disdrometer 
   ele float
        The elevation angle used to retrieve the polarimetric variables of a disdrometer 
   ScanPeriod float
        The scaning period of the radar in seconds. This parameter is defined in the 'loc' config file

COMPARE_POINT
""""""""""""""""""""""""""""""
description
   Plots in the same graph 2 time series of data . One time series is a point measurement of radar data while the other is from a co-located instrument (rain gauge or disdrometer)
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_timeseries_products.py#L381>`_
parameters
   dpi int
        The pixel density of the plot. Default 72 vmin, 
   vmax float
        The limits of the Y-axis. If none they will be obtained from the Py-ART config file. 
   sensor str
        The sensor type. Can be 'rgage' or 'disdro' 
   sensorid str
        The sensor ID. 
   location str
        A string identifying the location of the disdrometer 
   freq float
        The frequency used to retrieve the polarimetric variables of a disdrometer 
   ele float
        The elevation angle used to retrieve the polarimetric variables of a disdrometer

COMPARE_TIME_AVG
""""""""""""""""""""""""""""""
description
   Creates a scatter plot of average radar data versus average sensor data.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_timeseries_products.py#L566>`_
parameters
   dpi int
        The pixel density of the plot. Default 72 
   sensor str
        The sensor type. Can be 'rgage' or 'disdro' 
   sensorid str
        The sensor ID. 
   location str
        A string identifying the location of the disdrometer 
   freq float
        The frequency used to retrieve the polarimetric variables of a disdrometer 
   ele float
        The elevation angle used to retrieve the polarimetric variables of a disdrometer 
   cum_time float
        Data accumulation time [s]. Default 3600. 
   base_time float
        Starting moment of the accumulation [s from midnight]. Default 0.

PLOT_AND_WRITE
""""""""""""""""""""""""""""""
description
   Writes and plots a trajectory time series.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_timeseries_products.py#L674>`_
parameters
   ymax float
        The minimum and maximum value of the Y-axis. If none it will be obtained from the Py-ART config file.

PLOT_AND_WRITE_POINT
""""""""""""""""""""""""""""""
description
   Plots and writes a time series of radar data at a particular point
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_timeseries_products.py#L214>`_
parameters
   dpi int
        The pixel density of the plot. Default 72 vmin, 
   vmax float
        The limits of the Y-axis. If none they will be obtained from the Py-ART config file.

PLOT_CUMULATIVE_POINT
""""""""""""""""""""""""""""""
description
   Plots a time series of radar data accumulation at a particular point.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_timeseries_products.py#L301>`_
parameters
   dpi int
        The pixel density of the plot. Default 72 vmin, 
   vmax float
        The limits of the Y-axis. If none they will be obtained from the Py-ART config file. 
   ScanPeriod float
        The scaning period of the radar in seconds. This parameter is defined in the 'loc' config file

PLOT_HIST
""""""""""""""""""""""""""""""
description
   plots and writes a histogram of all the data gathered during the trajectory processing
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_timeseries_products.py#L715>`_
parameters
   step float
        or None The quantization step of the data. If None it will be obtained from the Py-ART config file

TRAJ_CAPPI_IMAGE
""""""""""""""""""""""""""""""
description
   Creates a CAPPI image with the trajectory position overplot on it.
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_timeseries_products.py#L748>`_
parameters
   color_ref str
        The meaning of the color code with which the trajectory is plotted. Can be 'None', 'altitude' (the absolute altitude), 'rel_altitude' (altitude relative to the CAPPI altitude), 'time' (trajectory time respect of the start of the radar scan leading to the CAPPI) 
   altitude float
        The CAPPI altitude [m] 
   wfunc str
        Function used in the gridding of the radar data. The function types are defined in pyart.map.grid_from_radars. Default 'NEAREST' 
   res float
        The CAPPI resolution [m]. Default 500.

WRITE_MULTIPLE_POINTS
""""""""""""""""""""""""""""""
description
   Writes a csv file containing data from multiple points obtained from a radar object
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_timeseries_products.py#L198>`_
parameters

TRAJ_ONLY
-----------------------------
TRAJ_MAP
""""""""""""""""""""""""""""""
description
   Plots the trajectory on a lat-lon map with the altitude color coded
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_traj_products.py#L136>`_
parameters

TRAJ_PLOT
""""""""""""""""""""""""""""""
description
   Plots time series of the trajectory respect to the radar elevation, azimuth or range
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_traj_products.py#L50>`_
parameters
   str
        The type of 
   parameter
       'EL', 'AZ', or 'RANGE'

VPR
-----------------------------
PLOT_VPR_THEO
""""""""""""""""""""""""""""""
description
   Plots and writes the retrieved theoretical VPR
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L932>`_
parameters
   dpi int
        The pixel density of the plot. Default 72

WRITE_VPR_THEO_PARAMS
""""""""""""""""""""""""""""""
description
   Writes the parameters of the theoretical VPR in a text file
 `[Source] <https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad/prod/process_product.py#L972>`_
parameters
