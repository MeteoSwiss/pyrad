The main configuration file is used to define the global settings, notably the paths to the different sources of data. The parameters of the main configuration file are described in Table 2.

Table 2: Configuration parameters of the main configuration file

====================  =======  =======================================================================================
Name                  Type     Description
====================  =======  =======================================================================================
name                  STRING   Name of the data processing. This name is used in the path of the saved products in the following manner:
                               ``<saveimgbasepath>/<name>/<YYYY-MM-DD>/<datasetname>/<prodname>/<outputname>``
datapath              STRING   Base directory of the rainbow raw data. This field must have a trailing '/'. The raw data files of a scan can be found using the following file path:
                               ``<datapath>/<scanname>/<YYYY-MM-DD>/<YYYYMMDDHHMMSS00datatype>.<ext>``
configpath            STRING   Base directory of the configuration files. This directory contains clutter maps, filter coefficients, antenna pattern, and the data processing configuration files.
locationConfigFile    STRING   File name (with full path) of the location configuration file. Described in Section 3.2.
productConfigFile     STRING   File name (with full path) of the product configuration file. Described in Section 4.
lastStateFile         STRING   File name (with full path) of the file containing the time of the last processed scan. Used in particular for real-time processing.
imgformat             STRING/  File format(s) of the images. The following formats are supported: eps, png, and jpg. If ``saveimg`` is set to 0, this field is not used.
                      STRARR   
saveimgbasepath       STRING   Base directory for the images to save. The directory structure looks as follows:
                               ``<saveimgbasepath>/<name>/<YYYY-MM-DD>/<datasetname>/<prodname>/<outputname>``
                               If ``saveimg`` is set to 0, this field is not used.
s3copypath			  STRING   OPTIONAL. Path to an S3 bucket. If provided all generated products will be written there as well using the same data structure. The format must be
							   https://bucket_name.endpoint.domain for example https://tests.fr-par-1.linodeobjects.com/. The S3 copy procedure will only work if the environment variables
							   AWS_KEY and AWS_SECRET are defined in the pyrad scope. AWS_KEY contains the S3 bucket AWS key and AWS_SECRET the associated secret.
loadbasepath          STRING   OPTIONAL. Base path of saved data. By default, this field is set to ``saveimgbasepath``.
loadname              STRING   OPTIONAL. Name of the saved data processing. Used for saved volume loading. By default, this field is set to ``name``.
dempath               STRING   OPTIONAL. Base directory of the Digital Elevation Model (DEM) files. Basically to load the radar visibility (Optional).
smnpath               STRING   OPTIONAL. Base directory of the SwissMetNet stations data. Used in the comparison between radar data and rain gauges (Optional).
disdropath            STRING   OPTIONAL. Base directory of the disdrometer data. Used in the comparison between radar data and disdrometers (Optional).
solarfluxpath         STRING   OPTIONAL. Base directory of the solar flux data. Used to plot the calibration bias based on sun monitoring (Optional).
cosmopath             STRING   OPTIONAL. Base directory of the COSMO data files.
====================  =======  =======================================================================================

