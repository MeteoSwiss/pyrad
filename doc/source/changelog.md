# Changelog

## v1.8.9
**Bug corrections:**
- [FIX: fix np.Inf for numpy 2](https://github.com/MeteoSwiss/pyrad/commit/f0d3bf32c7dbb3ca9c8eb13708957072af872dd5)
- [fix conflict in src/pyrad_proc/pyrad/flow/flow_aux.py](https://github.com/MeteoSwiss/pyrad/commit/64c6c3064996e6be7682b34b27df9772692b896b)
- [revert of](https://github.com/MeteoSwiss/pyrad/commit/3955a619d1145b5d307bdc24742ab54215ff8398) https://github.com/MeteoSwiss/pyrad/commit/41382257dc5578d5cb121de11b3766b70e4c1878 [that removed by mistake call to read_smn2](https://github.com/MeteoSwiss/pyrad/commit/3955a619d1145b5d307bdc24742ab54215ff8398)
- [[BUG] fix in process_monitoring_products.py to allow copy to S3 of VOL_HISTOGRAM imgs](https://github.com/MeteoSwiss/pyrad/commit/33581ced90e5c45705fcba99328d4e76be5b4c40)
- [added coercion of product names to list in flow_aux for s3 writing](https://github.com/MeteoSwiss/pyrad/commit/369048bff7cdc275d541feb6aaff5a13471a1ded)
- [added support of visibility_polar varname in process_intercomp (GECSX format)](https://github.com/MeteoSwiss/pyrad/commit/e6adbba15df72ab6c6a782d2d87687b0c1aea7d8)
- [fixed issue in process_traj.py that lead to some timesteps not being considered](https://github.com/MeteoSwiss/pyrad/commit/d83c14064993203f80b5f204d823cca10d88d3bb)
- [provided alternative to fcntl for windows](https://github.com/MeteoSwiss/pyrad/commit/1fbc1f0786726bf1b68637d45e2142396340b689)

**New additions**
- [added possibility to copy saved products to S3 storage](https://github.com/MeteoSwiss/pyrad/commit/2d54aca2896dd02ab2f19b0bdb251f1f8982daa9)
- [modified function merge_scans_mfcfradial in read_data_radar.py to allow ingesting CFRadial files containing individual parameters and scans](https://github.com/MeteoSwiss/pyrad/commit/3050ffe8590206c4fb7776cfcd10fc497449bb3b)
- [minor changes to be able to read better SkyEcho files](https://github.com/MeteoSwiss/pyrad/commit/49dab5b2ffb2d5e4ea86308431e446a981743baf)
- [added spectraImageConfig struct in loc cfg file to specify spectral image parameters](https://github.com/MeteoSwiss/pyrad/commit/5cc42995f2e0791f2cc5b25a18239113a927d495)
- [function get_datatype_skyecho made available](https://github.com/MeteoSwiss/pyrad/commit/f49814667748f45f3cf52642b70d4e48e5a33fd3)
- [Added the possibility to read rain gauge data with reader read_smn2](https://github.com/MeteoSwiss/pyrad/commit/2e89dea47483b95299ac0dd1b1c5a27762fe47b2)
- [Added possibility to use Gate Filter in LP phase processing. Added user-defined settings](https://github.com/MeteoSwiss/pyrad/commit/fbba560740ed62de7e6119e2d103a090f0caa39c)
- [added new SkyEcho variables](https://github.com/MeteoSwiss/pyrad/commit/a83d1859ead9c1e315e06fd2892215dac14ab34d)

## v1.8.8
**Bug corrections**

- [bug fix for crash when temp_ref is not defined](https://github.com/MeteoSwiss/pyrad/commit/8877fc18ec598a16aa472ea4f2706c8309861571)
- [fix deprecation warning in get_cmap](https://github.com/MeteoSwiss/pyrad/commit/ca965ce702d80dcd35d44f90f04e7931c8fd3f30)


## v1.8.7
**Bug fixes:**
- [bug correction in function to read MFCFRADIAL data](https://github.com/MeteoSwiss/pyrad/commit/dfacd71c301d3ee37e7c6c002b2b54bc347a2480)
- [fix bug due to missing fzl in some routines](https://github.com/MeteoSwiss/pyrad/commit/02aae44c72d5bb7c2e140a68504c8e00066309d0)
- [fix for missing box variable in process_roi2](https://github.com/MeteoSwiss/pyrad/commit/997db1f956aad6e348f1d92ae75882e7296fc22f)
- [Fix in process_grid_products.py for wrong warning message when bg DDA field is missing](https://github.com/MeteoSwiss/pyrad/commit/a8ee6c7d96da3932724f5e5775482ca5289f326d)
- [bug correction in read_knmi](https://github.com/MeteoSwiss/pyrad/commit/af5a2e2acfb7f6f4e0fd76f643860f04a55acb80)
- [bug correction in function process_turbulence in process_Doppler.py](https://github.com/MeteoSwiss/pyrad/commit/ed44ae9fa062df8a8c868c10abdf8081590c5d24)
- [fixed old pandas command in read_radiosounding](https://github.com/MeteoSwiss/pyrad/commit/a3ef19b556994bf1e2b792e360e1eb2dfdbd7908)

**New additions**
- [modified function merge_scans_mfcfradial in read_data_radar.py to allow ingesting CFRadial files containing individual parameters and scans](https://github.com/MeteoSwiss/pyrad/commit/3050ffe8590206c4fb7776cfcd10fc497449bb3b)
- [ingestion of KNMI H5 radar data files and minor adaptations](https://github.com/MeteoSwiss/pyrad/commit/ad8bba2367479ed54d2b4f84aafda34dcf42e62a)
- [add possibility to use sounding data for fzl estimation](https://github.com/MeteoSwiss/pyrad/commit/59311182cabd0a79352964b7547a75ec8374d733)
- [added function to mask all data not belonging to a ROI. Added rectangular ROI.](https://github.com/MeteoSwiss/pyrad/commit/2555e68ac7c3e67ab468dff0a243fd9e11c89232)
- [modified function read_knmi of read_data_sensor.py](https://github.com/MeteoSwiss/pyrad/commit/4f43bf9051bdfd082d1b11299b8866cb5a2d98e3)
- [added IGRA sounding reader + support for IGRA sounding in DDA wrapper](https://github.com/MeteoSwiss/pyrad/commit/fe98f900d30b67c3e9d50fa277f34428079e46c5)
- [modified function read_knmi of read_data_sensor.py](https://github.com/MeteoSwiss/pyrad/commit/4f43bf9051bdfd082d1b11299b8866cb5a2d98e3)
- [new functions to read parsivel disdrometer data and compare to radar data](https://github.com/MeteoSwiss/pyrad/commit/28bb465287017ca28e9f146df4aa8c83cad5bb50)
- [ingestion of data from SkyEcho proprietary netcdf files](https://github.com/MeteoSwiss/pyrad/commit/45d3cb2873d72f2cb1fa6e25449dcf98176b7c30)
- [adaptation of turbulence retrieval to accept lidar data](https://github.com/MeteoSwiss/pyrad/commit/76673e0837f8c323c32e365361c064a666397d15)

## v1.8.6

**Bug corrections**

- [fixed old pandas command in read_radiosounding](https://github.com/MeteoSwiss/pyrad/commit/a3ef19b556994bf1e2b792e360e1eb2dfdbd7908)

## v1.8.5

**Bug corrections**

- [added missing functions in __init__.py files](https://github.com/MeteoSwiss/pyrad/commit/b8404ca22550e331495e520446ebb0e61959ca76)

## v1.8.4

**Bug corrections**

- [fix of a bug where vmin and vmax were not applied to lat/lon grid profiles](https://github.com/MeteoSwiss/pyrad/commit/31cc9d2a3dd915a3dbe05f9aade9f7f5dbc50ecb) 
- [ fix in plots_grid.py for barbs DDA plots](https://github.com/MeteoSwiss/pyrad/commit/44bd2310fe683541d34f32c0762c9633203df73a)
- [fix for file locking issue on NFS systems for timeseries csv files](https://github.com/MeteoSwiss/pyrad/commit/95cd7e19059b4b6eb85cf67fb45287808e229fa8)
- [deprecation fixes for latest numpy version](https://github.com/MeteoSwiss/pyrad/commit/e5865f3f05c96762e8f0917018e12e7644e380dd)

**New additions**

- [added support for multi radar in rt processing](https://github.com/MeteoSwiss/pyrad/commit/861bb8b015d9417767216980450b779c1f909643)
- [Added new Product WINDSHEAR_LIDAR](https://github.com/MeteoSwiss/pyrad/commit/93745298fc42afa2c2c88d78af4a673cb5aa61b2)
- [made pyDDA wrapper compatible with new pyDDA version](https://github.com/MeteoSwiss/pyrad/commit/e1bac0f0b07f2b7e7e0f11a2e6f0acce22c6ce3d)
- [added GATEFILTER dataset which reproduces Py-ART gatefilter in pyrad + fixes in DDA processing](https://github.com/MeteoSwiss/pyrad/commit/90641235e496fe6bf3928f9502e1cb9d24c2054f)

## v1.8.3

**Bug corrections**

- [various fixes to be able to use GECSX products as visibility inputs](https://github.com/MeteoSwiss/pyrad/commit/ef4b0fbd99d1dadffff1f8f96aabd8dae2ec46ff)
- [fix in merge_cfradial2 for RADARV convention](https://github.com/MeteoSwiss/pyrad/commit/a540c0133246e6385c7919d2545ecf3c30f8a527)

**New additions**

- [added time tolerance when merging cfradial2 scans with RADARV convention](https://github.com/MeteoSwiss/pyrad/commit/3b60e3129cf100ae72e3854fbba4efd2d547bc54)
- [added possibility to suppress warnings in rt processing](https://github.com/MeteoSwiss/pyrad/commit/e48ef9e9e959ad0ae13699fe84ed4f264e48f7d6)

## v1.8.2

**New additions:**

- [added possibility to suppress warnings in rt processing](https://github.com/MeteoSwiss/pyrad/commit/e48ef9e9e959ad0ae13699fe84ed4f264e48f7d6)

## v1.8.1

**Bug fixes:**

- [fixed issue with list instead of dict in merge_scans_cfradial2](https://github.com/MeteoSwiss/pyrad/commit/a1155f0dbf3cf2d7eeae72014fd7d0312d6968fc)
- [bug correction in overplotting of grid data. Improvements in the reading of user-defined parameters](https://github.com/MeteoSwiss/pyrad/commit/cadff8eed43b9b4d77d9e5fa45343cf1da305313)

## v1.8

**Bug fixes:**
 - [correction of a bug in function process_vol2grid. Added xsecImageConfig](https://github.com/MeteoSwiss/pyrad/commit/e7eaac19f03d4e3ff6bb6fc7bf4afefb3a0250aa)
- [redesigned dem readers so they output proj read from file if possible](https://github.com/MeteoSwiss/pyrad/commit/6b3f78df7d4b68f74dbc06a6ec472b4e4233c2e3)


**New additions:**
- [added DDA process and DDA_LONGITUDE_SLICE DDA_LATITUDE_SLICE DDA_MAP products](https://github.com/MeteoSwiss/pyrad/commit/5bbb84d9f752a8466ce450bc9cc378366cad8e93) 
- [added possibility to specify selfconsistencypath + added poor support of NFS mounts in file locking](https://github.com/MeteoSwiss/pyrad/commit/f05f3d7a9772383c5bf8ba7edda8372239267c15) 
- [added exact_limits keyword in loc file to display exactly user-specified lon_lines and lat_lines](https://github.com/MeteoSwiss/pyrad/commit/1bd4b4789193894456a2b2c14a62b53ecc5cdbe7) 
- [added info on searched file directory if no file was found](https://github.com/MeteoSwiss/pyrad/commit/ce67680ee9faee6b466e16d95382543e42ff0285)
- [added RADARV path convention for windcube lidar data (CFRADIAL2)](https://github.com/MeteoSwiss/pyrad/commit/052ea9a228a44b9a7f7e48e364083264ff4e8808)
- [added script under ci that generates pyrad var name mappings references under doc/mappings](https://github.com/MeteoSwiss/pyrad/commit/837c2622da247ddd75aadab9bdf948afb01423c3) 
- [added possibility to specify nyquist velo in process_dealias_region_based](https://github.com/MeteoSwiss/pyrad/commit/522ab3e056b5891846da1827a0337efbc025d807) 
- [added the GECSX datatype to read pyrad GECSX input for visibility filter](https://github.com/MeteoSwiss/pyrad/commit/6fa2d8fd4c028dea0bc1eb6d8bef85a7a164a1e8)

## v1.7

**New additions:**
-  New polar product: _CROSS_SECTION_ in _RadarDisplay_ which can be used to display a cross-section of polar data between arbitrary coordinates
- Added possibility to overplot of a trajectory in a PPI map
-  Output data from regions of interest in a radar volume as kml and csv files
- Added some output on the location where the parsing of the config files fails (it is fails)


## v1.6.2

**New additions:**

* Added the possibility to create Cartesian ODIM outputs with the ODIMPYRADGRID keyword

## v1.6

**Bug fixes:**
 - changed the number of colors from cmap.N to len(boundaries) – 1 in get_cmap to avoid colorbar issues with discrete colorbars
 - fixed a bug that made it impossible to use any other method than main in flow_control.py 
 - set xlim and ylim of grid surface plots to the min/max of user defined lat/lon (this step was missing since the last release)
  - correction of a bug that did not allow import of pyrad when the metranet library was mssing
 
**New additions:**
- change in function process_point_measurement: When Truealt is false input keyword alt is not required and the output alt will be determined by the elevation angle
 - added new user-defined parameters in VPR correction scheme
-  modifications to allow to read multiple ODIM files containing parameters and scan
-  changed keyword add_lines to add_grid_lines to be consistent with Py-ART
-  new functions to accumulate gridded data and get data at multiple points in a grid
- added ODIMGRID keyword to be able to ingest gridded data from ODIM files
- added keyword MasterScanTimeTol that allows to combine ODIM files from different scans not having the same nominal time in file name

