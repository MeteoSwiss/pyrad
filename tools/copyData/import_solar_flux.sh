#!/bin/bash

# set permits
umask 0002

# import solar flux data from DRAO
wget -q --tries=2 --timeout=5 --spider ftp://ftp:ftp@ftp.geolab.nrcan.gc.ca
if [[ $? -eq 0 ]]; then
    cd /srn/analysis/solar_flux/
    rm fluxtable.txt
    wget ftp://ftp:ftp@ftp.geolab.nrcan.gc.ca/data/solar_flux/daily_flux_values/fluxtable.txt
else
   echo "Unable to retrieve solar flux. Server off-line"
fi