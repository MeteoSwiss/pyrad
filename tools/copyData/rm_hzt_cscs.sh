#!/bin/bash

# This script gets the HZT data from the CSCS
# repository and puts it in the right folder.
# To be run in CSCS

# set permits
umask 0002

# Config
data_destbase=/store/msrad/cosmo/cosmo1/HZT/
while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
        -d|--day)
        DAY="$2"
        OIFS=$IFS
        IFS=','
        read -r -a day_vec <<< "$DAY"
        IFS=$OIFS
        shift # past argument
        ;;
        -p|--dest_base)
        data_destbase="$2"
        shift # past argument        
        ;;    
    esac
    shift # past argument or value
done

nday=${#day_vec[@]}
for ((iday=0; iday<${nday}; iday++)); do    
    years=$(date --date "${day_vec[${iday}]}" +"%y")
    julday=$(date --date "${day_vec[${iday}]}" +"%j")

    echo "Processing day "${years}${julday}

    # create destination path for polar data    
    rm -rf ${data_destbase}${years}${julday}/
done
