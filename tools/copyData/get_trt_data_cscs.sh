#!/bin/bash

# This script gets the rad4alp data from the CSCS
# repository and puts it in the right folder.
# To be run in CSCS

# set permits
umask 0002

# Config
dateCmd="/bin/date"

rawdata_origbase=/store/msrad/radar/swiss/data/

# default variables
data_destbase=/store/msrad/radar/rad4alp/TRT/
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
    yearl=$(date --date "${day_vec[${iday}]}" +"%Y")
    years=$(date --date "${day_vec[${iday}]}" +"%y")
    julday=$(date --date "${day_vec[${iday}]}" +"%j")

    yearl_end=$(date -d "$(date --date "${day_vec[${iday}]}")+1 day" +"%Y")
    years_end=$(date -d "$(date --date "${day_vec[${iday}]}")+1 day" +"%y")
    julday_end=$(date -d "$(date --date "${day_vec[${iday}]}")+1 day" +"%j")

    echo "Processing day "${years}${julday}
    # transfer data from CSCS to destination folder, unzip it and remove zip file

    data_origpath=${rawdata_origbase}${yearl}/${years}${julday}/
    # check type of file that exists in the repository
    if [ `ls ${data_origpath}TRTC${years}${julday}.zip` ]
    then
        echo "File in repository"
    else
        echo "No file found in repository"
        continue
    fi

    data_destpath=${data_destbase}${years}${julday}/TRTC${years}${julday}/
    filebase=TRTC${years}${julday}

    # create destination path for polar data
    mkdir -p ${data_destpath}
    unzip -o ${data_origpath}${filebase}.zip *.trt -d ${data_destpath}
    chmod -R gu+rw ${data_destpath}
done
