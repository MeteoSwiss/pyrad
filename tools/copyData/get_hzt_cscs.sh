#!/bin/bash

# This script gets the HZT data from the CSCS
# repository and puts it in the right folder.
# To be run in CSCS

# set permits
umask 0002

# Config
dateCmd="/bin/date"

rawdata_origbase=/store/msrad/radar/swiss/data/
data_destbase=/store/msrad/cosmo/HZT/

hour_fcst_vec='all'
hour_run_vec='all'
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
        -t|--time)
        RUN="$2"
        OIFS=$IFS
        IFS=','
        read -r -a hour_run_vec <<< "$RUN"
        IFS=$OIFS
        shift # past argument
        ;;        
        -h|--hour)
        HOUR="$2"
        OIFS=$IFS
        IFS=','
        read -r -a hour_fcst_vec <<< "$HOUR"
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
nhour_fcst=${#hour_fcst_vec[@]}
nhour_run=${#hour_run_vec[@]}

for ((iday=0; iday<${nday}; iday++)); do
    yearl=$(date --date "${day_vec[${iday}]}" +"%Y")
    years=$(date --date "${day_vec[${iday}]}" +"%y")
    julday=$(date --date "${day_vec[${iday}]}" +"%j")
    
    echo "Processing day "${years}${julday}

    # create destination path for polar data and unzip files there
    data_origpath=${rawdata_origbase}${yearl}/${years}${julday}/
    data_destpath=${data_destbase}${years}${julday}/
    filebase=HZT${years}${julday}
    mkdir -p ${data_destpath}
    
    hour24=0
    if [ "${hour_fcst_vec}" = "all" ] && [ "${hour_run_vec}" = "all" ];then
        unzip -o ${data_origpath}${filebase}.zip -d ${data_destpath}
    else        
        for ((ihour_run=0; ihour_run<${nhour_run}; ihour_run++)); do
            if [ ${hour_run_vec[${ihour_run}]} = "all" ]; then
                hour_run=*
            else
                hour_run=${hour_run_vec[${ihour_run}]}
            fi
            if [ ${hour_run} = '24' ]; then
                hour24=1
                continue
            fi
            for ((ihour_fcst=0; ihour_fcst<${nhour_fcst}; ihour_fcst++)); do
                if [ ${hour_fcst_vec[${ihour_fcst}]} = "all" ]; then
                    hour_fcst=*
                else
                    hour_fcst="$((${hour_fcst_vec[${ihour_fcst}]} + 800))"
                fi
                unzip -o ${data_origpath}${filebase}.zip ${filebase}${hour_run}00*.${hour_fcst} -d ${data_destpath}
            done
        done
    fi
    chmod -R gu+rw ${data_destpath}
    
    # add file 00:00 UTC the next day
    if [ "${hour_run_vec}" = "all" ] || [ "$hour24" -eq 1 ];then
        yearl_end=$(date -d "$(date --date "${day_vec[${iday}]}")+1 day" +"%Y")
        years_end=$(date -d "$(date --date "${day_vec[${iday}]}")+1 day" +"%y")
        julday_end=$(date -d "$(date --date "${day_vec[${iday}]}")+1 day" +"%j")
        
        data_origpath=${rawdata_origbase}${yearl_end}/${years_end}${julday_end}/
        data_destpath=${data_destbase}${years_end}${julday_end}/
        filebase=HZT${years_end}${julday_end}
        mkdir -p ${data_destpath}
        if [ "${hour_fcst_vec}" = "all" ];then
            unzip -o ${data_origpath}${filebase}.zip ${filebase}0000*.* -d ${data_destpath}
        else
            for ((ihour=0; ihour<${nhour_fcst}; ihour++)); do
                hour="$((${hour_fcst_vec[${ihour}]} + 800))"
                unzip -o ${data_origpath}${filebase}.zip ${filebase}0000*.${hour} -d ${data_destpath}
            done
        fi
        chmod -R gu+rw ${data_destpath}
    fi
done
