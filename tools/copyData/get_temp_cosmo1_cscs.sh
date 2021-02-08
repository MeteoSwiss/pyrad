#!/bin/bash
# This script: 
# Gets COSMO-1 temperature file from CSCS repository and puts it in the right
# folder to be processed in CSCS 
# created 23.03.2017 (fvj)

# set permits
umask 0002

# Config
cosmobaseraw='/store/msrad/cosmo/TEMP/raw/'
cosmopathcscs='/store/s83/owm/COSMO-1/ORDERS/MDR/'

#module load idl
. /apps/albis/itt/idl/idl84/inst/idl84/bin/idl_setup.bash

export IDL_OTLDIR=/store/msrad/idl
export IDL_STARTUP=$IDL_OTLDIR/im_setup.idl;export IDL_STARTUP

export IDL_DIR=/apps/albis/itt/idl/idl84/inst/idl84/
export LD_LIBRARY_PATH=$HOME/malsgit/src/libDX50/lib:$LD_LIBRARY_PATH

export IDL_NO_XWINDOWS=1
export DISPLAY=

execpath=$HOME/malsgit/src/idl/cosmo/
idlpath=/apps/albis/itt/idl/idl84/inst/idl84/bin/idl

# defaults
radar="None"
res="L"
hour_run_vec='all'
while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in        
        -d|--day)
        DAY="$2"
        OIFS=$IFS
        IFS=','
        read -r -a date_vec <<< "$DAY"
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
        -r|--radar)
        radar="$2"
        shift # past argument
        ;;
        -e|--res)
        res="$2"
        shift # past argument
        ;;
        -p|--dest_base)
        cosmobaseraw="$2"
        shift # past argument        
        ;;   
    esac
    shift # past argument or value
done

nday=${#date_vec[@]}
nhour_run=${#hour_run_vec[@]}

# Log
date

for ((iday=0; iday<${nday}; iday++)); do
    datedir=$(date --date "${date_vec[${iday}]}" +"%Y-%m-%d")
    data_destpath=${cosmobaseraw}${datedir}
    filebase=cosmo-1_MDR_3D_${date_vec[${iday}]}
    mkdir -p ${data_destpath}    

    hour24=0
    if [ "${hour_run_vec}" = "all" ];then
        cp ${cosmopathcscs}${filebase}*.nc ${data_destpath}
    else        
        for ((ihour_run=0; ihour_run<${nhour_run}; ihour_run++)); do
            hour_run=${hour_run_vec[${ihour_run}]}
            if [ ${hour_run} = '24' ]; then
                hour24=1
                continue
            fi
            cp ${cosmopathcscs}${filebase}${hour_run}.nc ${data_destpath}            
            
            if [ $radar != "None" ]
            then                    
                years=$(date --date "${date_vec[${iday}]}" +"%y")
                julday=$(date --date "${date_vec[${iday}]}" +"%j")
                        
                temp_cosmo_radar=${radar}
                temp_cosmo_res=${res}
                temp_cosmo_year=$(date --date "${date_vec[${iday}]}" +"%y")
                temp_cosmo_day=$(date --date "${date_vec[${iday}]}" +"%j")
                temp_cosmo_hour=${hour_run}
    
                export temp_cosmo_radar
                export temp_cosmo_res
                export temp_cosmo_year
                export temp_cosmo_day
                export temp_cosmo_hour
    
                cd ${execpath}
                ${idlpath} -rt=temp_cosmo1_cscs.run
            fi
        done
    fi
    chmod -R gu+rw ${data_destpath}     

    # add file 00:00 UTC the next day
    if [ "${hour_run_vec}" = "all" ] || [ "$hour24" -eq 1 ];then
        datedir=$(date -d "$(date --date "${date_vec[${iday}]}")+1 day" +"%Y-%m-%d")
        day=$(date -d "$(date --date "${date_vec[${iday}]}")+1 day" +"%Y%m%d")
        data_destpath=${cosmobaseraw}${datedir}
        filebase=cosmo-1_MDR_3D_${day}00.nc
        cp ${cosmopathcscs}${filebase} ${data_destpath}
        chmod -R gu+rw ${data_destpath}
    fi
done

# Log
echo "All done!"
date
