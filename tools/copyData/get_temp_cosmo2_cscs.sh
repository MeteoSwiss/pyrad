#!/bin/bash
# get_temp_cosmo2_cscs.sh
# This script: 
# Gets COSMO-2 temperature file from CSCS and puts it radar coordinates
# created 20.06.2017 (fvj)

# set permits
umask 0002

# Config

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
        read -r -a hour_run <<< "$RUN"
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
    esac
    shift # past argument or value
done

nday=${#date_vec[@]}
nhour=${#hour_run[@]}

# Log
date

for ((iday=0; iday<${nday}; iday++)); do
    for ((ihour=0; ihour<${nhour}; ihour++)); do		
        if [ $radar != "None" ]
        then                    
            temp_cosmo_radar=${radar}
            temp_cosmo_res=${res}
            temp_cosmo_year=$(date --date "${date_vec[${iday}]}" +"%y")
            temp_cosmo_day=$(date --date "${date_vec[${iday}]}" +"%j")
            temp_cosmo_hour=${hour_run[${ihour}]}

            export temp_cosmo_radar
            export temp_cosmo_res
            export temp_cosmo_year
            export temp_cosmo_day
            export temp_cosmo_hour

            cd ${execpath}
            ${idlpath} -rt=temp_cosmo2_cscs.run
        fi
    done
done

# Log
echo "All done!"
date
