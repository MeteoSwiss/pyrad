#!/bin/bash
# This script: 
# Gets COSMO ISO0 data in HZT format from CSCS repository and puts it in the right
# folder to be processed in CSCS 
# created 23.03.2017 (fvj)

# set permits
umask 0002

#module load idl
. /apps/escha/itt/idl/idl85/inst/idl85/bin/idl_setup.bash

export IDL_OTLDIR=/store/msrad/idl
export IDL_STARTUP=$IDL_OTLDIR/im_setup.idl;export IDL_STARTUP

export IDL_DIR=/apps/escha/itt/idl/idl85/inst/idl85/
export LD_LIBRARY_PATH=$HOME/malsgit/src/libDX50/lib:$LD_LIBRARY_PATH

export IDL_NO_XWINDOWS=1
export DISPLAY=

execpath=$HOME/malsgit/src/idl/cosmo
idlpath=/apps/escha/itt/idl/idl85/inst/idl85/bin/idl

# defaults
nhour_fcst=7
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
        -f|--forecasted_hours)
        nhour_fcst="$2"
        shift # past argument
        ;;
    esac
    shift # past argument or value
done

nday=${#date_vec[@]}
nhour_run=${#hour_run[@]}

# Log
date

for ((iday=0; iday<${nday}; iday++)); do
    for ((ihour_run=0; ihour<${nhour_run}; ihour_run++)); do		
        for ((ihour_fcst=0; ihour<${nhour_fcst}; ihour_fcst++)); do
            iso0_cosmo_radar=${radar}
            iso0_cosmo_res=${res}
            iso0_cosmo_year=$(date --date "${date_vec[${iday}]}" +"%y")
            iso0_cosmo_day=$(date --date "${date_vec[${iday}]}" +"%j")
            iso0_cosmo_hour=${hour_run[${ihour_run}]}
            iso0_cosmo_hour_forecast=$(printf %02d $(( ${iso0_cosmo_hour} + ${ihour_fcst} )))

            export iso0_cosmo_radar
            export iso0_cosmo_res
            export iso0_cosmo_year
            export iso0_cosmo_day
            export iso0_cosmo_hour
            export iso0_cosmo_hour_forecast

            cd ${execpath}
            ${idlpath} -rt=iso0_cosmo_cscs.run
        done
    done
done

# Log
echo "All done!"
date
