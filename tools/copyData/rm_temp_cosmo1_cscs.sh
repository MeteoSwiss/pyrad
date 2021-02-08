#!/bin/bash
# This script: 
# Gets COSMO-1 temperature file from CSCS repository and puts it in the right
# folder to be processed in CSCS 
# created 23.03.2017 (fvj)

# set permits
umask 0002

# Config
dateCmd="/bin/date"
cosmobaseraw='/store/msrad/cosmo/cosmo1/TEMP/raw1/'

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
        -p|--dest_base)
        cosmobaseraw="$2"
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
    datedir=$(${dateCmd} --date "${date_vec[${iday}]}" +"%Y-%m-%d")

    for ((ihour=0; ihour<${nhour}; ihour++)); do		
        # remove cosmo file
        cosmoFileRaw=cosmo-1_MDR_3D_${date_vec[${iday}]}${hour_run[${ihour}]}.nc
        rm -f ${cosmobaseraw}${datedir}/${cosmoFileRaw}
		
        mkdir -p 
        cp ${cosmopathcscs}${cosmoFileRaw} ${datedir}
        chmod -R gu+rw ${datedir}
    done
    # remove day directory if empty
    if [ ! "$(ls -A ${cosmobaseraw}${datedir})" ]; then
        rm -rf ${cosmobaseraw}${datedir}
    fi    
done

# Log
echo "All done!"
date
