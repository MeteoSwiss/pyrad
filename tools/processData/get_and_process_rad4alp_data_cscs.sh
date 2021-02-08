#!/bin/bash
# This script controls the fetching of rad4alp and cosmo data from the
# archives and the data processing using pyrad
#
#
# Parameters:
# -d --day : days to process with format YYYYMMDD.
#            If more than one separate with a comma
# --start_time : start time of the processing each day
# --end_time : end time of the processing each day
# -c --cfgfile: config files used in the processing.
#               If more than one separate with a comma
# -p --cfgpath: path of the config files used in the processing.
#               If more than one separate with a comma
#               If more than one there must be as many as config files
# -r --radar: Name of the radar the data of which should be retrieved.
#             Can be A,D,L,P,W. If more than one separate with a comma.
#             There must be as many as config files
# -e --res: Resolution of the radar data. Can be H or L.
#           If more than one radar is retrieved separate with a comma.
#           There must be as many as config files
# -i --info: Info string to be added in the name of some output products
#            If more than one config file is used separate with a comma.
#            To be used there must be as many as config files
# -t --trajfile: Trajectory file. There must be as many as days to process
# --trajtype: Trajectory type. Either 'plane' or 'lightnint'. Default 'Plane'
# --flashnr: If type of trajectory is 'lightning' flash number the data of which
#            will be processed 0 means that all lightning data will be processed
# --get_data : if set, get radar data. Default 1
# --ele : radar elevations to retrieve in format 001 to 020
#         if "all" retrieve all the elevations
# --rm_data : if set, remove radar data after processing if data was fetched
#             from archive. Default 0
# --data_destbase : rad4alp data destination base path.
#                   Default /store/msrad/radar/rad4alp/rawdata/
# --get_cosmo : if set, get COSMO temperature data. Default 1
# --rm_cosmo : if set, remove COSMO temperature data after processing if data
#              was fetched from archive. Default 0
# --cosmo_destbase : COSMO temperature data destination base path.
#                    Default /store/msrad/cosmo/cosmo1/TEMP/raw1/
# --get_hzt : if set, get HZT data. Default 0
# --rm_hzt : if set, remove HZT data after processing if data was fecthed from
#            archive. Default 0
# -h|--hour : hour from run of iso0 to retrieved format 00 to 06
#             if "all" all forecast times will be retrieved
# --hzt_destbase : HZT data destination base path.
#                  Default /store/msrad/cosmo/cosmo1/HZT/
# --MP_DSET : Parallelize dataset generation. Default 0
# --MP_PROD : Parallelize product generation. Default 0
# --MP_PROF : Profile parallelization. Default 0
umask 0002;
export PATH="/store/msrad/utils/anaconda3-pyrad/bin/:$PATH"

pyradpath="$HOME/pyrad/src/pyrad_proc/scripts/"

# default variables
START_TIME="000001"
END_TIME="240000"

GET_DATA=1
RM_DATA=0
DATA_DESTBASE=/store/msrad/radar/rad4alp/rawdata/

GET_COSMO=1
RM_COSMO=0
COSMO_DESTBASE=/store/msrad/cosmo/TEMP/raw/

GET_HZT=0
RM_HZT=0
HZT_DESTBASE=/store/msrad/cosmo/HZT/

info_vec='None'
cfgpath_vec="$HOME/pyrad/config/processing/"
trajfile_vec='None'
TRAJTYPE='plane'
FLASHNR=0

ELE='all'
HOUR='all'

MP_DSET=0
MP_PROD=0
MP_PROF=0
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
        --start_time)
        START_TIME="$2"
        shift # past argument
        ;;
        --end_time)
        END_TIME="$2"
        shift # past argument
        ;;
        -c|--cfgfile)
        CFGFILE="$2"
        OIFS=$IFS
        IFS=','
        read -r -a cfgfile_vec <<< "$CFGFILE"
        IFS=$OIFS
        shift # past argument
        ;;
        -p|--cfgpath)
        CFGPATH="$2"
        OIFS=$IFS
        IFS=','
        read -r -a cfgpath_vec <<< "$CFGPATH"
        IFS=$OIFS
        shift # past argument
        ;;
        -t|--trajfile)
        TRAJFILE="$2"
        OIFS=$IFS
        IFS=','
        read -r -a trajfile_vec <<< "$TRAJFILE"
        IFS=$OIFS
        shift # past argument
        ;;
        --trajtype)
        TRAJTYPE="$2"
        shift # past argument
        ;;
        --flashnr)
        FLASHNR="$2"
        shift # past argument
        ;;
        -r|--radar)
        RADAR="$2"
        OIFS=$IFS
        IFS=','
        read -r -a radar_vec <<< "$RADAR"
        IFS=$OIFS
        shift # past argument
        ;;
        -e|--res)
        RES="$2"
        OIFS=$IFS
        IFS=','
        read -r -a res_vec <<< "$RES"
        IFS=$OIFS
        shift # past argument
        ;;
        -i|--info)
        INFO="$2"
        OIFS=$IFS
        IFS=','
        read -r -a info_vec <<< "$INFO"
        IFS=$OIFS
        shift # past argument
        ;;
        --get_data)
        GET_DATA="$2"
        shift # past argument
        ;;
        --ele)
        ELE="$2"
        shift # past argument
        ;;
        --rm_data)
        RM_DATA="$2"
        shift # past argument
        ;;
        --data_destbase)
        DATA_DESTBASE="$2"
        shift # past argument
        ;;
        --get_cosmo)
        GET_COSMO="$2"
        shift # past argument
        ;;
        --rm_cosmo)
        RM_COSMO="$2"
        shift # past argument
        ;;
        --cosmo_destbase)
        COSMO_DESTBASE="$2"
        shift # past argument
        ;;
        --get_hzt)
        GET_HZT="$2"
        shift # past argument
        ;;
        -h|--hour)
        HOUR="$2"
        shift # past argument
        ;;
        --rm_hzt)
        RM_HZT="$2"
        shift # past argument
        ;;
        --hzt_destbase)
        HZT_DESTBASE="$2"
        shift # past argument
        ;;
        --MP_DSET)
        MP_DSET="$2"
        shift # past argument
        ;;
        --MP_PROD)
        MP_PROD="$2"
        shift # past argument
        ;;
        --MP_PROF)
        MP_PROF="$2"
        shift # past argument
        ;;
    esac
    shift # past argument or value
done

nday=${#day_vec[@]}

ncfg=${#cfgfile_vec[@]}
nrad=${#radar_vec[@]}
nres=${#res_vec[@]}
ninfo=${#info_vec[@]}
ncfgpath=${#cfgpath_vec[@]}
ntrajfile=${#trajfile_vec[@]}

if [ "$ncfg" -ne "$nrad" ] || [ "$ncfg" -ne "$nres" ];then
    echo "In order to retrieve radar data the number of radars and resolutions should be the same as the number of config files"
    echo "nrad: "$nrad
    echo "ncfg: "$ncfg
    echo "nres: "$nres
    GET_DATA=0
    RM_DATA=0
fi

# get COSMO run time
if [ "$GET_COSMO" -eq 1 ] || [ "$GET_HZT" -eq 1 ]; then
    run_period=3
    start_hour_aux=$(echo ${START_TIME} | cut -c1-2)
    start_hour=$[(${start_hour_aux}/${run_period})*${run_period}]
    ((start_hour=start_hour-start_hour%run_period))
    end_hour=$(echo ${END_TIME} | cut -c1-2)
    vals=($(seq $start_hour $run_period $end_hour))
    nvals=${#vals[@]}
    RUN_TIME=$(printf %02d ${vals[0]})
    for ((ivals=1; ivals<${nvals}; ivals++)); do
        RUN_TIME="${RUN_TIME},$(printf %02d ${vals[${ivals}]})"
    done
fi

for ((iday=0; iday<${nday}; iday++)); do
    procday_start=`date +%s`

    DAY=${day_vec[${iday}]}
    for ((icfg=0; icfg<${ncfg}; icfg++)); do
        proccfg_start=`date +%s`

        CFGFILE=${cfgfile_vec[${icfg}]}
        IFS='.' read -r CFGFILE_BASE TERMINATION <<< "$CFGFILE"

        if [ "$ncfg" -eq "$nrad" ] || [ "$ncfg" -eq "$nres" ];then
            RADAR=${radar_vec[${icfg}]}
            RES=${res_vec[${icfg}]}
        fi
        if [ "$ncfg" -eq "$ninfo" ];then
            INFO=${info_vec[${icfg}]}
        else
            INFO='None'
        fi
        if [ "$ncfgpath" -eq "$ncfg" ];then
            CFGPATH=${cfgpath_vec[${icfg}]}
        else
            CFGPATH="$HOME/pyrad/config/processing/"
        fi

        if [ "$GET_DATA" -eq 1 ]; then
            echo $(date --utc)" ${CFGFILE} Getting data for radar ${RADAR} with resolution ${RES} and day ${DAY}"
            proc_start=`date +%s`

            # Copy data from CSCS repository into the folder structure to process it
            LOGFILE=$HOME/log/${DAY}_${CFGFILE_BASE}_get_rad4alp_data_CSCS.log
            bash $HOME/pyrad/tools/copyData/get_rad4alp_data_CSCS_2.sh -r $RADAR -e $RES -d $DAY --start_time ${START_TIME} --end_time ${END_TIME} -p $DATA_DESTBASE --ele $ELE >$LOGFILE 2>$LOGFILE

            proc_end=`date +%s`
            runtime=$((proc_end-proc_start))
            echo $(date --utc)" ${CFGFILE} Finished getting data for radar ${RADAR} with resolution ${RES} and day ${DAY}. Run time: ${runtime} s"
        fi

        if [ "$GET_COSMO" -eq 1 ]; then
            echo $(date --utc)" ${CFGFILE} Getting COSMO data for runs ${RUN_TIME} and day ${DAY}"
            proc_start=`date +%s`

            # Get COSMO temperature information files
            LOGFILE=$HOME/log/${DAY}_${CFGFILE_BASE}_temp_cosmo1_rad4alp.log
            bash $HOME/pyrad/tools/copyData/get_temp_cosmo1_cscs.sh -d $DAY -t ${RUN_TIME} -p $COSMO_DESTBASE >$LOGFILE 2>$LOGFILE

            proc_end=`date +%s`
            runtime=$((proc_end-proc_start))
            echo $(date --utc)" ${CFGFILE} Finished getting COSMO data for runs ${RUN_TIME} and day ${DAY}. Run time: ${runtime} s"
        fi

        if [ "$GET_HZT" -eq 1 ]; then
            echo $(date --utc)" Getting HZT data for day ${DAY}"
            proc_start=`date +%s`

            # copy HZT files to right folder
            LOGFILE=$HOME/log/${DAY}_${CFGFILE_BASE}_hzt_rad4alp.log
            bash $HOME/pyrad/tools/copyData/get_hzt_cscs.sh -d $DAY -t ${RUN_TIME} -h $HOUR -p $HZT_DESTBASE >$LOGFILE 2>$LOGFILE

            proc_end=`date +%s`
            runtime=$((proc_end-proc_start))
            echo $(date --utc)" ${CFGFILE} Finished getting HZT data for day ${DAY}. Run time: ${runtime} s"
        fi

        echo $(date --utc)" ${CFGFILE} Processing config file ${CFGFILE} on day ${DAY}"
        proc_start=`date +%s`

        # process data
        source activate pyrad

        LOGFILE=$HOME/log/${DAY}_${CFGFILE_BASE}.log
        cd ${pyradpath}

        if [ "${trajfile_vec}" = "None" ];then
            # No trajectory file used

            # check start and end time
            hour_end=$(echo ${END_TIME} | cut -c1-2)
            if [ "$hour_end" -ge 24 ]; then
                day_end=$(date -d "$(date --date "${DAY}")+1 day" +"%Y%m%d")
                hour_end=$[${hour_end}-24]

                hour_end=$(printf %02d ${hour_end})
                min_end=$(echo ${END_TIME} | cut -c3-4)
                sec_end=$(echo ${END_TIME} | cut -c5-6)

                time_end=${day_end}${hour_end}${min_end}${sec_end}
            else
                time_end=${DAY}${END_TIME}
            fi
            time_start=${DAY}${START_TIME}

            python -u main_process_data.py ${CFGFILE} --starttime ${time_start} --endtime ${time_end} -i ${INFO} --cfgpath ${CFGPATH} --MULTIPROCESSING_DSET ${MP_DSET} --MULTIPROCESSING_PROD ${MP_PROD} --PROFILE_MULTIPROCESSING ${MP_PROF} >$LOGFILE 2>&1
        else
            for ((itraj=0; itraj<${ntrajfile}; itraj++)); do
                python -u main_process_data.py ${CFGFILE} -i ${INFO} -t ${trajfile_vec[${itraj}]} --trajtype $TRAJTYPE --flashnr $FLASHNR --cfgpath ${CFGPATH} --MULTIPROCESSING_DSET ${MP_DSET} --MULTIPROCESSING_PROD ${MP_PROD} --PROFILE_MULTIPROCESSING ${MP_PROF} >$LOGFILE 2>&1
            done
        fi

        source deactivate

        proc_end=`date +%s`
        runtime=$((proc_end-proc_start))
        echo $(date --utc)" ${CFGFILE} Finished Processing config file ${CFGFILE} on day ${DAY}. Run time: ${runtime} s"

        if [ "$RM_DATA" -eq 1 ] && [ "$GET_DATA" -eq 1 ]; then
            echo $(date --utc)" ${CFGFILE} Removing data for radar ${RADAR} with resolution ${RES} and day ${DAY}"
            proc_start=`date +%s`

            LOGFILE=$HOME/log/${DAY}_${RADAR}${RES}_rm_rad4alp_data_CSCS.log
            bash $HOME/pyrad/tools/copyData/rm_rad4alp_data_CSCS_2.sh -r $RADAR -e $RES -d $DAY --start_time ${START_TIME} --end_time ${END_TIME} -p $DATA_DESTBASE --ele $ELE >$LOGFILE 2>$LOGFILE

            proc_end=`date +%s`
            runtime=$((proc_end-proc_start))
            echo $(date --utc)" ${CFGFILE} Finished removing data for radar ${RADAR} with resolution ${RES} and day ${DAY}. Run time: ${runtime} s"
        fi
        if [ "$RM_COSMO" -eq 1 ] && [ "$GET_COSMO" -eq 1 ]; then
            echo $(date --utc)" ${CFGFILE} Removing COSMO data for runs ${RUN_TIME} and day ${DAY}"
            proc_start=`date +%s`

            LOGFILE=$HOME/log/${DAY}_${RADAR}${RES}_rm_temp_cosmo1_cscs.log
            bash $HOME/pyrad/tools/copyData/rm_temp_cosmo1_cscs.sh -d $DAY -t ${RUN_TIME} -p $COSMO_DESTBASE >$LOGFILE 2>$LOGFILE

            proc_end=`date +%s`
            runtime=$((proc_end-proc_start))
            echo $(date --utc)" ${CFGFILE} Finished removing COSMO data for runs ${RUN_TIME} and day ${DAY}. Run time: ${runtime} s"
        fi
        if [ "$RM_HZT" -eq 1 ] && [ "$GET_HZT" -eq 1 ]; then
            echo $(date --utc)" ${CFGFILE} Removing HZT data for day ${DAY}"
            proc_start=`date +%s`

            LOGFILE=$HOME/log/${DAY}_${RADAR}${RES}_rm_hzt_cscs.log
            bash $HOME/pyrad/tools/copyData/rm_hzt_data_cscs.sh -d $DAY -p $HZT_DESTBASE >$LOGFILE 2>$LOGFILE

            proc_end=`date +%s`
            runtime=$((proc_end-proc_start))
            echo $(date --utc)" ${CFGFILE} Finished removing HZT data for day ${DAY}. Run time: ${runtime} s"
        fi

        proccfg_end=`date +%s`
        runtime=$((proccfg_end-proccfg_start))
        echo $(date --utc)" ${CFGFILE} ${DAY} Run time: ${runtime} s"
    done

    procday_end=`date +%s`
    runtime=$((procday_end-procday_start))
    echo $(date --utc)" ${DAY} Run time: ${runtime} s"
done
