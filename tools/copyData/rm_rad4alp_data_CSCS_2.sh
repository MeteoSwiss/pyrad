#!/bin/bash

# This script removes the rad4alp Cartesian data from the CSCS
# repository
# To be run in CSCS

# set permits
umask 0002

# Config
dateCmd="/bin/date"

# default variables
START_TIME="000001"
END_TIME="240000"
time_vec='all'

data_destbase=/store/msrad/radar/rad4alp/rawdata/
while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
        --prod)
        PROD="$2"
        OIFS=$IFS
        IFS=','
        read -r -a prod_vec <<< "$PROD"
        IFS=$OIFS
        shift # past argument
        ;;
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
        -p|--dest_base)
        data_destbase="$2"
        shift # past argument
        ;;
    esac
    shift # past argument or value
done

# get vector of times to retrieve

# set start time as multiple of 5 minutes
h_s=$(echo ${START_TIME} | cut -c1-2)
m_s=$(echo ${START_TIME} | cut -c3-4)
s_s=$(echo ${START_TIME} | cut -c5-6)

if [ $s_s -gt 0 ]; then
    s_s=00
    m_s="$(($m_s + 5))"
fi
m_s=$[(${m_s}/5)*5]
if [ $m_s -gt 60 ]; then
    h_s="$(($h_s + 1))"
    h_s="$(($m_s - 60))"
fi
st=$h_s:$m_s:$s_s
st_s=$(date --date $st +%s)
st_ref=$(date --date 00:00:00 +%s)
dt_st=$((st_ref-st_s))

# set end time as multiple of 5 minutes
hour24=0
h_e=$(echo ${END_TIME} | cut -c1-2)
m_e=$(echo ${END_TIME} | cut -c3-4)
s_e=$(echo ${END_TIME} | cut -c5-6)

if [ $s_e -gt 0 ]; then
    s_e=00
    m_e="$(($m_e + 5))"
fi
m_e=$[(${m_e}/5)*5]
m_e=$(printf %02d ${m_e})
if [ $m_e -gt 60 ]; then
    h_e="$(($h_e + 1))"
    m_e="$(($m_e - 60))"
fi
if [ $h_e = "24" ]; then
    h_e=23
    m_e=55
    s_e=00
    hour24=1
fi
et=$h_e:$m_e:$s_e
et_s=$(date --date $et +%s)
et_ref=$(date --date 23:55:00 +%s)
dt_et=$((et_ref-et_s))

# time span does not cover the whole day. Get time stamps to retrieve
if [ $dt_st -lt 0 ] || [ $dt_et -gt 0 ]; then
    time_vec=()
    i=0
    time_sample=$(date --date $st +%H:%M)
    time_vec[0]=$(date --date $time_sample +%H%M)
    time_sample_end=$(date --date $et +%H:%M)
    while [ $time_sample != $time_sample_end ]; do
        time_sample=$(date -d "$(date --date "$time_sample")+5 min" +"%H:%M")
        let i+=1
        time_vec[i]=$(date --date $time_sample +%H%M)
    done
fi

nday=${#day_vec[@]}
nprod=${#prod_vec[@]}
ntime=${#time_vec[@]}

for ((iday=0; iday<${nday}; iday++)); do
    yearl=$(date --date "${day_vec[${iday}]}" +"%Y")
    years=$(date --date "${day_vec[${iday}]}" +"%y")
    julday=$(date --date "${day_vec[${iday}]}" +"%j")

    yearl_end=$(date -d "$(date --date "${day_vec[${iday}]}")+1 day" +"%Y")
    years_end=$(date -d "$(date --date "${day_vec[${iday}]}")+1 day" +"%y")
    julday_end=$(date -d "$(date --date "${day_vec[${iday}]}")+1 day" +"%j")

    echo "Removing day "${years}${julday}

    for ((iprod=0; iprod<${nprod}; iprod++)); do
        prod=${prod_vec[${iprod}]}
        echo "Removing product "${prod}

        if [ "${res}" == "H" ]
        then
            if [ "${ele_vec}" = "all" ] && [ "${time_vec}" = "all" ];then
                # remove polar data from 00:05 to 23:55
                find ${data_destbase}${years}${julday}/M${res}${radar}${years}${julday}/ -type f -not -name M${res}${radar}${years}${julday}0000*.* -delete
                find ${data_destbase}${years}${julday}/P${res}${radar}${years}${julday}/ -type f -not -name P${res}${radar}${years}${julday}0000*.* -delete
            else
                for ((itime=0; itime<${ntime}; itime++)); do
                    if [ "${time_vec[${itime}]}" = "all" ]; then
                        time_rad=*
                    else
                        time_rad=${time_vec[${itime}]}
                    fi
                    if [ "${ele_vec}" = "all" ]; then
                        rm -f ${data_destbase}${years}${julday}/M${res}${radar}${years}${julday}/M${res}${radar}${years}${julday}${time_rad}*.*
                        rm -f ${data_destbase}${years}${julday}/P${res}${radar}${years}${julday}/P${res}${radar}${years}${julday}${time_rad}*.*
                    else
                        for ((iele=0; iele<${nele}; iele++)); do
                            ele=${ele_vec[${iele}]}
                            rm -f ${data_destbase}${years}${julday}/M${res}${radar}${years}${julday}/M${res}${radar}${years}${julday}${time_rad}*.${ele}
                            rm -f ${data_destbase}${years}${julday}/P${res}${radar}${years}${julday}/P${res}${radar}${years}${julday}${time_rad}*.${ele}
                        done
                    fi
                done
            fi

            # remove directory if empty
            if [ ! "$(ls -A ${data_destbase}${years}${julday}/M${res}${radar}${years}${julday})" ]; then
                rm -rf ${data_destbase}${years}${julday}/M${res}${radar}${years}${julday}
            fi
            if [ ! "$(ls -A ${data_destbase}${years}${julday}/P${res}${radar}${years}${julday})" ]; then
                rm -rf ${data_destbase}${years}${julday}/P${res}${radar}${years}${julday}
            fi

            # remove polar data at 00:00 the next day
            if [ "${time_vec}" = "all" ];then
                if [ "${ele_vec}" = "all" ];then
                    rm -f ${data_destbase}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end}0000*.*
                    rm -f ${data_destbase}${years_end}${julday_end}/P${res}${radar}${years_end}${julday_end}/P${res}${radar}${years_end}${julday_end}0000*.*
                else
                    for ((iele=0; iele<${nele}; iele++)); do
                        ele=${ele_vec[${iele}]}
                        rm -f ${data_destbase}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end}*.${ele}
                        rm -f ${data_destbase}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end}/P${res}${radar}${years_end}${julday_end}*.${ele}
                    done
                fi

                # remove directory if empty
                if [ ! "$(ls -A ${data_destbase}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end})" ]; then
                    rm -rf ${data_destbase}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end}
                fi
                if [ ! "$(ls -A ${data_destbase}${years_end}${julday_end}/P${res}${radar}${years_end}${julday_end})" ]; then
                    rm -rf ${data_destbase}${years_end}${julday_end}/P${res}${radar}${years_end}${julday_end}
                fi
            fi

            # check if there is L polar data directory. If not remove also status data
            path_ML=${data_destbase}${years}${julday}/ML${radar}${years}${julday}/
            path_PL=${data_destbase}${years}${julday}/PL${radar}${years}${julday}/
            if [ ! -d "${path_ML}" ] && [ ! -d "${path_PL}" ]; then
                if [ "${ele_vec}" = "all" ] && [ "${time_vec}" = "all" ];then
                    find ${data_destbase}${years}${julday}/ST${radar}${years}${julday}/ -type f -not -name ST${radar}${years}${julday}0000*.xml -delete
                fi
            fi
            # remove directory if empty
            if [ ! "$(ls -A ${data_destbase}${years}${julday}/ST${radar}${years}${julday})" ]; then
                rm -rf ${data_destbase}${years}${julday}/ST${radar}${years}${julday}
            fi

            if [ "${time_vec}" = "all" ];then
                # check if there is L polar data directory next day. If not remove also status data
                path_ML=${data_destbase}${years_end}${julday_end}/ML${radar}${years_end}${julday_end}/
                path_PL=${data_destbase}${years_end}${julday_end}/PL${radar}${years_end}${julday_end}/
                if [ ! -d "${path_ML}" ] && [ ! -d "${path_PL}" ]; then
                    rm -f ${data_destbase}${years_end}${julday_end}/ST${radar}${years_end}${julday_end}/ST${radar}${years_end}${julday_end}0000*.xml
                fi
                # remove directory if empty
                if [ ! "$(ls -A ${data_destbase}${years_end}${julday_end}/ST${radar}${years_end}${julday_end})" ]; then
                    rm -rf ${data_destbase}${years_end}${julday_end}/ST${radar}${years_end}${julday_end}
                fi
            fi
        else
            if [ "${ele_vec}" = "all" ] && [ "${time_vec}" = "all" ];then
                # remove polar data from 00:05 to 23:55
                find ${data_destbase}${years}${julday}/M${res}${radar}${years}${julday}/ -type f -not -name M${res}${radar}${years}${julday}0000*.* -delete
                find ${data_destbase}${years}${julday}/P${res}${radar}${years}${julday}/ -type f -not -name P${res}${radar}${years}${julday}0000*.* -delete
            else
                for ((itime=0; itime<${ntime}; itime++)); do
                    if [ "${time_vec[${itime}]}" = "all" ]; then
                        time_rad=*
                    else
                        time_rad=${time_vec[${itime}]}
                    fi
                    if [ "${ele_vec}" = "all" ]; then
                        rm -f ${data_destbase}${years}${julday}/M${res}${radar}${years}${julday}/M${res}${radar}${years}${julday}${time_rad}*.*
                        rm -f ${data_destbase}${years}${julday}/P${res}${radar}${years}${julday}/P${res}${radar}${years}${julday}${time_rad}*.*
                    else
                        for ((iele=0; iele<${nele}; iele++)); do
                            ele=${ele_vec[${iele}]}
                            rm -f ${data_destbase}${years}${julday}/M${res}${radar}${years}${julday}/M${res}${radar}${years}${julday}${time_rad}*.${ele}
                            rm -f ${data_destbase}${years}${julday}/P${res}${radar}${years}${julday}/P${res}${radar}${years}${julday}${time_rad}*.${ele}
                        done
                    fi
                done
            fi

            # remove directory if empty
            if [ ! "$(ls -A ${data_destbase}${years}${julday}/M${res}${radar}${years}${julday})" ]; then
                rm -rf ${data_destbase}${years}${julday}/M${res}${radar}${years}${julday}
            fi
            if [ ! "$(ls -A ${data_destbase}${years}${julday}/P${res}${radar}${years}${julday})" ]; then
                rm -rf ${data_destbase}${years}${julday}/P${res}${radar}${years}${julday}
            fi

            # remove polar data at 00:00 the next day
            if [ "${time_vec}" = "all" ];then
                if [ "${ele_vec}" = "all" ];then
                    rm -f ${data_destbase}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end}0000*.*
                    rm -f ${data_destbase}${years_end}${julday_end}/P${res}${radar}${years_end}${julday_end}/P${res}${radar}${years_end}${julday_end}0000*.*
                else
                    for ((iele=0; iele<${nele}; iele++)); do
                        ele=${ele_vec[${iele}]}
                        rm -f ${data_destbase}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end}*.${ele}
                        rm -f ${data_destbase}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end}/P${res}${radar}${years_end}${julday_end}*.${ele}
                    done
                fi
            fi
            # remove directory if empty
            if [ ! "$(ls -A ${data_destbase}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end})" ]; then
                rm -rf ${data_destbase}${years_end}${julday_end}/M${res}${radar}${years_end}${julday_end}
            fi
            if [ ! "$(ls -A ${data_destbase}${years_end}${julday_end}/P${res}${radar}${years_end}${julday_end})" ]; then
                rm -rf ${data_destbase}${years_end}${julday_end}/P${res}${radar}${years_end}${julday_end}
            fi

            # remove hydrometeor classification data
            if [ "${ele_vec}" = "all" ] && [ "${time_vec}" = "all" ];then
                find ${data_destbase}${years}${julday}/YM${radar}${years}${julday}/ -type f -not -name YM${radar}${years}${julday}0000*.* -delete
                rm -f ${data_destbase}${years_end}${julday_end}/YM${radar}${years_end}${julday_end}/YM${radar}${years_end}${julday_end}0000*.*
            else
                for ((itime=0; itime<${ntime}; itime++)); do
                    if [ "${time_vec[${itime}]}" = "all" ]; then
                        time_rad=*
                    else
                        time_rad=${time_vec[${itime}]}
                    fi
                    if [ "${ele_vec}" = "all" ]; then
                        rm -f ${data_destbase}${years}${julday}/YM${radar}${years}${julday}/YM${radar}${years}${julday}${time_rad}*.*
                    else
                        for ((iele=0; iele<${nele}; iele++)); do
                            ele="$(expr ${ele_vec[${iele}]} + 800)"
                            rm -f ${data_destbase}${years}${julday}/YM${radar}${years}${julday}/YM${radar}${years}${julday}${time_rad}*.${ele}
                        done
                    fi
                done
            fi

            # remove directory if empty
            if [ ! "$(ls -A ${data_destbase}${years}${julday}/YM${radar}${years}${julday})" ]; then
                rm -rf ${data_destbase}${years}${julday}/YM${radar}${years}${julday}
            fi
            if [ ! "$(ls -A ${data_destbase}${years_end}${julday_end}/YM${radar}${years_end}${julday_end})" ]; then
                rm -rf ${data_destbase}${years_end}${julday_end}/YM${radar}${years_end}${julday_end}
            fi

            # remove dealiased Doppler velocity data
            if [ "${ele_vec}" = "all" ] && [ "${time_vec}" = "all" ];then
                find ${data_destbase}${years}${julday}/DV${radar}${years}${julday}/ -type f -not -name DV${radar}${years}${julday}0000*.* -delete
                rm -f ${data_destbase}${years_end}${julday_end}/DV${radar}${years_end}${julday_end}/DV${radar}${years_end}${julday_end}0000*.*
            else
                for ((itime=0; itime<${ntime}; itime++)); do
                    if [ "${time_vec[${itime}]}" = "all" ]; then
                        time_rad=*
                    else
                        time_rad=${time_vec[${itime}]}
                    fi
                    if [ "${ele_vec}" = "all" ]; then
                        rm -f ${data_destbase}${years}${julday}/DV${radar}${years}${julday}/DV${radar}${years}${julday}${time_rad}*.*
                    else
                        for ((iele=0; iele<${nele}; iele++)); do
                            ele="$(expr ${ele_vec[${iele}]} + 800)"
                            rm -f ${data_destbase}${years}${julday}/DV${radar}${years}${julday}/DV${radar}${years}${julday}${time_rad}*.${ele}
                        done
                    fi
                done
            fi

            # remove directory if empty
            if [ ! "$(ls -A ${data_destbase}${years}${julday}/YM${radar}${years}${julday})" ]; then
                rm -rf ${data_destbase}${years}${julday}/DV${radar}${years}${julday}
            fi
            if [ ! "$(ls -A ${data_destbase}${years_end}${julday_end}/YM${radar}${years_end}${julday_end})" ]; then
                rm -rf ${data_destbase}${years_end}${julday_end}/DV${radar}${years_end}${julday_end}
            fi

            path_MH=${data_destbase}${years}${julday}/MH${radar}${years}${julday}/
            path_PH=${data_destbase}${years}${julday}/PH${radar}${years}${julday}/

            # check if there is L polar data directory. If not remove also status data
            if [ ! -d "${path_MH}" ] && [ ! -d "${path_PH}" ]; then
                find ${data_destbase}${years}${julday}/ST${radar}${years}${julday}/ -type f -not -name ST${radar}${years}${julday}0000*.xml -delete
            fi
            # remove directory if empty
            if [ ! "$(ls -A ${data_destbase}${years}${julday}/ST${radar}${years}${julday})" ]; then
                rm -rf ${data_destbase}${years}${julday}/ST${radar}${years}${julday}
            fi

            if [ "${time_vec}" = "all" ];then
                # check if there is L polar data directory next day. If not remove also status data
                path_ML=${data_destbase}${years_end}${julday_end}/ML${radar}${years_end}${julday_end}/
                path_PL=${data_destbase}${years_end}${julday_end}/PL${radar}${years_end}${julday_end}/
                if [ ! -d "${path_ML}" ] && [ ! -d "${path_PL}" ]; then
                    rm -f ${data_destbase}${years_end}${julday_end}/ST${res}${radar}${years_end}${julday_end}/ST${res}${radar}${years_end}${julday_end}0000*.xml
                fi
                # remove directory if empty
                if [ ! "$(ls -A ${data_destbase}${years_end}${julday_end}/ST${radar}${years_end}${julday_end})" ]; then
                    rm -rf ${data_destbase}${years_end}${julday_end}/ST${radar}${years_end}${julday_end}
                fi
            fi
        fi
    done

    # remove day directory if empty
    if [ ! "$(ls -A ${data_destbase}${years}${julday})" ]; then
        rm -rf ${data_destbase}${years}${julday}
    fi
    if [ ! "$(ls -A ${data_destbase}${years_end}${julday_end})" ]; then
        rm -rf ${data_destbase}${years_end}${julday_end}
    fi
done
