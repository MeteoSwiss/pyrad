#!/bin/bash

# This script gets the rad4alp data from the CSCS
# repository and puts it in the right folder.
# To be run in CSCS

# set permits
umask 0002

# Config
dateCmd="/bin/date"

phdata_origbase=/store/msrad/radar/polarHR/data/
rawdata_origbase=/store/msrad/radar/swiss/data/

file_type=M

# default variables
START_TIME="000001"
END_TIME="240000"
time_vec='all'

ele_vec='all'

data_destbase=/store/msrad/radar/rad4alp/rawdata/
while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
        -r|--radar)
        RADAR="$2"
        OIFS=$IFS
        IFS=','
        read -r -a radar_vec <<< "$RADAR"
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
        -e|--res)
        res="$2"
        shift # past argument
        ;;
        --ele)
        ELE="$2"
        OIFS=$IFS
        IFS=','
        read -r -a ele_vec <<< "$ELE"
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
nrad=${#radar_vec[@]}
nele=${#ele_vec[@]}
ntime=${#time_vec[@]}

for ((irad=0; irad<${nrad}; irad++)); do
    radar=${radar_vec[${irad}]}
    echo "Processing radar "${radar}

    for ((iday=0; iday<${nday}; iday++)); do
        yearl=$(date --date "${day_vec[${iday}]}" +"%Y")
        years=$(date --date "${day_vec[${iday}]}" +"%y")
        julday=$(date --date "${day_vec[${iday}]}" +"%j")

        yearl_end=$(date -d "$(date --date "${day_vec[${iday}]}")+1 day" +"%Y")
        years_end=$(date -d "$(date --date "${day_vec[${iday}]}")+1 day" +"%y")
        julday_end=$(date -d "$(date --date "${day_vec[${iday}]}")+1 day" +"%j")

        echo "Processing day "${years}${julday}
        # transfer polar data from CSCS to destination folder, unzip it and remove zip file
        if [ "${res}" == "H" ]
        then
            data_origpath=${phdata_origbase}${yearl}/${years}${julday}/
            # check type of file that exists in the repository
            if [ `ls ${data_origpath}M${res}${radar}${years}${julday}.zip` ]
            then
                file_type=M
            elif [ `ls ${data_origpath}P${res}${radar}${years}${julday}.zip` ]
            then
                file_type=P
            else
                echo "No file found in repository"
                continue
            fi
            data_destpath=${data_destbase}${years}${julday}/${file_type}${res}${radar}${years}${julday}/
            filebase=${file_type}${res}${radar}${years}${julday}

            # create destination path for polar data and unzip files there
            mkdir -p ${data_destpath}
            if [ "${ele_vec}" = "all" ] && [ "${time_vec}" = "all" ];then
                unzip -o ${data_origpath}${filebase}.zip -d ${data_destpath}
            else
                for ((itime=0; itime<${ntime}; itime++)); do
                    if [ "${time_vec[${itime}]}" = "all" ]; then
                        time_rad=*
                    else
                        time_rad=${time_vec[${itime}]}
                    fi
                    if [ "${ele_vec}" = "all" ]; then
                        unzip -o ${data_origpath}${filebase}.zip ${filebase}${time_rad}*.* -d ${data_destpath}
                    else
                        for ((iele=0; iele<${nele}; iele++)); do
                            ele=${ele_vec[${iele}]}
                            unzip -o ${data_origpath}${filebase}.zip ${filebase}${time_rad}*.${ele} -d ${data_destpath}
                        done
                    fi
                done
            fi
            chmod -R gu+rw ${data_destpath}

            # add file 00:00 UTC the next day
            if [ "${time_vec}" = "all" ] || [ "$hour24" -eq 1 ];then
                data_origpath=${phdata_origbase}${yearl_end}/${years_end}${julday_end}/
                # check type of file that exists in the repository
                if [ `ls ${data_origpath}M${res}${radar}${years_end}${julday_end}.zip` ]
                then
                    file_type=M
                elif [ `ls ${data_origpath}P${res}${radar}${years_end}${julday_end}.zip` ]
                then
                    file_type=P
                else
                    echo "No file found in repository"
                    continue
                fi
                data_destpath=${data_destbase}${years_end}${julday_end}/${file_type}${res}${radar}${years_end}${julday_end}/
                filebase=${file_type}${res}${radar}${years_end}${julday_end}

                # create destination path for polar data and unzip files there
                mkdir -p ${data_destpath}
                if [ "${ele_vec}" = "all" ];then
                    unzip -o ${data_origpath}${filebase}.zip ${filebase}0000*.* -d ${data_destpath}
                else
                    for ((iele=0; iele<${nele}; iele++)); do
                        ele=${ele_vec[${iele}]}
                        unzip -o ${data_origpath}${filebase}.zip ${filebase}0000*.${ele} -d ${data_destpath}
                    done
                fi
                chmod -R gu+rw ${data_destpath}
            fi
        else
            data_origpath=${rawdata_origbase}${yearl}/${years}${julday}/
            # check type of file that exists in the repository
            if [ `ls ${data_origpath}M${res}${radar}${years}${julday}.zip` ]
            then
                file_type=M
            elif [ `ls ${data_origpath}P${res}${radar}${years}${julday}.zip` ]
            then
                file_type=P
            else
                echo "No file found in repository"
                continue
            fi
            data_destpath=${data_destbase}${years}${julday}/${file_type}${res}${radar}${years}${julday}/
            filebase=${file_type}${res}${radar}${years}${julday}

            # create destination path for polar data
            mkdir -p ${data_destpath}
            if [ "${ele_vec}" = "all" ] && [ "${time_vec}" = "all" ];then
                unzip -o ${data_origpath}${filebase}.zip -d ${data_destpath}
            else
                for ((itime=0; itime<${ntime}; itime++)); do
                    if [ "${time_vec[${itime}]}" = "all" ]; then
                        time_rad=*
                    else
                        time_rad=${time_vec[${itime}]}
                    fi
                    if [ "${ele_vec}" = "all" ]; then
                        unzip -o ${data_origpath}${filebase}.zip ${filebase}${time_rad}*.* -d ${data_destpath}
                    else
                        for ((iele=0; iele<${nele}; iele++)); do
                            ele=${ele_vec[${iele}]}
                            unzip -o ${data_origpath}${filebase}.zip ${filebase}${time_rad}*.${ele} -d ${data_destpath}
                        done
                    fi

                done
            fi
            chmod -R gu+rw ${data_destpath}

            # create destination path for hydrometeor classification data
            data_destpath=${data_destbase}${years}${julday}/YM${radar}${years}${julday}/
            filebase=YM${radar}${years}${julday}
            mkdir -p ${data_destpath}
            if [ "${ele_vec}" = "all" ] && [ "${time_vec}" = "all" ];then
                unzip -o ${data_origpath}${filebase}.zip -d ${data_destpath}
            else
                for ((itime=0; itime<${ntime}; itime++)); do
                    if [ "${time_vec[${itime}]}" = "all" ]; then
                        time_rad=*
                    else
                        time_rad=${time_vec[${itime}]}
                    fi
                    if [ "${ele_vec}" = "all" ]; then
                        unzip -o ${data_origpath}${filebase}.zip ${filebase}${time_rad}*.* -d ${data_destpath}
                    else
                        for ((iele=0; iele<${nele}; iele++)); do
                            ele="$(expr ${ele_vec[${iele}]} + 800)"
                            unzip -o ${data_origpath}${filebase}.zip ${filebase}${time_rad}*.${ele} -d ${data_destpath}
                        done
                    fi
                done
            fi
            chmod -R gu+rw ${data_destpath}

            # create destination path for dealised Doppler velocity data
            data_destpath=${data_destbase}${years}${julday}/DV${radar}${years}${julday}/
            filebase=DV${radar}${years}${julday}
            mkdir -p ${data_destpath}
            if [ "${ele_vec}" = "all" ] && [ "${time_vec}" = "all" ];then
                unzip -o ${data_origpath}${filebase}.zip -d ${data_destpath}
            else
                for ((itime=0; itime<${ntime}; itime++)); do
                    if [ "${time_vec[${itime}]}" = "all" ]; then
                        time_rad=*
                    else
                        time_rad=${time_vec[${itime}]}
                    fi
                    if [ "${ele_vec}" = "all" ]; then
                        unzip -o ${data_origpath}${filebase}.zip ${filebase}${time_rad}*.* -d ${data_destpath}
                    else
                        for ((iele=0; iele<${nele}; iele++)); do
                            ele="$(expr ${ele_vec[${iele}]} + 800)"
                            unzip -o ${data_origpath}${filebase}.zip ${filebase}${time_rad}*.${ele} -d ${data_destpath}
                        done
                    fi
                done
            fi
            chmod -R gu+rw ${data_destpath}

            # add file 00:00 UTC the next day
            if [ "${time_vec}" = "all" ] || [ "$hour24" -eq 1 ];then
                data_origpath=${rawdata_origbase}${yearl_end}/${years_end}${julday_end}/
                # check type of file that exists in the repository
                if [ `ls ${data_origpath}M${res}${radar}${years_end}${julday_end}.zip` ]
                then
                    file_type=M
                elif [ `ls ${data_origpath}P${res}${radar}${years_end}${julday_end}.zip` ]
                then
                    file_type=P
                else
                    echo "No file found in repository"
                    continue
                fi
                data_destpath=${data_destbase}${years_end}${julday_end}/${file_type}${res}${radar}${years_end}${julday_end}/
                filebase=${file_type}${res}${radar}${years_end}${julday_end}

                # create destination path for polar data
                mkdir -p ${data_destpath}
                if [ "${ele_vec}" = "all" ];then
                    unzip -o ${data_origpath}${filebase}.zip ${filebase}0000*.* -d ${data_destpath}
                else
                    for ((iele=0; iele<${nele}; iele++)); do
                        ele=${ele_vec[${iele}]}
                        unzip -o ${data_origpath}${filebase}.zip ${filebase}0000*.${ele} -d ${data_destpath}
                    done
                fi
                chmod -R gu+rw ${data_destpath}

                # create destination path for hydrometeor classification data
                data_destpath=${data_destbase}${years_end}${julday_end}/YM${radar}${years_end}${julday_end}/
                filebase=YM${radar}${years_end}${julday_end}
                mkdir -p ${data_destpath}
                if [ "${ele_vec}" = "all" ];then
                    unzip -o ${data_origpath}${filebase}.zip ${filebase}0000*.* -d ${data_destpath}
                else
                    for ((iele=0; iele<${nele}; iele++)); do
                        ele="$(expr ${ele_vec[${iele}]} + 800)"
                        unzip -o ${data_origpath}${filebase}.zip ${filebase}0000*.${ele} -d ${data_destpath}
                    done
                fi
                chmod -R gu+rw ${data_destpath}

                # create destination path for dealised Doppler velocity data
                data_destpath=${data_destbase}${years_end}${julday_end}/DV${radar}${years_end}${julday_end}/
                filebase=DV${radar}${years_end}${julday_end}
                mkdir -p ${data_destpath}
                if [ "${ele_vec}" = "all" ];then
                    unzip -o ${data_origpath}${filebase}.zip ${filebase}0000*.* -d ${data_destpath}
                else
                    for ((iele=0; iele<${nele}; iele++)); do
                        ele="$(expr ${ele_vec[${iele}]} + 800)"
                        unzip -o ${data_origpath}${filebase}.zip ${filebase}0000*.${ele} -d ${data_destpath}
                    done
                fi
                chmod -R gu+rw ${data_destpath}
            fi
        fi

        # create destination path for status data
        data_origpath=${rawdata_origbase}${yearl}/${years}${julday}/
        data_destpath=${data_destbase}${years}${julday}/ST${radar}${years}${julday}/
        filebase=ST${radar}${years}${julday}
        mkdir -p ${data_destpath}
        if [ "${time_vec}" = "all" ];then
            unzip -o ${data_origpath}${filebase}.zip ${filebase}*.xml -d ${data_destpath}
        else
            for ((itime=0; itime<${ntime}; itime++)); do
                time_rad=${time_vec[${itime}]}
                for ((iele=0; iele<${nele}; iele++)); do
                    unzip -o ${data_origpath}${filebase}.zip ${filebase}${time_rad}*.xml -d ${data_destpath}
                done
            done
        fi
        chmod -R gu+rw ${data_destpath}

        # add file 00:00 UTC the next day
        if [ "${time_vec}" = "all" ] || [ "$hour24" -eq 1 ];then
            data_origpath=${rawdata_origbase}${yearl_end}/${years_end}${julday_end}/
            data_destpath=${data_destbase}${years_end}${julday_end}/ST${radar}${years_end}${julday_end}/
            filebase=ST${radar}${years_end}${julday_end}
            mkdir -p ${data_destpath}
            unzip -o ${data_origpath}${filebase}.zip ${filebase}0000*.xml -d ${data_destpath}
            chmod -R gu+rw ${data_destpath}
        fi
    done
done
