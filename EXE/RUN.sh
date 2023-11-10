#!/bin/bash

# set where you have placed MEDSLIK_II folder
HOME_MEDSLIK=HOME_MEDSLIK=${PWD%/EXE*}
MEDSLIK=$HOME_MEDSLIK/EXE

# generate update timelog file
start_time=$(date +%s.%N)
echo "`date +'%Y-%m-%d %H:%M:%S'` `date +%s` `hostname` start EXE" > ${MEDSLIK}/timelog.log

# launches model
source medslik_II_ens.sh 

# update timelog file
end_date=$(date +%s.%N)
exec_time=$(echo "(${end_date}-${start_time})" | bc)
echo "`date +'%Y-%m-%d %H:%M:%S'` `date +%s` `hostname` end EXE - Execution time: ${exec_time} sec. " >> ${MEDSLIK}/timelog.log

chmod -R 777 $MEDSLIK"/output/"$DIR_output

