#!/bin/bash

HOME_MEDSLIK=$PWD
MEDSLIK=$HOME_MEDSLIK/EXE
export NCARG_ROOT=$HOME_MEDSLIK/ncl
export NCARG_USRRESFILE=$HOME_MEDSLIK/.hluresfile
export PATH=$HOME_MEDSLIK/ncl/bin:$PATH




rm $MEDSLIK/status.txt
cd $MEDSLIK

source medslik_II.sh 
echo 'VISUALIZATION IS RUNNING'
cd medslik_plots
sh medslik_plots.sh $DIR_output

chmod -R 777 $MEDSLIK"/output/"$DIR_output


if [ -f $MEDSLIK"/output/"$DIR_output"/plots/index.html" ]; then
echo "FINISH" >>   $MEDSLIK/status.txt
else
echo "ERROR" >>   $MEDSLIK/status.txt
fi



