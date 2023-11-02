# Launcher for MEDSLIK output visualisation
# 
#-----------------------------------------------------------------------------------
#  Copyright (C) <2012>
#  This program was originally written
#  by INGV, Slava Lyubartsev, 2011.02.22 
#----------------------------------------------------------------------------------
#  The development of the MEDSLIK-II model is supported by a formal agreement
#  Memorandum of Agreement for the Operation and Continued Development of MEDSLIK-II
#  signed by the following institutions:
#  INGV - Istituto Nazionale di Geofisica e Vulcanologia
#  OC-UCY - Oceanography Center at the University of Cyprus
#  CNR-IAMC - Consiglio Nazionale delle Ricerche – Istituto per 
#  lo Studio dell’Ambiente Marino Costiero
#  CMCC - Centro Euro-Mediterraneo sui Cambiamenti Climatici
# 
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------------
HOME_MEDSLIK=$HOME/MEDSLIK_II_1.01
MEDSLIK=$HOME_MEDSLIK/EXE
MEDSLIK_DIR=$MEDSLIK/output
INPUT_DIR=$MEDSLIK_DIR"/"$1
BASE_DIR=`dirname $0`

if [ -z "$1" ] || [ ! -e $INPUT_DIR ] ; then
	echo
	echo "Usage: "`basename $0`" INPUT_DIR [ PLOT_DIR [ TASK_NCL ] ]"
	echo
	echo "MEDSLIK_DIR is defined here as "$MEDSLIK_DIR
	echo "It must be the MEDSLIK output directory, edit this shell script if necessary"
	echo
	echo "The mandatory parameter INPUT_DIR is a directory name"
	echo "This directory must exist in MEDSLIK_DIR directory"
	echo
	echo "If PLOT_DIR is defined, then plots will be placed into "
	echo "MEDSLIK_DIR/INPUT_DIR/plots/PLOT_DIR/"
	echo "Otherwise plots will be placed into MEDSLIK_DIR/OUTPUT_DIR/plots/"
	echo
	echo "File medslik_plots.ncl in MEDSLIK_DIR/INPUT_DIR/ controls the plotting procedure"
	echo "It will be created, if absent, and will contain default parameters"
	echo "Default settings are in "$BASE_DIR"/medslik_plots.ncl"
	echo
	echo "If TASK_FILE is defined, then it will be used to control plotting procedure"
	echo "Otherwise TASK_FILE is MEDSLIK_DIR/INPUT_DIR/medslik_plots.ncl"
	echo
	echo "Parameters PLOT_DIR and TASK_NCL can be convenient for the batch visualisation"
	if [ ! -e $DIR ] ; then
		echo 
		echo $DIR" directory does not exist!"
	fi
	exit
fi

cd $BASE_DIR  # all NCL scripts must be in the same directory as this shell script

export MEDSLIK_OUTPUT_DIR=$INPUT_DIR"/"  # pass the parameter to NCL scripts

echo
ncl -n main.ncl  # launch visualistion
