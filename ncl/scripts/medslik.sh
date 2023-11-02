#source paths.sh

MEDSLIK_OUTPUT_DIR=/home/medslik/EXE/output/

if [ "$1" != "" ]; then
	MEDSLIK_OUTPUT_DIR=$MEDSLIK_OUTPUT_DIR$1"/"
fi

if [ ! -e $MEDSLIK_OUTPUT_DIR ]; then
	echo "There is no directory "$MEDSLIK_OUTPUT_DIR
	exit
fi

export MEDSLIK_OUTPUT_DIR=$MEDSLIK_OUTPUT_DIR
cd `dirname $0`
ncl -n main.ncl 
