#!/bin/bash
# ------------------------------------------------------------
# This script runs Serious Game DRIFTERS test cases:
# simulations are performed with the 4th order Rung-Kutta 
# numerical scheme.
#
# (D. Bruciaferri)
# ------------------------------------------------------------

OUT_DIR="../OUTPUT/RUNGE_SG_drifter_test_case"
mkdir ${OUT_DIR}

for i in {1,3}; do

    if [ $i == 1 ]; then
       WAV_S="NONE"
       WAV="NONE"
       exp="cur"
       rm  medslik5.par
       cp medslik5_0.par medslik5.par
    fi
    if [ $i == 2 ]; then
       WAV_S="MEDSLIK-II JONSWAP PARAMETERIZATION"
       WAV="NONE"
       exp="cur_stok"
       rm  medslik5.par
       cp medslik5_1.par medslik5.par
    fi
    if [ $i == 3 ]; then
       WAV_S="MFS-WW3 6.5 km"
       WAV="MFS-WW3"
       exp="cur_wav"
       rm  medslik5.par
       cp medslik5_2.par medslik5.par
    fi

    echo "RUNGE-KUTTA DRIFTER EXPERIMENT:"
    echo "==============================="
    echo "OCEAN: MFS 6.5 km "
    echo "METEO: ECMWF 25 km"
    echo "WAVE:  "${WAV_S}
    echo ""

    OCE="MFS-COP01hm"
    MET="ECMWF025"

    for j in {1,4,6,7,8,9}; do

        rm medslik_inputfile.txt
        sim_name="test_case_drft_${j}_v1.02_"${exp}

        if [ $j == 1 ]; then
           dur="0024"
           hr="11"
           mn="38"
           lat="42.97053"
           lon="9.98382"
        fi

        if [ $j == 4 ]; then
           dur="0024"
           hr="11"
           mn="39"
           lat="42.96968"
           lon="9.98383"
        fi

        if [ $j == 6 ]; then
           dur="0048"
           hr="11"
           mn="39"
           lat="42.96907"
           lon="9.98350"
        fi 

        if [ $j == 7 ]; then
           dur="0024"
           hr="12"
           mn="08"
           lat="42.97157"
           lon="9.98055"
        fi
        if [ $j == 8 ]; then
           dur="0048"
           hr="12"
           mn="00"
           lat="42.96943"
           lon="9.98105"
        fi
        if [ $j == 9 ]; then
           dur="0024"
           hr="12"
           mn="04"
           lat="42.96938"
           lon="9.97975"
        fi
        
        INP_FILE="./medslik_inputfile.txt"

        cat > ${INP_FILE} << EOF

SIM_NAME=${sim_name}
RESTART=0000
TRACKMODE=1
OCEAN=${OCE}
WIND=${MET}
WAVE=${WAV}
SIM_LENGTH=${dur}
NUM_SCHEME=1
GRID_SIZE=150.0
OUTPUT_NAME=dr${j}_
STEP_OUTPUT=1

NSLICK=1
SAT_DATA=NO        
NAMEFILE_GML=      
CONTOURSLICK=NO   
OIL=API
OIL_TYPE=220.0
AGE=0

S1DD=17
S1MM=05
S1YY=14
S1HR=${hr}
S1MN=${mn}
S1durath=0
S1spllrt=0.
S1lon[1]=${lon}
S1lat[1]=${lat}
EOF

        ./medslik_II.sh 
        mv ../OUTPUT/*${sim_name} ${OUT_DIR}
    
    done

done
