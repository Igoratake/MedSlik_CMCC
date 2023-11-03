#!/bin/bash

#================================================================================
#  MEDSLIK-II_1.02                                                               |
#                                                                                |
#  Oil spill fate and transport numerical model                                  |
#--------------------------------------------------------------------------------|
#  medslik_II.sh                                                                 |
#                                                                                |
#  This script coordinates the model run                                         |
#--------------------------------------------------------------------------------|
#                                                                                |
#  Copyright (C) <2012>                                                          |
#                                                                                |
#  This program was originally written by Robin Lardner and George Zodiatis.     |
#                                                                                |
#  Subsequent additions and modifications have been made by Michela De Dominicis |
#  and Diego Bruciaferri.                                                        |
#                                                                                |
#  michela.dedominicis@ingv.it                        diego.bruciaferri@ingv.it  |
#                                                                                |
#--------------------------------------------------------------------------------|
#  The development of the MEDSLIK-II model is supported by a formal agreement    |
#  Memorandum of Agreement for the Operation and Continued Development of        |
#  MEDSLIK-II signed by the following institutions:                              |
#                                                                                |
#  INGV     - Istituto Nazionale di Geofisica e Vulcanologia                     |
#  OC-UCY   - Oceanography Center at the University of Cyprus                    |
#  CNR-IAMC - Consiglio Nazionale delle Ricerche – Istituto per                  |
#             lo Studio dell’Ambiente Marino Costiero                            |
#  CMCC     - Centro Euro-Mediterraneo sui Cambiamenti Climatici                 |
#--------------------------------------------------------------------------------|  
#  This program is free software: you can redistribute it and/or modify          |
#  it under the terms of the GNU General Public License as published by          |
#  the Free Software Foundation, either version 3 of the License, or             |
#  any later version.                                                            |
#                                                                                |
#  This program is distributed in the hope that it will be useful,               |
#  but WITHOUT ANY WARRANTY; without even the implied warranty of                |
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                 |
#  GNU General Public License for more details.                                  |
#  You should have received a copy of the GNU General Public License             |
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.         |
#================================================================================

#================================================================================
# 1. DECLARATION OF ENVIRONMENT VARIABLES OF THE MODEL
#================================================================================

HOME_MEDSLIK=${PWD%/RUN*}
HOME_DATA=${HOME_MEDSLIK}
EXE_FLDR=${HOME_MEDSLIK}"/EXE"

NCF_DATA=${HOME_DATA}"/DATA"
OCE_DATA=${NCF_DATA}"/OCE"
MET_DATA=${NCF_DATA}"/MET"
WAV_DATA=${NCF_DATA}"/WAV"

RUN_FLDR=${HOME_MEDSLIK}"/RUN"

INP_DATA=${RUN_FLDR}"/INP_DATA"
OCE_FLDR=${INP_DATA}"/OCE"
MET_FLDR=${INP_DATA}"/MET"
WAV_FLDR=${INP_DATA}"/WAV"

OUT_DATA=${HOME_MEDSLIK}"/OUTPUT"

PLOT_DIR=${HOME_MEDSLIK}"/PLOT"

INP_FILE="medslik_inputfile.txt"

if [ -f medslik5.inp ]; then
   rm medslik5.inp
fi
if [ -f medslik.tmp ]; then
   rm medslik.tmp
fi
if [ -f oil_file.txt ]; then
   rm oil_file.txt
fi
if [ -f medslik_plgs.inp ]; then
   rm medslik_plgs.inp
fi
if [ -f medslik_pnts.inp ]; then
   rm medslik_pnts.inp
fi

rm *.rso
rm flag*.tmp
rm tmp*.tmp
rm ${INP_DATA}/{OCE,MET,WAV}/*
rm files-{oce,met,wav}


source ${INP_FILE} 
python ${EXE_FLDR}/read_oil_data.py $OIL "$OIL_TYPE" 


#==============================================================================
# 2. POLYGON/SATELLITE OPTION
#==============================================================================

if [ "$CONTOURSLICK" == "YES" ] || [ "$SAT_DATA" == "YES" ]; then 
     isat=1
else  
     isat=0
fi 

#==============================================================================
# 3. CURRENTS, WINDS AND WAVES DATA
#
# ${OCEAN},${WIND},${WAVE} defined by ${INP_FILE}.  
#==============================================================================

# CURRENT----------------------- 

echo "OCEAN DATA CHOSEN FOR THE SIMULATION: $OCEAN"

# DAYLY DATA
if [ "$OCEAN" == "MFS-24dm" ]; then
    currents=14
    region=medf
    res_time_cu="dm"
fi
# HOURLY DATA
if [ "$OCEAN" == "MFS-01hm" ]; then
    currents=76
    region=medf
    res_time_cu="hm"
fi
if [ "$OCEAN" == "MFS-COP01hm" ]; then
    currents=77
    region=medf
    res_time_cu="hm"
fi


# WIND ----------------------------------

echo "METEO DATA CHOSEN FOR THE SIMULATION: $WIND"

if [ "$WIND" == "ECMWF025" ]; then
    wind=11
    res_time_wd="hi"
fi
if [ "$WIND" == "ECMWF05" ]; then
    wind=12
    res_time_wd="hi"
fi
if [ "$WIND" == "ECMWF0125" ]; then
    wind=11
    res_time_wd="hi"
fi

# WAVE -----------------------------------

echo "WAVE DATA CHOSEN FOR THE SIMULATION: $WAVE"

if [ "$WAVE" == "NONE" ]; then
    wave=000
fi
if [ "$WAVE" == "MFS-WW3" ]; then
    wave=101
    res_time_wv="hi"
fi


#==================================================================================
# 4. Date, time and oil spill positions defined by ${INP_FILE}
#==================================================================================

if [ $isat -eq 0 ] || [ "$CONTOURSLICK" == "YES" ]; then

   day=$S1DD                   
   month=$S1MM               
   year=$S1YY                
   hour=$S1HR                 
   minutes=$S1MN
fi           

echo "$day" " $month" "$year"

iage=$AGE

if [ $isat -eq 0 ]; then
   n=1
   echo "$NSLICK" '      Number of Total Spill Sources' >> medslik_pnts.inp
   while [ $n -le $NSLICK ]; do

         var_name_dur=$`echo '{S'$n'durath}'`
         apply_name_dur="Sdur=$var_name_dur"
         eval $apply_name_dur

         var_name_lon=$`echo '{S'$n'lon[1]}'`
         apply_name_lon="Slon=$var_name_lon"
         eval $apply_name_lon

         var_name_lat=$`echo '{S'$n'lat[1]}'`
         apply_name_lat="Slat=$var_name_lat"
         eval $apply_name_lat

         var_name_spl=$`echo '{S'$n'spllrt}'`
         apply_name_spl="Sspl=$var_name_spl"
         eval $apply_name_spl

         if [ $NSLICK -gt 1 ]; then
            var_name_day=$`echo '{S'$n'DD}'`
            apply_name_day="day=$var_name_day"
            eval $apply_name_day
     
            var_name_month=$`echo '{S'$n'MM}'`
            apply_name_month="month=$var_name_month"
            eval $apply_name_month

            var_name_hour=$`echo '{S'$n'HR}'`
            apply_name_hour="hour=$var_name_hour"
            eval $apply_name_hour

            var_name_minutes=$`echo '{S'$n'MN}'`
            apply_name_minutes="minutes=$var_name_minutes"
            eval $apply_name_minutes
         fi

         echo "$n" '      Spill Source' >> medslik_pnts.inp
         echo "$day" " $month" "$year" '      Date of Spill' >> medslik_pnts.inp
         echo "$hour$minutes" '      Hour of Spill' >> medslik_pnts.inp
         echo "$Sdur" '      Duration of Spill' >> medslik_pnts.inp
         echo "$Slat" '      Latitude of Spill' >> medslik_pnts.inp
         echo "$Slon" '      Longitude of Spill' >> medslik_pnts.inp
         echo "$Sspl" '      Spill Rate (units/hr)' >> medslik_pnts.inp
         n=`expr $n + 1`
   done
fi

if [ "$CONTOURSLICK" == "YES" ]; then

   Sdur=$S1durath
   Sspl=$S1spllrt

   echo "$day" " $month" "$year" '      Date of Spill' >> medslik_plgs.inp
   echo "$hour$minutes" '      Hour of Spill' >> medslik_plgs.inp
   echo "$Sdur" '      Duration of Spill' >> medslik_plgs.inp
   echo "$NSLICK" '      Number of Total Oil Slicks' >> medslik_plgs.inp   

   n=1
   p=0
   N=1

   while [ $n -le $NSLICK ]; do

         var_name_spl=$`echo '{S'$n'spllrt}'`
         apply_name_spl="Sspl=$var_name_spl"
         eval $apply_name_spl

         echo "$Sspl" 'Total oil of the' "$n" 'slick' >> medslik_plgs.inp
         n=`expr $n + 1`
   done

   while [ $N -le $NSLICK ]; do

         var_name=$`echo '{#S'$N'lon[*]}'`
         apply_count="element_count=$var_name"
         eval $apply_count

         index=1
         while [ "$index" -le "$element_count" ]; do
               var_name_lon=$`echo '{S'$N'lon[$index]}'`
               apply_name_lon="Slon[$index]=$var_name_lon"
               eval $apply_name_lon

               var_name_lat=$`echo '{S'$N'lat[$index]}'`
               apply_name_lat="Slat[$index]=$var_name_lat"
               eval $apply_name_lat

               echo ${Slat[$index]} ${Slon[$index]} >> initial0.txt
               index=`expr  $index + 1`
               p=`expr $p + 1`
         done
         N=`expr $N + 1`
         p=`expr $p + 1`
         echo ${Slat[1]} ${Slon[1]} >> initial0.txt
   done
   echo $p'        Number of data points' >> medslik_plgs.inp
   echo '  lat     lon' >> medslik_plgs.inp
   cat initial0.txt>>medslik_plgs.inp
   rm initial0.txt
fi

if [ "$SAT_DATA" == "YES" ]; then

    python ${EXE_FLDR}/ReadSatData.py $namefileGML $NSLICK

fi

#===============================================================================
# 6. SAVE INPUT DATA 
#===============================================================================
sim_length=$SIM_LENGTH

if [ $iage -eq 24 ]; then
   Day_24=20$year$month$day
   Day_0=`${EXE_FLDR}/jday $Day_24 -1`
   year=`echo $Day_0 |cut -c3-4`
   month=`echo $Day_0 |cut -c5-6`
   day=`echo $Day_0 |cut -c7-8`
   sim_length=`expr $sim_length + 24`  
fi

if [ $iage -eq 48 ]; then
   Day_24=20$year$month$day
   Day_0=`${EXE_FLDR}/jday $Day_24 -2`
   year=`echo $Day_0 |cut -c3-4`
   month=`echo $Day_0 |cut -c5-6`
   day=`echo $Day_0 |cut -c7-8`
   sim_length=`expr $sim_length + 48`
fi
# NUMBER OF FILES NEEDED --------------------------------------------------------- 

numfiles=`expr $sim_length / 24 + 1`
#int1=`expr $numfiles \* 24`
#int2=`expr $sim_length / 1`
int1=`expr $sim_length / 24`
#echo 'int1:' $int1
#int2=`expr $sim_length / 1`
int2=`expr $int1 \* 24`
#echo 'int2:' $int2
int3=`expr $sim_length / 1`
#echo 'sim_length:' $int3

if [ $currents = 76 ]; then
   if [ $hour -lt 13 ] || [ $int1 -lt $int2 ]; then
      numfiles=`expr $numfiles + 1`
   fi
fi
if [ $currents = 77 ] || [ $currents = 14 ]; then
   if [ $hour -ge 13 ] || [ $int2 -lt $int3 ]; then
        numfiles=`expr $numfiles + 1`
   fi
fi

numfiles_wd=$numfiles
numfiles_wv=$numfiles

if [ $iage -eq 24 ]; then
     numfiles=`expr $numfiles + 1`     
     numfiles_wd=`expr $numfiles_wd + 1`
     numfiles_wv=`expr $numfiles_wv + 1`
fi

# INITIAL DATE --------------------------------------------------------------------
    
FCStart1=20$year$month$day

FcStart1=$FCStart1 # currents
FcStart2=$FCStart1 # wind
FcStart3=$FCStart1 # waves

# MFS-24dm case
if [ $currents = 14 ] && [ $hour -ge 13 ]; then
   FcStart1=`${EXE_FLDR}/jday $FCStart1 +1`
fi
# MFS-01hm case
if [ $currents = 76 ]; then
   if [ $hour -lt 13 ]; then
      FcStart1=`${EXE_FLDR}/jday $FCStart1 -1`
   fi
fi
# MFS-COP01h case
if [ $currents = 77 ] && [ $hour -ge 13 ]; then
   FcStart1=`${EXE_FLDR}/jday $FCStart1 +1`
fi

# MFS-WW3 case
if [ $wave = 101 ] && [ $hour -ge 13 ]; then
   FcStart3=`${EXE_FLDR}/jday $FCStart1 +1`
fi


FcStart1=`echo $FcStart1 |cut -c1-8` # currents
FcStart2=`echo $FcStart2 |cut -c1-8` # wind
FcStart3=`echo $FcStart3 |cut -c1-8` # wave

# WRITING the medslik5.inp file -----------------------------------------------------

cp medslikYYYY.inp medslik0.inp

if [ $RESTART = 0000  ]; then
   RESTART="0  0"
else
   RESTART="1  $RESTART"
fi

sed -e "s/region/$region/"\
    -e "s/restart/$RESTART/"\
    -e "s/scheme/$NUM_SCHEME/"\
    -e "s/trackmode/$TRACKMODE/"\
    -e "s/NSlick/$NSLICK/"\
    -e "s/day/$day/"\
    -e "s/month/$month/"\
    -e "s/year/$year/"\
    -e "s/hour/$hour/"\
    -e "s/minutes/$minutes/"\
    -e "s/name/$OUTPUT_NAME/"\
    -e "s/length/$sim_length/"\
    -e "s/step/$STEP_OUTPUT/"\
    -e "s/currents/$currents/"\
    -e "s/wind/$wind/"\
    -e "s/wave/$wave/"\
    -e "s/age/$iage/"\
    -e "s/sat/$isat/"\
    -e "s/oilgrid/$GRID_SIZE/"\
    -e "s/files_number/$numfiles/"\
    medslik0.inp>medslik1.inp
rm medslik0.inp

sed -n '1,17 p' medslik1.inp >medslik5.inp
cat oil_file.txt>>medslik5.inp
sed -n '26,29 p' medslik1.inp >medslik0.inp
cat medslik0.inp>>medslik5.inp
rm medslik[01].inp

#====================================================================
# 7. AREA SELECTION (MIN/MAX Longitudes & Latitudes)
#====================================================================

cp medslikYYYY.tmp medslik0.tmp

sed -e "s/REGION/$region/"\
    -e "s/CURRENTS/$currents/"\
    -e "s/WIND/$wind/"\
    -e "s/WAVES/$wave/"\
     medslik0.tmp>medslik1.tmp
rm medslik0.tmp

if [ -f medslik_pnts.inp ]; then
   ${EXE_FLDR}/lat_lon.exe medslik_pnts.inp
else
   ${EXE_FLDR}/lat_lon.exe medslik_plgs.inp
fi
mv medslik1.tmp medslik.tmp

#====================================================================
# 8. PRE-PROCESSING OF CURRENTS & WIND FILES NEEDED FOR SIMULATION
#====================================================================

# DAYLY CURRENT FIELDS-----------------------------------------------

if [ $currents = 14 ]; then
   n=0
#   loop=`expr $numfiles + 1`
   loop=$numfiles
else
   n=0
   loop=$numfiles
fi

while [ $n != $loop ]; do

      n=`expr $n + 1`
      nn=`expr $n - 1`

#--------------------------------------------------------------------
#                          OCEAN-CURRENTS
#--------------------------------------------------------------------

      Datafc=`${EXE_FLDR}/jday $FcStart1 +$nn`
      DataFC=`echo $Datafc |cut -c3-8`
      DataFc=$DataFC
      echo ""
      echo '== CHECK CURRENTS == ' $DataFc
      echo ""
      DataFc_out=${DataFc}

      if [ $currents = 76 ]; then
         file_U=`ls ${OCE_DATA}/'MEDffE'${DataFc}'_01_grid_U'*`
         file_V=`ls ${OCE_DATA}/'MEDffE'${DataFc}'_01_grid_V'*`
         file_T=`ls ${OCE_DATA}/'MEDffE'${DataFc}'_01_grid_T'*`
#         cmd_U="ncks -O -d depth,0,9 $file_U $file_U"
#         cmd_V="ncks -O -d depth,0,9 $file_V $file_V"
#         cmd_T="ncks -O -d depth,0,9 $file_T $file_T"
         echo $res_time_cu > files-oce
      fi

      if [ $currents = 14 ] || [ $currents = 77 ]; then
         file_U=`ls ${OCE_DATA}/'20'${DataFc}'_'${res_time_cu}'-INGV--RFVL-'*`
         file_V=`ls ${OCE_DATA}/'20'${DataFc}'_'${res_time_cu}'-INGV--RFVL-'*`
         file_T=`ls ${OCE_DATA}/'20'${DataFc}'_'${res_time_cu}'-INGV--TEMP-'*`
         cmd_U="ncks -O -d depth,0,9 $file_U $file_U"
         cmd_V="ncks -O -d depth,0,9 $file_V $file_V"
         cmd_T="ncks -O -d depth,0,9 $file_T $file_T"
         echo $res_time_cu > files-oce
      fi

     
      if [ -n "${file_U}" ]; then
         if [ -n "${cmd_U}" ]; then   
            echo ${cmd_U}
            eval ${cmd_U}
         fi
         if [ -e ${OCE_DATA}/${DataFc}'_U.nc' ]; then
            echo ${file_U} ": File exists"
            echo "${DataFc}_U.nc 1" >> tmp1.tmp
         else
            ln -s ${file_U} ${OCE_DATA}/${DataFc}'_U.nc'
            echo ${file_U} ": File exists"
            echo "${DataFc}_U.nc 1" >> tmp1.tmp
         fi
      else
         echo "${DataFc}_U.nc 0" >> tmp1.tmp
         echo "For this run you need the U-CURRENTS file for the date "${DataFc}
      fi

      if [ -n "${file_V}" ]; then
         if [ -n "${cmd_V}" ]; then
            echo ${cmd_V}
            eval ${cmd_V}
         fi
         if [ -e ${OCE_DATA}/${DataFc}'_V.nc' ]; then
            echo ${file_V} ": File exists"
            echo "${DataFc}_V.nc 1" >> tmp1.tmp
         else
            ln -s ${file_V} ${OCE_DATA}/${DataFc}'_V.nc'
            echo ${file_V} ": File exists"
            echo "${DataFc}_V.nc 1" >> tmp1.tmp
         fi
      else
         echo "${DataFc}_V.nc 0" >> tmp1.tmp
         echo "For this run you need the V-CURRENTS file for the date "${DataFc}
      fi

      if [ -n "${file_T}" ]; then
         if [ -n "${cmd_T}" ]; then
            echo ${cmd_T}
            eval ${cmd_T}
         fi
         if [ -e ${OCE_DATA}/${DataFc}'_T.nc' ]; then
            echo ${file_T} ": File exists"
            echo "${DataFc}_T.nc 1" >> tmp1.tmp
         else
            ln -s ${file_T} ${OCE_DATA}/${DataFc}'_T.nc'
            echo ${file_T} ": File exists"
            echo "${DataFc}_T.nc 1" >> tmp1.tmp
         fi
      else
         echo "${DataFc}_T.nc 0" >> tmp1.tmp
         echo "For this run you need the T-CURRENTS file for the date "${DataFc}
      fi

#--------------------------------------------------------------------
#                          METEO-WINDS
#--------------------------------------------------------------------

      Datafc_wd=`${EXE_FLDR}/jday $FcStart2 +$nn`
      DataFC_wd=`echo $Datafc_wd |cut -c3-8`
      DataFc_wd=$DataFC_wd

      Data_wind="20"${DataFc_wd}
      File_wind="met_"${DataFc_wd}".met"

      echo ""
      echo '== CHECK WINDS == ' $Data_wind
      echo ""

      if [ $wind = 11 ] || [ $wind = 12 ]; then
         file_WIN=`ls ${MET_DATA}/${Data_wind}'_6hi-INGV-ECMWF-AM0'*'.nc'`
         echo "file_WIN =" $file_WIN
         echo $res_time_wd > files-met
      fi


      if [ -n "${file_WIN}" ]; then
         if [ -e ${MET_DATA}/${Data_wind}'_WIND.nc' ]; then
            echo "${File_wind} 1" >> tmp2.tmp
            echo ${file_WIN} ": File exists" 
         else
            ln -s ${file_WIN} ${MET_DATA}/${Data_wind}'_WIND.nc'
            echo "${File_wind} 1" >> tmp2.tmp
            echo ${file_WIN} ": File exists" 
         fi
      else
         echo "${File_wind} 0" >> tmp2.tmp
         echo "For this run you need the METEO file for the date " ${Data_wind}

      fi

#--------------------------------------------------------------------
#                          OCEAN-WAVES
#--------------------------------------------------------------------

      Datafc_wv=`${EXE_FLDR}/jday $FcStart3 +$nn`
      DataFC_wv=`echo $Datafc_wv |cut -c3-8`
      DataFc_wv=$DataFC_wv

      Data_wave="20"${DataFc_wv}
      File_wave="wav_"${DataFc_wv}".wav"

      echo ""
      echo "== CHECK WAVES ==" $Data_wave
      echo ""

      if [ $wave = 101 ]; then
         file_WAV=`ls ${WAV_DATA}'/'${Data_wave}'_hi-INGV--WAVE-'*`
         echo $res_time_wv > files-wav
      fi
     
      if [ $wave = 101 ]; then
         if [ -n "${file_WAV}" ]; then
            if [ -e ${WAV_DATA}/${Data_wave}'_WAVE.nc' ]; then
               echo "${Data_wave}_WAVE.nc 1" >> tmp3.tmp
               echo ${file_WAV} ": File exists"
            else 
               ln -s ${file_WAV} ${WAV_DATA}/${Data_wave}'_WAVE.nc'
               echo "${Data_wave}_WAVE.nc 1" >> tmp3.tmp
               echo ${file_WAV} ": File exists"
            fi 
         else
            echo "${Data_wave}_WAVE.nc 0" >> tmp3.tmp
            echo "For this run you need the WAVE file" ${Data_wave}.nc
         fi
      fi
done

#------------------------------------------------------------------------------------

echo $numfiles >> medslik5.inp
cat tmp1.tmp>>medslik5.inp
echo $numfiles_wd >> medslik5.inp
cat tmp2.tmp>>medslik5.inp
if [ $wave = 101 ]; then
   echo $numfiles_wv >> medslik5.inp
   cat tmp3.tmp>>medslik5.inp
fi

# OCEAN -------------------------

echo $numfiles >> medslik.tmp

if [ $currents = 76 ]; then
   numfiles_tmp=`expr $numfiles + 1`
else
   numfiles_tmp=$numfiles
fi

if [ $currents = 77 ]; then
   Dataaux=`${EXE_FLDR}/jday $FcStart1 -1`
   DataAUX=`echo $Dataaux |cut -c3-8`
   echo ${DataAUX}24 >> medslik.tmp
fi

n=0

while [ $n != $numfiles_tmp ]; do
      n=`expr $n + 1`
      nn=`expr $n - 1`
      Datafc=`${EXE_FLDR}/jday $FcStart1 +$nn`
      DataFC=`echo $Datafc |cut -c3-8`
      DataFc=$DataFC
      echo ${DataFc}24 >> medslik.tmp
done

# WIND --------------------------

echo $numfiles_wd >> medslik.tmp
n=0
while [ $n != $numfiles_wd ]; do
      n=`expr $n + 1`
      nn=`expr $n - 1`
      Datafc=`${EXE_FLDR}/jday $FcStart2 +$nn`
      echo ${Datafc} >> medslik.tmp
done

# WAVE --------------------------
if [ $wave = 101 ]; then
   echo $numfiles_wv >> medslik.tmp
   if [ $wave = 101 ]; then
      Dataaux=`${EXE_FLDR}/jday $FcStart3 -1`
      echo ${Dataaux} >> medslik.tmp
   fi

   n=0

   while [ $n != $numfiles_wv ]; do
         n=`expr $n + 1`
         nn=`expr $n - 1`
         Datafc=`${EXE_FLDR}/jday $FcStart3 +$nn`
         echo ${Datafc} >> medslik.tmp
   done
fi
echo " 0" >> medslik.tmp 
# -------------------------------

tr -d "\015" < medslik5.inp > medslik52.inp
cp medslik52.inp medslik5.inp
rm medslik52.inp
if [ -f medslik_multi.txt ]; then
   cat medslik_multi.txt>>medslik5.inp 
fi
tr -d "\015" < medslik.tmp > medslik2.tmp
cp medslik2.tmp medslik.tmp
rm medslik2.tmp


#====================================================================
# 8. EXTRACT CURRENTS, WIND AND SST DATA
#====================================================================

echo ''
echo '                               ********************************'
echo '                               *                              *'
echo '                               * READING CURRENTS, WIND and   *'
echo '                               *         WAVE DATA            *'
echo '                               *                              *'
echo '                               ********************************'
echo ''

${EXE_FLDR}/Extract_II.exe ${NCF_DATA}

#===================================================================
# 9. RUN 
#====================================================================

echo ''
echo '                                  ***************************'
echo '                                  *                         *'
echo '                                  * PLEASE WAIT: MEDSLIK-II *'
echo '                                  *  SIMULATION IS RUNNING  *'
echo '                                  *                         *'
echo '                                  ***************************'
echo ''

${EXE_FLDR}/medslik_II.exe

#!!!!!!!!!!!!! fino a qua
#####################################################################
#9. ARCHIVE OUTPUT FILES
#####################################################################

rm tmp*.tmp
model=MEDSLIKII
DIR_output=$model'_20'$year'_'$month'_'$day'_'$hour$minutes'_'$SIM_NAME

mkdir $OUT_DATA/$DIR_output
cp medslik5.inp  $OUT_DATA/$DIR_output
cp medslik5.par  $OUT_DATA/$DIR_output
cp medslik.tmp $OUT_DATA/$DIR_output

cp medslik_inputfile.txt $OUT_DATA/$DIR_output
cp files-{oce,met} $OUT_DATA/$DIR_output
if [ $wave -ne 000 ]; then
   cp files-wav $OUT_DATA/$DIR_output
fi
if [ -f medslik_plgs.inp ]; then
   cp medslik_plgs.inp $OUT_DATA/$DIR_output
fi
if [ -f medslik_pnts.inp ]; then
   cp medslik_pnts.inp $OUT_DATA/$DIR_output
fi
# agginugere medslik_plgs.inp medslik_pnts.inp

mv *.srf $OUT_DATA/$DIR_output
mv *.cst $OUT_DATA/$DIR_output
mv *.dsp $OUT_DATA/$DIR_output
mv *.fte $OUT_DATA/$DIR_output
mv *.trj $OUT_DATA/$DIR_output


cp -r $INP_DATA/OCE/ $OUT_DATA/$DIR_output/
cp -r $INP_DATA/MET/ $OUT_DATA/$DIR_output/
if [ $wave -ne 000 ]; then
   cp -r $INP_DATA/WAV/ $OUT_DATA/$DIR_output/
fi

#######################################################################
#10. CREATE THE MEDESS NETCDF OUTPUT FILE
#######################################################################

if [ $TRACKMODE -eq 0 ]; then

   ${EXE_FLDR}/create_output.exe $OUT_DATA/$DIR_output/ 
   mv coast.map $OUT_DATA/$DIR_output

#######################################################################
#11. CREATE SLICK OUTPUT PLOTS
#######################################################################

   if [ "$SIM_NAME" == "ALGERIA_testcase" ]; then

        python $PLOT_DIR/plot_slick_ALGERIA.py $OUT_DATA/$DIR_output/ ${OUTPUT_NAME}.nc ${wave}

   elif [ "$SIM_NAME" == "LEBANON_testcase" ]; then

        python $PLOT_DIR/plot_slick_LEBANON.py $OUT_DATA/$DIR_output/ ${OUTPUT_NAME}.nc ${wave}

   elif [ "$SIM_NAME" == "SRG_testcase_eul" ] || [ "$SIM_NAME" == "SRG_testcase_run" ]; then

        python $PLOT_DIR/plot_slick_SRG.py $OUT_DATA/$DIR_output/ ${OUTPUT_NAME}.nc ${wave}

   else
      
        python $PLOT_DIR/MEDSLIKII_plot_slick.py $OUT_DATA/$DIR_output/ ${OUTPUT_NAME}.nc ${wave}
       
   fi

   python $PLOT_DIR/MEDSLIKII_plot_oilfate.py $OUT_DATA/$DIR_output/
   mv *.png $OUT_DATA/$DIR_output
fi

if [ $TRACKMODE -eq 1 ]; then

   tst=`echo $SIM_NAME | cut -c1-14`

   if [ "$tst" == "test_case_drft" ]; then

      echo "TEST CASE DRIFTER EXPERIMENT"

   else
 
      python $PLOT_DIR/MEDSLIKII_plot_traj.py $OUT_DATA/$DIR_output/ ${OUTPUT_NAME}${year}${month}${day}_${hour}${minutes}_F.trj $SIM_LENGTH 
      mv trajectory.png $OUT_DATA/$DIR_output

   fi

fi
