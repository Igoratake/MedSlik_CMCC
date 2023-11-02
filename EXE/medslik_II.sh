#!/bin/sh
#-----------------------------------------------------------------------------------
#  MEDSLIK-II_1.01 
#  oil spill fate and transport model 
#-----------------------------------------------------------------------------------
#  medslik_II.sh
#  This script coordinates the model run
#-----------------------------------------------------------------------------------
#  Copyright (C) <2012>
#  This program was originally written
#  by Robin Lardner and George Zodiatis.
#  Subsequent additions and modifications
#  have been made by Michela De Dominicis. 
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


rm *.rso
rm medslik5.inp
rm medslik.tmp

rm output/*.srf
rm output/*.cst
rm output/*.dsp
rm output/*.tot
rm output/*.fte
rm flag*.tmp
rm tmp*.tmp

rm oil_file.txt
rm initial*.txt

HOME_MEDSLIK=$HOME/MEDSLIK_II_1.01
F_DATA=$HOME_MEDSLIK/DATA

source medslik_inputfile.txt
python read_oil_data.py $OIL "$OIL_TYPE" 


###############################################################################
#0.                               OPTIONS
###############################################################################
if [ "$ContourSlick" == "YES" ] || [ "$SAT_DATA" == "YES" ]; then 
isat=1; else  
isat=0
fi 

###############################################################################
#1.                               CURRENTS AND WINDS
###############################################################################
if [ "$MODEL" == "MFS" ]
then
currents=70
region=medf
fi
if [ "$MODEL" == "MFS24" ]
then
currents=10
region=medf
fi
if [ "$MODEL" == "AFS" ]
then
currents=72
region=adri
fi 
if [ "$MODEL" == "AFS24" ]
then
currents=11
region=adri
fi 
if [ "$MODEL" == "SCRM" ]
then
currents=71
region=sici
fi 
if [ "$MODEL" == "SCRM24" ]
then
currents=12
region=sici
fi 
if [ "$MODEL" == "TYRR" ]
then
currents=73
region=medf
fi 
if [ "$MODEL" == "TYRR24" ]
then
currents=13
region=medf
fi 
if [ "$MODEL" == "REL" ]
then
currents=74
region=medf
fi 
if [ "$MODEL" == "WME" ]
then
currents=75
region=medf
fi


if [ "$WIND" == "ECMWF05" ]
then
wind=06
fi
if [ "$WIND" == "ECMWF025" ]
then
wind=25
fi 


step_output=001       # DO NOT CHANGE IT! in hours, 3 character example: 001
output_name=out       # DO NOT CHANGE IT! 3 letters, example:out


if [ $isat -eq 0 ] || [ "$ContourSlick" == "YES" ]
then
restart=0
#hrestart=0
day=$day                   # 2 character example: 07
month=$month               # 2 character example: 08
year=$year                 # 2 character example: 08
hour=$hour                 # 2 character example: 09
minutes=$minutes           # 2 character example: 05
duration=$duration         # in hours, 4 character example: 0024
lat_degree=$lat_degree     # degrees, 2 character example: 06
lat_minutes=$lat_minutes   # minutes, example: 06.62
lon_degree=$lon_degree     # degrees, 2 character example: 06
lon_minutes=$lon_minutes   # minutes, example: 06.62
spillrate=$spillrate       # in tons/hours, example: 000055.00
spillvolume=$spillvolume   # in tons example: 000055.00
fi
iage=$age


###############################################################################
#2.                  READ CONTOUR from the input file
###############################################################################
if [ "$ContourSlick" == "YES" ]
then
echo '20'$year'/'$month'/'$day' '$hour':'$minutes >> initial.txt
p=0
N=1
while [ $N -le $NSlick ]
 do

 var_name=$`echo '{#S'$N'lon[*]}'`
 apply_count="element_count=$var_name"
 eval $apply_count

  index=1 
  while [ "$index" -le "$element_count" ]
  do  
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
  echo $p'        Number of data points' >> initial.txt
  echo '  lat     lon' >> initial.txt
  cat initial0.txt>>initial.txt
fi
rm initial0.txt
################################################################################
#3.                   READ INPUT DATA FROM SATELLITE DATA 
################################################################################
if [ "$SAT_DATA" == "YES" ]
then

python source/ReadSatData_EMSA.py $namefileGML $N_OS


myFile="medslik_sat.inp"
count=0
for var in `cat $myFile` ; do
count=`expr $count + 1`
   if [ $count -eq 2 ]; then  day=$var; fi
   if [ $count -eq 4 ]; then  month=$var; fi
   if [ $count -eq 6 ]; then  year=$var; fi
   if [ $count -eq 8 ]; then  hour=$var; fi
   if [ $count -eq 10 ]; then  minutes=$var; fi
   if [ $count -eq 12 ]; then  lat_degree=$var; fi
   if [ $count -eq 14 ]; then  lat_minutes=$var; fi
   if [ $count -eq 16 ]; then  lon_degree=$var; fi
   if [ $count -eq 18 ]; then  lon_minutes=$var; fi
   if [ $count -eq 20 ]; then  spillrate=$var; fi
   if [ $count -eq 22 ]; then  duration=$var; fi
   if [ $count -eq 24 ]; then  output_name=$var; fi
   if [ $count -eq 26 ]; then  step_output=$var; fi

done
fi
################################################################################
#4. SAVE INPUT DATA 
################################################################################

     if [ $iage -eq 24 ]
     then
        Day_24=20$year$month$day
        Day_0=`./jday $Day_24 -1`
        year=`echo $Day_0 |cut -c3-4`
        month=`echo $Day_0 |cut -c5-6`
        day=`echo $Day_0 |cut -c7-8`
        sim_length=`expr $sim_length + 24`  
     fi
     if [ $iage -eq 48 ]
     then
        Day_24=20$year$month$day
        Day_0=`./jday $Day_24 -2`
        year=`echo $Day_0 |cut -c3-4`
        month=`echo $Day_0 |cut -c5-6`
        day=`echo $Day_0 |cut -c7-8`
        sim_length=`expr $sim_length + 48`
     fi


multiple=01 #number of sumperimposed spills
numfiles=`expr $sim_length / 24 + 1`
int1=`expr $numfiles \* 24`
int2=`expr $sim_length / 1`


if [ $currents = 70 ] || [ $currents = 72 ] || [ $currents = 10 ] 
then
     if [ $hour -lt 13 ] || [ $int1 -lt $int2 ]
 
    then 
       numfiles=`expr $numfiles + 1`
     fi
fi

     if [ $iage -eq 24 ]
     then
        numfiles=`expr $numfiles + 1`     
     fi
     
     
cp source/medslikYYYY.inp medslik0.inp
sed -e "s/giorno/$day/"\
    -e "s/mese/$month/"\
    -e "s/anno/$year/"\
    -e "s/ora/$hour/"\
    -e "s/minuti/$minutes/"\
    -e "s/durata/$duration/"\
    -e "s/lat_gradi/$lat_degree/"\
    -e "s/lat_primi/$lat_minutes/"\
    -e "s/lon_gradi/$lon_degree/"\
    -e "s/lon_primi/$lon_minutes/"\
    -e "s/nome/$output_name/"\
    -e "s/lunghezza/$sim_length/"\
    -e "s/step/$step_output/"\
    -e "s/eta/$iage/"\
    -e "s/sat/$isat/"\
    -e "s/portata/$spillrate/"\
    -e "s/regione/$region/"\
    -e "s/correnti/$currents/"\
    -e "s/vento/$wind/"\
    -e "s/multi/$multiple/"\
    -e "s/riinizio/$restart/"\
    -e "s/griglia/$grid_size/"\
    -e "s/numero_files/$numfiles/"\
    medslik0.inp>medslik1.inp
rm medslik0.inp

sed -n '1,18 p' medslik1.inp >medslik5.inp
cat oil_file.txt>>medslik5.inp
sed -n '27,31 p' medslik1.inp >medslik0.inp
cat medslik0.inp>>medslik5.inp
rm medslik[01].inp

#####################################################################
#5. PRE-PROCESSING OF CURRENTS & WIND FILES NEEDED FOR SIMULATION
#####################################################################

FCStart1=20$year$month$day
if [ $currents = 70 ] || [ $currents = 72 ] || [ $currents = 10 ] 
then
 if [ $hour -lt 13 ]
 then
 FcStart1=`./jday $FCStart1 -1`
 else
 FcStart1=$FCStart1
fi
else 
FcStart1=$FCStart1
fi

FcStart1=`echo $FcStart1 |cut -c1-8`

if [ $currents = 70 ]
then
dir='O1h'
U='vozocrtx'
V='vomecrty'
T='votemper'
model='MFS'
pre_nameFC='MEDffE'
pre_nameAN='MEDaaE'
pre_nameSI='MEDssE'
tail_name='_01_grid'
fi

if  [ $currents = 71 ]
then
dir='S1h'
U='U'
V='V'
T='TEM'
model='SICILY'
pre_name='SCRMfsE'
tail_name='_01_grid'
fi

if  [ $currents = 72 ]
then
dir='A1h'
U='vozocrtx'
V='vomecrty'
T='votemper'
model='AFS'
pre_name='ADRIfaE'
tail_name='_01_grid'
fi

if  [ $currents = 73 ]
then
dir='T1h'
U='vozocrtx'
V='vomecrty'
T='votemper'
model='TYRRHENIAN'
pre_name='TIRRfpE'
tail_name='_01_grid'
fi

if  [ $currents = 74 ]
then
dir='H3k'
U='vozocrtx'
V='vomecrty'
T='votemper'
model='RELOCATABLE'
pre_name='HOPS3km_H_'
tail_name=''
fi


if  [ $currents = 75 ]
then
dir='WME'
U='U'
V='V'
T='TEM'
model='WEST_MED'
pre_name='WMEDfs'
tail_name='_01_grid'
fi

if  [ $currents = 10 ]
then 
dir='OPA'
U='vozocrtx'
V='vomecrty'
T='votemper'
model='MFS'
pre_name='MEDaaE'
tail_name='_24_grid'
fi

if  [ $currents = 11 ]
then
dir='A24'
U='vozocrtx'
V='vomecrty'
T='votemper'
model='AFS'
pre_name='ADRIssE'
tail_name='_24_grid'
fi

if  [ $currents = 12 ]
then
dir='S24'
U='U'
V='V'
T='TEM'
model='SICILY'
pre_name='SCRMaE'
tail_name='_24_grid'
fi

if  [ $currents = 13 ]
then
dir='T24'
U='vozocrtx'
V='vomecrty'
T='votemper'
model='TYRRHENIAN'
pre_name='TIRRaaE'
tail_name='_24_grid'
fi

fcst_data=fcst_data
FD=$F_DATA/$fcst_data/$dir
rm tmp*.tmp

n=0
while [ $n != $numfiles ]; do
n=`expr $n + 1`
nn=`expr $n - 1`
Datafc=`./jday $FcStart1 +$nn`
DataFC=`echo $Datafc |cut -c3-8`
DataFc=$DataFC
echo 'CHECK CURRENTS' $DataFc
DataFc_out=${DataFc}

if [ $currents = 70 ] & [ -f $FD/$pre_nameFC${DataFc}$tail_name'_U.nc' ] & [ -f $FD/$pre_nameFC${DataFc}$tail_name'_V.nc' ] & [ -f $FD/$pre_nameFC${DataFc}$tail_name'_T.nc' ]
then
pre_name=$pre_nameFC
fi
if [ $currents = 70 ] & [ -f $FD/$pre_nameSI${DataFc}$tail_name'_U.nc' ] & [ -f $FD/$pre_nameSI${DataFc}$tail_name'_V.nc' ] & [ -f $FD/$pre_nameSI${DataFc}$tail_name'_T.nc' ]
then
pre_name=$pre_nameSI
fi
if [ $currents = 70 ] & [ -f $FD/$pre_nameAN${DataFc}$tail_name'_U.nc' ] & [ -f $FD/$pre_nameAN${DataFc}$tail_name'_V.nc' ] & [ -f $FD/$pre_nameAN${DataFc}$tail_name'_T.nc' ]
then
pre_name=$pre_nameAN
fi

if [ -f $FD/$pre_name${DataFc}$tail_name'_U.nc' ]
then
   ln -s $FD/$pre_name${DataFc}$tail_name'_U.nc' $FD/${DataFc}'_U.nc'
   echo "${DataFc}_U.nc 1" >> tmp1.tmp
else
   echo "${DataFc}_U.nc 0" >> tmp1.tmp
   if [ $currents = 70 ]
   then
   echo "For this run you need the file:" $pre_nameFC${DataFc}$tail_name'_U.nc' "or" $pre_nameAN${DataFc}$tail_name'_U.nc' "or" $pre_nameSI${DataFc}$tail_name'_U.nc'
   else
   echo "For this run you need the file" $pre_name${DataFc}$tail_name'_U.nc'
   fi
fi
if [ -f $FD/$pre_name${DataFc}$tail_name'_V.nc' ]
then
   ln -s $FD/$pre_name${DataFc}$tail_name'_V.nc' $FD/${DataFc}'_V.nc'
   echo "${DataFc}_V.nc 1" >> tmp1.tmp
   
else
   echo "${DataFc}_V.nc 0" >> tmp1.tmp
   if [ $currents = 70 ]
   then
   echo "For this run you need the file:" $pre_nameFC${DataFc}$tail_name'_V.nc' "or" $pre_nameAN${DataFc}$tail_name'_V.nc' "or" $pre_nameSI${DataFc}$tail_name'_V.nc' 
   else
   echo "For this run you need the file" $pre_name${DataFc}$tail_name'_V.nc'
   fi
fi

if [ -f $FD/$pre_name${DataFc}$tail_name'_T.nc' ]
then
   ln -s $FD/$pre_name${DataFc}$tail_name'_T.nc' $FD/${DataFc}'_T.nc'
   echo "${DataFc}_T.nc 1" >> tmp1.tmp
   
else
   echo "${DataFc}_T.nc 0" >> tmp1.tmp
   if [ $currents = 70 ]
   then
   echo "For this run you need the file:" $pre_nameFC${DataFc}$tail_name'_T.nc' "or" $pre_nameAN${DataFc}$tail_name'_T.nc' "or" $pre_nameSI${DataFc}$tail_name'_T.nc' 
   else
   echo "For this run you need the file" $pre_name${DataFc}$tail_name'_T.nc'
   fi
fi

if [ $wind != 01 ]
then

if [ $wind = 06 ]
then
dir_wind='ECM'
fi

if [ $wind = 25 ]
then
dir_wind='E25'
fi

FD_wind_in=$F_DATA/$fcst_data/$dir_wind
FD_wind=$F_DATA/$fcst_data/$dir_wind

if [ $currents = 70 ] || [ $currents = 72 ] || [ $currents = 10 ] 
then
if [ $hour -lt 13 ]
then
Datafc=`./jday $Datafc +1`
fi
fi
DataFC=`echo $Datafc |cut -c3-8`
DataFc=$DataFC


Data_wind="20"${DataFc}
File_wind="ecm_"${DataFc}".ecm"
echo 'CHECK WINDS' $Data_wind
if [ -f $FD_wind_in/${Data_wind}.nc ]
then
   echo "${File_wind} 1" >> tmp2.tmp
   echo $FD_wind_in/${Data_wind}.nc ": File exists" 
else
   echo "${File_wind} 0" >> tmp2.tmp
   echo "For this run you need the file" ${Data_wind}.nc

fi


fi


done


cat tmp1.tmp>>medslik5.inp
echo $numfiles >> medslik5.inp

cat tmp2.tmp>>medslik5.inp
cat medslik_multi.tmp>>medslik5.inp 



#############################################################
#6. AREA SELECTION (MIN/MAX Longitudes & Latitudes)
#############################################################
cp source/medslikYYYY.tmp medslik0.tmp
sed -e "s/regione/$region/"\
    -e "s/correnti/$currents/"\
    -e "s/vento/$wind/"\
     medslik0.tmp>medslik1.tmp
rm medslik0.tmp


./lat_lon.exe
mv medslik1.tmp medslik.tmp

echo $numfiles >> medslik.tmp

numfiles_tmp=`expr $numfiles + 1`
n=0
while [ $n != $numfiles_tmp ]; do
n=`expr $n + 1`
nn=`expr $n - 1`
Datafc=`./jday $FcStart1 +$nn`
DataFC=`echo $Datafc |cut -c3-8`
DataFc=$DataFC
echo ${DataFc}24 >> medslik.tmp
done

if [ $currents = 70 ] || [ $currents = 72 ] || [ $currents = 10 ] 
then

if [ $hour -ge 13 ] 
then
n=0
numfiles_tmp=`expr $numfiles + 1`
echo $numfiles_tmp >> medslik.tmp 
while [ $n != $numfiles_tmp ]; do
n=`expr $n + 1`
nn=`expr $n - 1`
Datafc=`./jday $FcStart1 +$nn` 
echo ${Datafc} >> medslik.tmp
done

else 
n=0
numfiles_tmp=`expr $numfiles`
echo $numfiles_tmp >> medslik.tmp 
while [ $n != $numfiles_tmp ]; do
n=`expr $n + 1`
Datafc=`./jday $FcStart1 +$n` 
echo ${Datafc} >> medslik.tmp
done
fi
fi


echo " 0" >> medslik.tmp 


tr -d "\015" < medslik5.inp > medslik52.inp
cp medslik52.inp medslik5.inp
rm medslik52.inp

tr -d "\015" < medslik.tmp > medslik2.tmp
cp medslik2.tmp medslik.tmp
rm medslik2.tmp
#####################################################################
#7. EXTRACT CURRENTS, WIND AND SST DATA
#####################################################################
echo 'READING CURRENTS & WIND DATA'
./Extract_II.exe $F_DATA
######################################################################
#8. RUN 
#####################################################################
echo 'PLEASE WAIT: SIMULATION IS RUNNING'
./medslik_II.exe

#####################################################################
#9. ARCHIVE OUTPUT FILES
#####################################################################
rm tmp*.tmp
DIR_output=$model'_20'$year'_'$month'_'$day'_'$hour$minutes'_'$SIM_NAME

mkdir output/$DIR_output
cp medslik5.inp  output/$DIR_output
cp medslik5.par  output/$DIR_output
cp medslik.tmp output/$DIR_output
cp initial.txt output/$DIR_output
cp medslik_inputfile.txt output/$DIR_output
mv output/*.srf output/$DIR_output
mv output/*.cst output/$DIR_output
mv output/*.dsp output/$DIR_output
mv output/*.tot output/$DIR_output
mv output/*.fte output/$DIR_output

mv obs* output/$DIR_output

if  [ $wind = 25 ]
then
mkdir output/$DIR_output/E25
mv fcst_data/E25/*.ecm output/$DIR_output/E25
fi

if  [ $wind = 06 ]
then
mkdir output/$DIR_output/ECM
mv fcst_data/ECM/*.ecm output/$DIR_output/ECM
fi


if  [ $currents = 10 ]
then
mkdir output/$DIR_output/OPA
mv fcst_data/OPA/*.opa output/$DIR_output/OPA
fi
if  [ $currents = 70 ]
then
mkdir output/$DIR_output/O1h
mv fcst_data/O1h/*.opa output/$DIR_output/O1h
fi
if  [ $currents = 71 ]
then
mkdir output/$DIR_output/S1h
mv fcst_data/S1h/*.sic output/$DIR_output/S1h
fi
if  [ $currents = 72 ]
then
mkdir output/$DIR_output/A1h
mv fcst_data/A1h/*.adr output/$DIR_output/A1h
fi
if  [ $currents = 73 ]
then
mkdir output/$DIR_output/T1h
mv fcst_data/T1h/*.tyr output/$DIR_output/T1h
fi
if  [ $currents = 74 ]
then
mkdir output/$DIR_output/H3k
mv fcst_data/H3k/*.rel output/$DIR_output/H3k
fi
if  [ $currents = 75 ]
then
mkdir output/$DIR_output/WME
mv fcst_data/WME/*.wme output/$DIR_output/WME
fi


