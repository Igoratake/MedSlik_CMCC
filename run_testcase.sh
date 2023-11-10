#!/bin/bash

chmod +x requirements.sh
sudo ./requirements.sh

TESTCASENAME=$1
echo "Your testcase is ${TESTCASENAME}".
mkdir test_case
cd test_case

case $TESTCASENAME in

  paria)
    LINK="https://www.dropbox.com/sh/50ak5xarrc3rqns/AABRjgilPbaeBSHGgdiDp10wa/paria_casestudy.zip?raw=1"
    REGION=vnzl
    ;;

  lebanon)
    LINK="https://www.dropbox.com/sh/50ak5xarrc3rqns/AADjdU_5PisJFJbiAV_GvlIMa/lebanon_casestudy.zip?raw=1"
    REGION=lebn
    ;;

  *)
    echo -n "unknown"
    ;;
esac

wget LINK
tar -zxvf "${TESTCASENAME}_casestudy.zip"
cp ${TESTCASENAME}_casestudy/conf.json ../
cp ${TESTCASENAME}_casestudy/METOCE_INP/ORIGINAL/OCE/ ../METOCE_INP/ORIGINAL/OCE/
cp ${TESTCASENAME}_casestudy/METOCE_INP/ORIGINAL/MET/ ../METOCE_INP/ORIGINAL/MET/
cp ${TESTCASENAME}_casestudy/config1.txt ../RUN/config1.txt
cp ${TESTCASENAME}_casestudy/config2.txt ../RUN/config2.txt
cd ..
sed "s/INSTALLATION_FOLDER/PWD%\/MEDSLIK_II_2.01*/g" conf.json > temp_file && mv temp_file conf.json
python3 ./METOCE_INP/PREPROC_SCRIPTS/preproc_currents_mdk2.py ./conf.json
python3 ./METOCE_INP/PREPROC_SCRIPTS/preproc_winds_mdk2.py ./conf.json
python3 ./DTM_INP/PREPROC_SCRIPTS/preproc_gebco_mdk2.py ./conf.json
python3 ./DTM_INP/PREPROC_SCRIPTS/preproc_gshhs_mdk2.py ./conf.json
cd RUN
sed "82s/.*/${REGION}/" medslik_II.sh > temp_file && mv temp_file medslik_II.sh
./MODEL_SRC/compile.sh
./RUN.sh
cd ..
python3 ./PLOT/oil_track_mdk2.py ./conf.json
python3 ./PLOT/oil_beached_mdk2.py ./conf.json
cd OUT