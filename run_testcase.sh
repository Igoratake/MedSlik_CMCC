#!/bin/bash

chmod +x requirements.sh
sudo ./requirements.sh

TESTCASENAME=$1
echo "Your testcase is ${TESTCASENAME}".
mkdir test_case
cd test_case
wget "http://www.medslik-ii.org/data/cases/${TESTCASENAME}_test_case.tar.gz"
tar -zxvf "${TESTCASENAME}_test_case.tar.gz"
cd "${TESTCASENAME}_test_case"
cp OCE/* ../../DATA/OCE/
cp MET/* ../../DATA/MET/
gunzip ../../DATA/*/*
cp SET-UP/medslik_inputfile.txt ../../RUN/medslik_inputfile.txt
cp SET-UP/medslik_inputfile_euler.txt ../../RUN/medslik_inputfile.txt
cp SET-UP/medslik5.par ../../RUN/
cd ../../SRC
./compile.sh
cd ../RUN
./medslik_II.sh
cd ../OUTPUT