#!/bin/bash

chmod +x requirements.sh
sudo ./requirements.sh

mkdir test_case
cd test_case
wget "http://www.medslik-ii.org/data/v2.0/cases/paria_casestudy.tar.gz"
tar -zxvf "paria_casestudy.tar.gz"
mkdir -p ../DATA/fcst_data/H3k
cp Users/stage/Downloads/paria_casestudy/oce_files/* ../DATA/fcst_data/H3k/
mkdir -p ../DATA/fcst_data/SK1
cp Users/stage/Downloads/paria_casestudy/met_files/* ../DATA/fcst_data/SK1/
cp Users/stage/Downloads/paria_casestudy/bnc_files/* ../EXE/data/
cp Users/stage/Downloads/paria_casestudy/xp_files/* ../EXE/
cd ../EXE
sh source/compile.sh
./RUN.sh
cd output