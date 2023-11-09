#!/bin/bash

# Install libraries and programs
sudo apt-get install gcc -y
sudo apt-get install gfortran -y
sudo apt-get install netcdf-bin -y
sudo apt-get install libnetcdf-dev -y
sudo apt-get install libnetcdff-dev -y
# Create conda env with ncl
conda create -n mdk1.01 -c conda-forge ncl
conda activate mdk1.01
# Download forecast data: currents (O1h) and winds (ECM)
cd DATA/fcst_data
wget http://www.medslik-ii.org/data/O1h.tar.gz
tar -zxvf O1h.tar.gz O1h
rm O1h.tar.gz
wget http://www.medslik-ii.org/data/ECM.tar.gz
tar -zxvf ECM.tar.gz ECM
rm ECM.tar.gz
cd ../..
# Copy Algeria test case setup
cd EXE
cp test_cases/TEST_ALGERIA/medslik_inputfile.txt .
cp test_cases/TEST_ALGERIA/observation_0808071050.txt .
# Compile and run
sh source/compile.sh
./RUN.sh
# Go to the oputput folder and list files
cd output/MFS_2008_08_06_0951_TEST_ALGERIA_V1.01
ls -l