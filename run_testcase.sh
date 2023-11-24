#!/bin/bash

# Install libraries and programs.
# chmod +x requirements.sh
# sudo ./requirements.sh
# # Download forecast data: currents (O1h) and winds (ECM).
# cd DATA/fcst_data
# wget http://www.medslik-ii.org/data/O1h.tar.gz
# tar -zxvf O1h.tar.gz O1h
# rm O1h.tar.gz
# wget http://www.medslik-ii.org/data/ECM.tar.gz
# tar -zxvf ECM.tar.gz ECM
# rm ECM.tar.gz
# cd ../..
# Copy Algeria test case setup.
cp EXE/test_cases/TEST_ALGERIA/medslik_inputfile.txt ./EXE
cp EXE/test_cases/TEST_ALGERIA/observation_0808071050.txt ./EXE
# Compile and run.
sh EXE/source/compile.sh
./EXE/RUN.sh
# Go to the oputput folder and list files.
cd output/MFS_2008_08_06_0951_TEST_ALGERIA_V1.01
ls -l