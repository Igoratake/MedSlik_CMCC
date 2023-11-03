# Medslik-II version 1.02 07/2018

This version currently has all the files, pre and post processing scripts as in its original releases.

Please bear in minf that some of the scripts might be outdated and need older versions in order to be able to run.

Modifying files is adviced to use this version properly.

Down below, the instructions to use it are displayed as it was intended to be.

# Quick-start guide

The currently supported architectures is Linux.
Software Requirements:
1. FORTRAN 90/95 compiler (gfortran is fully compatible):
to Install gfortran
compiler:
sudo apt-get install gfortran
2. NetCDF library.
sudo apt-get install netcdf-bin
sudo apt-get install libnetcdf-dev
Installation
1. Download the model and put the tarball in your home directory.
2. Extract the contents of the tarball:
tar -xvzf MEDSLIK_II_v2.0.tar.gz MEDSLIK_II_v2.0
3. Open /home/user/MEDSLIK_II_v2.0/EXE/source/compile.sh and give the correct path for your netCDF installation to HOME_MEDSLIK.

4. Now go to /home/user/MEDSLIK_II_v2.0/EXE and compile the code:
sh source/compile.sh
5. The model should be compiled now. Two last steps:
open /home/user/MEDSLIK_II_v2.0/EXE/medslik_II_ens.sh and inform the correct MEDSLIK installation folder (/home/user/MEDSLIK_II_v2.0) to the HOME_MEDSLIK variable
repeat the same procedure for /home/user/MEDSLIK_II_v2.0/EXE/RUN.sh
6. You should be now ready to run the test case. Start by extracting the file paria_casestudy.tar.gz

7. Copy the input files to your MEDSLIK-II installation
copy the contents from the oce_files folder into
/home/user/MEDSLIK_II_v2/DATA/fcst/H3k
copy the contents from the met_files folder into
/home/user/MEDSLIK_II_v2/DATA/fcst/SK1
copy the contents from the bnc_files folder into
/home/user/MEDSLIK_II_v2/EXE/data
copy the contents from the xp_files folder into
/home/user/MEDSLIK_II_v2/EXE
8. You are ready to launch your simulation:
go to /home/user/MEDSLIK_II_v2/EXE/
launch the model: ./RUN.sh
9. The MEDSLIK-II outputs can be found in /home/user/MEDSLIK_II_v2/EXE/outputs