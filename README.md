## Medslik-II version 1.01 10/2012

This version currently has all the files, pre and post processing scripts as in its original releases.

Please bear in minf that some of the scripts might be outdated and need older versions in order to be able to run.

Modifying files is adviced to use this version properly.

Down below, the instructions to use it are displayed as it was intended to be.

# Quick-start guide
The currently supported architectures is Linux (tested on Ubuntu 10.04 LTS and 11.04 LTS).
Software Requirements:
1. FORTRAN 90/95 compiler (gfortran is fully compatible):
to Install gfortran
compiler:
sudo apt-get install gfortran
2. NetCDF library.
sudo apt-get install netcdf-bin
sudo apt-get install libnetcdf-dev
3. NCL (provided within the model for the visualization).

Installation
1. Download the model and put the tarball in your home directory.
The file downloaded will have a different name according to the version.

2. Uncompress and extract the contents of the tarball.
tar –xvzf MEDSLIK_II_1.01.tar.gz MEDSLIK_II_1.01

This operation will create and populate the directories MEDSLIK_II_1.01/EXE containing the source
files and script files and executables for running and visualizing the test case.
The MEDSLIK_II_1.01 system is composed of six main parts, i.e.:

source code (directory source);
input data files (directory data);
output data files (directory output);
script files to compile and execute in a Linux operative system;
visualization software (directory medslik_plots);
Test Case set-up (directory test_cases).
3. Enter MEDSLIK_II_1.01 and compile the code. For example:
cd $HOME/MEDSLIK_II_1.01/EXE
sh source/compile.sh

4. At this stage, you need the sample current and wind files.
Pick up the sample current file Currents
Pick up the sample wind file Wind

5. Unzip and place currents in $HOME/MEDSLIK_II_1.01/DATA/fcst_data/O1h and winds in $HOME/MEDSLIK_II_1.01/DATA/fcst_data/ECM.

6. Copy the Algeria test case input file
cd $HOME/MEDSLIK_II_1.01/EXE
cp test_cases/TEST_ALGERIA/medslik_inputfile.txt .
cp test_cases/TEST_ALGERIA/observation_0808071050.txt .

7. Now you are ready to run the code. Just type
./RUN.sh

The simulation will start and after few minutes in the directory
$HOME/MEDSLIK_II_1.01/EXE/output/MFS_2008_08_06_0951_TEST_ALGERIA/plots
you will find the pictures of the initial position of the slick and its predicted position every 6 hours.
You can also visualize the pictures by opening the file index.html (in the diretory plots) with a web browser.