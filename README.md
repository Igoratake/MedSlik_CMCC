# Medslik-II version 2.01 06/2020

This version currently has all the files, pre and post processing scripts as in its original releases.

Please bear in minf that some of the scripts might be outdated and need older versions in order to be able to run.

Modifying files is adviced to use this version properly.

Down below, the instructions to use it are displayed as it was intended to be.

# Quick-start guide

Software Requirements
The list of software requirements is given below with a link to the source web page (the command lines reported refer if you are using a Debian GNU/Linux as Operating System)

1. The fortran compiler gfortran
sudo apt-get update
sudo apt-get install gcc gfortran
NOTE: gfortran may be pre-installed in your linux systems. Check the version using 'gfortran --version'
2. The NetCDF libraries (>v4.2)
The easiest way is to install the precompiled NetCDF library (NetCDF C and Fortran compiler)
sudo apt-get install libnetcdf-dev
sudo apt-get install netcdf-bin
sudo apt-get install libnetcdf-dev
sudo apt-get install libnetcdff-dev
NOTE: This approach can cause compatibility problem with your Fortran compiler. The problem is that for a library to be compatible with your Fortran compiler it has to be compiled with the same compiler and with the same flags you will be using and that may be not the case when you get a precompiled library from a repository.
A safer way is to install NetCDF from the source code. The Fortran netCDF library need to be built after the NetCDF-C library is built and installed.
Download the NetCDF source code from unidata website
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.4.1.1.tar.gz
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.4.4.tar.gz
Download also the supporting libraries (HDF5, zlib, szip)
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/hdf5-1.8.13.tar.gz
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/zlib-1.2.8.tar.gz
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/szip-2.1.tar.gz
Start to install the szip library. Untar the file:
tar -zxvf szip-2.1.tar.gz
Configure and build the szip library in the directory /usr/local/szip/szip-2.1
cd szip-2.1
./configure --prefix=/usr/local/szip/szip-2.1
make
make check
make install
Then install the zlib library. Untar the file:
tar -zxvf zlib-1.2.8.tar.gz
Configure and build the zlib library in the directory /usr/local/zlib/zlib-1.2.8
cd zlib-1.2.8
./configure --prefix=/usr/local/zlib/ zlib-1.2.8
make
make test
make install
Then install the hdf5 library. Untar the file:
tar -zxvf hdf5-1.8.13.tar.gz
Configure and build the hdf5 library in the directory /usr/local/hdf5/hdf5-1.8.13
cd hdf5-1.8.13
./configure --enable-fortran --with-zlib=/usr/local/zlib/zlib-1.2.8
--with-szip=/usr/local/szip/szip-2.1 --prefix=/usr/local/hdf5/hdf5-1.8.13
make
make test
make install
Now you can install the netcdf-library. First Netcdf-C library. Untar the file:
tar -zxvf netcdf-4.4.1.1.tar.gz
tar -zxvf netcdf-fortran-4.4.4.tar.gz
Configure and build the NetCDF-C library in the directory /usr/local/netcdf/netcdf-4.4.1.1
cd netcdf-4.4.1.1
export LDFLAGS='-L/usr/local/zlib/zlib-1.2.8/lib -L/usr/local/szip/szip-2.1/lib
-L/usr/local/hdf5/hdf5-1.8.13/lib'
export CPPFLAGS='-I/usr/local/zlib/zlib-1.2.8/include
-I/usr/local/szip/szip-2.1/include -I/usr/local/hdf5/hdf5-1.8.13/include'
./configure --prefix=/usr/local/netcdf/netcdf-4.4.1.1 --disable-dap
--disable-dap-remote-tests --enable-netcdf4 --enable-shared
make
make check
make install
You need to add the path of the netcdf library to etc/ld.so.conf and update shared library cache. In this way our library can be correctly found
echo "/usr/local/netcdf/netcdf-4.4.1.1/lib" > /etc/ld.so.conf.d/netcdf-4.4.1.1.conf
sudo ldconfig
Now you can install the netcdf-fortran library. Untar the file:
tar -zxvf netcdf-fortran-4.4.4.tar.gz
Configure and install the NetCDF-fortran libraries
cd ../netcdf-fortran-4.4.4
export LDFLAGS='-L/usr/local/netcdf/netcdf-4.4.1.1/lib'
export CPPFLAGS='-I/usr/local/netcdf/netcdf-4.4.1.1/include'
./configure --prefix=/usr/local/netcdf/netcdf-4.4.1.1 --enable-shared
make
make check
make install
Finally add the /bin directory to the PATH environment variable. Insert the following command to the ~/.bashrc
export PATH=$PATH:/usr/local/netcdf/netcdf-4.4.1.1/bin
Check the version using ‘nc-config --version’

3. The Climate Data Operator (CDO) software

sudo apt-get install cdo
4. The Python 3.x. The easiest and common way to install Python is to install Miniconda distribution software

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
Follow the prompts on the installer screens. Miniconda contains also the command-line tool called conda which is a package manager you can use to install additional python packages. Test your installation by run the command conda list. A list of installed packages appears if it has been installed correctly.

5. The additional python modules netCDF4, matplotlib, scipy, basemap, basemap-data-hires.

conda install -c conda-forge netcdf4
conda install -c conda-forge matplotlib
conda install -c conda-forge scipy
conda install -c conda-forge basemap
conda install -c conda-forge basemap-data-hires

Installation of model code
1. After downloading the MEDSLIK-II v2.01 code, move and extract the zip compressed archive in the chosen installation path (${INSTALLATION_FOLDER}/).
2. The source code is located inside the folder ${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/RUN/MODEL_SRC/. Before executing compile.sh, open this file with any editor and give the correct path to the ${INSTALLATION_FOLDER} variable.
3. Once you are sure that the file has execute permission, just run the compilation process with the following command:

./compile.sh
4. Two last steps:
Open /${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/RUN/medslik_II.sh with any editor and give the correct path to the ${INSTALLATION_FOLDER} variable.
Repeat the same procedure for /${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/RUN/RUN.sh.
Running the Paria test case
1. After downloading the paria_casestudy.zip, extract the compressed archive in the path (${PARIA_PATH}/).

2. Copy the configuration file conf.json located in the directory ${PARIA_PATH}/ paria_casestudy/) in the MEDSLIK-II 2.01 root directory /${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/

cp ${PARIA_PATH}/paria_casestudy/conf.json ${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/
Before using conf.json, open this file with any editor and give the correct path to the ${INSTALLATION_FOLDER} variable.

3. Copy the ocean current input files located in the directory
${PARIA_PATH}/paria_casestudy/METOCE_INP/ORIGINAL/OCE/
in the directory
${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/METOCE_INP/ORIGINAL/OCE/
and execute the preprocessing for the currents using the script preproc_currents_mdk2.py from command line:

python preproc_currents_mdk2.py conf.json
4. Execute the preprocessing for the bathymetry using the script preproc_gebco_mdk2.py from command line:

python preproc_gebco_mdk2.py conf.json
5. Execute the preprocessing for the coastline using the script preproc_gshhs_mdk2.py from command line:

python preproc_gshhs_mdk2.py conf.json
6. Copy the atmospheric input files located in the directory
${PARIA_PATH}/paria_casestudy/METOCE_INP/ORIGINAL/MET/
in the directory
${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/METOCE_INP/ORIGINAL/MET/
and execute the preprocessing for the atmospheric data using the script preproc_winds_mdk2.py from command line:

python preproc_winds_mdk2.py conf.json
7. Copy config1.txt and config2.txt located in the directory
${PARIA_PATH}/paria_casestudy/RUN/
in the directory
${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/RUN/
8. Execute the RUN.sh bash script that is located in
${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/RUN/
Check the medslik log file located in this directory.
9. The output of the paria test case can be found in
${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/OUT/
10. Visualize the output using the script oil_track_mdk2.py and oil_beached_mdk2.py located in
${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/PLOT/.
python oil_track_mdk2.py conf.json
python oil_beached_mdk2.py conf.json
Running the Lebanon test case
1. After downloading the lebanon_casestudy.zip, extract the compressed archive in the path (${LEBANON_PATH}/).

2. Copy the configuration file conf.json located in the directory ${LEBANON_PATH}/ lebanon_casestudy/) in the MEDSLIK-II 2.01 root directory /${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/

cp ${LEBANON_PATH}/lebanon_casestudy/conf.json ${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/
Before using conf.json, open this file with any editor and give the correct path to the ${INSTALLATION_FOLDER} variable.

3. Copy the ocean current input files located in the directory
${LEBANON_PATH}/lebanon_casestudy/METOCE_INP/ORIGINAL/OCE/
in the directory
${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/METOCE_INP/ORIGINAL/OCE/
and execute the preprocessing for the currents using the script preproc_currents_mdk2.py from command line:

python preproc_currents_mdk2.py conf.json
4. Execute the preprocessing for the bathymetry using the script preproc_gebco_mdk2.py from command line:

python preproc_gebco_mdk2.py conf.json
5. Execute the preprocessing for the coastline using the script preproc_gshhs_mdk2.py from command line:

python preproc_gshhs_mdk2.py conf.json
6. Copy the atmospheric input files located in the directory
${LEBANON_PATH}/lebanon_casestudy/METOCE_INP/ORIGINAL/MET/
in the directory
${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/METOCE_INP/ORIGINAL/MET/
and execute the preprocessing for the atmospheric data using the script preproc_winds_mdk2.py from command line:

python preproc_winds_mdk2.py conf.json
7. Copy config1.txt and config2.txt located in the directory
${LEBANON_PATH}/lebanon_casestudy/RUN/
in the directory
${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/RUN/
8. Execute the RUN.sh bash script that is located in
${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/RUN/
Check the medslik log file located in this directory.
9. The output of the Lebanon test case can be found in
${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/OUT/
10. Visualize the output using the script oil_track_mdk2.py and oil_beached_mdk2.py located in
${INSTALLATION_FOLDER}/MEDSLIK_II_2.01/PLOT/.
python oil_track_mdk2.py conf.json
python oil_beached_mdk2.py conf.json