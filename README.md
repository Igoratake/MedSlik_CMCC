# Medslik-II version 1.02 09/2015

This version currently has all the files, pre and post processing scripts as in its original releases.

Please bear in minf that some of the scripts might be outdated and need older versions in order to be able to run.

Modifying files is adviced to use this version properly.

Down below, the instructions to use it are displayed as it was intended to be.

# Quick-start guide

Medslik-II version 1.02 is written in FORTRAN-77/90, with in no machine- dependent elements, so that it can be installed without modifications on most platforms. The architecture currently supported is Linux (tested on Ubuntu 10.04 LTS, Lubuntu 14.04, Debian 7 and Centos 5).

The following external libraries are needed in order to run Medslik-II version 1.02:

1.MEDSLIK-II requires NetCDF 4.3 fortran libraries or later. Documentation on how to install properly NetCDF4 FORTRAN libraries can be found at http://www.unidata.ucar.edu/netcdf
2.Documentation on how to install properly PyNGL can be found at
https: //www.pyngl.ucar.edu/Download/
3.Documentation on how to install properly PyNIO can be found at https://www.pyngl.ucar.edu/Download/
A GNU Fortran compiler is needed.

## INSTALLATION

After downloading the version 1.02 code, it is needed to move and extract the gzip compressed tar archive in the chosen installation path:

mv /${DOWNLOAD_FOLDER}/MEDSLIK_II_1.02.tar.gz /${INSTALLATION_FOLDER}
cd /${INSTALLATION_FOLDER}
tar -xvzf MEDSLIK_II_1.02.tar.gz

The variables ${DOWNLOAD_FOLDER} and ${INSTALLATION_FOLDER} have to be replaced by the user specific path of the DOWNLOAD and user chosen INSTALLATION folders, respectively. The source code is located inside /${INSTALLATION_FOLDER}/MEDSLIK_II_1.02/SRC directory:

To compile properly the model FORTRAN code, first it is needed to open the compile.sh file with any editor and change the bash variable NETCDF_DIR adding the exact absolute path of the NetCDF libraries installation directory. Once you are sure that the file has execute permission, just run the build process with the following command:


./compile.sh

A detailed description of the model and all the steps needed to run a simulation is given in the
User manual