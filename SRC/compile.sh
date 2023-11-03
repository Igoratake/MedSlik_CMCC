#!/bin/bash

HOME_MEDSLIK=${PWD%/SRC*}
DIR_SRC=${HOME_MEDSLIK}/SRC
DIR_EXE=${HOME_MEDSLIK}/EXE

NETCDF_DIR=

gfortran  -I${NETCDF_DIR}include -L${NETCDF_DIR}lib  $DIR_SRC/Extract_II.F90 -lnetcdf -lnetcdff -o $DIR_EXE/Extract_II.exe

gfortran -o $DIR_EXE/jday $DIR_SRC/jday.f 

gfortran -g -o $DIR_EXE/medslik_II.exe $DIR_SRC/medslik_II.for 

gfortran -o $DIR_EXE/lat_lon.exe $DIR_SRC/lat_lon.F90

gfortran -g -fbounds-check -I${NETCDF_DIR}include -o $DIR_EXE/create_output.exe $DIR_SRC/module_phymath.F90 $DIR_SRC/module_interpolation.F90 $DIR_SRC/module_netcdf.F90 $DIR_SRC/create_netcdf.F90 -L${NETCDF_DIR}lib -lnetcdf -lnetcdff
