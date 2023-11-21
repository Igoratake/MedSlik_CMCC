# set folders
DIR_EXE=-I"/Users/iatake/Dropbox (CMCC)/Work/MEDSLIK-II and Pyslick/Medslik-II/MEDSLIK_II_3.01/RUN"
DIR_SRC=$DIR_EXE/MODEL_SRC

gfortran -o $DIR_EXE/jday $DIR_SRC/jday.f
gfortran  -I/usr/include -L/usr/lib -fno-align-commons $DIR_SRC/medslik_II.for -lnetcdf -lnetcdff -o $DIR_EXE/medslik_II.exe
gfortran -o $DIR_EXE/lat_lon.exe $DIR_SRC/lat_lon.for
