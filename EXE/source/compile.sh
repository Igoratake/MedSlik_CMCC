DIR_EXE=/scratch/surf/MEDSLIK_II_v2/EXE
DIR_SRC=/scratch/surf/MEDSLIK_II_v2/EXE/source

gfortran  -I/usr/local/netcdf/current/include -L/usr/local/netcdf/current/lib  $DIR_SRC/Extract_II.for -lnetcdf -lnetcdff -o $DIR_EXE/Extract_II.exe
gfortran -o $DIR_EXE/jday $DIR_SRC/jday.f 
gfortran  -I/usr/local/netcdf/current/include -L/usr/local/netcdf/current/lib  $DIR_SRC/medslik_II.for -lnetcdf -lnetcdff -o $DIR_EXE/medslik_II.exe
gfortran -o $DIR_EXE/lat_lon.exe $DIR_SRC/lat_lon.for
