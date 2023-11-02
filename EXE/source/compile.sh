DIR_EXE=$HOME/MEDSLIK_II_1.01/EXE
DIR_SRC=$HOME/MEDSLIK_II_1.01/EXE/source

gfortran  -I/usr/include -L/usr/lib  $DIR_SRC/Extract_II.for -lnetcdf -lnetcdff -o $DIR_EXE/Extract_II.exe
gfortran -o $DIR_EXE/jday $DIR_SRC/jday.f 
gfortran -o $DIR_EXE/medslik_II.exe $DIR_SRC/medslik_II.for 
gfortran -o $DIR_EXE/lat_lon.exe $DIR_SRC/lat_lon.for
