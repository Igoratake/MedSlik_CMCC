DIR_EXE=$PWD
DIR_SRC=$PWD/MODEL_SRC

NETCDF_DIR=$HOME/.local/netcdf/netcdf-c-4.9.2
gfortran  -I$NETCDF_DIR/lib  $DIR_SRC/Extract_II.for -lnetcdf -lnetcdff -o $DIR_EXE/Extract_II.exe
gfortran -o $DIR_EXE/jday $DIR_SRC/jday.f
gfortran  -I$NETCDF_DIR/lib -fno-align-commons $DIR_SRC/medslik_II.for -lnetcdf -lnetcdff -o $DIR_EXE/medslik_II.exe
gfortran -o $DIR_EXE/lat_lon.exe $DIR_SRC/lat_lon.for
