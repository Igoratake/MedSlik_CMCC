# load  modules
module load  gcc-12.2.0/12.2.0   intel-2021.6.0/2021.6.0  intel-2021.6.0/libszip/2.1.1-tvhyi  intel-2021.6.0/hdf5/1.13.3-xwdun intel-2021.6.0/netcdf-c/4.9.0-cjqig intel-2021.6.0/impi-2021.6.0/netcdf-cxx-threadsafe/4.2-p2ziq intel-2021.6.0/impi-2021.6.0/netcdf-fortran-threadsafe/4.6.0-75oow   intel-2021.6.0/netcdf-c-threadsafe/4.9.0-25h5k intel-2021.6.0/udunits/2.2.28-5obkm intel-2021.6.0/cdo-threadsafe/2.1.1-lyjsw intel-2021.6.0/magics/4.9.3-jrpbm intel-2021.6.0/ncview/2.1.8-sds5t intel-2021.6.0/nco/5.0.6-jp6y4 intel-2021.6.0/eccodes/2.25.0-jwx44 intel-2021.6.0/netcdf-fortran/4.6.0-5vpiq

# read config file
# source ../mdk2.conf

# set folders
DIR_EXE=.
DIR_SRC=$DIR_EXE/MODEL_SRC
NETCDF=/juno/opt/spacks/0.20.0/opt/spack/linux-rhel8-icelake/intel-2021.6.0/netcdf-fortran

# compile jday
echo "gfortran -o $DIR_EXE/jday $DIR_SRC/jday.f"
gfortran -o $DIR_EXE/jday $DIR_SRC/jday.f

# compile medslik
echo "gfortran -I$NETCDF/4.6.0-5vpiqkqbf5m6prlbpl4gie5rdixzybks/include   -L$NETCDF/4.6.0-5vpiqkqbf5m6prlbpl4gie5rdixzybks/lib   $DIR_SRC/medslik_II.for -lnetcdf -lnetcdff -o $DIR_EXE/medslik_II.exe"

gfortran -I$NETCDF/4.6.0-5vpiqkqbf5m6prlbpl4gie5rdixzybks/include   -L$NETCDF/4.6.0-5vpiqkqbf5m6prlbpl4gie5rdixzybks/lib   $DIR_SRC/medslik_II.for -lnetcdf -lnetcdff -fallow-argument-mismatch -o $DIR_EXE/medslik_II.exe

# compile lat_lon
echo "gfortran -o $DIR_EXE/lat_lon.exe $DIR_SRC/lat_lon.for"
gfortran -o $DIR_EXE/lat_lon.exe $DIR_SRC/lat_lon.for
