# INSTALLING REQUIREMENTS for MEDSLIK-II v1.02

1. Install fortran compiler from apt-get.
```
sudo apt-get install gcc
sudo apt-get install gfortran
```
2. Download NetCDF libraries. You can do that using apt-get.
```
sudo apt-get install netcdf-bin
sudo apt-get install libnetcdf-dev
sudo apt-get install libnetcdff-dev
```
If you have issues installing netcdf libraries this way, you can install them from source, following the next steps.

3. Create a .local directory within your home.
```
mkdir -p $HOME/.local
cd $HOME/.local
```
4. Create folders for libraries.
```
mkdir szip
mkdir zlib
mkdir hdf5
mkdir netcdf
```
5. Download required libraries.
```
wget -P ./szip https://support.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz
wget -P ./zlib https://zlib.net/current/zlib-1.3.tar.gz
wget -P ./hdf5 https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.11/src/hdf5-1.10.11.tar.gz
wget -P ./netcdf https://downloads.unidata.ucar.edu/netcdf-c/4.9.2/netcdf-c-4.9.2.tar.gz
wget -P ./netcdf https://downloads.unidata.ucar.edu/netcdf-fortran/4.6.1/netcdf-fortran-4.6.1.tar.gz
```
6. Install szip.
```
cd szip
tar -zxvf szip-2.1.1.tar.gz
cd szip-2.1.1
./configure --prefix=$HOME/.local/szip/szip-2.1.1
make
make check
make install
cd ..
```
7. Instal zlib.
```
cd zlib
tar -zxvf zlib-1.3.tar.gz
cd zlib-1.3
./configure --prefix=$HONE/.local/zlib/zlib-1.3
make
make test
make install
cd ..
```
8. Install hdf5.
```
cd hdf5
tar -zxvf hdf5-1.10.11.tar.gz
cd hdf5-1.10.11
./configure --enable-fortran --with-zlib=$HOME/.local/zlib/zlib-1.3
--with-szip=$HOME/.local/szip/szip-2.1.1 --prefix=$HOME/.local/hdf5/hdf5-1.10.11
make
make test
make install
cd ..
```
9. Install netcdf-c.
```
cd netcdf
tar -zxvf netcdf-c-4.9.2.tar.gz
cd netcdf-c-4.9.2
export LDFLAGS='-L$HOME/.local/zlib/zlib-1.3/lib -L$HOME/.local/szip/szip-2.1.1/lib -L$HOME/.local/hdf5/hdf5-1.10.11/lib'
export CPPFLAGS='-I$HOME/.local/zlib/zlib-1.3/include -I$HOME/.local/szip/szip-2.1.1/include -I$HOME/.local/hdf5/hdf5-1.10.11/include'
./configure --prefix=$HOME/.local/netcdf/netcdf-c-4.9.2 --disable-dap --disable-dap-remote-tests --enable-netcdf4 --enable-shared
make
make check
make install
cd ..
```
10. You need to add the path of the netcdf library to etc/ld.so.conf and update shared library cache. In this way our library can be correctly found.
```
sudo echo "$HOME/.local/netcdf/netcdf-c-4.9.2/lib" > /etc/ld.so.conf.d/netcdf-c-4.9.2.conf
sudo ldconfig
```
11. Install netcdf-fortran.
```
cd netcdf
tar -zxvf netcdf-fortran-4.6.1.tar.gz
cd netcdf-fortran-4.6.1
export LDFLAGS='-L$HOME/.local/netcdf/netcdf-c-4.9.2/lib'
export CPPFLAGS='-I$HOME/.local/netcdf/netcdf-c-4.9.2/include'
./configure --prefix=$HOME/.local/netcdf/netcdf-c-4.9.2 --enable-shared
cmake .
make
make check
make install
```
12. You need to install conda or miniconda with python3.x, following the instructions on conda [website](https://docs.conda.io/projects/miniconda/en/latest/). Then, create a new conda environment with PyNGL and PyNIO libraries.
```
conda create --name mdk1.02 --channel conda-forge pynio pyngl
```
You can now run medslik following the steps in the [explanatory file](https://github.com/Igoratake/Medslik-II/tree/medslik_II_1_02/doc/running_medslik.md).
