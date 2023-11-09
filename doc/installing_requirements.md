# INSTALLING REQUIREMENTS for MEDSLIK-II v1.01

1. You need to install conda or miniconda with python3.x, following the instructions on conda [website](https://docs.conda.io/projects/miniconda/en/latest/). Then, create a new conda environment with ncl library.
```
conda create -n mdk1.01 -c conda-forge ncl
```
2. Install fortran compiler.
```
sudo apt-get install gcc -y
sudo apt-get install gfortran -y
```
3. Install netcdf libraries
```
sudo apt-get install netcdf-bin -y
sudo apt-get install libnetcdf-dev -y
sudo apt-get install libnetcdff-dev -y
```

You can now run medslik following the steps in the [explanatory file](https://github.com/Igoratake/Medslik-II/tree/medslik_II_1_01/doc/running_medslik.md).