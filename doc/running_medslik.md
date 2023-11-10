# RUNNING MEDSLIK-II v2.00

If not done yet, you can install the needed packages following the instructions in the [explanatory file](https://github.com/Igoratake/Medslik-II/tree/medslik_II_2_00/doc/installing_requirements.md).

1. Activate your conda environment and enter the MEDSLIK_II folder.
```
conda activate mdk2.00
cd MEDSLIK_II_v2
```
2. Create a folder for test cases.
```
mkdir test_case
cd test_case
```
3. Download test case and extract the tarball.
```
wget "http://www.medslik-ii.org/data/v2.0/cases/paria_casestudy.tar.gz"
tar -zxvf "paria_casestudy.tar.gz"
```
4. Copy the input files within the DATA folder and unzip them. Create data folder if it does not exist.
```
mkdir -p ../DATA/fcst_data/H3k
cp Users/stage/Downloads/paria_casestudy/oce_files/* ../DATA/fcst_data/H3k/
mkdir -p ../DATA/fcst_data/SK1
cp Users/stage/Downloads/paria_casestudy/met_files/* ../DATA/fcst_data/SK1/
cp Users/stage/Downloads/paria_casestudy/bnc_files/* ../EXE/data/
cp Users/stage/Downloads/paria_casestudy/xp_files/* ../EXE/
```
5. Compile and run.
```
cd ../EXE
sh source/compile.sh
./RUN.sh
cd ../output
```
8. You may found the output files within the specific folder.
```
cd output
```

For details about running the test cases and more, please see medslik [manual](https://github.com/Igoratake/Medslik-II/blob/medslik_II_2_00/doc/MEDSLIK_II_v2.0_user_manual.pdf).
