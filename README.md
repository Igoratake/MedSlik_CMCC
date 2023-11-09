# [MEDSLIK-II](http://www.medslik-ii.org/index.html) v1.01 (10/2012)

MEDSLIK-II software can be run on GNU/Linux Operative Systems.
> You can download MEDSLIK-II source code from the [website](http://www.medslik-ii.org/users/login.php), after registering.

## GET STARTED!
1. If not done yet, you need to install conda or miniconda, following the instructions on conda [website](https://docs.conda.io/projects/miniconda/en/latest/), with python3.x.

2.  Move the tarball to your home directory. Extract the tarball contents and enter the main folder.
```
mv MEDSLIK_II_1.01.tar.gz $HOME
cd $HOME
tar –zxvf MEDSLIK_II_1.01.tar.gz MEDSLIK_II_1.01
cd MEDSLIK_II_1.01
```
3. You can install the required libraries and programs through the requirements bash script.
```
chmod +x requirements.sh
sudo ./requirements.sh
```
4. Install ncl within a new conda environment and activate it.
```
conda create -n mdk1.01 -c conda-forge ncl
conda activate mdk1.01
```
5. Download the _sample currents_ within the forecast data folder.
```
cd DATA/fcst_data
wget http://www.medslik-ii.org/data/O1h.tar.gz
tar -zxvf O1h.tar.gz O1h
rm O1h.tar.gz
```
6. Download the _sample wind_ within the forecast data folder.
```
wget http://www.medslik-ii.org/data/ECM.tar.gz
tar -zxvf ECM.tar.gz ECM
rm ECM.tar.gz
cd ../..
```
7. Copy the _Algeria_ test case input file in the execution folder.
```
cd EXE
cp test_cases/TEST_ALGERIA/medslik_inputfile.txt .
cp test_cases/TEST_ALGERIA/observation_0808071050.txt .
```
8. Compile the source code and execute.
```
sh source/compile.sh
./RUN.sh
```
9. You can now find the output files in the relative output folder and list them.
```
cd output/MFS_2008_08_06_0951_TEST_ALGERIA_V1.01
ls -l 
```

# Enjoy!
