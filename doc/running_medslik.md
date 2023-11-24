# RUNNING MEDSLIK-II v1.01

If not done yet, you can install the needed packages following the instructions in the [explanatory file](https://github.com/Igoratake/Medslik-II/tree/medslik_II_1_01/doc/installing_requirements.md).

1. Activate your conda environment and enter the MEDSLIK_II folder.
```
conda activate mdk1.01
cd MEDSLIK-II_1.01
```
2. Download the _sample currents_ within the forecast data folder.
```
cd DATA/fcst_data
wget http://www.medslik-ii.org/data/O1h.tar.gz
tar -zxvf O1h.tar.gz O1h
rm O1h.tar.gz
```
3. Download the _sample wind_ within the forecast data folder.
```
wget http://www.medslik-ii.org/data/ECM.tar.gz
tar -zxvf ECM.tar.gz ECM
rm ECM.tar.gz
cd ../..
```
4. Copy the _Algeria_ test case input file in the execution folder.
```
cp EXE/test_cases/TEST_ALGERIA/medslik_inputfile.txt ./EXE
cp EXE/test_cases/TEST_ALGERIA/observation_0808071050.txt ./EXE
```
5. Compile the source code and execute.
```
sh EXE/source/compile.sh
./EXE/RUN.sh
```
6. You can now find the output files in the relative output folder and list them.
```
cd output/MFS_2008_08_06_0951_TEST_ALGERIA_V1.01
ls -l 
```