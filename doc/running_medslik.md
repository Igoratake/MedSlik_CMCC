# RUNNING MEDSLIK-II v2.01

If not done yet, you can install the needed packages following the instructions in the [explanatory file](https://github.com/Igoratake/Medslik-II/tree/medslik_II_1_02/doc/installing_requirements.md).

1. Activate your conda environment and enter the MEDSLIK_II folder.
```
conda activate mdk2.01
cd MEDSLIK_II_2.01
```
2. Choose the name of your test case. The argument "name" should be one of the following words: lebanon, paria.
```
TESTCASENAME=name
```
3. Create a folder for test cases.
```
mkdir test_case
cd test_case
```
4. Download one of the test case dataset (choose one of the following two commands).
```
wget https://www.dropbox.com/sh/50ak5xarrc3rqns/AADjdU_5PisJFJbiAV_GvlIMa/lebanon_casestudy.zip?raw=1 (LEBANON test case)
wget https://www.dropbox.com/sh/50ak5xarrc3rqns/AABRjgilPbaeBSHGgdiDp10wa/paria_casestudy.zip?raw=1 (PARIA test case)
```
5. Extract the folder, copy the data and config files.
```
tar -zxvf "${TESTCASENAME}_casestudy.zip"
cp ${TESTCASENAME}_casestudy/conf.json ../
cp ${TESTCASENAME}_casestudy/METOCE_INP/ORIGINAL/OCE/ ../METOCE_INP/ORIGINAL/OCE/
cp ${TESTCASENAME}_casestudy/METOCE_INP/ORIGINAL/MET/ ../METOCE_INP/ORIGINAL/MET/
cp ${TESTCASENAME}_casestudy/config1.txt ../RUN/config1.txt
cp ${TESTCASENAME}_casestudy/config2.txt ../RUN/config2.txt
cd ..
```
6. In file _conf.json_, you need to substitute ${INSTALLATION_FOLDER} with path to the MEDSLIK_II_2.01 folder. You can use the following command.
```
sed "s/INSTALLATION_FOLDER/PWD%\/MEDSLIK_II_2.01*/g" conf.json > temp_file && mv temp_file conf.json
```
7. Execute the preprocessing.
```
python3 ./METOCE_INP/PREPROC_SCRIPTS/preproc_currents_mdk2.py ./conf.json
python3 ./METOCE_INP/PREPROC_SCRIPTS/preproc_winds_mdk2.py ./conf.json
python3 ./DTM_INP/PREPROC_SCRIPTS/preproc_gebco_mdk2.py ./conf.json
python3 ./DTM_INP/PREPROC_SCRIPTS/preproc_gshhs_mdk2.py ./conf.json
```
8. Change name of region in RUN/medslik_II.sh file.
```
cd RUN
sed "82s/.*/${REGION}/" medslik_II.sh > temp_file && mv temp_file medslik_II.sh
```
9. Compile and run.
```
./MODEL_SRC/compile.sh
./RUN.sh
```
10. Run the postprocessing for visualizing results.
```
cd ..
python3 ./PLOT/oil_track_mdk2.py ./conf.json
python3 ./PLOT/oil_beached_mdk2.py ./conf.json
```
11. You may found the output results in the specific folder
```
cd OUT
ls -l 
```
For details about running the test cases and more, please see medslik [manual](https://github.com/Igoratake/Medslik-II/blob/medslik_II_2_01/doc/Manual_MEDSLIK_II_v2.01.pdf).
