# RUNNING MEDSLIK-II v1.02

If not done yet, you can install the needed packages following the instructions in the [explanatory file](https://github.com/Igoratake/Medslik-II/tree/medslik_II_1_02/doc/installing_requirements.md).

1. Activate your conda environment and enter the MEDSLIK_II folder.
```
conda activate mdk1.02
cd MEDSLIK-II_1.02
```
2. Choose the name of your test case. The argument "name" should be one of the following words: Algeria, Lebanon, Serious_Game (please, notice the capital letters and the underscore).
```
TESTCASE=name
```
3. Create a folder for test cases.
```
mkdir test_case
cd test_case
```
4. Download test case and extract the tarball.
```
wget "http://www.medslik-ii.org/data/cases/${TESTCASENAME}_test_case.tar.gz"
tar -zxvf "${TESTCASENAME}_test_case.tar.gz"
cd "${TESTCASENAME}_test_case"
```
5. Copy the input files within the DATA folder and unzip them.
```
cp OCE/* ../../DATA/OCE/
cp MET/* ../../DATA/MET/
gunzip ../../DATA/*/*
```
6. Copy the set-up files in the RUN folder.
```
cp SET-UP/medslik_inputfile.txt ../../RUN/medslik_inputfile.txt
cp SET-UP/medslik_inputfile_euler.txt ../../RUN/medslik_inputfile.txt
cp SET-UP/medslik5.par ../../RUN/
```
7. Compile and run.
```
cd ../../SRC
./compile.sh
cd ../RUN
./medslik_II.sh
```
8. You may found the output files within the specific folder.
```
cd ../OUTPUT
```

For details about running the test cases and more, please see medslik [manual](https://github.com/Igoratake/Medslik-II/blob/medslik_II_1_02/doc/MEDSLIK_II_v1.02_user_manual.pdf).
