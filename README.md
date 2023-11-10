# [MEDSLIK-II](http://www.medslik-ii.org/index.html) v1.02 (09/2015)

MEDSLIK-II software can be run on GNU/Linux Operative Systems.
> You can download MEDSLIK-II source code from the [website](http://www.medslik-ii.org/users/login.php), after registering.

## GET STARTED in 3 STEPS!
1. You need to install conda or miniconda with python3.x, following the instructions on conda [website](https://docs.conda.io/projects/miniconda/en/latest/). Then, create a new conda environment with PyNGL and PyNIO library.
```
conda create --name mdk1.02 --channel conda-forge pynio pyngl
conda activate mdk1.02
```
2. Extract the tarball contents and enter the main software folder.
```
tar –zxvf MEDSLIK_II_1.02.tar.gz MEDSLIK_II_1.02
cd MEDSLIK_II_1.02
```
3. You can now run the test case for Medslik II (v1.02). You need to replace "name" with one of the following words: Algeria, Lebanon, Serious_Game (please, notice the capital letters and the underscore).
```
chmod +x run_testcase.sh
sudo ./run_testcase.sh name
```

For more details about running MEDSLIK-II 1.02 and installing software requirements, please see the [documentation folder](https://github.com/Igoratake/Medslik-II/tree/medslik_II_1_02/doc/).

# Enjoy!