# [MEDSLIK-II](http://www.medslik-ii.org/index.html) v2.01 (06/2020)

MEDSLIK-II software can be run on GNU/Linux Operative Systems.
> You can download MEDSLIK-II source code from the [website](http://www.medslik-ii.org/users/login.php), after registering.

## GET STARTED in 3 STEPS!
1. You need to install conda or miniconda with python3.x, following the instructions on conda [website](https://docs.conda.io/projects/miniconda/en/latest/). Then, create a new conda environment and install pip requirements.
```
conda create --name mdk2.01
conda activate mdk2.01
conda install pip
```
2. Extract the tarball contents and enter the main software folder. Install pip requirements.
```
tar –zxvf MEDSLIK_II_2.01.zip MEDSLIK_II_2.01
cd MEDSLIK_II_2.01
pip install -r pip_requirements.txt
```
3. You can now run the test case for Medslik II (v2.01). You need to substitute the argument "name" with one of the followig words: lebanon, paria.
```
chmod +x run_testcase.sh
sudo ./run_testcase.sh name
```

For more details about running MEDSLIK-II 1.02 and installing software requirements, please see the [documentation folder](https://github.com/Igoratake/Medslik-II/tree/medslik_II_2_01/doc/).

# Enjoy!
