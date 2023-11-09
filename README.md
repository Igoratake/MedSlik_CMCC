# [MEDSLIK-II](http://www.medslik-ii.org/index.html) v1.01 (10/2012)

MEDSLIK-II software can be run on GNU/Linux Operative Systems.
> You can download MEDSLIK-II source code from the [website](http://www.medslik-ii.org/users/login.php), after registering.

## GET STARTED in 3 STEPS!
1. You need to install conda or miniconda with python3.x, following the instructions on conda [website](https://docs.conda.io/projects/miniconda/en/latest/). Then, create a new conda environment with ncl library.
```
conda create -n mdk1.01 -c conda-forge ncl
conda activate mdk1.01
```
2.  Move the tarball to your home directory. Extract the tarball contents and enter the main folder.
```
mv MEDSLIK_II_1.01.tar.gz $HOME
cd $HOME
tar –zxvf MEDSLIK_II_1.01.tar.gz MEDSLIK_II_1.01
cd MEDSLIK_II_1.01
```

3. You can now run the Algeria test case for Medslik II (v1.01)
```
chmod +x run_testcase.sh
sudo ./run_testcase.sh
```

For more details about running MEDSLIK-II 1.01 and installing software requirements, please see the [documentation folder](https://github.com/Igoratake/Medslik-II/tree/medslik_II_1_01/doc/).

# Enjoy!
