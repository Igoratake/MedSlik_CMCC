#geo libs
import xarray as xr
import geopandas as gpd
import cartopy.crs as ccrs
from shapely.geometry import Point

#numerical and plotting libs
import numpy as np
import matplotlib.pyplot as plt

#system libs
import os
import sys
import time
import datetime
import subprocess
import threading
from glob import glob

# Functions outside this script
from functions.medslik_utils import *
from scripts import *

simdir         = 'cases/'
simname        = 'paria'
sim_date       = '24/04/2017' ### Simulation start day  - format DD/MM/YYYY (string)
sim_hour       = '13:00'      ### Simulation start hour - format HH:mm (string)

 # Simulation dates
dt_sim = validate_date(sim_date + ' ' + sim_hour)

if isinstance(dt_sim,str):
    raise ValueError('Wrong date format.')

if __name__ == '__main__':


    #sending curr asc files files to TEMP dir to execute medlisk
    subprocess.run([f'cp {simdir}{simname}/oce_files/*.mrc MEDSLIK_II_3.01/RUN/TEMP/OCE/'],shell=True)
    
    #sending wind asc files files to TEMP dir to execute medlisk
    subprocess.run([f'cp {simdir}{simname}/met_files/*.eri MEDSLIK_II_3.01/RUN/TEMP/MET/'],shell=True)

    # Process for Bathymetry and Coastline Files
    grid_string = glob(f'{simdir}{simname}/oce_files/*.nc')[0] 

    # Bathymetry for gebco 2023
    subprocess.run([f'{sys.executable}', 'scripts/pre_processing/preproc_gebco_mdk2.py', 
                    'data/gebco/GEBCO_2023.nc',
                    grid_string, f'{simdir}{simname}/bnc_files/'])

    # # gshhs in intermediate resolution
    subprocess.run([f'{sys.executable}', 'scripts/pre_processing/preproc_gshhs_mdk2.py', 
                    'data/gshhs/GSHHS_shp/i/GSHHS_i_L1.shp',
                    grid_string, f'{simdir}{simname}/bnc_files/'])

    # prepare medslik_II.for and config1.txt
    print('Preparing configuration files... ')

    # get dimensions from ncfiles
    my_o = xr.open_dataset(glob(f'cases/{simname}/oce_files/*nc')[0]).isel(deptht=0,time_counter=0).votemper.values.shape

    # my_o = xr.open_dataset(glob(f'cases/{simname}/oce_files/*nc')[0]).isel(depth=0,time=0).votemper.values.shape
    my_w = xr.open_dataset(glob(f'cases/{simname}/met_files/*nc')[0]).isel(time=0).U10M.values.shape
    # my_w = my_o

    nmax = np.max([np.max(my_o),np.max(my_w)])
    imx_o = my_o[1]
    jmx_o = my_o[0]
    imx_w = my_w[1]
    jmx_w = my_w[0]
    
    # modify medslik_ii
    print('...medslik_ii.for...')

    med_for = f'cases/{simname}/xp_files/medslik_II.for'

    subprocess.run([f'cp scripts/templates/medslik_II_template.for {med_for}'],shell=True)

    # Replacing NMAX in medslik fortran with a python function
    search_and_replace(med_for, 'NMAX', str(nmax))

    # copy METOCEAN files to MEDSLIK-II installation
    subprocess.run([f'cp {simdir}{simname}/oce_files/* MEDSLIK_II_3.01/METOCE_INP/PREPROC/OCE/'],shell=True)
    subprocess.run([f'cp {simdir}{simname}/oce_files/* MEDSLIK_II_3.01/RUN/TEMP/OCE/'],shell=True)

    subprocess.run([f'cp {simdir}{simname}/met_files/* MEDSLIK_II_3.01/METOCE_INP/PREPROC/MET/'],shell=True)
    subprocess.run([f'cp {simdir}{simname}/met_files/* MEDSLIK_II_3.01/RUN/TEMP/MET/'],shell=True)
    # copy bnc files
    subprocess.run([f'cp {simdir}{simname}/bnc_files/* MEDSLIK_II_3.01/DTM_INP/'],shell=True)
    # copy Extract and config files
    subprocess.run([f'cp {simdir}{simname}/xp_files/medslik_II.for MEDSLIK_II_3.01/RUN/MODEL_SRC/'],shell=True)
    subprocess.run([f'cp {simdir}{simname}/xp_files/config1.txt MEDSLIK_II_3.01/RUN/'],shell=True)
    subprocess.run([f'cp {simdir}{simname}/xp_files/config2.txt MEDSLIK_II_3.01/RUN/'],shell=True)

    # Compile and start running
    subprocess.run([f'cd MEDSLIK_II_3.01/RUN/; sh MODEL_SRC/compile.sh; ./RUN.sh'],shell=True,check=True)

    subprocess.run([f'cp -r MEDSLIK_II_3.01/OUT/MDK_SIM_{dt_sim.year}_{str(dt_sim.month).zfill(2)}_{str(dt_sim.day).zfill(2)}*/ {simdir}{simname}/out_files/'],shell=True)
    subprocess.run([f'rm -rf {simdir}{simname}/out_files/MET {simdir}{simname}/out_files/OCE'],shell=True)