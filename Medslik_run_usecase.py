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
simname        = 'baniyas_medslik_v3'
separate_slicks = False         ### If the sim needs different slicks to have different properties, different sims have to be run

def run_medslik_sim(separate_slicks = separate_slicks):

    if separate_slicks == False:
        # copy METOCEAN files to MEDSLIK-II installation
        subprocess.run([f'cp {simdir}{simname}/oce_files/*.mrc MEDSLIK_II_3.01/RUN/TEMP/OCE/'],shell=True)
        subprocess.run([f'cp {simdir}{simname}/met_files/*.eri MEDSLIK_II_3.01/RUN/TEMP/MET/'],shell=True)
        # copy bnc files
        subprocess.run([f'cp {simdir}{simname}/bnc_files/* MEDSLIK_II_3.01/DTM_INP/'],shell=True)
        # copy Extract and config files
        subprocess.run([f'cp {simdir}{simname}/xp_files/medslik_II.for MEDSLIK_II_3.01/RUN/MODEL_SRC/'],shell=True)
        subprocess.run([f'cp {simdir}{simname}/xp_files/config2.txt MEDSLIK_II_3.01/RUN/'],shell=True)
        subprocess.run([f'cp {simdir}{simname}/xp_files/config1.txt MEDSLIK_II_3.01/RUN/'],shell=True)
        # Compile and start running
        subprocess.run([f'cd MEDSLIK_II_3.01/RUN/; sh MODEL_SRC/compile.sh; ./RUN.sh'],shell=True,check=True)        
    
    else:
        slicks = glob(f'{simdir}{simname}/xp_files/*/')
        for i in range(0,len(slicks)):
            subprocess.run([f'cp {simdir}{simname}/oce_files/*.mrc MEDSLIK_II_3.01/RUN/TEMP/OCE/'],shell=True)
            subprocess.run([f'cp {simdir}{simname}/met_files/*.eri MEDSLIK_II_3.01/RUN/TEMP/MET/'],shell=True)
            # copy bnc files
            subprocess.run([f'cp {simdir}{simname}/bnc_files/* MEDSLIK_II_3.01/DTM_INP/'],shell=True)
            # copy Extract and config files
            subprocess.run([f'cp {simdir}{simname}/xp_files/medslik_II.for MEDSLIK_II_3.01/RUN/MODEL_SRC/'],shell=True)
            subprocess.run([f'cp {simdir}{simname}/xp_files/config2.txt MEDSLIK_II_3.01/RUN/'],shell=True)
            subprocess.run([f'cp {simdir}{simname}/xp_files/slick{i+1}/config1.txt MEDSLIK_II_3.01/RUN/'],shell=True)
            # Compile and start running
            subprocess.run([f'cd MEDSLIK_II_3.01/RUN/; sh MODEL_SRC/compile.sh; ./RUN.sh'],shell=True,check=True)
    
    #Send files to case dir and remove temp files
    subprocess.run([f'cp -r MEDSLIK_II_3.01/OUT/MDK_SIM_*{simname}_*/ {simdir}{simname}/out_files/'],shell=True)
    subprocess.run([f'rm -rf {simdir}{simname}/out_files/MET {simdir}{simname}/out_files/OCE'],shell=True)

if __name__ == '__main__':

    run_medslik_sim(separate_slicks=separate_slicks)