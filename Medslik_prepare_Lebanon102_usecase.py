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

#paralell code
import multiprocessing

# Functions outside this script
from functions.medslik_utils import *
from scripts import *

simdir         = 'cases/'
simname        = 'lebanon102'
sim_date       = '13/07/2006' ### Simulation start day  - format DD/MM/YYYY (string)
sim_hour       = '08:00'      ### Simulation start hour - format HH:mm (string)
longitude       = 35.16667        ### Longitude of Simulation spill location - format Decimal degrees (float)
latitude        = 33.38333          ### Latitude of Simulation spill  - format Decimal degrees (float)
sim_lenght      = 660              ### Length of the simulation - format hours (int)
spill_duration  = 144              ### Duration of the spill - format hours (int)
oil_api         = 20               ### Oil API - format (float)
oil_volume      = 130.35*144           ### Volume of oil in tons - format (float) 
use_satellite   = False            ### Usage of Satellite imagery to model multiple slicks - True/False
use_slk_contour = False            ### Usage of slicks countours - True/False
separate_slicks = False         ### If the sim needs different slicks to have different properties, different sims have to be run
number_slick    = 1                ### Number of slicks to be simulated - format (int)

if separate_slicks:
    #Using this will create different config files according to the length of the information list
    #Algeria
    # s_volume = [240,225,215]
    # s_rate   = [0,0,0]
    #Elba
    s_volume = [34,1,6,8,3,6,5]
    s_rate   = [0,0,0,0,0,0,0]

# Obtaining spill rate from oil volume and spill duration
if spill_duration != 0:
    spill_rate = oil_volume/spill_duration
else:
    spill_rate = oil_volume

def write_config_files(separate_slicks=False,s_volume=None,s_rate=None,s_num=None):

     # # modify config_1.txt
    print('...config1.txt...')

    # Iterating through slicks or doing for single simulation
    if separate_slicks == False:
        config_file = f'cases/{simname}/xp_files/config1.txt'
    
    else:
        config_file = f'cases/{simname}/xp_files/slick{s_num+1}/config1.txt'

    subprocess.run([f'cp scripts/templates/config1_template_0.txt {config_file}'],shell=True)

    #adding spill Name - Add slick number if separate slicks
    if separate_slicks == False:
        search_and_replace(config_file, 'RUNNAME', simname)
    else:
        search_and_replace(config_file, 'RUNNAME', simname+f'_slick{s_num+1}')

    #adding spill date and hour information
    search_and_replace(config_file, 'DD', f'{dt_sim.day:02d}')
    search_and_replace(config_file, 'MM', f'{dt_sim.month:02d}')
    search_and_replace(config_file, 'YY', f'{dt_sim.year-2000:02d}')
    search_and_replace(config_file, 'c_Hour', f'{dt_sim.hour:02d}')
    search_and_replace(config_file, 'c_minute', f'{dt_sim.minute:02d}')

    #adding simulation length
    search_and_replace(config_file, 'SIMLENGTH', f'{sim_lenght:04d}') 

    #  adding spill coordinates - dd for degrees and mm for minutes
    # Latitude
    dd = int(latitude)
    mm = (float(latitude)-dd)*60
    search_and_replace(config_file, 'LATd', f'{dd:02d}')
    search_and_replace(config_file, 'LATm', f'{mm:.3f}')          
    
    # Longitude
    dd = int(longitude)
    mm = (float(longitude)-dd)*60
    search_and_replace(config_file, 'LONd', f'{dd:02d}')
    search_and_replace(config_file, 'LONm', f'{mm:.3f}')

    # spill duration
    if separate_slicks == False:
        search_and_replace(config_file, 'SDUR', f'{spill_duration:04d}')
    else:
        search_and_replace(config_file, 'SDUR', f'{s_rate:04d}')

    # spill volume
    if separate_slicks == False:
        search_and_replace(config_file, 'SRATE', f'{spill_rate:08.2f}')
    else:
        search_and_replace(config_file, 'SRATE', f'{s_volume:08.2f}')

    # oil characteristics
    search_and_replace(config_file, 'APIOT', f'{oil_api}') 

    #number of slicks
    search_and_replace(config_file, 'N_SLICK', f'{number_slick}')

    #slick countour
    if use_slk_contour == True:
        slik = 'YES'

        if separate_slicks == False:
            with open(f'cases/{simname}/xp_files/slick_countour.txt', 'r') as file1:
                content = file1.read()
            with open(config_file, 'a') as file2:
            # Append the contents of the first file to config file
                file2.write(content)
        else:
            with open(f'cases/{simname}/xp_files/slick{s_num+1}/slick_countour.txt', 'r') as file1:
                content = file1.read()
            with open(config_file, 'a') as file2:
            # Append the contents of the first file to config file
                file2.write(content)
    else:
        slik = 'NO'  

    #Writing that will use slick countor
    search_and_replace(config_file, 'SSLICK', f'{slik}')


if __name__ == '__main__':
    __spec__ = None

     # Simulation dates
    dt_sim = validate_date(sim_date + ' ' + sim_hour)

    if isinstance(dt_sim,str):
        raise ValueError('Wrong date format.')
    
    #Initial date is always one day prior due to Medslik-II interpolation
    dtini = dt_sim - datetime.timedelta(days=1)
         
    #End date, by safety margin is two days after the sim start + sim length
    dtend = dt_sim + datetime.timedelta(days = (sim_lenght/24) + 2) 

    #List of dates between initial and end date
    date_list = pd.date_range(dtini, dtend, freq='D')
    
    ##### SEA / OCEAN ##### 

    #opening all files in the directory and concatenating them automatically through open_mfdataset
    concat = xr.open_mfdataset(f'cases/{simname}/oce_files/*.nc',combine='nested')
    concat = concat.drop_duplicates(dim="time", keep="last")

    #Interpolating the values in time, transforming it from daily to hourly values
    concat = concat.resample(time="1H").interpolate("linear")

    #iterating at each hour to generate the .mrc files
    parallel_processing_mrc(concat,simname=simname)          
        
    ##### WINDS ##### 

    concat = xr.open_mfdataset(f'cases/{simname}/met_files/*.nc',combine='nested')
    concat = concat.drop_duplicates(dim="time", keep="first")
    concat = concat.resample(time="1H").interpolate("linear")

    #iterating at each hour to generate the .eri files
    for date in date_list:            
        
        #Call write eri function located in medslik.utils file
        write_eri(concat,date,simname=simname)
    
    print('Met State variables written')

    ##### BATHYMETRY AND COASTLINE ##### 

    # Process for Bathymetry and Coastline Files
    grid_string = glob(f'{simdir}{simname}/oce_files/*.nc')[0] 

    # Bathymetry for gebco 2023
    subprocess.run([f'{sys.executable}', 'scripts/pre_processing/preproc_gebco_mdk2.py', 
                    'data/gebco/GEBCO_2023.nc',
                    grid_string, f'{simdir}{simname}/bnc_files/'])

    # # gshhs in intermediate resolution
    subprocess.run([f'{sys.executable}', 'scripts/pre_processing/preproc_gshhs_mdk2.py', 
                    'data/gshhs/GSHHS_shp/f/GSHHS_f_L1.shp',
                    grid_string, f'{simdir}{simname}/bnc_files/'])
    
    ##### CONFIG FILES ##### 

    # prepare medslik_II.for and config1.txt
    print('Preparing configuration files... ')

     # get dimensions from ncfiles
    #currents
    my_o = xr.open_dataset(glob(f'cases/{simname}/oce_files/*nc')[0])

    variables_to_rename = {'depthu': 'depth', 'depthv': 'depth', 'deptht': 'depth',
                           'votemper':'thetao',
                           'vozocrtx':'uo',
                           'vomecrty':'vo',
                           'time_counter':'time'}
    
    # Rename variables only if they exist in the dataset
    for old_name, new_name in variables_to_rename.items():
        if old_name in my_o.variables:
            my_o = my_o.rename({old_name: new_name})

    for var in ['uo','vo','thetao']:
        try:
            my_o = my_o.isel(time=0,depth=0)[var].values.shape
            found_variable = True
        except:
            continue
        if not found_variable:
            raise ValueError("Check the air state variables available in your dataset")

    #wind
    for var in ['U10M','u10']:
        try:
            my_w = xr.open_dataset(glob(f'cases/{simname}/met_files/*nc')[0]).isel(time=0)[var].values.shape
            found_variable = True
        except:
            continue
        if not found_variable:
            raise ValueError("Check the air state variables available in your dataset")

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

    if separate_slicks:
        for i, (vol, dur) in enumerate(zip(s_volume,s_rate)):
            write_config_files(separate_slicks=True,s_volume=vol,s_rate=dur,s_num=i)

    else:
        write_config_files()
