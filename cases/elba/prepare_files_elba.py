'''
This script is used to slice and transform the original files found on MEDSLIK-II website (http://www.medslik-ii.org/index.html)

The intent is to reduce the size of the files, so one can run the siulation without having a large content of ocean and meteorological data.

After this transformation, original pre processing scripts will not work on these slice data. Please use the ones found on the container

This script works exclusively with Serious Game - Elba Case for version 1.02, which can be found on the following link if you have been subscribed

http://www.medslik-ii.org/data/cases/Algeria_test_case.tar.gz
'''

#Necessary Libraries
import numpy as np
import xarray as xr
from glob import glob

#Update these paths according to where the original files have been downloaded/saved
path_origin = '/Users/iatake/Downloads/Simulation results/V1.02/original/Serious_Game_test_case/'
#Update these paths according to where the cropped files will be saved
path_destin = '/Users/iatake/Dropbox (CMCC)/Work/MEDSLIK-II and Pyslick/Medslik-II/cases/elba/'

#Characteristics of Algeria Oil spill to facilitate the data slicing
lon_min, lon_max = 9.5,10.5
lat_min, lat_max = 42.5,43.5

#Dictionary containing names with the possibility to rename
variables_to_rename = {'depthu': 'depth', 'depthv': 'depth', 'deptht': 'depth',
                        'nav_lat':'lat','nav_lon':'lon', 'y':'lat','x':'lon',
                        'votemper':'thetao',
                        'vozocrtx':'uo',
                        'vomecrty':'vo',
                        'time_counter':'time'}

### SEA STATE FILES ###

for day in ['17','18','19','20','21']:

    #Selecting all the netcdf files from the chosen day
    dts = glob(path_origin + f'OCE/201405{day}*.nc')

    #opening netcdf files available on the same day
    ds = xr.open_mfdataset(dts)

    # Rename variables only if they exist in the dataset
    for old_name, new_name in variables_to_rename.items():
        if old_name in ds.variables:
            ds = ds.rename({old_name: new_name})

    #Trimming into smaller regional box
    ds = ds.sel(lon = slice(lon_min,lon_max),
                lat = slice(lat_min,lat_max))

    #Selecting only 4 layers
    ds = ds.sel(depth=[0,10,30,120],method='nearest')
    #Modifying labels to simplfy drop in temperature columns
    ds['depth'] = [0,10,30,120]

    #saves the daily current or temperature netcdf in the case dir
    ds.to_netcdf(path_destin + f'oce_files/Sea_Med_slice_201405{day}_Elba.nc')

### MET STATE FILES ###
    
#Selecting all the netcdf files from the chosen day
for day in ['17','18','19','20','21']:

    met = xr.open_mfdataset(path_origin + f'MET/201405{day}*.nc')

    # Rename variables only if they exist in the dataset
    for old_name, new_name in variables_to_rename.items():
        if old_name in met.variables:
            met = met.rename({old_name: new_name})

    met = met[['U10M','V10M']]

    #Trimming into smaller regional box
    met = met.sel(lon = slice(lon_min,lon_max),
                  lat = slice(lat_max,lat_min))
    
    #saving the sampled netcdf in cases directory
    met.to_netcdf(path_destin + f'met_files/Wind_Med_slice_{day}_Elba.nc')