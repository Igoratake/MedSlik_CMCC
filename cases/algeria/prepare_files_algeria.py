'''
This script is used to slice and transform the original files found on MEDSLIK-II website (http://www.medslik-ii.org/index.html)

The intent is to reduce the size of the files, so one can run the siulation without having a large content of ocean and meteorological data.

After this transformation, original pre processing scripts will not work on these slice data. Please use the ones found on the container

This script works exclusively with Lebanon Case for version 1.02, which can be found on the following link if you have been subscribed

http://www.medslik-ii.org/data/cases/Algeria_test_case.tar.gz
'''

#Necessary Libraries
import numpy as np
import xarray as xr
from glob import glob

#Update these paths according to where the original files have been downloaded/saved
path_origin = '/Users/iatake/Downloads/Simulation results/V1.02/original/Lebanon_test_case/'
#Update these paths according to where the cropped files will be saved
path_destin = '/Users/iatake/Dropbox (CMCC)/Work/MEDSLIK-II and Pyslick/Medslik-II/cases/lebanon/'


#Characteristics of Algeria Oil spill to facilitate the data slicing
lon_min, lon_max = 4,8
lat_min, lat_max = 37,39

#Dictionary containing names with the possibility to rename
variables_to_rename = {'depthu': 'depth', 'depthv': 'depth', 'deptht': 'depth',
                        'nav_lat':'lat','nav_lon':'lon', 'y':'lat','x':'lon',
                        'votemper':'thetao',
                        'vozocrtx':'uo',
                        'vomecrty':'vo',
                        'time_counter':'time'}

### SEA STATE FILES ###

'''
Unfortunately in Algeria test case, the current files have different lat and lon values along the days

Therefore it is difficult to concatenate current files over days. 

In this code, it will be enforced to use the lat and lon values from the first dataset

'''

ds_base = xr.open_dataset(path_origin + f'OCE/MEDffE080806_01_grid_V.nc')
latt = ds_base.nav_lat.values
lonn = ds_base.nav_lon.values

for day in ['05','06','07','08']:

    #Selecting all the netcdf files from the chosen day
    d1 = glob(path_origin + f'OCE/MEDffE0808{day}_*.nc')

    #opening each file separately, since they have different variable names
    ds1 = xr.open_dataset(d1[0])
    ds2 = xr.open_dataset(d1[1])
    ds3 = xr.open_dataset(d1[2])

    # Rename variables only if they exist in the dataset
    for old_name, new_name in variables_to_rename.items():
        if old_name in ds1.variables:
            ds1 = ds1.rename({old_name: new_name})
        if old_name in ds2.variables:
            ds2 = ds2.rename({old_name: new_name})
        if old_name in ds3.variables:
            ds3 = ds3.rename({old_name: new_name})

    #Mergin the three datasets to crop on just one file
    ds = xr.merge([ds1,ds2,ds3],compat='override')

    ds.lat.values = latt
    ds.lon.values = lonn

    #In this case the dataset have coordinates as variables. 
    #Obtaining unique values from lat and lon variables. 
    unique_lat = np.unique(ds.lat)
    unique_lon = np.unique(ds.lon)

    #Storing these dimensions as the unique values
    ds['x'] = np.sort(unique_lon)
    ds['y'] = np.sort(unique_lat)

    #dropping meshgrid variables
    ds = ds.drop(['lon','lat'])

    #renaming x and y dimensions
    ds = ds.rename({'x':'lon','y':'lat'})

    #Trimming into smaller regional box
    ds = ds.sel(lon = slice(lon_min,lon_max),
                lat = slice(lat_min,lat_max))

    #Selecting only 4 layers
    ds = ds.sel(depth=[0,10,30,120],method='nearest')
    #Modifying labels to simplfy drop in temperature columns
    ds['depth'] = [0,10,30,120]

    #saves the daily current or temperature netcdf in the case dir
    ds.to_netcdf(path_destin + f'oce_files/Curr_Med_slice_200808{day}_Algeria.nc')

### MET STATE FILES ###
    
#Selecting all the netcdf files from the chosen day
for day in ['06','07','08']:

    met = xr.open_dataset(path_origin + f'MET/200808{day}_6hi-INGV-ECMWF-AM050-MED.nc')

    # Rename variables only if they exist in the dataset
    for old_name, new_name in variables_to_rename.items():
        if old_name in met.variables:
            met = met.rename({old_name: new_name})

    met = met[['U10M','V10M']]

    #Trimming into smaller regional box
    met = met.sel(lon = slice(lon_min,lon_max),
                  lat = slice(lat_max,lat_min))
    
    #saving the sampled netcdf in cases directory
    met.to_netcdf(path_destin + f'met_files/Wind_Med_slice_200808{day}_Algeria.nc')