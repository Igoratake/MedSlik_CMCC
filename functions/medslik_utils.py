#geo libs
import xarray as xr
import geopandas as gpd
import cartopy.crs as ccrs
from shapely.geometry import Point
# from folium.plugins import HeatMap
# from streamlit_folium import st_folium,folium_static

#numerical and plotting libs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import plotly.figure_factory as ff

#system libs
import os
import sys
import time
import datetime
import subprocess
import threading
from glob import glob
import multiprocessing

'''

This script contains several functions that can be used across MEDSLIK-II modeling.

'''

def check_land(lon,lat):

    '''
    This script receives a lon and lat value and  checks if the position is within land or sea
    It uses a shapefile of world boundaries currently based in geopandas database

    In case a position is on land:
        Script returns 0, since it will not be possible to simulate oil spill on land
    
    Otherwise returns 1 

    '''

    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

    point = Point(lon, lat)

    is_within_land = world.geometry.contains(point).any()
    
    if is_within_land:
        sea=0
    else:
        sea=1

    return sea

def validate_date(date):

    '''
    Convert a date in string format to datetime, checking if the date provided is valid

    It also checks if date is in the future, blocking the user to procced

    '''

    try: 
        dt=datetime.datetime.strptime(date,'%d/%m/%Y %H:%M')
        if dt > datetime.datetime.today():
            dt = 'Date provided is in the future. No data will be available'               
        
    except:
        dt = 'Wrong date format'

    return dt

def search_and_replace(file_path, search_word, replace_word):
    with open(file_path, 'r') as file:
        file_contents = file.read()

        updated_contents = file_contents.replace(search_word, replace_word)

    with open(file_path, 'w') as file:
        file.write(updated_contents)

def write_cds(key):

    with open('~.cdsapirc_test') as f:
        f.write('url: https://cds.climate.copernicus.eu/api/v2\n')
        f.write(f'key: {key}\n')
        f.write('verify: 0')
        f.close()

def write_mrc(ds,simname=None):

    '''
    This function write the currents and sea surface temperature files for Medslik-II.

    Data has to be passed in hourly format and in netcdf format

    '''

#iterating at each hour to generate the .mrc files
    for i in range(0,len(ds.time)):            
        #select the i time instant
        rec = ds.isel(time=i)
        #get the datetime values from that instant
        try:
            dt = pd.to_datetime(rec.time.values)
        except:
            try:
                dt = datetime.datetime.strptime(str(rec.time.values),'%Y-%m-%d %H:%M:%S')
            except:
                raise ValueError('Datetime from the dataset in on an unknown format')
        
        #transforms it from xarray datset to pandas dataframe to facilitate the processes of adjusting values
        df = rec.to_dataframe().reset_index()
        df = df.fillna(0)
        df = df.drop(['time'],axis=1)
        #pivoting it in order to create the same pattern in .mrc files
        df = df.pivot(index = ['lat','lon'],columns='depth',values = ['thetao','uo','vo']).reset_index()
        #join colum information to 
        df.columns = [pair[0]+str(pair[1]).split('.')[0] for pair in df.columns]
        #dropping temperature columns
        df = df.drop(['thetao10','thetao30','thetao120'],axis=1)
        #sort first by latitude and then by longitude
        df = df.sort_values(['lon','lat'])
        df.columns = ['lat','lon','SST','u_srf','u_10m','u_30m','u_120m','v_srf','v_10m','v_30m','v_120m']
        df = df[['lat','lon','SST','u_srf','v_srf','u_10m','v_10m','u_30m','v_30m','u_120m','v_120m']]

        #making sure that .mrc files in 0 hour are written as 24
        #this code also makes sure that the file is written correctly even in the transition of months
        if dt.hour == 0:
            hour = 24
            day = (dt - datetime.timedelta(hours=1)).day
        else:
            hour = dt.hour
            day = dt.day

        #writing the current files
        with open(f'cases/{simname}/oce_files/merc{dt.year-2000:02d}{dt.month:02d}{day:02d}{hour:02d}.mrc', 'w') as f:
            f.write(f"Ocean forecast data for {day:02d}/{dt.month:02d}/{dt.year} {hour:02d}:00\n")
            f.write("Subregion of the Global Ocean:\n")
            f.write(f"{df.lon.min():02.2f}  {df.lon.max():02.2f}  {df.lat.min():02.2f} {df.lat.max():02.2f}   {len(rec.lon)}   {len(rec.lat)}   Geog. limits\n")
            f.write(f"{len(df)}   0.0\n")
            f.write("lat        lon        SST        u_srf      v_srf      u_10m      v_10m       u_30m      v_30m      u_120m     v_120m\n")
            
            for index, row in df.iterrows():
                f.write(f"{row['lat']:<10.4f}    {row['lon']:<10.4f}    {row['SST']:<10.4f}     {row['u_srf']:<10.4f}    {row['v_srf']:<10.4f}     {row['u_10m']:<10.4f}    {row['v_10m']:<10.4f}     {row['u_30m']:<10.4f}    {row['v_30m']:<10.4f}     {row['u_120m']:<10.4f}    {row['v_120m']:<10.4f}\n")

    print('Sea State variables written')

def write_eri(ds,date,simname=None):

    '''
    This function write the wind velocity as 10 m files for Medslik-II.

    Data has to be passed in hourly format and in netcdf format

    '''

    #iterating at each hour to generate the .eri files     
    try:
        date1 = f'{date.year}-{date.month:02d}-{date.day:02d} 00:00'
        date2 = f'{date.year}-{date.month:02d}-{date.day:02d} 23:00'
        met = ds.sel(time = slice(date1,date2))
        #getting date from the netcdf
        try:
            dt = pd.to_datetime(met.time[0].values)
        except:
            dt = datetime.datetime.strptime(str(met.time[0].values),'%Y-%m-%d %H:%M:%S')
        df = met.to_dataframe().reset_index()
        df = df.fillna(9999)
        df = df.pivot(index=['lat','lon'],columns='time',values=['U10M','V10M']).reset_index()
        df.columns = [pair[0]+str(pair[1]) for pair in df.columns]
        
        df = df.rename({'lonNaT':'lon','latNaT':'lat'},axis=1)
        df = df.sort_values(['lon','lat'], ascending=[True,False])
        
        #writing the wind files
        with open(f'cases/{simname}/met_files/erai{dt.year-2000:02d}{dt.month:02d}{dt.day:02d}.eri', 'w') as file:
            file.write(f" 10m winds forecast data for {dt.day:02d}/{dt.month:02d}/{dt.year}\n")
            file.write(" Subregion of the Global Ocean with limits:\n")
            file.write(f"  {df.lon.min():02.5f}  {df.lon.max():02.5f}  {df.lat.max():02.5f}  {df.lat.min():02.5f}   {len(met.lon)}   {len(met.lat)}   Geog. limits\n")
            file.write(f"   {len(df)}   0.0\n")
            file.write("                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  lat        lon        u00        v00        u01       v01      u02      v02     u03     v03     u04     v04     u05     v05     u06     v06\n")
            
            for index, row in df.iterrows():
                file.write(f"   {row['lat']: .4f}   {row['lon']: .4f}")
                
                for h in range(1,25):
                    try:
                        file.write(f"    {row.iloc[1+h]: .4f}    {row.iloc[25+h]: .4f}")
                    except:
                        file.write(f"    {0: .4f}    {0: .4f}")
                                
                    if h == 24:
                        file.write('\n')     

    except:
        print(f'{date} has no files in met directory')
        pass

def process_mrc(i, concat,simname=None):
    
    rec = concat.isel(time=i)
    #get the datetime values from that instant
    try:
        dt = pd.to_datetime(rec.time.values)
    except:
        try:
            dt = datetime.datetime.strptime(str(rec.time.values),'%Y-%m-%d %H:%M:%S')
        except:
            raise ValueError('Datetime from the dataset in on an unknown format')
    
    #transforms it from xarray datset to pandas dataframe to facilitate the processes of adjusting values
    df = rec.to_dataframe().reset_index()
    df = df.fillna(0)
    df = df.drop(['time'],axis=1)
    #pivoting it in order to create the same pattern in .mrc files
    df = df.pivot(index = ['lat','lon'],columns='depth',values = ['thetao','uo','vo']).reset_index()
    #join colum information to 
    df.columns = [pair[0]+str(pair[1]).split('.')[0] for pair in df.columns]
    #dropping temperature columns
    df = df.drop(['thetao10','thetao30','thetao120'],axis=1)
    #sort first by latitude and then by longitude
    df = df.sort_values(['lon','lat'])
    df.columns = ['lat','lon','SST','u_srf','u_10m','u_30m','u_120m','v_srf','v_10m','v_30m','v_120m']
    df = df[['lat','lon','SST','u_srf','v_srf','u_10m','v_10m','u_30m','v_30m','u_120m','v_120m']]

    #making sure that .mrc files in 0 hour are written as 24
    #this code also makes sure that the file is written correctly even in the transition of months
    if dt.hour == 0:
        hour = 24
        day = (dt - datetime.timedelta(hours=1)).day
        month = (dt - datetime.timedelta(hours=1)).month
    else:
        hour = dt.hour
        day = dt.day
        month = dt.month

    #writing the current files
    with open(f'cases/{simname}/oce_files/merc{dt.year-2000:02d}{month:02d}{day:02d}{hour:02d}.mrc', 'w') as f:
        f.write(f"Ocean forecast data for {day:02d}/{month:02d}/{dt.year} {hour:02d}:00\n")
        f.write("Subregion of the Global Ocean:\n")
        f.write(f"{df.lon.min():02.2f}  {df.lon.max():02.2f}  {df.lat.min():02.2f} {df.lat.max():02.2f}   {len(rec.lon)}   {len(rec.lat)}   Geog. limits\n")
        f.write(f"{len(df)}   0.0\n")
        f.write("lat        lon        SST        u_srf      v_srf      u_10m      v_10m       u_30m      v_30m      u_120m     v_120m\n")
        
        for index, row in df.iterrows():
            f.write(f"{row['lat']:<10.4f}    {row['lon']:<10.4f}    {row['SST']:<10.4f}     {row['u_srf']:<10.4f}    {row['v_srf']:<10.4f}     {row['u_10m']:<10.4f}    {row['v_10m']:<10.4f}     {row['u_30m']:<10.4f}    {row['v_30m']:<10.4f}     {row['u_120m']:<10.4f}    {row['v_120m']:<10.4f}\n")

# Define a function to process multiple time steps in parallel
def parallel_processing_mrc(concat,simname=None):
    num_time_steps = len(concat.time)
    pool = multiprocessing.Pool()  # Create a pool of workers
    results = [pool.apply_async(process_mrc, args=(i, concat,simname)) for i in range(num_time_steps)]
    pool.close()  # Close the pool, no more tasks can be submitted
    pool.join()   # Wait for all worker processes to finish