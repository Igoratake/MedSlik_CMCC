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
simname        = 'elba'
sim_date       = '17/05/2014' ### Simulation start day  - format DD/MM/YYYY (string)
sim_hour       = '05:38'      ### Simulation start hour - format HH:mm (string)
longitude       = 10        ### Longitude of Simulation spill location - format Decimal degrees (float)
latitude        = 42.975          ### Latitude of Simulation spill  - format Decimal degrees (float)
sim_lenght      = 24              ### Length of the simulation - format hours (int)
spill_duration  = 0              ### Duration of the spill - format hours (int)
oil_api         = 22.30               ### Oil API - format (float)
oil_volume      = 63           ### Volume of oil in tons - format (float) 
use_satellite   = False            ### Usage of Satellite imagery to model multiple slicks - True/False
use_slk_contour = True            ### Usage of slicks countours - True/False
separate_slicks = True         ### If the sim needs different slicks to have different properties, different sims have to be run
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

     # Simulation dates
    dt_sim = validate_date(sim_date + ' ' + sim_hour)

    if isinstance(dt_sim,str):
        raise ValueError('Wrong date format.')
    
    #Initial date is always one day prior due to Medslik-II interpolation
    dtini = dt_sim - datetime.timedelta(days=1)
         
    #End date, by safety margin is two days after the sim start + sim length
    dtend = dt_sim + datetime.timedelta(days= (sim_lenght/24) + 2) 

    #List of dates between initial and end date
    date_list = pd.date_range(dtini, dtend, freq='D')
    
     ##### SEA / OCEAN ##### 

    #opening all files in the directory and concatenating them automatically through open_mfdataset
    concat = xr.open_mfdataset(f'cases/{simname}/oce_files/*.nc',combine='nested')
    concat = concat.drop_duplicates(dim="time", keep="last")

    #Interpolating the values in time, transforming it from daily to hourly values
    concat = concat.resample(time="1H").interpolate("linear")

    #iterating at each hour to generate the .mrc files
    for i in range(0,len(concat.time)):            
        #select the i time instant
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

     ##### WINDS ##### 

    concat = xr.open_mfdataset(f'cases/{simname}/met_files/*.nc',combine='nested')
    concat = concat.drop_duplicates(dim="time", keep="first")
    concat = concat.resample(time="1H").interpolate("linear")

    #iterating at each hour to generate the .eri files
    for date in date_list:            
        
        try:
            date1 = f'{date.year}-{date.month:02d}-{date.day:02d} 01:00'
            incremental = date + datetime.timedelta(days=1)
            date2 = f'{incremental.year}-{incremental.month:02d}-{incremental.day:02d} 01:00'
            met = concat.sel(time = slice(date1,date2))
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
                            file.write(f"    {row[1+h]: .4f}    {row[26+h]: .4f}")
                        except:
                            file.write(f"    {0: .4f}    {0: .4f}")
                                    
                        if h == 24:
                            file.write('\n')
        except:
            print(f'{date} has no files in met directory')
            pass

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
