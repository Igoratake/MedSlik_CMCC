#geo libs
import xarray as xr
import geopandas as gpd
import cartopy.crs as ccrs
from shapely.geometry import Point
from folium.plugins import HeatMap
from streamlit_folium import st_folium,folium_static

#numerical and plotting libs
import numpy as np
import matplotlib.pyplot as plt
import plotly.figure_factory as ff

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

# ================= DATA FORMAT INPUTS ================= #
type            = 'analysis' ### other options: reanalysis
time_res        = 'day'      ### other options: 1h
process_files   = True       ### False if data is already present

# ================= MEDSLIK MODEL INPUTS ================= #
simdir          = 'cases/'
simname         = 'taranto_manual' ### Simulation name - format (string)
sim_date        = '23/11/2023'     ### Simulation start day  - format DD/MM/YYYY (string)
sim_hour        = '13:00'          ### Simulation start hour - format HH:mm (string)
longitude       = 40.09278          ### Longitude of Simulation spill location - format Decimal degrees (float)
latitude        = 17.21558          ### Latitude of Simulation spill  - format Decimal degrees (float)
sim_lenght      = 24              ### Length of the simulation - format hours (int)
spill_duration  = 0              ### Duration of the spill - format hours (int)
oil_api         = 26.96               ### Oil API - format (float)
oil_volume      = 1000           ### Volume of oil in tons - format (float) 
use_satellite   = False            ### Usage of Satellite imagery to model multiple slicks - True/False
use_slk_contour = False            ### Usage of slicks countours - True/False
number_slick    = 1                ### Number of slicks to be simulated - format (int)

# ================= ENVIRONMENTAL FORCINGS  ================= #
#Select either a fixed delta for lat and lon or specify the bounding box
loc                 = 'med'
delta               = 0.0              ### Standard delta to collect an area around lat and discharge point - format degrees (float)
if delta == 0:
    lat_min         = 38.2927818
    lat_max         = 41.8927803
    lon_min         = 14.8626471
    lon_max         = 19.5685139
else:
    lat_min         = latitude - delta
    lat_max         = latitude + delta
    lon_min         = longitude - delta
    lon_max         = longitude + delta

# Obtaining spill rate from oil volume and spill duration
if spill_duration != 0:
    spill_rate = oil_volume/spill_duration
else:
    spill_rate = oil_volume

if __name__ == '__main__':
   
    ############# INITIAL CHECKING AND DOWNLOAD #############
    # Check if coordinates are on sea
    sea =  check_land(longitude,latitude)    

    if sea==0:
        raise ValueError('Coordinates values are on land. Try another location')
    
    # Simulation dates
    dt_sim = validate_date(sim_date + ' ' + sim_hour)

    if isinstance(dt_sim,str):
        raise ValueError('Wrong date format.')
    
    # if loc == None:
    #     if 30.37 < float(latitude) < 45.7 and -17.25 < float(longitude) < 35.9: 
    #         print('Coordinates lie within Mediterranean Sea')
    #         loc = 'med'
    #     else:
    #         loc = 'glob'
    
    if type == 'analysis':
        if loc == 'med' and time_res == 'day':
            if dt_sim.year < 2021:
                raise ValueError('Data not available for analysis')
            else:
                path = '/data/inputs/METOCEAN/historical/model/ocean/CMS/MedSea/CMCC/analysis/day/'
            
        elif loc == 'med' and time_res == 'hour':
            if dt_sim.year < 2023 and dt_sim.month < 5:
                raise ValueError('Data not available for analysis')
            else:
                raise ValueError('Pre process not yet implemented')
                # path = '/data/inputs/METOCEAN/historical/model/ocean/CMS/MedSea/CMCC/analysis/1h/'
           
        elif loc == 'glob':
            if dt_sim.year < 2012:
                raise ValueError('Data not available for analysis')
            else:
                raise ValueError('Pre process not yet implemented')
                # path = '/data/inputs/METOCEAN/historical/model/ocean/CMS/GlobOce/MERCATOR/analysis/day/'
            
    elif type == 'reanalysis':
        if loc == 'med' and time_res == 'hour':
            raise ValueError('Temperature Data not available for this time resolution')

        if loc == 'med' and time_res == 'day':
            if dt_sim.year < 2021:
                raise ValueError('Data not available for reanalysis')
            else:
                path = '/data/inputs/METOCEAN/historical/model/ocean/CMS/MedSea/CMCC/reanalysis/day/'
    
        elif loc == 'glob':
            if dt_sim.year > 2010:
                raise ValueError('Data not available for reanalysis')
            else:
                path = '/data/inputs/METOCEAN/historical/model/ocean/CMS/GlobOce/CMCC/reanalysis/day/'            
    
    else:
       raise ValueError('Mispelled type of data')     
        
    ############# PRE PROCESSING DATA AND FILES #############
    #Creating all directories    
    subprocess.run(f'mkdir {simdir}{simname}/',shell=True)
    subprocess.run(f'mkdir {simdir}{simname}/oce_files',shell=True)
    subprocess.run(f'mkdir {simdir}{simname}/met_files',shell=True)
    subprocess.run(f'mkdir {simdir}{simname}/bnc_files',shell=True)
    subprocess.run(f'mkdir {simdir}{simname}/xp_files',shell=True)
    subprocess.run(f'mkdir {simdir}{simname}/out_files',shell=True)
    subprocess.run(f'mkdir {simdir}{simname}/detections',shell=True)

    #Initial date is always one day prior due to Medslik-II interpolation
    dtini = dt_sim - datetime.timedelta(days=1)
         
    #End date, by safety margin is two days after the sim start + sim length
    dtend = dt_sim + datetime.timedelta(days= (sim_lenght/24) + 2) 

    #List of dates between initial and end date
    date_list = pd.date_range(dtini, dtend, freq='D')

    if process_files == True:
    ############################ PRE PROCCESING SEA / OCEAN AND WIND FILES ############################
    
    ##### SEA / OCEAN ##### 
        for date in date_list:
            day = str(date.day).zfill(2)
            month = str(date.month).zfill(2)
            year = str(date.year)

            #reading all the daily files in respective directory and storing it into case dir
            for var,varr in zip(['TEMP','RFVL'],['temp','curr']):
                filename = f'{year}{month}{day}*{var}*.nc'

                path_f = os.path.join(path,year,month,filename)
                #searches for the file in the data dir
                file = glob(path_f)[0]

                ds = xr.open_dataset(file)  

                #Trimming into smaller regional box
                ds = ds.sel(lon = slice(lon_min,lon_max),
                            lat = slice(lat_min,lat_max))
                
                #Selecting only 4 layers
                tot = ds.sel(depth=[0,10,30,120],method='nearest')
                #Modifying labels to simplfy drop in temperature columns
                tot['depth'] = [0,10,30,120]

                #saves the daily current or temperature netcdf in the case dir
                tot.to_netcdf(f'cases/{simname}/oce_files/{varr}_{year}{month}{day}_{simname}.nc')
            
        #opening all files in the directory and concatenating them automatically through open_mfdataset
        concat = xr.open_mfdataset(f'cases/{simname}/oce_files/*.nc')

        #Interpolating the values in time, transforming it from daily to hourly values
        concat = concat.resample(time="1H").interpolate("linear")

        #iterating at each hour to generate the .mrc files
        for i in range(0,len(concat.time)):            
            #select the i time instant
            rec = concat.isel(time=i)
            #get the datetime values from that instant
            dt = pd.to_datetime(rec.time.values)
            
            #transforms it from xarray datset to pandas dataframe to facilitate the processes of adjusting values
            df = rec.to_dataframe().reset_index()
            df = df.fillna(0)
            df = df.drop(['time','bottomT'],axis=1)
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
                f.write(f"{df.lon.min():02.2f}  {df.lon.max():02.2f}  {df.lat.min():02.2f} {df.lat.max():02.2f}   {len(tot.lon)}   {len(tot.lat)}   Geog. limits\n")
                f.write(f"{len(df)}   0.0\n")
                f.write("lat        lon        SST        u_srf      v_srf      u_10m      v_10m       u_30m      v_30m      u_120m     v_120m\n")
                
                for index, row in df.iterrows():
                    f.write(f"{row['lat']:<10.4f}    {row['lon']:<10.4f}    {row['SST']:<10.4f}     {row['u_srf']:<10.4f}    {row['v_srf']:<10.4f}     {row['u_10m']:<10.4f}    {row['v_10m']:<10.4f}     {row['u_30m']:<10.4f}    {row['v_30m']:<10.4f}     {row['u_120m']:<10.4f}    {row['v_120m']:<10.4f}\n")

            
        ##### WIND #####
        for date in date_list:
            day = str(date.day).zfill(2)
            month = str(date.month).zfill(2)
            year = str(date.year)
            #Searching the path based on dates
            path = f'/data/inputs/METOCEAN/historical/model/atmos/ECMWF/ERA5/reanalysis/1h/netcdf/{year}/{month}/'
            file = path + f'{year}{month}{day}.nc'

            met = xr.open_dataset(file)

            #Trimming into smaller regional box
            met = met.sel(lon = slice(lon_min,lon_max),
                          lat = slice(lat_max,lat_min))
            
            #saving the sampled netcdf in cases directory
            met.to_netcdf(f'cases/{simname}/met_files/{varr}_{year}{month}{day}_{simname}.nc')

            #getting date from the netcdf
            dt = pd.to_datetime(met.time[0].values)
            
            met = met.isel(time=slice(0,-1))
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
                        dt_str = str(dt)
                        u_str = 'U10M' + dt_str.replace(' 00',f' {h:02d}')
                        v_str = 'V10M' + dt_str.replace(' 00',f' {h:02d}')
                        try:
                            file.write(f"    {row[u_str]: .4f}    {row[v_str]: .4f}")
                        except:
                            file.write(f"    {0: .4f}    {0: .4f}")
                                    
                        if h == 24:
                            file.write('\n')

    #sending curr asc files files to TEMP dir to execute medlisk
    subprocess.run([f'cp {simdir}{simname}/oce_files/*.mrc MEDSLIK_II_3.01/RUN/TEMP/OCE/'],shell=True)
    
    #sending wind asc files files to TEMP dir to execute medlisk
    subprocess.run([f'cp {simdir}{simname}/met_files/*.eri MEDSLIK_II_3.01/RUN/TEMP/MET/'],shell=True)

    ############################ GENERATING BATHYMETRY AND COASTLINE ############################

    # Process for Bathymetry and Coastline Files
    grid_string = glob(f'{simdir}{simname}/oce_files/*.nc')[0] 

    # Bathymetry for gebco 2023
    subprocess.run([f'{sys.executable}', 'scripts/pre_processing/preproc_gebco_mdk2.py', 
                    'data/gebco/GEBCO_2023.nc',
                    grid_string, f'{simdir}{simname}/bnc_files/'])

    # gshhs in intermediate resolution
    subprocess.run([f'{sys.executable}', 'scripts/pre_processing/preproc_gshhs_mdk2.py', 
                    'data/gshhs/h/GSHHS_h_L1.shp',
                    grid_string, f'{simdir}{simname}/bnc_files/'])

    ############################ PREPARING FILES TO RUN THE SIMULATION ############################

    # prepare medslik_II.for and config1.txt
    print('Preparing configuration files... ')

    # get dimensions from ncfiles
    for var in ['votemper','thetao','uo','vo']:
        try:
            my_o = xr.open_dataset(glob(f'cases/{simname}/oce_files/*nc')[0]).isel(depth=0,time=0)[var].values.shape
            found_variable = True
        except:
            continue
        if not found_variable:
            raise ValueError("Check the sea state variables available in your dataset")
        found_variable = False

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

    # # modify config_1.txt
    print('...config1.txt...')

    config = "LC_CTYPE=C && LANG=C"
    config_file = f'cases/{simname}/xp_files/config1.txt'

    subprocess.run([f'cp scripts/templates/config1_template_0.txt {config_file}'],shell=True)

    #adding spill date
    subprocess.run([f"{config} sed -i s#RUNNAME#{simname}#g {config_file}"],shell=True)
    subprocess.run([f"{config} sed -i s/DD/{dt_sim.day:02d}/ {config_file}"],shell=True)
    subprocess.run([f"{config} sed -i s/MM/{dt_sim.month:02d}/ {config_file}"],shell=True)
    subprocess.run([f"{config} sed -i s/YY/{dt_sim.year-2000:02d}/ {config_file}"],shell=True)
    subprocess.run([f"{config} sed -i s/c_Hour/{dt_sim.hour:02d}/ {config_file}"],shell=True)
    subprocess.run([f"{config} sed -i s/c_minute/{dt_sim.minute:02d}/ {config_file}"],shell=True)   

    #adding simulation length
    subprocess.run([f"{config} sed -i s/SIMLENGTH/{sim_lenght:04d}/ {config_file}"],shell=True)       

    #  adding spill coordinates - dd for degrees and mm for minutes
    # Latitude
    dd = int(latitude)
    mm = (float(latitude)-dd)*60          
    subprocess.run([f"{config} sed -i s/LATd/{dd:02d}/ {config_file}"],shell=True)
    subprocess.run([f"{config} sed -i s/LATm/{mm:.3f}/ {config_file}"],shell=True)
    
    # Longitude
    dd = int(longitude)
    mm = (float(longitude)-dd)*60
    subprocess.run([f"{config} sed -i s/LONd/{dd:02d}/ {config_file}"],shell=True)
    subprocess.run([f"{config} sed -i s/LONm/{mm:.3f}/ {config_file}"],shell=True)

    # spill duration
    subprocess.run([f"{config} sed -i s/SDUR/{spill_duration:04d}/ {config_file}"],shell=True)

    # spill volume
    subprocess.run([f"{config} sed -i s/SRATE/{spill_rate:08.2f}/ {config_file}"],shell=True)

    # oil characteristics
    subprocess.run([f"{config} sed -i s/APIOT/{oil_api}/ {config_file}"],shell=True)

    #number of slicks
    subprocess.run([f"{config} sed -i s/N_SLICK/{number_slick}/ {config_file}"],shell=True)

    #slick countour
    if use_slk_contour == True:
        slik = 'YES'
        with open('MEDSLIK_II_3.01/RUN/slick_countour.txt', 'r') as file1:
            content = file1.read()
        subprocess.run([f'cp MEDSLIK_II_3.01/RUN/slick_countour.txt {simdir}{simname}/xp_files/'],shell=True)
        with open(config_file, 'a') as file2:
        # Append the contents of the first file to config file
            file2.write(content)
    else:
        slik = 'NO'  

    #Writing that will use slick countor
    subprocess.run([f"{config} sed -i s/SSLICK/{slik}/ {config_file}"],shell=True)

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

    ############################ RUN THE SIMULATION ############################

    # Compile and start running
    subprocess.run([f'cd MEDSLIK_II_3.01/RUN/; sh MODEL_SRC/compile.sh; ./RUN.sh'],shell=True,check=True)

    ############################ SEND OUTPUT FILES TO CASE DIR   ############################

    subprocess.run([f'cp -r MEDSLIK_II_3.01/OUT/MDK_SIM_{dt_sim.year:02d}_{str(dt_sim.month).zfill(2)}_{str(dt_sim.day).zfill(2)}_{simname}/ {simdir}{simname}/out_files/'],shell=True)
    subprocess.run([f'rm -rf {simdir}{simname}/out_files/MET {simdir}{simname}/out_files/OCE'],shell=True)