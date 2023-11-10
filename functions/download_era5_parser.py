import numpy as np
import pdb
import datetime
import os
import cdsapi
import argparse
import subprocess
import time

def get_era5(xmin,xmax,ymin,ymax,start_date,end_date,output_path):
    server = cdsapi.Client()

    start_date = datetime.datetime.strptime(start_date,'%Y-%m-%d')
    end_date = datetime.datetime.strptime(end_date,'%Y-%m-%d')  

    days = (end_date-start_date).days

    outputs = []

    for i in range(0,days):
  
        date = start_date + datetime.timedelta(days=i)

        outputname = output_path + f'era5_winds10_{str(date.year)+str(date.month).zfill(2)+str(date.day).zfill(2)}.nc'
       
        server.retrieve(
        'reanalysis-era5-single-levels',
        {
          'product_type': 'reanalysis',
          'format': 'netcdf',
          'variable': [
              '10m_u_component_of_wind', '10m_v_component_of_wind',
          ],
          'year' :  str(date.year),
          'month':  str(date.month).zfill(2),
          'day'  :  str(date.day).zfill(2),
          'time' : [
              '00:00', '01:00', '02:00',
              '03:00', '04:00', '05:00',
              '06:00', '07:00', '08:00',
              '09:00', '10:00', '11:00',
              '12:00', '13:00', '14:00',
              '15:00', '16:00', '17:00',
              '18:00', '19:00', '20:00',
              '21:00', '22:00', '23:00',
          ],
          'area': [
              ymax, xmin, ymin,
              xmax,
          ],
        },
        outputname)
        
        time.sleep(2)
        subprocess.run(f'ncrename -O -d longitude,lon -d latitude,lat -v longitude,lon -v latitude,lat -v u10,U10M -v v10,V10M {outputname}',shell=True,check=True)


if __name__ == '__main__':

    # Script to download daily ERA-5 files

    parser = argparse.ArgumentParser(description='Download wind fields for a specific area and time window.')
    parser.add_argument('lat', type=float, help='Latitude value')
    parser.add_argument('lon', type=float, help='Longitude value')
    parser.add_argument('delta_latlon', type=float, help='Delta latitude and longitude value')
    parser.add_argument('date_min', type=str, help='Start date in yyyy-mm-dd format')
    parser.add_argument('date_max', type=str, help='End date in yyyy-mm-dd format')
    parser.add_argument('output_path', type=str, default='./', help='Output path (default: current directory)')
    args = parser.parse_args()

    #Set your area of interest
    xmin = np.floor(float(args.lon - args.delta_latlon))
    xmax = np.ceil(float(args.lon + args.delta_latlon))
    ymin = np.floor(float(args.lat - args.delta_latlon))
    ymax = np.ceil(float(args.lat + args.delta_latlon))

    # Set your period of interest
    start_date=args.date_min
    end_date=args.date_max

    print('********************************************')
    print('PREPARING ERA5 WIND DATA - MEDSLIK II FORMAT')
    print('Start date :' + start_date)
    print('End date :' + end_date)
    print('********************************************')


    get_era5(xmin,xmax,ymin,ymax,start_date,end_date,args.output_path)
