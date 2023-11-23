import motuclient
import argparse
import os
import subprocess
import datetime
import pandas as pd

class MotuOptions:
    def __init__(self, attrs: dict):
        super(MotuOptions, self).__setattr__("attrs", attrs)

    def __setattr__(self, k, v):
        self.attrs[k] = v

    def __getattr__(self, k):
        try:
            return self.attrs[k]
        except KeyError:
            return None  


def download(date_start, date_end, 
                    depth_min, depth_max,
                    lat, lon, delta_latlon,
                    service_id,product_id,var,
                    output_path,output_name,
                    user,psd):


    data_request_options_dict_manual = {
        "service_id": service_id,
        "product_id": product_id,
        "date_min": date_start,
        "date_max": date_end,
        "longitude_min": lon-delta_latlon,
        "longitude_max": lon+delta_latlon,
        "latitude_min": lat-delta_latlon,
        "latitude_max": lat+delta_latlon,
        "depth_min": depth_min,
        "depth_max": depth_max,
        "variable": var,
        # "motu": "https://nrt.cmems-du.eu/motu-web/Motu",
        "motu": "http://my.cmems-du.eu/motu-web/Motu",
        "out_dir": output_path,
        "out_name": output_name,
        "auth_mode": "cas",
        "user": user,
        "pwd": psd
        }
    
    return data_request_options_dict_manual

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Download Mercator")

    # Define command-line arguments
    parser.add_argument('date_start', type=str, help='Start date of the simulation')
    parser.add_argument('date_end', type=str, help='End date of the simulation')
    parser.add_argument('depth_min', type=float, help='Minimum depth value')
    parser.add_argument('depth_max', type=float, help='Maximum depth value')
    parser.add_argument('lat', type=float, help='Latitude value')
    parser.add_argument('lon', type=float, help='Longitude value')
    parser.add_argument('delta_latlon', type=float, help='Delta latitude and longitude value')
    parser.add_argument('down', type=str, help='Indication of global or local data')
    parser.add_argument('output_path', type=str, default='.', help='Output path (default: current directory)')
    parser.add_argument('user', type=str, default='.', help='Copernicus Username')
    parser.add_argument('psd', type=str, default='.', help='Copernicus Password')

    parser

    # Parse the command-line arguments
    args = parser.parse_args()

    outputs = []

    #get datetime of sim and compare with limits to mercator datasets
    conf_date = args.date_start.split('T')[0]

    if args.down == 'global':
        loc='glo'
        if pd.to_datetime(conf_date) < pd.to_datetime('2020-11-21') == False:
            service_id = 'GLOBAL_ANALYSISFORECAST_PHY_001_024-TDS'
            product_id_curr = 'cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m'
            product_id_temp = 'cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m'
        else:
            service_id = 'GLOBAL_MULTIYEAR_PHY_001_030-TDS'
            product_id_curr = 'cmems_mod_glo_phy_my_0.083_P1D-m'
            product_id_temp = 'cmems_mod_glo_phy_my_0.083_P1D-m'

    else:
        loc='med'
        if pd.to_datetime(conf_date) < pd.to_datetime('2020-11-21') == False:
            service_id = 'MEDSEA_ANALYSISFORECAST_PHY_006_013-TDS'
            product_id_curr = 'cmems_mod_med_phy-cur_anfc_4.2km_P1D-m'
            product_id_temp = 'cmems_mod_med_phy-tem_anfc_4.2km_P1D-m'
        else:
            service_id = 'MEDSEA_MULTIYEAR_PHY_006_004-TDS'
            product_id_curr = 'med-cmcc-cur-rean-d'
            product_id_temp = 'med-cmcc-tem-rean-d'

    for var in ['curr','temp']:

        if var == 'curr':
            product_id = product_id_curr
            varr = ['uo','vo']
        elif var == 'temp':
            product_id = product_id_temp
            varr = ['thetao']
        
        outputname = f'{var}_tempo.nc'
        outputs.append(outputname)

        request = download(args.date_start, args.date_end, 
                        args.depth_min, args.depth_max, 
                        args.lat, args.lon, args.delta_latlon, 
                        service_id,product_id,varr,
                        args.output_path,outputname,
                        args.user,args.psd)

        motuclient.motu_api.execute_request(MotuOptions(request))

#   real start date
    real_date = datetime.datetime.strptime(args.date_start[0:10],'%Y-%m-%d')
    real_date = real_date + datetime.timedelta(days=1) 

    subprocess.run(f'cdo merge {args.output_path}{outputs[0]} {args.output_path}{outputs[1]}\
                    {args.output_path}CMEMS_{loc}_{real_date.year}{real_date.month}{real_date.day}_raw.nc',shell=True,check=True)

    subprocess.run(f'rm {args.output_path}{outputs[0]} {args.output_path}{outputs[1]}',shell=True,check=True)

