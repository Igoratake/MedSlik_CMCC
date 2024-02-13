import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point, Polygon

import matplotlib
norm = matplotlib.colors.Normalize()

from shapely.geometry import box, Polygon

land = gpd.read_file('data/gshhs/GSHHS_shp/f/GSHHS_f_L1.shp')

if __name__ == '__main__':

    path = '/Medslik-II/cases/paria_remake/out_files/MDK_SIM_2017_04_24_1300_paria_remake/'
    filename = path + 'spill_properties.nc'
    ds = xr.open_dataset(filename)

    lon_min,lon_max = -64.1,-63.5
    lat_min,lat_max = 10.6,11

    rec = land.cx[lon_min:lon_max, lat_min:lat_max]

    for t in range(0,48):
        dss = ds.isel(time=t)
        
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.set_facecolor('#ADD8E6')
        rec.plot(ax=ax,color="#FFFDD0", edgecolor='black', zorder = 1000)
        ax.scatter(dss.longitude,dss.latitude,dss.non_evaporative_volume,color='r')
        
        plt.xlim(lon_min,lon_max)
        plt.ylim(lat_min,lat_max)

        plt.xlabel('Longitude')
        plt.ylabel('Latitude')

        plt.title(f'time={t}')

        plt.grid()

        plt.savefig(f'{path}/figures/fig_test_{t:03d}.png',dpi=50)