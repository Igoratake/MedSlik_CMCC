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

    path = '/Medslik-II/cases/elba/out_files/MDK_SIM_2014_05_17_0538_elba_slick1/'
    filename = path + 'spill_properties.nc'
    ds1 = xr.open_dataset(filename)

    path = '/Medslik-II/cases/elba/out_files/MDK_SIM_2014_05_17_0538_elba_slick2/'
    filename = path + 'spill_properties.nc'
    ds2 = xr.open_dataset(filename)

    path = '/Medslik-II/cases/elba/out_files/MDK_SIM_2014_05_17_0538_elba_slick3/'
    filename = path + 'spill_properties.nc'
    ds3 = xr.open_dataset(filename)

    path = '/Medslik-II/cases/elba/out_files/MDK_SIM_2014_05_17_0538_elba_slick4/'
    filename = path + 'spill_properties.nc'
    ds4 = xr.open_dataset(filename)

    path = '/Medslik-II/cases/elba/out_files/MDK_SIM_2014_05_17_0538_elba_slick5/'
    filename = path + 'spill_properties.nc'
    ds5 = xr.open_dataset(filename)

    path = '/Medslik-II/cases/elba/out_files/MDK_SIM_2014_05_17_0538_elba_slick6/'
    filename = path + 'spill_properties.nc'
    ds6 = xr.open_dataset(filename)

    path = '/Medslik-II/cases/elba/out_files/MDK_SIM_2014_05_17_0538_elba_slick7/'
    filename = path + 'spill_properties.nc'
    ds7 = xr.open_dataset(filename)

    lon_min, lon_max = 9.5,10.5
    lat_min, lat_max = 42.5,43.5

    rec = land.cx[lon_min:lon_max, lat_min:lat_max]

    for t in range(0,24):
        dss1 = ds1.isel(time=t)
        dss2 = ds2.isel(time=t)
        dss3 = ds3.isel(time=t)
        dss4 = ds4.isel(time=t)
        dss5 = ds5.isel(time=t)
        dss6 = ds6.isel(time=t)
        dss7 = ds7.isel(time=t)
        
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.set_facecolor('#ADD8E6')
        rec.plot(ax=ax,color="#FFFDD0", edgecolor='black', zorder = 1000)
        ax.scatter(dss1.longitude,dss1.latitude,dss1.non_evaporative_volume,color='r')
        ax.scatter(dss2.longitude,dss3.latitude,dss2.non_evaporative_volume,color='r')
        ax.scatter(dss3.longitude,dss3.latitude,dss3.non_evaporative_volume,color='r')
        ax.scatter(dss4.longitude,dss4.latitude,dss4.non_evaporative_volume,color='r')
        ax.scatter(dss5.longitude,dss5.latitude,dss5.non_evaporative_volume,color='r')
        ax.scatter(dss6.longitude,dss6.latitude,dss6.non_evaporative_volume,color='r')
        ax.scatter(dss7.longitude,dss7.latitude,dss7.non_evaporative_volume,color='r')
        
        plt.xlim(lon_min,lon_max)
        plt.ylim(lat_min,lat_max)

        plt.xlabel('Longitude')
        plt.ylabel('Latitude')

        plt.title(f'time={t}')

        plt.grid()

        plt.savefig(f'/Medslik-II/cases/elba/out_files/figures/fig_elba_{t:03d}.png',dpi=50)
