import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point, Polygon

import matplotlib
norm = matplotlib.colors.Normalize()

from shapely.geometry import box, Polygon

land = gpd.read_file('/Users/iatake/Dropbox (CMCC)/Work/MEDSLIK-II and Pyslick/Medslik-II/data/gshhs/GSHHS_shp/f/GSHHS_f_L1.shp')

if __name__ == '__main__':

    path = '/Users/iatake/Dropbox (CMCC)/Work/MEDSLIK-II and Pyslick/Medslik-II_results/ORIGINAL_DATA/LEBANON201/'

    ds_particles = xr.open_dataset(path + 'spill_properties.nc')
    ds_concentration = xr.open_dataset(path + 'spill_properties_concentration.nc')

    lon_min,lon_max = 35,36
    lat_min,lat_max = 33.6,34.6

    rec = land.cx[lon_min:lon_max, lat_min:lat_max]

    for t in [0,12,24,48,76,92,96,120,200,300,400,439]:
        ds_p = ds_particles.isel(time=t)
        ds_c = ds_concentration.isel(time=t)

        for plot in ['particle','concentration']:
        
            fig, ax = plt.subplots(figsize=(10, 8))
            ax.set_facecolor('#ADD8E6')
            
            # Ploting coastline
            rec.plot(ax=ax,color="#FFFDD0", edgecolor='black', zorder = 1000,aspect=1)

            if plot == 'particle':
                ax.scatter(ds_p.longitude,ds_p.latitude,ds_p.non_evaporative_volume,color='r')

            else:
                ds_c = xr.where(ds_c.concentration > 0.000000000001, ds_c.concentration,np.nan)
                ds_c.plot(ax=ax)
            
            plt.xlim(lon_min,lon_max)
            plt.ylim(lat_min,lat_max)

            plt.xlabel('Longitude')
            plt.ylabel('Latitude')

            plt.title(f'time={t}')

            plt.grid()

            plt.savefig(path + f'figures/fig_lebanon201_{t:03d}_{plot}.png',dpi=100)

            plt.close()