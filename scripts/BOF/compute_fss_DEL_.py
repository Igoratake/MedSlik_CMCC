import pysteps
from shapely.geometry import Polygon
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import pandas as pd
import pdb
import shapely
import sys
import time
from scipy.interpolate import NearestNDInterpolator
import os
import re
from datetime import  *
from matplotlib import colors as c
from mpl_toolkits.basemap import Basemap
import geopandas as gpd
import argparse
from scipy.ndimage.filters import uniform_filter


parser = argparse.ArgumentParser()
parser.add_argument('--forecastfolder', '-f', help='path to your medslik ii forecast')
parser.add_argument('--observationfolder', '-o', help='path to spill shapefile')
parser.add_argument('--outputfolder', '-u', help='output path')
args = parser.parse_args()

"""
Script para calculo do Fracticional Skill Score em previsoes de deriva de petroleo
"""
def erre_sbl_2021(a,b,c,d):
    return (a+b)/(a+b+c+d)
def esse_sbl_2021(a,b,c,d):
    return (a+c)/(a+b+c+d)
def fbias_sbl_2021(a,b,c,d):
    return (a+b)/(a+c)

def barrelToTonnes(oil_density):
    rbm3 = 0.158987
    return 1 / (rbm3 * (oil_density / 1000))
  
def get_surface_parcels(fname,time_index):

    f = netCDF4.Dataset(fname)

    lat = np.longdouble(f.variables['latitude'][time_index,:])
    lon = np.longdouble(f.variables['longitude'][time_index,:])
    evaporative_volume = f.variables['evaporative_volume'][time_index,:]
    non_evaporative_volume = f.variables['non_evaporative_volume'][time_index,:]
    particle_status = f.variables['particle_status'][time_index, :]
    oil_density = f.variables['non_evaporative_volume'].oil_density
    barrel2tonnes = barrelToTonnes(oil_density)   

    floatingParticles=np.logical_and(particle_status > 0, particle_status <= 2).nonzero()[0]

    lons_f=lon[floatingParticles]
    lats_f=lat[floatingParticles]
    surf_volume = (evaporative_volume[floatingParticles] + non_evaporative_volume[floatingParticles]) / barrel2tonnes

    return lons_f,lats_f,surf_volume

def make_poly_grid(xmin,xmax,ymin,ymax,cell_size,crs):

    cols = list(np.arange(xmin, xmax+cell_size, cell_size))
    rows = list(np.arange(ymin, ymax+cell_size, cell_size))

    polygons = []
    for x in cols[:-1]:
        for y in rows[:-1]:
            polygons.append(Polygon([(x,y), (x+cell_size, y), (x+cell_size, y+cell_size), (x, y+cell_size)]))
            
    cell = gpd.GeoDataFrame({'geometry':polygons}).set_crs('EPSG:4326',allow_override=True)

    return cell
 
def get_obs_date(observation_shp):
    observation_df = gpd.read_file(observation_shp)
    observation_date_string = observation_df.IDENTIFIER[0]
    year=int(observation_date_string[0:4])
    month=int(observation_date_string[4:6])
    day=int(observation_date_string[6:8])
    hour=int(observation_date_string[9:11])
    minute=int(observation_date_string[11:13])
    
    obs_date = date(year,month,day).toordinal()
    obs_date = obs_date + hour/24.
    return obs_date
    
def get_mdksim_id(simulation_folder):
    
    # get sim _ info
    with open(simulation_folder + '/config1.txt',encoding='iso-8859-1') as f:
        lines = f.readlines()
    ss = lines[0]
    rr = re.search('=(.*)\n',ss)
    sim_name = rr.group(1)
    return sim_name
    
def get_mdksim_date(simulation_folder):
    
    # get sim _ info
    with open(simulation_folder + '/config1.txt',encoding='iso-8859-1') as f:
        lines = f.readlines()
    # sim_length
    ss = lines[3]
    rr = re.search('=(.*)\n',ss)
    sim_length = int(lines[1][11:15])
    # sim_day
    ss = lines[5]
    rr = re.search('=(.*)\n',ss)
    dd = int(rr.group(1))
    # sim_month
    ss = lines[6]
    rr = re.search('=(.*)\n',ss)
    mm = int(rr.group(1))
    # sim_year
    ss = lines[7]
    rr = re.search('=(.*)\n',ss)
    yy = rr.group(1)
    yy = '20' + yy
    year=int(yy)
    # sim_hora
    ss = lines[8]
    rr = re.search('=(.*)\n',ss)
    hh = int(rr.group(1))
    # sim_minute
    ss = lines[9]
    rr = re.search('=(.*)\n',ss)
    mmin = int(rr.group(1))

    # set the date when the simulation started
    iStartDay = date(year,mm,dd).toordinal()
    iStartHour = hh
    iStartMinute = mmin

    # set time steps of interest (hours by default -- Python counting starts from 0).
    # It may be a single number e.g. [146] or a list of numbers e.g. np.arange(0,15)
    # outputs can be every 6h, for instance, by changing the steps in np.arange to 6,
    # for instance.
    time_line = np.arange(0,sim_length,1)
    real_time = time_line/24. + (iStartHour+1.)/24. + iStartDay
    
    return real_time
    
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
           
def main():   
    # filter out beached, dispersed, sunk and unreleased parcels
    # lons_f = longitude coordinates of each oil parcel
    # lats_f = latitude coordinates of each oil parcel    
    # surface_volumes = total volumes (evaporative and non evaporative) for each parcel found on the surface
    lons_f,lats_f,surface_volumes=get_surface_parcels(fname,time_index)

    # create dataframe with the location and volume of surface parcels
    df_gpan = pd.DataFrame({"lat":lats_f,
                   "lon":lons_f,
                   "vol":surface_volumes})
                   
    # create geodataframe from pandas dataframe                 
    # parcels_gdf = gpd.GeoDataFrame(df_gpan, 
    #             geometry=gpd.points_from_xy(df_gpan.lon, df_gpan.lat),
    #             crs=crs) 
    parcels_gdf = gpd.GeoDataFrame(df_gpan, 
                geometry=gpd.points_from_xy(df_gpan.lon, df_gpan.lat),crs=crs).set_crs('EPSG:4326',allow_override=True)
    parcels_gdf = parcels_gdf.drop(columns=['lon', 'lat']) # remove useless variables
                   
    # load satellite detection shapefile
    observation_df = gpd.read_file(observation_shp)
    # observation_gdf = gpd.GeoDataFrame(observation_df[['geometry']],crs=crs)

    observation_gdf = gpd.GeoDataFrame(observation_df[['geometry']]).set_crs('EPSG:4326',allow_override=True)

    
    # generate comparison grid
    # comparison grid is a common grid covering observed and modelled spills
    # the grid resolution has been set, for now, at 150m 
    obs_bounds = observation_gdf.bounds # first, getting observation bounds


    lonmin=np.min(lons_f)
    latmin=np.min(lats_f)
    lonmax=np.max(lons_f)
    latmax=np.max(lats_f)
       
    if obs_bounds.minx[0] < lonmin:
        lonmin=obs_bounds.minx[0]
    if obs_bounds.maxx[0] > lonmax:
        lonmax=obs_bounds.maxx[0]
    if obs_bounds.miny[0] < latmin:
        latmin=obs_bounds.miny[0]
    if obs_bounds.maxy[0] > latmax:
        latmax=obs_bounds.maxy[0]
    
    # create grid polygon consisting of multiple polygons describing squared grid cells    
    output_frame = make_poly_grid(lonmin,lonmax,latmin,latmax,grid_resolution,crs) 
    # (AUGUSTO SEPP ----- MODIFIED)
    output_frame.set_crs('EPSG:4326',allow_override=True)
    # create numpy-compatible grid 
    # by first getting the central coordinates of each grid cell
    #grid_centroid = np.asarray(output_frame['geometry'].centroid) 
    grid_centroid = output_frame['geometry'].to_crs(crs).centroid.to_crs(output_frame.crs)
    # then running through grid points (indices) and creating an output matrix "event_set"
    # containing:
    # column 0: cell index
    # column 1: longitude coordinate of grid cell centroid
    # column 2: latitude coordinate of grid cell centroid  
    # column 3: simulated oil presence/absence
    # column 4: observed oil presence/absence
    # column 5: intersection between simulated and observed oil
    counter=0
    event_set=np.zeros((len(grid_centroid),6))*np.nan
    
    for ii in grid_centroid:
        event_set[counter,0]=counter
        event_set[counter,1]=ii.coords[0][0]
        event_set[counter,2]=ii.coords[0][1]   
        counter=counter+1 
    # and save grid polygon
    # output_frame.to_file('grid.shp')    
    
    ##################################################################
    # group your parcels into visualization grid
    gridded_parcels = gpd.sjoin(parcels_gdf,output_frame, how='left', op='within').set_crs('EPSG:4326',allow_override=True)
    # aggregate volumes to grid cells with dissolve 
    gridded_parcels['cell_total_volume']=1
    cell_total_volume = gridded_parcels.dissolve(by="index_right", aggfunc="count")   
    # place them into "cell" array with their associated coordinates   
    output_frame.loc[cell_total_volume.index, 'cell_total_volume'] = cell_total_volume.cell_total_volume.values    
    modelled_spill = output_frame[output_frame['cell_total_volume'] > 0]
    
    for ii in range(0,len(modelled_spill)):
        counter = modelled_spill.iloc[ii].name
        event_set[counter,3]=modelled_spill.iloc[ii].cell_total_volume/modelled_spill.iloc[ii].cell_total_volume
    ##################################################################
    # fit satellite observation (spill shape) into visualization grid
    gridded_observation = gpd.sjoin(left_df=output_frame, right_df=observation_gdf[['geometry']], how='inner').set_crs('EPSG:4326',allow_override=True)
    for ii in range(0,len(gridded_observation)):
        counter = gridded_observation.iloc[ii].name
        event_set[counter,4]=1      
    # gridded_observation.to_file('observation.shp')   
    ##################################################################
    # find areas where model and observations coincide on oil detection
    model_and_obs= gpd.sjoin(left_df=output_frame[output_frame['cell_total_volume'] > 0], right_df=observation_gdf[['geometry']], how='inner')
    for ii in range(0,len(model_and_obs)):
        counter = model_and_obs.iloc[ii].name
        event_set[counter,5]=1
    ##################################################################

    # compute your Fractional Skill Score 
    X,Y=np.meshgrid(np.unique(event_set[:,1]),np.unique(event_set[:,2]))
    interp_model=NearestNDInterpolator(list(zip(event_set[:,1],event_set[:,2])),event_set[:,3])
    interp_observation=NearestNDInterpolator(list(zip(event_set[:,1],event_set[:,2])),event_set[:,4])
    interp_intersection=NearestNDInterpolator(list(zip(event_set[:,1],event_set[:,2])),event_set[:,5])
    
    array_model = interp_model(X, Y)
    array_observation = interp_observation(X, Y)   
    array_intersection = interp_intersection(X, Y)    
    
    array_model[np.isnan(array_model)] = 0
    array_observation[np.isnan(array_observation)] = 0

    array_union = array_model+2*array_observation
    array_union[array_union==0]=np.nan
    
    # set search window sizes  (in number of pixes)
    horizontal_scales = range(1,150,2) # horizontal_scales = range(1,11,2) per aumentare la scala
    fss_output=np.zeros((len(horizontal_scales),2))
    
    cc = 0
    for hh in horizontal_scales:         
        fss_=pysteps.verification.spatialscores.fss(array_model, array_observation, 1, hh)
        fss_output[cc,0]=hh
        fss_output[cc,1]=fss_
        cc=cc+1

    

    plt.figure()
    cMap = c.ListedColormap(['y','b','m'])
    m = Basemap(llcrnrlon=lonmin,llcrnrlat=latmin,\
                urcrnrlon=lonmax,urcrnrlat=latmax,\
                rsphere=(6378137.00,6356752.3142),\
                resolution='i',projection='merc',\
                lat_0=(latmax + latmin)/2.,\
                lon_0=(lonmax + lonmin)/2.,epsg=4326)
    x_map,y_map=m(X,Y)
    m.pcolor(x_map,y_map,array_union,cmap=cMap)
    m.drawcoastlines()
    m.fillcontinents(alpha=1,zorder=3)
    m.drawmeridians(np.arange(lonmin,lonmax,(lonmax-lonmin)/4.), labels=[0,0,0,1],color='white',linewidth=0.03, fontsize = 5) # draw parallels
    m.drawparallels(np.arange(latmin,latmax,(latmax-latmin)/4.),labels=[1,0,0,0],color='white',linewidth=0.03, fontsize = 5) # draw meridians
    plt.savefig(output_folder + '/overlay_' + xp_identifier + '.png',dpi=600,bbox_inches='tight')

    np.savetxt(output_folder + '/fss_' + xp_identifier + '.txt',fss_output)
    np.savetxt(output_folder + '/event_set_' + xp_identifier + '.txt',event_set)
    
    plt.figure()    
    plt.plot(fss_output[:,0]*verif_grid_resolution,fss_output[:,1],'.-')
    plt.savefig(output_folder + '/agg_fss_' + xp_identifier + '.png',dpi=600,bbox_inches='tight')

    plt.close('all')


# SCRIPT FOR FRACTIONAL SKILL SCORE CALCULATION APPLIED TO OIL SPILL FORECASTS
# start by informing your simulation folder (where you placed your MDK-II outputs)
# you may consider, later, to include a loop through a "mother" folder...
simulation_folder = args.forecastfolder 

# set the folder where you placed your observations
observation_shp = args.observationfolder #'/Users/asepp/work/imagine/observations/20210824-0343-SYR-PL-B-01-S1/20210824-0343-SYR-PL-B-01-S1.shp'

# set folder where script results will be placed
'''
if os.path.isdir(args.outputfolder):
    output_folder = args.outputfolder
else:
    os.mkdir(args.outputfolder)
    output_folder = args.outputfolder
'''

output_folder = args.outputfolder
#print('>>>>>>>>>>>>>>>>>>>>>',output_folder, '+')

# set the spatial resolution for the verification grid (in km)
verif_grid_resolution = .15

# set mapping projection
crs = "+proj=merc +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"

grid_resolution = np.longdouble(verif_grid_resolution)/110.
area = (grid_resolution*110.)**2
fillval = 0
fname = simulation_folder + '/spill_properties.nc'
nc_read = netCDF4.Dataset(fname, 'r')

simulation_date = get_mdksim_date(simulation_folder)
observation_date = get_obs_date(observation_shp)
time_index=find_nearest(simulation_date, observation_date)
sim_identifier=get_mdksim_id(simulation_folder)
xp_identifier = sim_identifier + '_' + '%02d' % (time_index+1) + 'h_' 


# for time_index in output_timesteps:
print(':::::::::: TIME_INDEX :: ' + str(time_index))
print(xp_identifier)


main()
