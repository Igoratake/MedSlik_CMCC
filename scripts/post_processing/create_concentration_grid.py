import numpy as np
import pandas as pd
import xarray as xr

from glob import glob

def create_concentration_dataset(resolution = 150, multiple_slick=False,filepath=''):

    '''
    This script creates a gridded netcdf from point output netcdf from Medslik-II

    The file is usually called as spill_properties.nc

    In order to create the netcdf file on a grid, some premises are made:

    1 - Bounding box: a bounding box is defined from all time steps in the spill properties netcdf, therefore collecting 

        maximum and mininum latitudes and longitudes in all particles in the given output.

    2 - Resolution: Resolution can be modified, standard is 150 meters to respect oil grid in medslik.for.

        IMPORTANT! This code creates area values for each grid cell considering this value, so do not use degrees ads unit
    
    3 - Multiple slick flag allows the user to construct the netcdf from multiple slicks, combining all outputs on a single
        point netcdf and the gridded as well.
    
    '''

    #Opens all separate slicks if needed and group on a single netcdf
    if multiple_slick == True:
        files = glob(filepath + '*/spill_properties.nc')
        ds = xr.open_mfdataset(files,combine='nested',concat_dim='parcel_id')
    #Get the single file for regular simulations
    else:
        ds = xr.open_mfdataset(filepath + 'spill_properties.nc')

    #get the oil density in kg/m3 from the netcdf
    oil_density = ds.non_evaporative_volume.oil_density
        
    # Calculate latitude and longitude resolution fro the given output
    lat_min = ds.latitude.values.min()
    lon_min = ds.longitude.values.min()
    lat_max = ds.latitude.values.max()
    lon_max = ds.longitude.values.max()

    # Obtain the resolution on degrees for latitude and longitude given the mean latitude value
    lat_mean = (lat_min + lat_max) / 2
    
    lat_resolution_degree = resolution / (111320)
    lon_resolution_degree = resolution / (111320 * np.cos(np.pi*(lat_mean)/180))

    # Create latitude and longitude arrays
    lats_array = np.arange(lat_min, lat_max, lat_resolution_degree)
    lons_array = np.arange(lon_min, lon_max, lon_resolution_degree)

    # Create concentration grid - Later to be the base for the netcdf
    concentration_grid = np.zeros((len(ds.time), len(lats_array), len(lons_array)))

    #Loop over time and create the concentrations for each grid element
    for t in range(0,len(ds.time)):

        print(f'timestep {t}')

        # Select a single timesetp
        rec = ds.isel(time=t)

        #Obtain the coordinates for each point and its repective volume converted to tonnes
        lats = rec.latitude.values
        lons = rec.longitude.values

        #Get the total amount of volume on water surface
        volumes = (rec.non_evaporative_volume.values + rec.evaporative_volume.values)
        #obtain the mass of each particle by multiplying by the oil density
        mass = volumes * oil_density  

        # Create a Pandas DataFrame with information for each particle
        df = pd.DataFrame({'latitude': lats, 'longitude': lons, 'mass': mass})

        # Keep only paticles that have a mass bigger than 0
        particle_data = df[df.mass > 0]

        # Assign particles to grid cells -  Core of the solution with np.digitize
        particle_data['lat_bin'] = np.digitize(particle_data['latitude'].values, lats_array)
        particle_data['lon_bin'] = np.digitize(particle_data['longitude'].values, lons_array)

        # Aggregate masses within each grid cell and normalize by grid area
        aggregated = particle_data.groupby(['lat_bin', 'lon_bin']).agg({'mass': 'sum'})

        #Reset index
        aggregated = aggregated.reset_index()

        #Transform from mass on cell to concentration, dividing by the resolution area
        aggregated['concentration'] = aggregated['mass']/(lat_resolution_degree*lon_resolution_degree)

        #insert each grid concentration on the matrix
        for x in range(0,len(aggregated)):
            concentration_grid[t,int(aggregated.iloc[x].lat_bin)-1, int(aggregated.iloc[x].lon_bin)-1] = aggregated.iloc[x].concentration

    # Create xarray Dataset from the array
    concentration_dataset = xr.Dataset({'concentration': (['time', 'lat', 'lon'], concentration_grid)},
                                        coords={'time': ds.time.values, 'lat': lats_array, 'lon': lons_array})
    
    #Saves merged netcdf if needed
    if multiple_slick == True:
        ds.to_netcdf(filepath + 'spill_properties_merged.nc')

    #Saves concentration grid netcdf
    concentration_dataset.to_netcdf(filepath + 'spill_properties_concentration.nc')
