#!/usr/bin/python
# Modified by Marco Seracini - Released 08 Nov 2023

# Application:
# to extract subsets from GEBCO 30'' files and transform it into MEDSLIK-II
# compatible .bath files

import netCDF4
import numpy as np
from scipy.interpolate import RectBivariateSpline
import os
import sys
import json


# Functions
# open and crop GEBCO
def load_interp_gebco(filename, x_mod, y_mod):

	fGEBCO = netCDF4.Dataset(filename)
	
	cols, rows = fGEBCO.variables["dimension"]
	iDs = fGEBCO.variables["spacing"][0]

	# Coordinates for GEBCO corners - LAT,LON
	iNECorner=np.array([89.+(59./60)+45./(60*60), -179.-(59./60)-45./(60*60)])
	iSWCorner=np.array([-90.,180.])

	iLatitude=np.arange(iNECorner[0],iSWCorner[0],-iDs)
	iLongitude=np.arange(iNECorner[1],iSWCorner[1],iDs)

	iLatitudeMin= np.min(y_mod)
	iLatitudeMax= np.max(y_mod)
	iLongitudeMin= np.min(x_mod)
	iLongitudeMax= np.max(x_mod)

    
    # Crop to area of interest
	iLonIndex=np.argwhere((iLongitude>=iLongitudeMin) & (iLongitude<=iLongitudeMax))
	iLatIndex=np.argwhere((iLatitude>=iLatitudeMin) & (iLatitude<=iLatitudeMax))

	Offset=iLatIndex[0]*cols+iLonIndex[0]

	rowsNumber=np.count_nonzero(iLatIndex)
	colsNumber=np.count_nonzero(iLonIndex)

	x_crop, y_crop = np.meshgrid(iLongitude[iLonIndex],iLatitude[iLatIndex])

	z_crop = np.empty(shape=(np.max(iLatIndex)-np.min(iLatIndex)+1,np.max(iLonIndex)-np.min(iLonIndex)+1))



	for i in range(rowsNumber):
		z_crop[i,:] = fGEBCO.variables["z"][int(Offset+i*cols):int(Offset+i*cols+colsNumber)]


	# Generate interpolator
	y_crop = np.flipud(y_crop)
	x_crop = np.flipud(x_crop)
	z_crop = np.flipud(z_crop)

	z_int = RectBivariateSpline(y_crop[:,1],x_crop[1,:],z_crop)

	# Interpolate
	z_proc = z_int(y_mod,x_mod)

     # Fix orientation 

	z_proc = np.flipud(z_proc)
    
	# Convert bathymetry to MDK-II
	mdk_z=[]
	mdk_x=[]
	mdk_y=[]


	r_,c_ = z_proc.shape
	for i in range(0,r_,1):
		for j in range(0,c_,1):
			mdk_z.append(z_proc[i,j])
			mdk_x.append(x_mod[j])
			mdk_y.append(y_mod[i])

	mdk_z = np.array(mdk_z)
	land_mask = np.where(mdk_z >= 0)
	mdk_z[land_mask]=-9999
	mdk_z = -1.*(mdk_z)
	return mdk_x, mdk_y, mdk_z, c_, r_

# extract grid corners from ocean fields
def oce_grid(filename):
	fOCM = netCDF4.Dataset(filename)
	x_mod = fOCM.variables['nav_lon'][:]
	y_mod = fOCM.variables['nav_lat'][:]
	return x_mod, y_mod

def open_mercator(sFileName):
	# load file and geo coordinates
	fM=netCDF4.Dataset(sFileName,"r")
	y=fM.variables["nav_lat"][:]
	x=fM.variables["nav_lon"][:]
	iT=fM.variables["vomecrtx"][:] # water temperature in deg C
	return x,y,iT

################################################################################
# USER INPUTS
################################################################################

# - gebco_filename:set path for your original gebco netCDF file
# - output_dir: define where outputs will be placed (MEDSLIK-adapted nc files)
# - oce_dir: the bathymetry fields will be interpolated to your hydrodynamic grid
#   therefore, let us know where you've placed your current files
gebco_filename = None
output_dir= None
oce_dir = None

# parse configuration file, if any
if len(sys.argv) > 1:

        configFile = sys.argv[1]
        print("Parsing configuration file %s" % configFile)
        
        with open(configFile) as f:
                try:
                        config = json.load(f)
                        gebco_filename = config["preproc_gebco"]["gebco_filename"]
                        output_dir = config["preproc_gebco"]["output_dir"]
                        oce_dir = config["preproc_gebco"]["oce_dir"]
                except KeyError:
                        print("Missing entry in configuration file! Abort!")
                        sys.exit()
                        
else:
        print("You need to provide the config file!")
        sys.exit()



################################################################################
# USER INPUTS - OVER!
################################################################################
# From here onwards, the script should do everything pretty much automatic
# bugs/errors are expected and, in case you do find one,
# feel free to send us comments/corrections.





# open an ocean forecast file
###oce_filename = (oce_dir + "/" + os.listdir(oce_dir)[2])
oce_filename = (oce_dir + "/MDK_ocean_170423_T.nc")
x_mod, y_mod = oce_grid(oce_filename)


# load GEBCO file and extract bathymetry
mdk_x, mdk_y, mdk_z, c_, r_ = load_interp_gebco(gebco_filename, x_mod, y_mod)

# write .bath file
BathFile=open(output_dir + "vnzl_.bath", "w")
BathFile.write("MEDSLIK-II compatible bathymetry file. Degraded resolution based on GEBCO 30''\n")
BathFile.write("%-7.4f %-7.4f %-7.4f %-7.4f \n" % (np.min(mdk_x),np.max(mdk_x),np.min(mdk_y),np.max(mdk_y)))
BathFile.write("%d %d \n" % (c_,r_))
np.savetxt(BathFile,mdk_z,fmt="%04.0f")
BathFile.close()

