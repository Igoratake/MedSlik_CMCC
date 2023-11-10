#!/usr/bin/python

# Application:
# to transform ERA-Interim reanalysis daily files into MEDSLIK-II
# compatible netCDF files

from netCDF4 import Dataset
import numpy as np
from datetime import  *


def open_era_subset(sFileName,grid_corners):
        fM=Dataset(sFileName,"r",format="NETCDF4")
        iLatitude=fM.variables["latitude"][:]
        iLongitude=fM.variables["longitude"][:]
        
        iLongitude[iLongitude>180]=iLongitude[iLongitude>180]-360.

        # cropping before extraction
        xmin = grid_corners[0]
        xmax = grid_corners[1]
        ymin = grid_corners[2]
        ymax = grid_corners[3]

        xwindow = np.logical_and(iLongitude>=xmin, iLongitude<=xmax)
        ywindow = np.logical_and(iLatitude>=ymin, iLatitude<=ymax)

        x = iLongitude[xwindow]
        y = iLatitude[ywindow]

        # Extraction and conversion of variables
        iU=fM.variables["u10"][:,ywindow,xwindow]
        iV=fM.variables["v10"][:,ywindow,xwindow]
        iTime=fM.variables["time"][:]  # hours since 1900-01-01 00:00:00
        iTime=date(1900,01,1).toordinal() + iTime/24. # time converted to standard pythonic time
        return x,y,iU,iV,iTime

def interp_t(iOutputUd, iOutputRefHours, iInputRefHours):
	# prepare storing matrixes
	iShapeO = np.shape(iOutputUd)
	iOutputU = np.zeros((len(iOutputRefHours),iShapeO[1],iShapeO[2]))

	#interpolating
	for ii in range(0,iShapeO[1]):
		for jj in range(0,iShapeO[2]):
				iOutputU[:,ii,jj]=np.interp(iOutputRefHours,iInputRefHours,iOutputUd[:,ii,jj])

	return iOutputU


################################################################################
# USER INPUTS
################################################################################
# Define where your CMEMS-GLO files are placed
iInputFolder = '/Users/asepp/work/data/lebanon/era/'

# Define where outputs will be placed (MEDSLIK-adapted nc files)
iStorageDirectory= "/Users/asepp/work/data/lebanon/met_files/"

# Define start and end data of the pre-processing
# It will depend on your input data, of course, and on your simulation setup
iStartDate=date(2006,7,11).toordinal()
iEndDate=date(2006,8,13).toordinal()

# Requirements for temporal interpolation
# reference hours
iOutputRefHours=np.arange(0,24) # basically, you're saying you want hourly outputs
iInputRefHours=[0,6,12,18,24] # and that your ERA files have 6h resolution

# Now, set the geographic limits of your preprocessing (you may have a large CMEMS file)
xE = 37. #easternmost coordinate
xW = 33. #westernmost coordinate
yN = 36.  #northernmost coordinate
yS = 32.   #southernmost coordinate
grid_corners = [xW,xE,yS,yN]

################################################################################
# USER INPUTS - OVER!
################################################################################
# From here onwards, the script should do everything pretty much automatic
# bugs/errors are expected and in case you unfortunate enough to find out one,
# feel free to send us comments/corrections.

for i in range(iStartDate,iEndDate+1):

	iDate1=date.fromordinal(i)
	iDate2=date.fromordinal(i+1)

	print('Generating ' + iDate1.strftime('%Y') + iDate1.strftime('%m') + iDate1.strftime('%d') + '.nc')
	
	# open era files 
	x1,y1,u1,v1,t1 = open_era_subset(iInputFolder + iDate1.strftime('%Y') + iDate1.strftime('%m') + iDate1.strftime('%d') + '.nc',grid_corners)
	x2,y2,u2,v2,t2 = open_era_subset(iInputFolder + iDate2.strftime('%Y') + iDate2.strftime('%m') + iDate2.strftime('%d') + '.nc',grid_corners)

	# Temporal interpolation
	iOutputUd = np.concatenate((u1,u2),axis=0)[0:5,:,:]
	iOutputVd = np.concatenate((v1,v2),axis=0)[0:5,:,:]

	iOutputU = interp_t(iOutputUd, iOutputRefHours, iInputRefHours)
	iOutputV = interp_t(iOutputVd, iOutputRefHours, iInputRefHours)

	# Masking
	iOutputU=np.where(iOutputU < -100,9999,iOutputU)
	iOutputV=np.where(iOutputV < -100,9999,iOutputV)

	#Generates NC file
	sConvertedFilename=iStorageDirectory + iDate1.strftime('%Y') + iDate1.strftime('%m') + iDate1.strftime('%d') + '.nc'
	f = Dataset(sConvertedFilename,"w","NETCDF3_CLASSIC")
	f.history = 'ERA Interim 10m winds -- Adapted to Medslik-II'
	f.createDimension('time',None)
	f.createDimension('lat', len(y1))
	f.createDimension('lon', len(x1))
	time = f.createVariable('time', 'd', ('time',))

	# Geo coordinates
	yCoord= f.createVariable('lat', 'd', ('lat',))
	xCoord= f.createVariable('lon', 'd', ('lon',))
	# Zonal velocity component
	uCoord= f.createVariable('U10M', 'd',('time','lat','lon'))
	vCoord= f.createVariable('V10M', 'd',('time','lat','lon'))
	time[:] = i + np.arange(1./24,24./24,1./24)
	#Filling the sausage
	yCoord[:]=y1
	xCoord[:]=x1
	#Filling the sausage chap. II
	uCoord[:]=iOutputU
	vCoord[:]=iOutputV
	#Adding units
	yCoord.units = 'degrees_north'
	xCoord.units = 'degrees_east'
	xCoord.lon_min = np.min(x1)
	yCoord.lat_min = np.min(y1)
	time.units = 'pythonic days'
	uCoord.units='m*s-1'
	uCoord.missing_value=9999
	vCoord.units='m*s-1'
	vCoord.missing_value=9999
	f.close()
