#!/usr/bin/python

# Application: Paria testcase
# to transform CMEMS-Global (MERCATOR) daily current files into MEDSLIK-II
# compatible netCDF files

from netCDF4 import Dataset
import numpy as np
from datetime import  *
import pdb

import os
import sys
import json

# some functions which will be used by the main program
def open_mercator_subset(sFileName,grid_corners):

	print('=== sFilename:' + sFileName)

	# load file and geo coordinates
	fM=Dataset(sFileName,"r")
	iLatitude=fM.variables["latitude"][:]
	iLongitude=fM.variables["longitude"][:]

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
	iU=fM.variables["uo"][:,:,ywindow,xwindow] # u component of water velocity in m/s
	iV=fM.variables["vo"][:,:,ywindow,xwindow] # v component
	iT=fM.variables["thetao"][:,:,ywindow,xwindow] # water temperature in deg C
	iTime=fM.variables["time"][:]  # hours since 1950-01-01 00:00:00
	iTime=date(1950,1,1).toordinal() + iTime/24 # time converted to standard pythonic time
	iDepth=fM.variables["depth"][:]		 # in m
	return x,y,iU,iV,iT,iTime,iDepth

def interp_z(iU):
	iShapeO=np.shape(iU)
	# daily resolution and four depths - 0,10,30,120m
	iOutputUd=np.zeros((iShapeO[0],4,iShapeO[2],iShapeO[3]))	# u(time,depth,y,x)
	# surface
	iOutputUd[:,0,:,:]=iU[:,0,:,:]
	# 10 m
	iOutputUd[:,1,:,:]=iU[:,7,:,:] + (iU[:,8,:,:]-iU[:,7,:,:])*((10. - iDepth[7])/(iDepth[8] - iDepth[7]))
	# 30 m
	iOutputUd[:,2,:,:]=iU[:,14,:,:] + (iU[:,15,:,:]-iU[:,14,:,:])*((30. - iDepth[14])/(iDepth[15] - iDepth[14]))
	# 120 m
	iOutputUd[:,3,:,:]=iU[:,22,:,:] + (iU[:,23,:,:]-iU[:,22,:,:])*((120. - iDepth[22])/(iDepth[23] - iDepth[22]))
	return np.squeeze(iOutputUd)

def interp_t(iOutputUd, iHRTimeLine, iLRTimeLine):
		# prepare storing matrixes
		iShapeO = np.shape(iOutputUd)
		iOutputU=np.zeros((24,4,iShapeO[2],iShapeO[3]))-32767.0	# u(time,depth,y,x)
		#interpolating
		for ii in range(0,iShapeO[2]):
			for jj in range(0,iShapeO[3]):
				for dd in range(0,4):
					if iOutputUd[0,dd,ii,jj]>-20:
						iOutputU[:,dd,ii,jj]=np.interp(iHRTimeLine,iLRTimeLine,iOutputUd[:,dd,ii,jj])
						#Extrapolating data
						if (dd!=0 and any(iOutputU[:,dd,ii,jj]<-2) and all(iOutputU[:,dd-1,ii,jj]>-2)):
							iOutputU[:,dd,ii,jj]=iOutputU[:,dd-1,ii,jj]
		return iOutputU

################################################################################
# USER INPUTS
################################################################################

# the following parameters will be read from a configuration file:
# - iInputFolder: where your CMEMS-GLO files are placed
# - iStorageDirectory: where outputs will be placed (MEDSLIK-adapted nc files)
# - iStartDate and iEndDate: start and end data of the pre-processing
# - grid_corners: the geographic limits of your preprocessing (you may have a large CMEMS file)
iInputFolder = None
iStorageDirectory = None
iStartDate = None
iEndDate = None
iOutputRefHours = None
iInputRefHours = None
xE = xW = yN = yS = None
grid_corners = None

# parse configuration file, if any
if len(sys.argv) > 1:

	configFile = sys.argv[1]
	print("Parsing configuration file %s" % configFile)
	
	with open(configFile) as f:
		try:
			config = json.load(f)

			# read folders
			iInputFolder = config["preproc_currents"]["input_folder"]
			iStorageDirectory = config["preproc_currents"]["output_folder"]

			# read dates
			start_date = config["preproc_currents"]["start_date"]
			end_date = config["preproc_currents"]["end_date"]
			iStartDate=date(start_date["year"], start_date["month"], start_date["day"]).toordinal()
			iEndDate=date(end_date["year"], end_date["month"], end_date["day"]).toordinal()

			# read hours
			orh = config["preproc_currents"]["output_ref_hours"]
			irh = config["preproc_currents"]["input_ref_hours"]
			iOutputRefHours=np.arange(orh["start"], orh["end"])
			iInputRefHours=irh

			# read coords
			gc =  config["preproc_currents"]["grid_corners"]
			grid_corners = [ gc["xW"], gc["xE"], gc["yS"], gc["yN"] ]

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
# bugs/errors are expected and in case you unfortunate enough to find out one,
# feel free to send us comments/corrections.

for i in range(iStartDate,iEndDate,1):

	# Keep track of what you are doing
	iDateY=date.fromordinal(int(i-1))
	iDate=date.fromordinal(int(i))
	iDateT=date.fromordinal(int(i+1))
	print('Generating ' +  'HOPS' + iDate.strftime('%y') + iDate.strftime('%m') + iDate.strftime('%d') + '_.nc')

	# Open MERCATOR files and interpolate in z for preset dephts 0,10,30,120m
	# day before
	if i == iStartDate:
		#mercatorpsy4v3r1_gl12_mean_YYYYMMDD_RYYYYMMDD.nc MERCATOR production filename
		fn = 'mercatorpsy4v3r1_gl12_mean_20' + iDateY.strftime('%y') + iDateY.strftime('%m') + iDateY.strftime('%d') + '_R20170503.nc'
		if not os.path.exists(iInputFolder + fn):
			print("File not found %s" % fn)
			print("Input folder %s contains:" % iInputFolder)
			for f in os.listdir(iInputFolder):
				print(" - %s" % f)
			fn = input("Please select one of the file listed above to continue: ")
		
		iLongitude,iLatitude,iU_Y,iV_Y,iT_Y,iTime_Y,iDepth = open_mercator_subset(iInputFolder + fn,grid_corners)
		#iLongitude,iLatitude,iU_Y,iV_Y,iT_Y,iTime_Y,iDepth = open_mercator_subset(iInputFolder + 'mercator_paria2_20' + iDateY.strftime('%y') + iDateY.strftime('%m') + iDateY.strftime('%d') + '.nc',grid_corners)
		iU_Y = (interp_z(iU_Y))
		iV_Y = (interp_z(iV_Y))
		iT_Y = (interp_z(iT_Y))
		# day of interest
		fn = 'mercatorpsy4v3r1_gl12_mean_20' + iDate.strftime('%y') + iDate.strftime('%m') + iDate.strftime('%d') + '_R20170503.nc'
		if not os.path.exists(iInputFolder + fn):
			print("File not found %s" % fn)
			print("Input folder %s contains:" % iInputFolder)
			for f in os.listdir(iInputFolder):
				print(" - %s" % f)
			fn = input("Please select one of the file listed above to continue: ")

		iLongitude,iLatitude,iU_P,iV_P,iT_P,iTime,iDepth = open_mercator_subset(iInputFolder + fn,grid_corners)
		#iLongitude,iLatitude,iU_P,iV_P,iT_P,iTime,iDepth = open_mercator_subset(iInputFolder + 'mercator_paria2_20' + iDate.strftime('%y') + iDate.strftime('%m') + iDate.strftime('%d') + '.nc',grid_corners)
		iU_P = (interp_z(iU_P))
		iV_P = (interp_z(iV_P))
		iT_P = (interp_z(iT_P))
		#  day after
		fn = 'mercatorpsy4v3r1_gl12_mean_20' + iDateT.strftime('%y') + iDateT.strftime('%m') + iDateT.strftime('%d') + '_R20170503.nc'
		if not os.path.exists(iInputFolder + fn):
			print("File not found %s" % fn)
			print("Input folder %s contains:" % iInputFolder)
			for f in os.listdir(iInputFolder):
				print(" - %s" % f)
			fn = input("Please select one of the file listed above to continue: ")
		
		iLongitude,iLatitude,iU_T,iV_T,iT_T,iTime_T,iDepth = open_mercator_subset(iInputFolder + fn,grid_corners)
		#iLongitude,iLatitude,iU_T,iV_T,iT_T,iTime_T,iDepth = open_mercator_subset(iInputFolder + 'mercator_paria2_20' + iDateT.strftime('%y') + iDateT.strftime('%m') + iDateT.strftime('%d') + '.nc',grid_corners)
		iU_T = (interp_z(iU_T))
		iV_T = (interp_z(iV_T))
		iT_T = (interp_z(iT_T))

	else:
		# Present becomes yesterday -- oh life...
		iU_Y = iU_P
		iV_Y = iV_P
		iT_Y = iT_P
		# Tomorrow becomes today -- ...time flies.
		iU_P = iU_T
		iV_P = iV_T
		iT_P = iT_T

		fn = 'mercatorpsy4v3r1_gl12_mean_20' + iDateT.strftime('%y') + iDateT.strftime('%m') + iDateT.strftime('%d') + '_R20170503.nc'
		if not os.path.exists(iInputFolder + fn):
			print("File not found %s" % fn)
			print("Input folder %s contains:" % iInputFolder)
			for f in os.listdir(iInputFolder):
				print(" - %s" % f)
			fn = input("Please select one of the file listed above to continue: ")
		
		iLongitude,iLatitude,iU_T,iV_T,iT_T,iTime_T,iDepth = open_mercator_subset(iInputFolder + fn,grid_corners)
		#iLongitude,iLatitude,iU_T,iV_T,iT_T,iTime_T,iDepth = open_mercator_subset(iInputFolder + 'mercator_paria2_20' + iDateT.strftime('%y') + iDateT.strftime('%m') + iDateT.strftime('%d') + '.nc',grid_corners)
		iU_T = (interp_z(iU_T))
		iV_T = (interp_z(iV_T))
		iT_T = (interp_z(iT_T))

	# Put the three fields (t - 1, t, t + 1) into a single matrix to interpolate
	iOutputUd = np.stack((iU_Y,iU_P,iU_T),axis=0)
	iOutputVd = np.stack((iV_Y,iV_P,iV_T),axis=0)
	iOutputTd = np.stack((iT_Y,iT_P,iT_T),axis=0)

	# Temporal interpolation:
	iOutputU = interp_t(iOutputUd, iOutputRefHours, iInputRefHours)
	iOutputV = interp_t(iOutputVd, iOutputRefHours, iInputRefHours)
	iOutputT = interp_t(iOutputTd, iOutputRefHours, iInputRefHours)

	# Masking
	iOutputU=np.where(iOutputU < -3,9999.,iOutputU)
	iOutputV=np.where(iOutputV < -3,9999.,iOutputV)
	iOutputT=np.where(iOutputT < -4,9999.,iOutputT)

	iOutputU=np.where(iOutputU > 3,9999.,iOutputU)
	iOutputV=np.where(iOutputV > 3,9999.,iOutputV)
	iOutputT=np.where(iOutputT > 40,9999.,iOutputT)


	#Generates NC files
	# U component
	
	sConvertedFilename=iStorageDirectory + 'MDK_ocean_' + iDate.strftime('%y') + iDate.strftime('%m') + iDate.strftime('%d') + '_U.nc'
	
	print('=== sConvertedFilename: ' + sConvertedFilename)
	
	f = Dataset(sConvertedFilename, "w", "NETCDF3_CLASSIC")
	f.history = 'MERCATOR Forecast - CMEMS -- Adapted to Medslik-II'
	f.createDimension('time_counter',None)
	f.createDimension('y', len(iLatitude))
	f.createDimension('x', len(iLongitude))
	f.createDimension('depthu',4)
	time = f.createVariable('time_counter', 'd', ('time_counter',))
	# Geo coordinates
	yCoord= f.createVariable('nav_lat', 'd', ('y'))
	xCoord= f.createVariable('nav_lon', 'd', ('x'))
	# Zonal velocity component
	uCoord= f.createVariable('vozocrtx', 'd',('time_counter','depthu','y','x'))
	time[:] = i + np.arange(1./24,25./24,1./24)
	#Filling in variables
	yCoord[:]=iLatitude
	xCoord[:]=iLongitude
	uCoord[:]=iOutputU
	#Adding units
	yCoord.units = 'degrees_north'
	xCoord.units = 'degrees_east'
	time.units = 'julian days'
	uCoord.units='m*s-1'
	uCoord.missing_value=9999.
	f.close()
	
	# V component
	sConvertedFilename=iStorageDirectory + 'MDK_ocean_' + iDate.strftime('%y') + iDate.strftime('%m') + iDate.strftime('%d') + '_V.nc'
	
	print('=== sConvertedFilename: ' + sConvertedFilename)
	
	f = Dataset(sConvertedFilename, "w", "NETCDF3_CLASSIC")
	f.history = 'MERCATOR Forecast - CMEMS -- Adapted to Medslik-II'
	f.createDimension('time_counter',None)
	f.createDimension('y', len(iLatitude))
	f.createDimension('x', len(iLongitude))
	f.createDimension('depthv',4)
	time = f.createVariable('time_counter', 'd', ('time_counter',))
	# Geo coordinates
	yCoord= f.createVariable('nav_lat', 'd', ('y'))
	xCoord= f.createVariable('nav_lon', 'd', ('x'))
	# Zonal velocity component
	vCoord= f.createVariable('vomecrty', 'd',('time_counter','depthv','y','x'))
	time[:] = i + np.arange(1./24,25./24,1./24)
	#Filling in variables
	yCoord[:]=iLatitude
	xCoord[:]=iLongitude
	vCoord[:]=iOutputV
	#Adding units
	yCoord.units = 'degrees_north'
	xCoord.units = 'degrees_east'
	time.units = 'julian days'
	vCoord.units='m*s-1'
	vCoord.missing_value=9999.
	f.close()

	# Temperature
	sConvertedFilename=iStorageDirectory + 'MDK_ocean_' + iDate.strftime('%y') + iDate.strftime('%m') + iDate.strftime('%d') + '_T.nc'
	
	print('=== sConvertedFilename: ' + sConvertedFilename)
	
	f = Dataset(sConvertedFilename, "w", "NETCDF3_CLASSIC")
	f.history = 'MERCATOR Forecast - CMEMS -- Adapted to Medslik-II v2.0'
	f.createDimension('time_counter',None)
	f.createDimension('y', len(iLatitude))
	f.createDimension('x', len(iLongitude))
	f.createDimension('deptht',4)
	time = f.createVariable('time_counter', 'd', ('time_counter',))
	# Geo coordinates
	yCoord= f.createVariable('nav_lat', 'd', ('y'))
	xCoord= f.createVariable('nav_lon', 'd', ('x'))
	# Zonal velocity component
	tCoord= f.createVariable('votemper', 'd',('time_counter','deptht','y','x'))
	time[:] = i + np.arange(1./24,25./24,1./24)
	#Filling in variables
	yCoord[:]=iLatitude
	xCoord[:]=iLongitude
	tCoord[:]=iOutputT
	#Adding units
	yCoord.units = 'degrees_north'
	xCoord.units = 'degrees_east'
	time.units = 'julian days'
	tCoord.units='degree Celsius'
	tCoord.missing_value=9999.
	f.close()
