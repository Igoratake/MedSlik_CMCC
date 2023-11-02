;;;;;;;;;;;;;;;;; currents subgrid ;;;;;;;;;;;;;;;;

grid_lona = (floor(min_map_lon/MODEL_LOND)-2)*MODEL_LOND
grid_lata = (floor(min_map_lat/MODEL_LATD)-2)*MODEL_LATD
grid_lonb = ( ceil(max_map_lon/MODEL_LOND)+2)*MODEL_LOND
grid_latb = ( ceil(max_map_lat/MODEL_LATD)+2)*MODEL_LATD
grid_lonn = round((grid_lonb-grid_lona)/MODEL_LOND,3)+1
grid_latn = round((grid_latb-grid_lata)/MODEL_LATD,3)+1

if .not.isvar("VECTOR_LAT_STRIDE") then
	VECTOR_LAT_STRIDE = max((/1,(grid_latn)/10/))
end if
if .not.isvar("VECTOR_LON_STRIDE") then
	VECTOR_LON_STRIDE = max((/1,(grid_lonn)/15/))
end if

if True then ;VERBOSE then
	print("VECTOR_LON_STRIDE = "+VECTOR_LAT_STRIDE)
	print("VECTOR_LAT_STRIDE = "+VECTOR_LON_STRIDE)
	print("")
end if


if VERBOSE then 
	print("currents subgrid OK")
	print("latitude :  min = "+sprintf("%7.4f",grid_lata)+"  max = "+sprintf("%7.4f",grid_latb)+"  n = "+grid_latn)
	print("longitude:  min = "+sprintf("%7.4f",grid_lona)+"  max = "+sprintf("%7.4f",grid_lonb)+"  n = "+grid_lonn)
end if

ii = round((MODEL_LON-grid_lona)/MODEL_LOND,3)
jj = round((MODEL_LAT-grid_lata)/MODEL_LATD,3)
fg = (ii.ge.0).and.(ii.lt.grid_lonn).and.(jj.ge.0).and.(jj.lt.grid_latn)
grid_ind = ind(fg)
if ismissing(grid_ind(0)) then
	ngrid = 0
	print("no currents values inside subgrid")
	print("try to expand map boundaries")
	print("**** procedure aborted ****")
	exit
else
	ngrid = dimsizes(grid_ind)
	grid_i = ii(grid_ind)
	grid_j = jj(grid_ind)
	if VERBOSE then 
		print("number of currents values inside = "+ngrid)
		print("")
	end if
end if
delete(ii)
delete(jj)
delete(fg)

grid_lons = fspan(grid_lona,grid_lonb,grid_lonn)
grid_lats = fspan(grid_lata,grid_latb,grid_latn)
grid_lons@units = "degrees_east"
grid_lats@units = "degrees_north"

;;;;;;;;;;;;;;;;; wind subgrid ;;;;;;;;;;;;;;;;

if WIND then
	
	wind_lona = (floor(min_map_lon/WIND_LOND)-1)*WIND_LOND
	wind_lata = (floor(min_map_lat/WIND_LATD)-1)*WIND_LATD
	wind_lonb = ( ceil(max_map_lon/WIND_LOND)+1)*WIND_LOND
	wind_latb = ( ceil(max_map_lat/WIND_LATD)+1)*WIND_LATD
	wind_lonn = round((wind_lonb-wind_lona)/WIND_LOND,3)+1
	wind_latn = round((wind_latb-wind_lata)/WIND_LATD,3)+1
	
	if VERBOSE then 
		print("wind subgrid OK")
		print("latitude :  min = "+sprintf("%7.4f",wind_lata)+"  max = "+sprintf("%7.4f",wind_latb)+"  n = "+wind_latn)
		print("longitude:  min = "+sprintf("%7.4f",wind_lona)+"  max = "+sprintf("%7.4f",wind_lonb)+"  n = "+wind_lonn)
	end if
	
	ii = round((WIND_LON-wind_lona)/WIND_LOND,3)
	jj = round((WIND_LAT-wind_lata)/WIND_LATD,3)
	fg = (ii.ge.0).and.(ii.lt.wind_lonn).and.(jj.ge.0).and.(jj.lt.wind_latn)
	wind_ind = ind(fg)
	if ismissing(wind_ind(0)) then
		nwind = 0
		print("no wind values inside subgrid")
		print("try to expand map boundaries")
		print("**** procedure aborted ****")
		exit
	else
		nwind = dimsizes(wind_ind)
		wind_i = ii(wind_ind)
		wind_j = jj(wind_ind)
		if VERBOSE then 
			print("number of wind values inside = "+nwind)
			print("")
		end if
	end if
	delete(ii)
	delete(jj)
	delete(fg)
	
	wind_lons = fspan(wind_lona,wind_lonb,wind_lonn)
	wind_lats = fspan(wind_lata,wind_latb,wind_latn)
	wind_lons@units = "degrees_east"
	wind_lats@units = "degrees_north"


end if



