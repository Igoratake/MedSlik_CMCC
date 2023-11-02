_log("start reading wind files")
print("")

_UV = new((/WIND_FNK,4,2,nwind/),float)
inds = wind_ind+5
do ifile = 0,WIND_FNK-1
	lines = asciiread(WIND_FN(ifile),-1,"string")
	tmp = stringtofloat(_get_substr(lines(inds),(/0,11,22,33,44,55,66,77,88,99/), \
																						(/10,21,32,43,54,65,76,87,98,109/)))
	do i = 0,3
		k = i+2
		_UV(ifile,i,:,:) = tmp(k:k+1,:)
	end do
end do
delete(tmp)
delete(inds)
delete(lines)

_uv = new((/ntime,2,nwind/),float)
do itime = 0,ntime-1
	ja = WIND_FNAK(itime)
	jb = WIND_FNBK(itime)
	if (ja.lt.0).or.(jb.lt.0) then
		continue
	end if
	ia = WIND_IA(itime)
	ib = WIND_IB(itime)
;print(""+ja+" "+jb+" "+ia+" "+ib)
	w = WIND_WEIGHTS(itime)
	_uv(itime,:,:) = _UV(ja,ia,:,:)*w+_UV(jb,ib,:,:)*(1.-w)
end do
delete(_UV)
delete(ia)
delete(ib)
delete(w)

wind_uvw = new((/ntime,3,wind_latn,wind_lonn/),float)
do i = 0,nwind-1
	wind_uvw(:,0:1,wind_j(i),wind_i(i)) = _uv(:,:,i)
end do
delete(_uv)

wind_uvw(:,varW,:,:) = sqrt(wind_uvw(:,varU,:,:)^2+wind_uvw(:,varV,:,:)^2)

wind_uvw!2 = "lat"
wind_uvw!3 = "lon"
wind_uvw&lat = wind_lats
wind_uvw&lon = wind_lons

wind_model_uvw = new((/ntime,3,grid_latn,grid_lonn/),float)
wind_model_uvw = linint2(wind_lons,wind_lats,wind_uvw,False,grid_lons,grid_lats,0)
wind_model_uvw!2 = "lat"
wind_model_uvw!3 = "lon"
wind_model_uvw&lat = grid_lats
wind_model_uvw&lon = grid_lons

tmp = avg(wind_model_uvw(:,varW,:,:))
wind_main_scale   = _get_vector_scale(tmp,WIND_SCALE_GRID)
;print("avg="+tmp+" max="+max(wind_model_uvw(:,varW,:,:))+" "+wind_main_scale)
delete(tmp)

if MAIN_MAPS.and.MAIN_MAPS_WIND.and.WIND_CENTER_ONLY then

	wind_center_uvw = new((/ntime,3/),float)
	do itime = 0,ntime-1
		wind_center_uvw(itime,:) = linint2_points(grid_lons,grid_lats,wind_model_uvw(itime,:,:,:),False, \
												oil_center_lon(itime),oil_center_lat(itime),0)
	end do
;	print(""+wind_center_uvw(:,0)+TAB+wind_center_uvw(:,1)+TAB+wind_center_uvw(:,2))
	tmp = avg(wind_center_uvw(:,varW))
	wind_center_scale = _get_vector_scale(tmp,WIND_SCALE_CENTER)
;	print(wind_center_scale)
	delete(tmp)

end if

function _get_wind(itime)
local w,lat,lon,lats,lons,d
begin
	if WIND_CENTER_ONLY then
		d = 0.01
		lon = oil_center_lon(itime)
		lat = oil_center_lat(itime)
		lons = (/lon-d,lon,lon+d/)
		lats = (/lat-d,lat,lat+d/)
		lons@units = "degrees_east"
		lats@units = "degrees_north"
		w = new((/3,3,3/),float)
		w!1 = "lat"
		w!2 = "lon"
		w&lat = lats
		w&lon = lons
		w(:,1,1) = wind_center_uvw(itime,:)
		return w
	else
		return wind_model_uvw(itime,:,::VECTOR_LAT_STRIDE,::VECTOR_LON_STRIDE)
	end if
end



