_log("start reading currents files")
print("")

_UV = new((/MODEL_FNK,2,ngrid/),float)
inds = grid_ind+5
do ifile = 0,MODEL_FNK-1
	lines = asciiread(MODEL_FN(ifile),-1,"string")
	_UV(ifile,:,:) = stringtofloat(_get_substr(lines(inds),(/33,44/),(/43,54/)))
end do
delete(inds)
delete(lines)

_uv = new((/ntime,2,ngrid/),float)
do itime = 0,ntime-1
	ia = MODEL_FNAK(itime)
	ib = MODEL_FNBK(itime)
	if (ia.lt.0).or.(ib.lt.0) then
		continue
	end if
	w = MODEL_WEIGHTS(itime)
	_uv(itime,:,:) = _UV(ia,:,:)*w+_UV(ib,:,:)*(1.-w)
end do
delete(_UV)
delete(ia)
delete(ib)
delete(w)

model_uvw = new((/ntime,3,grid_latn,grid_lonn/),float)
do i = 0,ngrid-1
	model_uvw(:,0:1,grid_j(i),grid_i(i)) = _uv(:,:,i)
end do
delete(_uv)


;do itime = 0,ntime-1
;	print(num(ismissing(model_uvw(itime,varU,:,:))))
;end do

model_uvw(:,varW,:,:) = sqrt(model_uvw(:,varU,:,:)^2+model_uvw(:,varV,:,:)^2)

model_uvw!2 = "lat"
model_uvw!3 = "lon"
model_uvw&lat = grid_lats
model_uvw&lon = grid_lons

tmp = avg(model_uvw(:,varW,:,:))
model_main_scale   = _get_vector_scale(tmp,CURRENTS_SCALE_GRID)
model_map_scale    = _get_vector_scale(tmp,CURRENTS_SCALE_MAP)
;print("avg="+tmp+" max="+max(model_uvw(:,varW,:,:))+" "+model_main_scale+" "+model_map_scale)
delete(tmp)

if MAIN_MAPS.and.MAIN_MAPS_CURRENTS.and.CURRENTS_CENTER_ONLY then

	model_center_uvw = new((/ntime,3/),float)
	do itime = 0,ntime-1
		model_center_uvw(itime,:) = linint2_points(grid_lons,grid_lats,model_uvw(itime,:,:,:),False, \
												oil_center_lon(itime),oil_center_lat(itime),0)
	end do
;	print(""+model_center_uvw(:,0)+TAB+model_center_uvw(:,1)+TAB+model_center_uvw(:,2))
	tmp = avg(model_center_uvw(:,varW))
	model_center_scale = _get_vector_scale(tmp,CURRENTS_SCALE_CENTER)
	delete(tmp)

end if

function _get_model(itime)
local w,lat,lon,lats,lons,d
begin
	if CURRENTS_CENTER_ONLY then
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
		w(:,1,1) = model_center_uvw(itime,:)
		return w
	else
		return model_uvw(itime,:,::VECTOR_LAT_STRIDE,::VECTOR_LON_STRIDE)
	end if
end





