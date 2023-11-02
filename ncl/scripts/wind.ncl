
wind_fn = systemfunc(" ls "+DIR_E25+"*.ecm")
NWIND = dimsizes(wind_fn)
n = strlen(wind_fn(0))
tmp = _get_substr(wind_fn,n-10,n-5)
WIND_YYYYMMDD = 20000000+stringtointeger(tmp(0,:))
delete(tmp)

tmp = asciiread(wind_fn(0),-1,"string")
tmp2 = str_split(tmp(2)," ")
WIND_LONA = stringtofloat(tmp2(0))
WIND_LONB = stringtofloat(tmp2(1))
WIND_LATA = stringtofloat(tmp2(2))
WIND_LATB = stringtofloat(tmp2(3))
WIND_LONN = stringtointeger(tmp2(4))
WIND_LATN = stringtointeger(tmp2(5))
delete(tmp2)
tmp2 = str_split(tmp(3)," ")
WIND_N = stringtointeger(tmp2(0))
delete(tmp2)
tmp2 = _get_nums(tmp(5:WIND_N+4),(/11,11/))
WIND_LAT = tmp2(0,:)
WIND_LON = tmp2(1,:)
delete(tmp2)
delete(tmp)
_log("directory E25 OK")

WIND_JULIAN = new(NWIND*4,double)
WIND_U = new((/NWIND*4,WIND_N/),float)
WIND_V = new((/NWIND*4,WIND_N/),float)
do i = 0,NWIND-1
	tmp = asciiread(wind_fn(i),-1,"string")
	tmp2 = _get_nums(tmp(5:WIND_N+4),(/11,11,11,11,11,11,11,11,11,11/))
	delete(tmp)
	yyyymmdd = WIND_YYYYMMDD(i)
	yyyy = yyyymmdd/10000
	mm = yyyymmdd/100%100
	dd = yyyymmdd%100
	do j = 0,3
	k = i*4+j
		WIND_JULIAN(k) = ut_inv_calendar(yyyy,mm,dd,j*6,0,0,time_units,0)
		WIND_U(k,:) = tmp2(2+2*j,:)
		WIND_V(k,:) = tmp2(3+2*j,:)
	end do
	delete(tmp2)
end do
_log("files E25/*.ecm OK")

wind_lon = fspan(WIND_LONA,WIND_LONB,WIND_LONN)
wind_lon@units    = "degrees_east"
wind_lat = fspan(WIND_LATA,WIND_LATB,WIND_LATN)
wind_lat@units    = "degrees_north"

;wind_dx = (WIND_LONB-WIND_LONA)/(WIND_LONN-1)
;wind_dy = (WIND_LATB-WIND_LATA)/(WIND_LATN-1)
;wind_i = round((WIND_LON-WIND_LONA)/wind_dx,3)
;wind_j = round((WIND_LAT-WIND_LATA)/wind_dy,3)
;;;print(""+wind_dx+" "+wind_dy)

wind_u6 = new((/NWIND*4,WIND_LATN,WIND_LONN/),float)
wind_v6 = new((/NWIND*4,WIND_LATN,WIND_LONN/),float)
do i = 0,NWIND-1
	do j = 0,3
		k = 4*i+j
		wind_u6(k,:,:) = transpose(onedtond(WIND_U(k,:),(/WIND_LONN,WIND_LATN/)))
		wind_v6(k,:,:) = transpose(onedtond(WIND_V(k,:),(/WIND_LONN,WIND_LATN/)))
	end do
end do

wind_u = new((/NTIME+1,WIND_LATN,WIND_LONN/),float)
wind_v = new((/NTIME+1,WIND_LATN,WIND_LONN/),float)
do itime = 0,NTIME
	julo = JULIANS(itime)
	julw = WIND_JULIAN(0)
	k = doubletointeger(floor(4*(julo-julw)))
	julw = WIND_JULIAN(k)
	wb = doubletofloat((julo-julw)*4)
	wa = doubletofloat(1.-wb)
;	print(""+k+" "+WIND_YYYYMMDD(k/4)+" "+(k%4)+" "+julo+" "+julw+" "+wa)
	wind_u(itime,:,:) = wind_u6(k,:,:)*wa+wind_u6(k+1,:,:)*wb
	wind_v(itime,:,:) = wind_v6(k,:,:)*wa+wind_v6(k+1,:,:)*wb
end do

wind_center_u = new(NTIME+1,float)
wind_center_v = new(NTIME+1,float)
do itime = 0,NTIME
	if itime.eq.0 then
		x = LON
		y = LAT
	else
		x = oil_center_lon(i-1)
		y = oil_center_lat(i-1)
	end if
	wind_center_u(itime) = linint2_points(wind_lon,wind_lat,wind_u(itime,:,:),False,x,y,0)
	wind_center_v(itime) = linint2_points(wind_lon,wind_lat,wind_v(itime,:,:),False,x,y,0)
end do

wind_model_u = new((/NTIME+1,MODEL_LATN,MODEL_LONN/),float)
wind_model_v = new((/NTIME+1,MODEL_LATN,MODEL_LONN/),float)
wind_model_u = linint2(wind_lon,wind_lat,wind_u,False,model_lon,model_lat,0)
wind_model_v = linint2(wind_lon,wind_lat,wind_v,False,model_lon,model_lat,0)
wind_model_u!1 = "lat"
wind_model_u!2 = "lon"
wind_model_u&lat = model_lat
wind_model_u&lon = model_lon
wind_model_v!1 = "lat"
wind_model_v!2 = "lon"
wind_model_v&lat = model_lat
wind_model_v&lon = model_lon


;print(""+WIND_LONN+" "+min(wind_i)+" "+max(wind_i))
;print(""+WIND_LATN+" "+min(wind_j)+" "+max(wind_j))
;do i = 0,WIND_N-1
;	wind_u(:,wind_j(i),wind_i(i)) = WIND_U(:,i)
;	wind_v(:,wind_j(i),wind_i(i)) = WIND_V(:,i)
;end do
;do i = 0,WIND_N-1
;	do j = 1,4
;		wa = 6.-j
;		wb = j
;		wind_u(6*i+j,) = WIND_U(:,i)
;		wind_v(6*i+j,:,:) = WIND_V(:,i)
;	end do
;end do

procedure _process_wind(itime:integer)
begin
	print(""+itime)
end

_log("wind grid OK")

