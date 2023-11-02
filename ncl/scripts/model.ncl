
model_fn = systemfunc(" ls "+MODEL_DIR+"*."+MODEL_EXT)
NMODEL = dimsizes(model_fn)
n = strlen(model_fn(0))
tmp = _get_substr(model_fn,n-12,n-5)
MODEL_YYYYMMDDHH = 2000000000+stringtointeger(tmp(0,:))
delete(tmp)
yyyymmddhh = MODEL_YYYYMMDDHH(0)
yyyy = yyyymmddhh/1000000
mm = yyyymmddhh/10000%100
dd = yyyymmddhh/100%100
hh = yyyymmddhh%100
;;;;;print(""+yyyymmddhh+" "+yyyy+" "+mm+" "+dd+" "+hh)
MODEL_TIMEA = ut_inv_calendar(yyyy,mm,dd,hh,30,0,time_units,0)
yyyymmddhh = MODEL_YYYYMMDDHH(NMODEL-1)
yyyy = yyyymmddhh/1000000
mm = yyyymmddhh/10000%100
dd = yyyymmddhh/100%100
hh = yyyymmddhh%100
MODEL_TIMEB = ut_inv_calendar(yyyy,mm,dd,hh,30,0,time_units,0)
;;;;print(""+MODEL_TIMEA+" "+MODEL_TIMEB)
MODEL_JULIANS = fspan(MODEL_TIMEA,MODEL_TIMEB,NMODEL)
MODEL_IA = doubletointeger(floor((TIMEA-MODEL_TIMEA-0.)*24))-1
MODEL_IB = doubletointeger(ceil ((TIMEB-MODEL_TIMEA+0.)*24))+1
MODEL_IN = MODEL_IB-MODEL_IA+1
_log("directory "+model_dir+" OK")
;print(""+MODEL_IA+" "+MODEL_IB+" "+MODEL_YYYYMMDDHH(MODEL_IA)+" "+MODEL_YYYYMMDDHH(MODEL_IB))

tmp = asciiread(model_fn(0),-1,"string")
tmp2 = str_split(tmp(2)," ")
MODEL_LONA = stringtofloat(tmp2(0))
MODEL_LONB = stringtofloat(tmp2(1))
MODEL_LATA = stringtofloat(tmp2(2))
MODEL_LATB = stringtofloat(tmp2(3))
MODEL_LONN = stringtointeger(tmp2(4))
MODEL_LATN = stringtointeger(tmp2(5))
delete(tmp2)
tmp2 = str_split(tmp(3)," ")
MODEL_N = stringtointeger(tmp2(0))
delete(tmp2)
tmp2 = _get_nums(tmp(5:MODEL_N+4),(/11,11/))
MODEL_LAT = tmp2(0,:)
MODEL_LON = tmp2(1,:)
delete(tmp2)
delete(tmp)

_log("model metadata OK")

MODEL_U = new((/MODEL_IN,MODEL_N/),float)
MODEL_V = new((/MODEL_IN,MODEL_N/),float)
do i = MODEL_IA,MODEL_IB
	k = i-MODEL_IA
	tmp = asciiread(model_fn(i),-1,"string")
;	tmp2 = _get_nums(tmp(5:MODEL_N+4),(/11,11,11,11,11/))
;	MODEL_U(k,:) = tmp2(3,:)
;	MODEL_V(k,:) = tmp2(4,:)
	MODEL_U(k,:) = stringtofloat(_get_substr(tmp(5:MODEL_N+4),33,43))
	MODEL_V(k,:) = stringtofloat(_get_substr(tmp(5:MODEL_N+4),44,54))
	delete(tmp)
;	delete(tmp2)
end do
_log("files "+model_dir+"/*."+MODEL_EXT+" OK")

model_lon = fspan(MODEL_LONA,MODEL_LONB,MODEL_LONN)
model_lon@units    = "degrees_east"
model_lat = fspan(MODEL_LATA,MODEL_LATB,MODEL_LATN)
model_lat@units    = "degrees_north"
model_dx = (MODEL_LONB-MODEL_LONA)/(MODEL_LONN-1)
model_dy = (MODEL_LATB-MODEL_LATA)/(MODEL_LATN-1)
model_i = round((MODEL_LON-MODEL_LONA)/model_dx,3)
model_j = round((MODEL_LAT-MODEL_LATA)/model_dy,3)
model_u1 = new((/MODEL_IN,MODEL_LATN,MODEL_LONN/),float)
model_v1 = new((/MODEL_IN,MODEL_LATN,MODEL_LONN/),float)
print(""+MODEL_LONA+TAB+MODEL_LONB+TAB+MODEL_LONN+TAB+min(model_i)+TAB+max(model_i))
print(""+MODEL_LATA+TAB+MODEL_LATB+TAB+MODEL_LATN+TAB+min(model_j)+TAB+max(model_j))
do i = 0,MODEL_N-1
	model_u1(:,model_j(i),model_i(i)) = MODEL_U(:,i)
	model_v1(:,model_j(i),model_i(i)) = MODEL_V(:,i)
end do
;print(""+avg(MODEL_U)+" "+avg(MODEL_V))
;print(""+avg(model_u1)+" "+avg(model_v1))
model_ia = round(floor((min_map_lon-MODEL_LONA)/model_dx),3)-1
model_ib = round(ceil((max_map_lon-MODEL_LONA)/model_dx),3)+1
model_ja = round(floor((min_map_lat-MODEL_LATA)/model_dy),3)-1
model_jb = round(ceil((max_map_lat-MODEL_LATA)/model_dy),3)+1
print(""+model_ia+TAB+model_ib+TAB+model_ja+TAB+model_jb)
print(""+model_lon(model_ia)+TAB+model_lon(model_ib)+TAB+ \
				 model_lat(model_ja)+TAB+model_lat(model_jb))

_log("model grid OK")

model_u = new((/NTIME+1,MODEL_LATN,MODEL_LONN/),float)
model_v = new((/NTIME+1,MODEL_LATN,MODEL_LONN/),float)

;print(MODEL_YYYYMMDDHH+" "+MODEL_JULIANS)

do itime = 0,NTIME
	julo = JULIANS(itime)
	julm = MODEL_JULIANS(MODEL_IA)
	k = doubletointeger(floor(24*(julo-julm)))
	julm = MODEL_JULIANS(k+MODEL_IA)
	wb = doubletofloat(julo-julm)*24
	wa = doubletofloat(1.-wb)
	model_u(itime,:,:) = model_u1(k,:,:)*wa+model_u1(k+1,:,:)*wb
	model_v(itime,:,:) = model_v1(k,:,:)*wa+model_v1(k+1,:,:)*wb
;print(""+itime+" "+k+" "+MODEL_YYYYMMDDHH(k+MODEL_IA)+" "+julo+" "+julm+" "+wa+" "+wb+ \
;" "+avg(model_u(itime,:,:))+" "+avg(model_v(itime,:,:)))
end do

model_u!1 = "lat"
model_u!2 = "lon"
model_u&lat = model_lat
model_u&lon = model_lon
model_v!1 = "lat"
model_v!2 = "lon"
model_v&lat = model_lat
model_v&lon = model_lon



