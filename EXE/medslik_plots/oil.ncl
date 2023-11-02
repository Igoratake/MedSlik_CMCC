;;;;;;;;;;;;;;;; SRF processing ;;;;;;;;;;;;;;;;;

oil_n    = new(ntime,integer)
oil_center_lat = new(ntime,float)
oil_center_lon = new(ntime,float)
_oil_val = new((/ntime,MAX_OIL_N/),float)
_oil_lon = new((/ntime,MAX_OIL_N/),float)
_oil_lat = new((/ntime,MAX_OIL_N/),float)

_log("start reading oil files")
print("")

do i = 0,ntime-1
	if times(i).gt.0 then
		fn = DIR+PREFIX+MAP_FN(times(i))+SRF
		lines = asciiread(fn,-1,"string")
		m = dimsizes(lines)
		n = stringtointeger(_get_num(lines(7),0))
		oil_n(i) = n
		items = _get_nums(lines(9:n+8),(/9,9,20/))
		_oil_lat(i,0:n-1) = items(0,:)
		_oil_lon(i,0:n-1) = items(1,:)
		_oil_val(i,0:n-1) = items(2,:)
		delete(items)
		items = _get_nums(lines(n+12),(/12,12/))
		oil_center_lat(i) = items(0,0)
		oil_center_lon(i) = items(1,0)
		delete(items)
		delete(lines)
	else
		oil_center_lat(i) = LAT
		oil_center_lon(i) = LON
	end if
end do

REAL_OIL_N = max(oil_n)
if .not.ismissing(REAL_OIL_N) then
	oil_lon = new((/ntime,REAL_OIL_N/),float)
	oil_lat = new((/ntime,REAL_OIL_N/),float)
	oil_val = new((/ntime,REAL_OIL_N/),float)
	oil_lat = _oil_lat(:,0:REAL_OIL_N-1)
	oil_lon = _oil_lon(:,0:REAL_OIL_N-1)
	oil_val = _oil_val(:,0:REAL_OIL_N-1)
end if
;;;;delete(_oil_lat)
;;;;delete(_oil_lon)
;;;;delete(_oil_val)

if VERBOSE then
	print(PREFIX+"NNNN.srf files OK")
	print("maximum number of oil points = "+REAL_OIL_N)
	print("")
end if

;;;;;;;;;;;;;;;; CST processing ;;;;;;;;;;;;;;;;;

oil_cst_n = new(ntime,integer)
oil_cst_n = 0
_oil_cst      = new((/ntime,MAX_OIL_N/),float)
_oil_cst_lona = new((/ntime,MAX_OIL_N/),float)
_oil_cst_lata = new((/ntime,MAX_OIL_N/),float)
_oil_cst_lonb = new((/ntime,MAX_OIL_N/),float)
_oil_cst_latb = new((/ntime,MAX_OIL_N/),float)

do i = 0,ntime-1
	if times(i).gt.0 then
		fn = DIR+PREFIX+MAP_FN(times(i))+CST
		lines = asciiread(fn,-1,"string")
		m = dimsizes(lines)
		n = stringtointeger(_get_num(lines(7),0))
		oil_cst_n(i) = n
		if n.gt.0 then
			items = _get_nums(lines(9:n+8),(/11,11,11,11,20/))
			_oil_cst_lata(i,0:n-1) = items(0,:)
			_oil_cst_lona(i,0:n-1) = items(1,:)
			_oil_cst_latb(i,0:n-1) = items(2,:)
			_oil_cst_lonb(i,0:n-1) = items(3,:)
			_oil_cst(i,0:n-1) = items(4,:)
			delete(items)
		end if
		delete(lines)
	end if
end do

REAL_CST_N = max(oil_cst_n)
if REAL_CST_N.gt.0 then
	oil_cst_lona = new((/ntime,REAL_CST_N/),float)
	oil_cst_lonb = new((/ntime,REAL_CST_N/),float)
	oil_cst_lata = new((/ntime,REAL_CST_N/),float)
	oil_cst_latb = new((/ntime,REAL_CST_N/),float)
	oil_cst = new((/ntime,REAL_CST_N/),float)
	oil_cst_lona = _oil_cst_lona(:,0:REAL_CST_N-1)
	oil_cst_lonb = _oil_cst_lonb(:,0:REAL_CST_N-1)
	oil_cst_lata = _oil_cst_lata(:,0:REAL_CST_N-1)
	oil_cst_latb = _oil_cst_latb(:,0:REAL_CST_N-1)
	oil_cst = _oil_cst(:,0:REAL_CST_N-1)
end if
delete(_oil_cst_lona)
delete(_oil_cst_lonb)
delete(_oil_cst_lata)
delete(_oil_cst_latb)
delete(_oil_cst)

;VERBOSE = True


if (REAL_CST_N.gt.0) then
	print(PREFIX+"NNNN.cst files OK")
	print("total number of coastal points = "+num(.not.ismissing(oil_cst)))
	print("maximum number of coastal points = "+REAL_CST_N)
	print("min = "+min(oil_cst)+TAB+" max = "+max(oil_cst))
	print("")
end if



;;;;;;;;;;;;;;;; map boundaries  ;;;;;;;;;;;;;;;;;

min_oil_lat = min(oil_lat)
max_oil_lat = max(oil_lat)
min_oil_lon = min(oil_lon)
max_oil_lon = max(oil_lon)

if NOBS.gt.0 then
	max_oil_lat = max((/max_oil_lat,OBS_MAX_LAT/))
	max_oil_lon = max((/max_oil_lon,OBS_MAX_LON/))
	min_oil_lat = min((/min_oil_lat,OBS_MIN_LAT/))
	min_oil_lon = min((/min_oil_lon,OBS_MIN_LON/))
end if

d = 60./LON_CELL
min_map_lon = floor(d*min_oil_lon)/d-(LON_MARGIN/60.)
max_map_lon =  ceil(d*max_oil_lon)/d+(LON_MARGIN/60.)

d = 60./LAT_CELL
min_map_lat = floor(d*min_oil_lat)/d-(LAT_MARGIN/60.)
max_map_lat =  ceil(d*max_oil_lat)/d+(LAT_MARGIN/60.)

d = 60./LON_CELL
min_map_lon = floor(d*(min_oil_lon-(LON_MARGIN/60.)))/d
max_map_lon = floor(d*(max_oil_lon+(LON_MARGIN/60.)))/d

d = 60./LAT_CELL
min_map_lat = floor(d*(min_oil_lat-(LAT_MARGIN/60.)))/d
max_map_lat = floor(d*(max_oil_lat+(LAT_MARGIN/60.)))/d

delete(d)

if VERBOSE then
	print("map boundaries OK")
	print("data latitude :  min = "+sprintf("%7.4f",min_oil_lat)+ \
				"  max = "+sprintf("%7.4f",max_oil_lat)+"  range = "+sprintf("%6.4f",max_oil_lat-min_oil_lat))
	print("data longitude:  min = "+sprintf("%7.4f",min_oil_lon)+ \
				"  max = "+sprintf("%7.4f",max_oil_lon)+"  range = "+sprintf("%6.4f",max_oil_lon-min_oil_lon))
	print("auto latitude :  min = "+sprintf("%7.4f",min_map_lat)+ \
				"  max = "+sprintf("%7.4f",max_map_lat)+"  range = "+sprintf("%6.4f",max_map_lat-min_map_lat))
	print("auto longitude:  min = "+sprintf("%7.4f",min_map_lon)+ \
				"  max = "+sprintf("%7.4f",max_map_lon)+"  range = "+sprintf("%6.4f",max_map_lon-min_map_lon))
end if

if isvar("MIN_LAT") then
	min_map_lat = MIN_LAT
	if VERBOSE then
		print("manual  latitude min = "+sprintf("%7.4f",min_map_lat))
	end if
end if
if isvar("MAX_LAT") then
	max_map_lat = MAX_LAT
	if VERBOSE then
		print("manual  latitude max = "+sprintf("%7.4f",max_map_lat))
	end if
end if
if isvar("MIN_LON") then
	min_map_lon = MIN_LON
	if VERBOSE then
		print("manual longitude min = "+sprintf("%7.4f",min_map_lon))
	end if
end if
if isvar("MAX_LON") then
	max_map_lon = MAX_LON
	if VERBOSE then
		print("manual longitude max = "+sprintf("%7.4f",max_map_lon))
	end if
end if

if VERBOSE then
	print("")
end if

;;;;;;;;;;;;;;;; oil grid  ;;;;;;;;;;;;;;;;;

if MAIN_MAPS then

	oil_dy = PIXEL/degree2meter
;	oil_dx = oil_dy/cos(degree2radian*(min_map_lat+max_map_lat)/2)
	oil_dx = oil_dy/cos(degree2radian*LAT)

if False then
	oil_lonn = round((max_oil_lon-min_oil_lon)/oil_dx,3)+1
	oil_lona = min_oil_lon
	oil_lonb = max_oil_lon
	oil_lons = fspan(oil_lona,oil_lonb,oil_lonn)
	oil_lons@units    = "degrees_east"

	oil_latn = round((max_oil_lat-min_oil_lat)/oil_dy,3)+1
	oil_lata = min_oil_lat
	oil_latb = max_oil_lat
	oil_lats = fspan(oil_lata,oil_latb,oil_latn)
	oil_lats@units    = "degrees_north"
else
	oil_lonn = floattointeger(ceil((max_oil_lon-min_oil_lon)/oil_dx+1))
	oil_lona = min_oil_lon
	oil_lonb = max_oil_lon
;	oil_lonn = oil_lonn+1
	oil_lons = fspan(oil_lona,oil_lonb,oil_lonn)
	oil_lons@units    = "degrees_east"

	oil_latn = floattointeger(ceil((max_oil_lat-min_oil_lat)/oil_dy+1))
	oil_lata = min_oil_lat
	oil_latb = min_oil_lat+oil_dy*(oil_latn-0)
	oil_latb = max_oil_lat
;	oil_latn = oil_latn+1
	oil_lats = fspan(oil_lata,oil_latb,oil_latn)
	oil_lats@units    = "degrees_north"
end if


;print(oil_lats)
;print(oil_lons)
;print(oil_dy)
;print(oil_dx)
;exit



	if VERBOSE then
		print("oil map grid OK")
		print("latitude  :  step = "+sprintf("%7.5f",oil_dy)+"  n = "+oil_latn)
		print("longitude :  step = "+sprintf("%7.5f",oil_dx)+"  n = "+oil_lonn)
		print("")
	end if

	oil_map = new((/ntime,oil_latn,oil_lonn/),float)
	oil_map!1 = "lat"
	oil_map!2 = "lon"
	oil_map&lat = oil_lats
	oil_map&lon = oil_lons

	do itime = 0,ntime-1
		if times(itime).gt.0 then
			n = oil_n(itime)
			vals = oil_val(itime,0:n-1)
			lons = oil_lon(itime,0:n-1)
			lats = oil_lat(itime,0:n-1)
	;		ii = round((lons-oil_lona)/oil_dx,3)
	;		jj = round((lats-oil_lata)/oil_dy,3)
			ii = floattointeger((lons-oil_lona)/oil_dx)
			jj = floattointeger((lats-oil_lata)/oil_dy)
			do i = 0,n-1
				oil_map(itime,jj(i),ii(i)) = vals(i)
			end do
			delete(vals)
			delete(lons)
			delete(lats)
			delete(ii)
			delete(jj)
		end if
	end do
end if


;_log("min_oil="+min(oil_val)+" avg_oil="+avg(oil_val)+" max_oil="+max(oil_val))

pa = 0.001
pb = 0.001
tmp = ndtooned(oil_val)
m = dimsizes(tmp)
qsort(tmp)
n = num(.not.ismissing(tmp))
;print("all="+m+" data="+n)
;print("min="+tmp(m-n)+" max="+tmp(m-1))
k = round(m-n+pa*n,3)
va = tmp(k)
k = round(m-1-pb*n,3)
vb = tmp(k)
delete(tmp)
;print("va="+va+" vb="+vb)

ia = -5
ib =  5

if (log10(vb)-log10(va)).gt.4 then
	m = 3*(ib-ia+1)
	tmp = new(m,float)
	do i = ia,ib
		k = i
		j = 3*(i-ia)
		tmp(j+0) = 10.^k*1.
		tmp(j+1) = 10.^k*2.
		tmp(j+2) = 10.^k*5.
	end do
else
	m = 6*(ib-ia+1)
	tmp = new(m,float)
	do i = ia,ib
		k = i
		j = 6*(i-ia)
		tmp(j+0) = 10.^k*1.0
		tmp(j+1) = 10.^k*1.5
		tmp(j+2) = 10.^k*2.0
		tmp(j+3) = 10.^k*3.0
		tmp(j+4) = 10.^k*5.0
		tmp(j+5) = 10.^k*7.0
	end do
end if
ka = 0
kb = m-1
do i = 0,m-1
	if va.lt.tmp(i) then
		ka = i-0
		break
	end if
end do
do i = 1,m
	if vb.gt.tmp(m-i) then
		kb = m-i+0
		break
	end if
end do

;print("ka="+ka+" ka="+kb)
oil_scale  = tmp(ka:kb)
oil_labels = ""+oil_scale
;print(oil_scale)
;print("min="+tmp(ka)+" max="+tmp(kb)+" n="+dimsizes(oil_scale))
delete(tmp)

;

if COAST_MAP.and.(REAL_CST_N.gt.0) then

	pa = 0.001
	pb = 0.001
	tmp = ndtooned(oil_cst)
	m = dimsizes(tmp)
	qsort(tmp)
	n = num(.not.ismissing(tmp))
	;print("all="+m+" data="+n)
	;print("min="+tmp(m-n)+" max="+tmp(m-1))
	k = round(m-n+pa*n,3)
	va_cst = tmp(k)
	k = round(m-1-pb*n,3)
	vb_cst = tmp(k)
	delete(tmp)
	print("va_cst="+va_cst+" vb_cst="+vb_cst)

	ia = -5
	ib =  5
	if (log10(vb_cst)-log10(va_cst)).gt.4 then
		m = 3*(ib-ia+1)
		tmp = new(m,float)
		do i = ia,ib
			k = i
			j = 3*(i-ia)
			tmp(j+0) = 10.^k*1.
			tmp(j+1) = 10.^k*2.
			tmp(j+2) = 10.^k*5.
		end do
	else
		m = 6*(ib-ia+1)
		tmp = new(m,float)
		do i = ia,ib
			k = i
			j = 6*(i-ia)
			tmp(j+0) = 10.^k*1.0
			tmp(j+1) = 10.^k*1.5
			tmp(j+2) = 10.^k*2.0
			tmp(j+3) = 10.^k*3.0
			tmp(j+4) = 10.^k*5.0
			tmp(j+5) = 10.^k*7.0
		end do
	end if
	ka = 0
	kb = m-1
	do i = 0,m-1
		if va_cst.lt.tmp(i) then
			ka = i-0
			break
		end if
	end do
	do i = 1,m
		if vb_cst.gt.tmp(m-i) then
			kb = m-i+0
			break
		end if
	end do

	print("ka="+ka+" ka="+kb)
	cst_scale  = tmp(ka:kb)
	cst_labels = ""+oil_scale
	print(cst_scale)
	print("min="+tmp(ka)+" max="+tmp(kb)+" n="+dimsizes(cst_scale))
	delete(tmp)

;	exit

end if
































if False then

;dssetp("dmx",10*oil_dx)
;dssetp("erf","/dev/null")

new_oil_lona = min_oil_lon-oil_dx
new_oil_lonb = max_oil_lon+oil_dx
new_oil_lonn = round((new_oil_lonb-new_oil_lona)/oil_dx*4,3)+1
new_oil_lons = fspan(new_oil_lona,new_oil_lonb,new_oil_lonn)
new_oil_lons@units    = "degrees_east"

new_oil_lata = min_oil_lat
new_oil_latb = max_oil_lat
new_oil_latn = round((new_oil_latb-new_oil_lata)/oil_dy*4,3)+1
new_oil_lats = fspan(new_oil_lata,new_oil_latb,new_oil_latn)
new_oil_lats@units    = "degrees_north"


function _get_oil(itime:integer)
local n,lons,lats,vals,minv,tmp
begin
	n = oil_n(itime)
	vals = oil_val(itime,0:n-1)
	lons = oil_lon(itime,0:n-1)
	lats = oil_lat(itime,0:n-1)
	minv = 2*min(vals)
	tmp = dsgrid2s(lons,lats,vals,new_oil_lons,new_oil_lats)
printVarSummary(tmp)
printVarSummary(new_oil_lons)
printVarSummary(new_oil_lats)
	tmp!0 = "lat"
	tmp!1 = "lon"
	tmp&lat = new_oil_lats
	tmp&lon = new_oil_lons
	tmp@_FillValue = -999.
	tmp = where(tmp.lt.minv,tmp@_FillValue,tmp)
	return tmp
end

function _get_oil1(itime:integer)
local n,lons,lats,vals,minv,err,tmp
begin
	err = NhlGetErrorObjectId()
;  setvalues err
;    "errLevel" : "Error"
;  end setvalues
	tmp = linint2(oil_lons,oil_lats,oil_map(itime,:,:),False,new_oil_lons,new_oil_lats,0)
	tmp!0 = "lat"
	tmp!1 = "lon"
	tmp&lat = new_oil_lats
	tmp&lon = new_oil_lons
	return tmp
end


_log("")
tmp11 = _get_oil(4)
_log("")
printMinMax(oil_map(4,:,:),True)
print(sum(oil_map(4,:,:)))
print(num(ismissing(oil_map(4,:,:))))

printMinMax(tmp11,True)
printVarSummary(tmp11)
print(num(ismissing(tmp11)))
print(avg(tmp11))
print("")
exit

end if
