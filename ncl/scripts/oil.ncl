fn = DIR+"medslik*.inp"
tmp = systemfunc("ls "+DIR+"medslik*.inp")
fn = tmp(0)
delete(tmp)
tmp = asciiread(fn,-1,"string")
tmp2 = str_split(tmp(3)," ")
DAY   = stringtointeger(tmp2(0))
MONTH = stringtointeger(tmp2(1))
YEAR  = stringtointeger(tmp2(2))
delete(tmp2)
tmp2 = stringtointeger(_get_num(tmp(4),0))
HOUR   = tmp2/100
MINUTE = tmp2%100
delete(tmp2)
TIMEA = ut_inv_calendar(YEAR,MONTH,DAY,HOUR,MINUTE,0,time_units,0)
NTIME = stringtointeger(_get_num(tmp(9),0))
PIXEL = stringtofloat(_get_num(tmp(27),0))
;PIXEL = 500
delete(tmp)
TIMEB = TIMEA+(NTIME/24.)
;;;;;print(""+_julian2yyyymmddhhnn(TIMEA)+" - "+_julian2yyyymmddhhnn(TIMEB))
JULIANS = fspan(TIMEA,TIMEB,NTIME+1)
JULIANS@units = time_units
_log("file medslik5.inp OK")

fn = DIR+"medslik.fte"
_check_file_exists(fn)
tmp = asciiread(fn,-1,"string")
tmp2 = str_split(tmp(1)," ")
LAT   = stringtofloat(tmp2(2))
LON   = stringtofloat(tmp2(4))
delete(tmp2)
NTIME = stringtointeger(_get_num(tmp(4),3))
FTE = _get_nums(tmp(8:7+NTIME*2),(/9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9/))
;FTE = readAsciiTable(fn,17,"float",7)
delete(tmp)
_log("file medslik.fte OK")

oil     = new((/NTIME,MAXN/),float)
oil_n   = new(NTIME,integer)
oil_lon = new((/NTIME,MAXN/),float)
oil_lat = new((/NTIME,MAXN/),float)

oil_uv_n   = new(NTIME,float)
oil_uv_lon = new((/NTIME,MAXN/),float)
oil_uv_lat = new((/NTIME,MAXN/),float)
oil_u = new((/NTIME,MAXN/),float)
oil_v = new((/NTIME,MAXN/),float)
oil_center_lon = new(NTIME,float)
oil_center_lat = new(NTIME,float)
oil_center_u = new(NTIME,float)
oil_center_v = new(NTIME,float)

if False then
oil_tot     = new((/NTIME,MAXN/),float)
oil_tot_n   = new(NTIME,float)
oil_tot_lon = new((/NTIME,MAXN/),float)
oil_tot_lat = new((/NTIME,MAXN/),float)
end if

oil_slick_n = new(NTIME,float)
oil_slick_lon = new((/NTIME,MAXN/),float)
oil_slick_lat = new((/NTIME,MAXN/),float)

oil_cst_n = new(NTIME,float)
oil_cst      = new((/NTIME,MAXN/),float)
oil_cst_lona = new((/NTIME,MAXN/),float)
oil_cst_lata = new((/NTIME,MAXN/),float)
oil_cst_lonb = new((/NTIME,MAXN/),float)
oil_cst_latb = new((/NTIME,MAXN/),float)

;print("nnnn"+TAB+"srf"+TAB+"uv"+TAB+"tot"+TAB+"slick"+TAB+"coast")

do i = 0,NTIME-1

	iiii = sprinti("%0.4i",i+1)
	
	fn = DIR+"out"+iiii+".srf"
	tmp = asciiread(fn,-1,"string")
	m = dimsizes(tmp)
	n = stringtointeger(_get_num(tmp(7),0))
	oil_n(i) = n
;print("out"+iiii+".srf "+n)
;_GET_NUMS_DEBUG = (i.eq.9)
	tmp2 = _get_nums(tmp(9:n+8),(/9,9,20/))
	oil_lat(i,0:n-1) = tmp2(0,:)
	oil_lon(i,0:n-1) = tmp2(1,:)
	oil(i,0:n-1) = tmp2(2,:)
	delete(tmp2)
	
	nskipa = 9+n+3
	nuv = m-nskipa
	tmp2 = _get_nums(tmp(nskipa:m-1),(/12,12,9,9/))
	delete(tmp)
	oil_uv_n(i) = nuv-1
	oil_center_lat(i) = tmp2(0,0)
	oil_center_lon(i) = tmp2(1,0)
	oil_center_u(i) = tmp2(2,0)
	oil_center_v(i) = tmp2(3,0)
	oil_uv_lat(i,0:nuv-2) = tmp2(0,1:nuv-1)
	oil_uv_lon(i,0:nuv-2) = tmp2(1,1:nuv-1)
	oil_u(i,0:nuv-2) = tmp2(2,1:nuv-1)
	oil_v(i,0:nuv-2) = tmp2(3,1:nuv-1)
	delete(tmp2)
	
	fn = DIR+"out"+iiii+".tot"
	tmp = asciiread(fn,-1,"string")
	m = dimsizes(tmp)
if False then
	n = stringtointeger(_get_num(tmp(7),0))
	oil_tot_n(i) = n
	tmp2 = _get_nums(tmp(9:n+8),(/9,9,20/))
	oil_tot_lat(i,0:n-1) = tmp2(0,:)
	oil_tot_lon(i,0:n-1) = tmp2(1,:)
	oil_tot(i,0:n-1) = tmp2(2,:)
	delete(tmp2)

	nskipa = 9+n+5
end if
	nskipa = ind(tmp.eq." Oil slick contour")+2
	nslick = m-nskipa
	tmp2 = _get_nums(tmp(nskipa:m-1),(/9,9/))
	delete(tmp)
	oil_slick_n(i) = nslick
	oil_slick_lat(i,0:nslick-1) = tmp2(0,:)
	oil_slick_lon(i,0:nslick-1) = tmp2(1,:)
	delete(tmp2)
	
	fn = DIR+"out"+iiii+".cst"
;print(""+fn)
	tmp = asciiread(fn,-1,"string")
	m = dimsizes(tmp)
	n = stringtointeger(_get_num(tmp(7),0))
	oil_cst_n(i) = n
	if n.gt.0 then
		tmp2 = _get_nums(tmp(9:n+8),(/11,11,11,11,20/))
		oil_cst_lata(i,0:n-1) = tmp2(0,:)
		oil_cst_lona(i,0:n-1) = tmp2(1,:)
		oil_cst_latb(i,0:n-1) = tmp2(2,:)
		oil_cst_lonb(i,0:n-1) = tmp2(3,:)
		oil_cst(i,0:n-1) = tmp2(4,:)
		delete(tmp2)
	end if
	delete(tmp)

;print(iiii+TAB+oil_n(i)+TAB+oil_uv_n(i)+TAB+oil_tot_n(i)+TAB+oil_slick_n(i)+TAB+oil_cst_n(i))

end do

oil = where(ismissing(oil),1e20,oil)
oil@_FillValue = 1e20
oil_map_spacing = 0.0
oil_map_stride = 5
oil_map_a = 0.
oil_map_b = 100.
_log("min oil="+min(oil)+" max oil="+max(oil))

pa = 0.001
pb = 0.999
tmp = ndtooned(oil)
qsort(tmp)
n = num(.not.ismissing(tmp))
k = round(pa*n,3)
va = tmp(k)
k = round(pb*n,3)
vb = tmp(k)
delete(tmp)
print("va="+va+" vb="+vb)


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
print("ka="+ka+" ka="+kb)
oil_scale = tmp(ka:kb)
;print(oil_scale)
print("min="+tmp(ka)+" max="+tmp(kb)+" n="+dimsizes(oil_scale))
delete(tmp)

;exit


min_lat = min(oil_lat)
max_lat = max(oil_lat)
min_lon = min(oil_lon)
max_lon = max(oil_lon)
d=1./12.
n = 12
d = 1./12
min_map_lon = floor(n*min_lon)/n-d
min_map_lat = floor(n*min_lat)/n-d
max_map_lon = ceil(n*max_lon)/n+d
max_map_lat = ceil(n*max_lat)/n+d
;print(""+min_lat+TAB+max_lat+TAB+min_lon+TAB+max_lon)
;print(""+min_map_lat+TAB+max_map_lat+TAB+min_map_lon+TAB+max_map_lon)
DY = PIXEL/degree2meter
DX = DY/cos(degree2radian*(min_map_lat+max_map_lat)/2)
;print(""+DY+TAB+DX+TAB+ \ 
;round(((max_map_lat-min_map_lat)/DY+1),3)+TAB+ \
;round(((max_map_lon-min_map_lon)/DX+1),3)+TAB+ \
;round(((max_lat-min_lat)/DY+1),3)+TAB+ \
;round(((max_lon-min_lon)/DX+1),3))

oil_map_lonn = floattointeger(ceil((max_lon-min_lon)/DX+1))
oil_map_lona = min_lon
oil_map_lonb = min_lon+DX*oil_map_lonn
oil_map_lonn= oil_map_lonn+1
oil_map_lon = fspan(oil_map_lona,oil_map_lonb,oil_map_lonn)
oil_map_lon@units    = "degrees_east"

oil_map_latn = floattointeger(ceil((max_lat-min_lat)/DY+1))
oil_map_lata = min_lat
oil_map_latb = min_lat+DY*oil_map_latn
oil_map_latn = oil_map_latn+1
oil_map_lat = fspan(oil_map_lata,oil_map_latb,oil_map_latn)
oil_map_lat@units    = "degrees_north"

;print(oil_map_lon)
;print(oil_map_lat)
;print(oil_map_lona)
;print(oil_map_lonb)
;print(oil_map_lata)
;print(oil_map_latb)

oil_map = new((/oil_map_latn,oil_map_lonn/),float)
oil_map!0 = "lat"
oil_map!1 = "lon"
oil_map&lat = oil_map_lat
oil_map&lon = oil_map_lon
oil_map@_FillValue = 1e20

_log("files out*.* OK")

procedure _process_oil(itime:integer)
local i,n,data,lons,lats,ii,jj
begin
	n = oil_n(itime)
	data = oil(itime,0:n-1)
	lons = oil_lon(itime,0:n-1)
	lats = oil_lat(itime,0:n-1)
	ii = round((lons-oil_map_lona)/DX,3)
	jj = round((lats-oil_map_lata)/DY,3)
;printMinMax(ii,True)
;printMinMax(jj,True)
;printVarSummary(oil_map)
	oil_map = oil_map@_FillValue
	do i = 0,n-1
		oil_map(jj(i),ii(i)) = (/data(i)/)
	end do
end
