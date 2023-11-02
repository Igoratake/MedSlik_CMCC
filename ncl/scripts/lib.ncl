;================================================;
; Some useful routines and constants 
;================================================;

; put time to any output for benchmarking
procedure _log(s:string)
begin
  print(systemfunc("date +%T")+" "+s)
end

START=0.

function _seconds_to_timestr(sec:float)
local str,h,m,s,n
begin
  n = floattointeger(sec)
  h = n/3600
  m = n/60 - h*60
  s = n%60
  str = ""
  if h.gt.0 then
    if h.eq.1 then
      str = str+h+" hour "
    else 
      str = str+h+" hours "  
    end if
  end if      
  if m.gt.0 then
    if m.eq.1 then
      str = str+m+" minute "
    else 
      str = str+m+" minutes "  
    end if
  end if  
  str = str+s+" seconds"
  return str     
end

procedure _log_start()
local ss,ch,h,m,s
begin
  ss = systemfunc("date +%T'.'%N")
  ch = stringtochar(ss)
  ss = chartostring(ch(0:9))
  h = stringtointeger(chartostring(ch(0:1)))
  m = stringtointeger(chartostring(ch(3:4)))
  s = stringtofloat(chartostring(ch(6:9)))
  START = s+60*(m+60*h)
  _log("start")
end

procedure _log_finish()
local ss,ch,h,m,s,FINISH
begin
  ss = systemfunc("date +%T'.'%N")
  ch = stringtochar(ss)
  ss = chartostring(ch(0:9))
  h = stringtointeger(chartostring(ch(0:1)))
  m = stringtointeger(chartostring(ch(3:4)))
  s = stringtofloat(chartostring(ch(6:9)))
  FINISH = s+60*(m+60*h)
  _log("finish")
  print("Elapsed time: "+_seconds_to_timestr(FINISH-START))
end

_log_start

;================================================;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
;================================================;

setfileoption("bin","ReadByteOrder" ,"BigEndian")   ; all data files
setfileoption("bin","WriteByteOrder","BigEndian")   ; all data files

setfileoption("nc","MissingToFillValue",False)

LBR = chartostring(inttochar((/10,0/)))  ; line break
TAB = chartostring(inttochar((/ 9,0/)))  ; tabulation

months_full   = (/"January","February","March","April","May","June","July", \
                  "August","September", "October","November","December"/)
months_MMM    = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", "Oct","Nov","Dec"/)
months_mmm    = (/"jan","feb","mar","apr","may","jun","jul","aug","sep", "oct","nov","dec"/)
months_M      = (/"J",  "F",  "M",  "A",  "M",  "J",  "J",  "A",  "S",   "O",  "N",  "D"/)
month_days    =  (/31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./)
month_days366 =  (/31.,29.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./)
month_nums    =  (/ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11 /)

pi            = 3.14159
degree2radian = pi/180.
mile          = 1852.    ; nautical mile in meters
degree2meter  = 60.*mile

time_units    = "days after 1999-12-31 00:00:00"

today_year    = stringtointeger(systemfunc("date +%Y")) ; YYYY
today_month   = stringtointeger(systemfunc("date +%m")) ; MM
today_day     = stringtointeger(systemfunc("date +%d")) ; DD
today_weekday = stringtointeger(systemfunc("date +%w")) ; "3" means wednesday
today_julian  = ut_inv_calendar(today_year,today_month,today_day,0,0,0,time_units,0)  

;==============================

; integer y,m,d -> double julian day
function _get_julian_y_m_d(y:integer, m:integer, d:integer)
begin
  return ut_inv_calendar(y,m,d,0,0,0,time_units,0)
end

;==============================

; integer yyyymmdd -> double julian day
function _get_julian_yyyymmdd(yyyymmdd:integer)
local y, m, d, z
begin
  z = yyyymmdd
  z = 0
  d =  yyyymmdd % 100
  m = (yyyymmdd % 10000 - d) / 100
  y = (yyyymmdd - 100*m - d) / 10000
  return ut_inv_calendar(y,m,d,z,z,z,time_units,0)
end
;==============================

; integer yyyymmdd -> double julian day
function _yyyymmdd2julian(yyyymmdd:integer)
local y, m, d, z
begin
  z = yyyymmdd
  z = 0
  d =  yyyymmdd % 100
  m = (yyyymmdd % 10000 - d) / 100
  y = (yyyymmdd - 100*m - d) / 10000
  return ut_inv_calendar(y,m,d,z,z,z,time_units,0)
end

;==============================

; string YYYYMMDD -> double julian day
function _get_julian_YYYYMMDD(YYYYMMDD:string)
begin
  return _get_julian_yyyymmdd(stringtointeger(YYYYMMDD))
end

;==============================

; double julian -> integer yyyymmddhhnn
function _julian2yyyymmddhhnn(jul:double)
local j,tmp,tmp2,str
begin
  j = jul
  j@units = time_units
  tmp2 = ut_calendar(j,0)
	tmp = floattointeger(tmp2(0,:))
  str = ""+tmp(0)+"."+sprinti("%0.2i",tmp(1))+"."+sprinti("%0.2i",tmp(2))+" "+ \
          sprinti("%0.2i",tmp(3))+":"+sprinti("%0.2i",tmp(4))
  return str
end

;==============================

; double julian -> integer yyyymmdd
function _julian2yyyymmdd(jul:double)
local j
begin
  j = jul
  j@units = time_units
  return ut_calendar(j, -2)
end

;==============================

; double julian -> string DDMMYYYY
function _julian2DDMMYYYY(jul:double)
local z, y, m, d
begin
  z = _julian2yyyymmdd(jul)
  y =  z / 10000
  m = (z % 10000) / 100
  d =  z % 100
  return sprinti("%0.2i",d)+"/"+sprinti("%0.2i",m)+"/"+y
end

;==============================

; double julian -> string DDMMM
function _julian2DDMMM(jul:double)
local z, m, d
begin
  z = _julian2yyyymmdd(jul)
  m = (z % 10000) / 100
  d =  z % 100
  return sprinti("%0.2i",d)+"-"+months_MMM(m)
end

;==============================

; double julian -> string DDMM
function _julian2DDMM(jul:double)
local z, m, d
begin
  z = _julian2yyyymmdd(jul)
  m = (z % 10000) / 100
  d =  z % 100
  return sprinti("%0.2i",d)+"/"+ sprinti("%0.2i",m)
end

;==============================

; double julian -> string YYMMDD
function _julian2YYMMDD(jul:double)
local i
begin
  i = _julian2yyyymmdd(jul)
  if i.gt.20000000 then
    i = i - 20000000
  else
    i = i - 19000000
  end if
  return sprinti("%0.6i", i)
end

;==============================

procedure _abort(str:string)
begin
  _log("ERROR: "+str)
  _log_finish
  exit
end

procedure _check(true:logical,comment:string)
begin
  if .not.true then
    _abort(comment)
  end if
end

procedure _check_file_exists(fn:string)
begin
  _check(isfilepresent(fn),"no file: "+fn)
end

procedure _check_dir_exists(dirname:string)
begin
  if .not.isfilepresent(dirname) then
    system("mkdir "+dirname)
    print("directory created: "+dirname)
  end if  
  _check_file_exists(dirname)
end

function _str2int(str:string,comment:string)
local i
begin
  i = -99999
  begin
    i = stringtointeger(str)
  end  
  _check(i.ne.-99999,str+" is not integer")
  return i
end

;==============================

function _get_num(str:string,i:integer)
local s
begin
	s = str_split(str," ")
	return s(i)
end

_GET_NUMS_DEBUG = False

function _get_nums(lines:string,n:integer)
local i,k,m,chr,ia,ib,out,tmp,tmp2
begin
	m = dimsizes(lines)
	chr = stringtocharacter(lines)
	k = dimsizes(n)
	out = new((/k,m/),float)
	ia = new(k,integer)
	ib = new(k,integer)
	ia(0) = 0
	ib(0) = n(0)-1
	do i = 1,k-1
		ia(i) = ib(i-1)+1
		ib(i) = ia(i)+n(i)-1
	end do
if _GET_NUMS_DEBUG then
	do i = 0,m-1
		print(lines(i)+" "+i)
		tmp = stringtocharacter(lines(i))
		tmp2 = stringtofloat(charactertostring(tmp(ia(1):ib(1))))
	end do
end if
	do i = 0,k-1
;print("k="+k+" i="+i+" ia="+ia(i))
		if m.gt.1 then
			out(i,:) = stringtofloat(charactertostring(chr(:,ia(i):ib(i))))
		else
			out(i,:) = stringtofloat(charactertostring(chr(ia(i):ib(i))))
		end if
	end do
	return out
end

function _get_substr(lines:string,ia:integer,ib:integer)
local i,k,m,chr,out
begin
	m = dimsizes(lines)
	chr = stringtocharacter(lines)
	k = dimsizes(ia)
	out = new((/k,m/),string)
	do i = 0,k-1
		out(i,:) = charactertostring(chr(:,ia(i):ib(i)))
	end do
	return out
end

;==============================


procedure _get_scale(data:numeric,pa:float,pb:float,nstep:integer,a:float,b:float,spacing:float,stride:integer)
local v,q,k,n,va,vb,r,l10,flr,step
begin
  v = ndtooned(data)
;  v = where(ismissing(v),1e10,v)
  qsort(v)
  n = num(.not.ismissing(v))
  k = round(pa*n,3)
  va = v(k)
  k = round(pb*n,3)
  vb = v(k)
  delete(v)
  r = (vb-va)/nstep
  l10 = log(10.)
  lr = log10(r)
  flr = floor(lr)
  lr = lr-flr
  flr = exp(flr*l10)
  stride = 5
;print(""+dimsizes(data)+TAB+n+TAB+va+TAB+vb)
  if lr.gt.0.8 then
    spacing = flr*10.
  else
  if lr.gt.0.5 then
    spacing = flr*5.
    stride = 4
  else    
  if lr.gt.0.15 then
    spacing = flr*2.
  else    
    spacing = flr
  end if
  end if
  end if
  step = spacing*stride
  a = step*floor(va/step)
  b = step*ceil(vb/step)
end

procedure _get_log_scale(data:numeric,pa:float,pb:float,nstep:integer,a:float,b:float)
local tmp,v,q,k,n,va,vb,r,l10,flr,step
begin
	tmp = log10(data)
  v = ndtooned(tmp)
  delete(tmp)
  qsort(v)
  n = num(.not.ismissing(v))
  k = round(pa*n,3)
  va = v(k)
  k = round(pb*n,3)
  vb = v(k)
  delete(v)
  r = (vb-va)/nstep
  l10 = log(10.)
  lr = log10(r)
  flr = floor(lr)
  lr = lr-flr
  flr = exp(flr*l10)
  stride = 5
;print(""+dimsizes(data)+TAB+n+TAB+va+TAB+vb)
  if lr.gt.0.8 then
    spacing = flr*10.
  else
  if lr.gt.0.5 then
    spacing = flr*5.
    stride = 4
  else    
  if lr.gt.0.15 then
    spacing = flr*2.
  else    
    spacing = flr
  end if
  end if
  end if
  step = spacing*stride
  a = step*floor(va/step)
  b = step*ceil(vb/step)
end

