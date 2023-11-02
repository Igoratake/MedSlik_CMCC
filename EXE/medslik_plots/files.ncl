;;;;;;;;;;;;;;;; check extra files ;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;; parsing medslik5.inp ;;;;;;;;;;;;;;;;;

fn = DIR+INP
_check_file_exists(fn)
lines = asciiread(fn,-1,"string")

items = str_split(lines(3)," ")
DAY   = stringtointeger(items(0))
MONTH = stringtointeger(items(1))
YEAR  = stringtointeger(items(2))
delete(items)

items = stringtointeger(_get_num(lines(4),0))
HOUR   = items/100
MINUTE = items%100
delete(items)

TIMEA = ut_inv_calendar(YEAR,MONTH,DAY,HOUR,MINUTE,0,time_units,0)

items = str_split(lines(6)," ")
LAT = stringtofloat(items(0))+stringtofloat(items(1))/60.
delete(items)

items = str_split(lines(7)," ")
LON = stringtofloat(items(0))+stringtofloat(items(1))/60.
delete(items)

PREFIX   = str_right_strip(_get_substr(lines(8),0,10))

NTIME = stringtointeger(_get_num(lines(9),0))

PIXEL = stringtofloat(_get_num(lines(29),0))

delete(lines)

TIMEB = TIMEA+(NTIME/24.)
JULIANS = fspan(TIMEA,TIMEB,NTIME+1)
JULIANS@units = time_units

TITLE = _julian2DDMMYYYY(TIMEA)+" "+sprinti("%0.2i",HOUR)+":"+sprinti("%0.2i",MINUTE)+" UTC"

print("oil files available = "+NTIME)
print("zero hour 0000 = "+_julian2yyyymmddhhnn(TIMEA))
print("last hour "+sprinti("%0.4i",NTIME)+" = "+_julian2yyyymmddhhnn(JULIANS(1)))
print("")

;_log("file medslik5.inp OK")

;;;;;;;;;;;;;;;; parsing medslik.fte ;;;;;;;;;;;;;;;;;

if False then

fn = DIR+FTE
_check_file_exists(fn)
lines = asciiread(fn,-1,"string")
;tmp2 = str_split(tmp(1)," ")
;LAT   = stringtofloat(tmp2(2))
;LON   = stringtofloat(tmp2(4))
;delete(tmp2)
NTIME = stringtointeger(_get_num(lines(4),3))
FTE_DATA = _get_nums(lines(8:7+NTIME*2),(/9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9/))
delete(lines)
_log("file medslik.fte OK")
;print("")

end if

;;;;;;;;;;;;;;;; check oil files exist ;;;;;;;;;;;;;;;;;

MAP_FG = new(NTIME+1,logical)
MAP_FG = False
MAP_FN = new(NTIME+1,string)
MAP_FN = sprinti("%0.4i",ispan(0,NTIME,1))
do i = HOUR_BEGIN,NTIME,HOUR_STEP
	if i.eq.0 then
		MAP_FN(i) = "initial"
	end if
;	MAP_FG(i) = ((HOUR_END.le.0).or.(i.le.HOUR_END)).and.(i.gt.0)
	MAP_FG(i) = ((HOUR_END.le.0).or.(i.le.HOUR_END)).and.(i.ge.0)
	if MAP_FG(i) then
		if i.gt.0 then
			_check_file_exists(DIR+PREFIX+MAP_FN(i)+SRF)
			_check_file_exists(DIR+PREFIX+MAP_FN(i)+CST)
		;	_check_file_exists(DIR+PREFIX+MAP_FN(i)+DSP)
		end if
	end if
end do

ntime = num(MAP_FG)
times = new(ntime,integer)
times = ind(MAP_FG)

if ntime.eq.0 then
	print("no hours to process")
	print("**** procedure aborted ****")
	exit
end if

if (ntime.eq.1).and.(times(0).eq.0) then
	print("no hours to process")
	print("**** procedure aborted ****")
	exit
end if

;print("oil files available = "+ntime)
;print("oil files available = "+NTIME)
;print("latitude  = "+LAT)
;print("longitude = "+LON)
;print("zero hour 0000 = "+_julian2yyyymmddhhnn(TIMEA))
;print("hours    = "+NTIME)
;print("last hour "+sprinti("%0.4i",NTIME)+" = "+_julian2yyyymmddhhnn(JULIANS(1)))
;print("hour"+sprinti("%0.4i",NTIME)+" = "+_julian2yyyymmddhhnn(JULIANS(NTIME)))
;print("")

;str = "hours to process: from "+HOUR_BEGIN+" to "
;print("")

;;;;;;;;;;;;;;;; check currents files exist ;;;;;;;;;;;;;;;;;

if CURRENTS then

	juld = 24.d
	MODEL_FNA = new(ntime,string)
	MODEL_FNB = new(ntime,string)
	MODEL_WEIGHTS = new(ntime,float)
	jula = floor(JULIANS(times)*juld)/juld
	julb = jula+(1.d/24.d)
	jula@units = time_units
	julb@units = time_units
	MODEL_WEIGHTS = doubletofloat((julb-JULIANS(times))*juld)
;do i=0,NTIME
;print(""+_julian2yymmddhh(JULIANS(i))+TAB+_julian2yymmddhh(jula(i)) \
;+TAB+_julian2yymmddhh(julb(i))+TAB+MODEL_WEIGHTS(i))
;end do
	MODEL_FNA = MODEL_DIR+_julian2yymmddhh(jula)+MODEL_EXT
	MODEL_FNB = MODEL_DIR+_julian2yymmddhh(julb)+MODEL_EXT
;do i=0,ntime-1
;print(MAP_FN(times(i))+TAB+_julian2yymmddhh(JULIANS(times(i)))+TAB+_julian2yymmddhh(jula(i)) \
;+TAB+_julian2yymmddhh(julb(i))+TAB+MODEL_WEIGHTS(i))
;end do
;print(MODEL_FNA+TAB+_julian2yyyymmddhhnn(JULIANS(times)))
	delete(jula)
	delete(julb)
	delete(juld)	
	
	MODEL_FNAK = new(ntime,integer)
	MODEL_FNBK = new(ntime,integer)
	MODEL_FN   = new(2*ntime,string)
	MODEL_FNK = 0
	
	itime = 0
	
	function _model_files(fn:string)
	local inds,k
	begin
		inds = ind(MODEL_FN.eq.fn)
		if ismissing(inds(0)) then
			if .not.isfilepresent(fn) then
;				print("no currents file for hour"+MAP_FN(times(itime)))
;				print(""+_julian2yyyymmddhhnn(JULIANS(times(itime))))
;				print(""+fn)
;				print("")
				return -1
;;				print("**** procedure aborted ****")
;;				exit
			end if
			MODEL_FN(MODEL_FNK) = fn
			k = MODEL_FNK
			MODEL_FNK = MODEL_FNK+1
		else 
			k = inds(0)
		end if
		return k
	end
	
	NO_CURRENTS_DATA = False
	
	do itime = 0,ntime-1
		MODEL_FNAK(itime) = _model_files(MODEL_FNA(itime))
		MODEL_FNBK(itime) = _model_files(MODEL_FNB(itime))
		if (MODEL_FNAK(itime).lt.0).or.(MODEL_FNBK(itime).lt.0) then
			NO_CURRENTS_DATA = True
			MAP_FG(times(itime)) = False
			print("no model data for hour "+MAP_FN(times(itime)))
;			print("")
		end if
	end do
		
	if NO_CURRENTS_DATA then
		ntime = num(MAP_FG)
		tmp = times(0:ntime-1)
		delete(times)
		times = tmp
		delete(tmp)
		print("")
		print("last hour to process is changed to "+MAP_FN(times(ntime-1)))
		print("")
	end if

	if ntime.eq.0 then
		print("no hours to process")
		print("**** procedure aborted ****")
		exit
	end if

	lines = asciiread(MODEL_FNA(0),-1,"string")
	items = str_split(lines(2)," ")
	MODEL_LONA = stringtofloat(items(0))
	MODEL_LONB = stringtofloat(items(1))
	MODEL_LATA = stringtofloat(items(2))
	MODEL_LATB = stringtofloat(items(3))
	MODEL_LONN = stringtointeger(items(4))
	MODEL_LATN = stringtointeger(items(5))
	MODEL_LOND = (MODEL_LONB-MODEL_LONA)/(MODEL_LONN-1)
	MODEL_LATD = (MODEL_LATB-MODEL_LATA)/(MODEL_LATN-1)
	delete(items)
	items = str_split(lines(3)," ")
	MODEL_N = stringtointeger(items(0))
	delete(items)
	items = _get_nums(lines(5:MODEL_N+4),(/11,11/))
	MODEL_LAT = items(0,:)
	MODEL_LON = items(1,:)
	delete(items)
	delete(lines)
	
	if VERBOSE then
		print("currents files to process = "+MODEL_FNK)
		print("latitude :  min = "+sprintf("%7.4f",MODEL_LATA)+"  max = "+sprintf("%7.4f",MODEL_LATB)+ \
		"  range = "+sprintf("%6.4f",MODEL_LATB-MODEL_LATA)+ \
		"  step = "+sprintf("%6.4f",MODEL_LATD)+"  n = "+MODEL_LATN)
		print("longitude:  min = "+sprintf("%7.4f",MODEL_LONA)+"  max = "+sprintf("%7.4f",MODEL_LONB)+ \
		"  range = "+sprintf("%6.4f",MODEL_LONB-MODEL_LONA)+ \
		"  step = "+sprintf("%6.4f",MODEL_LOND)+"  n = "+MODEL_LONN)
		print("")
	end if
	
end if

;;;;;;;;;;;;;;;; check wind files exist ;;;;;;;;;;;;;;;;;

;do i=0,NTIME
;print(""+i+TAB+MAP_FG(i))
;end do

if WIND then

	juld = 4.d
	WIND_FNA = new(ntime,string)
	WIND_FNB = new(ntime,string)
	WIND_IA = new(ntime,integer)
	WIND_IB = new(ntime,integer)
	WIND_WEIGHTS = new(ntime,float)
	jula = floor(JULIANS(times)*juld)/juld
	julb = jula+(1.d/4.d)
	jula@units = time_units
	julb@units = time_units
	WIND_WEIGHTS = doubletofloat((julb-JULIANS(times))*juld)
	WIND_IA = ut_calendar(jula,-3)%100/6
	WIND_IB = ut_calendar(julb,-3)%100/6
	WIND_FNA = WIND_DIR+ECM+"_"+_julian2YYMMDD(jula)+"."+ECM
	WIND_FNB = WIND_DIR+ECM+"_"+_julian2YYMMDD(julb)+"."+ECM
;print(""+_julian2YYMMDD(jula)+TAB+WIND_IA+TAB+_julian2YYMMDD(julb)+TAB+WIND_IB+TAB+ \
;					WIND_WEIGHTS+TAB+_julian2yymmddhh(JULIANS(times)))
	delete(jula)
	delete(julb)
	delete(juld)
	
	WIND_FNAK = new(ntime,integer)
	WIND_FNBK = new(ntime,integer)
	WIND_FN   = new(2*ntime,string)
	WIND_FNK = 0
	
	function _wind_files(fn:string)
	local inds,k
	begin
;print(fn+TAB+times(itime)+TAB+MAP_FN(times(itime)))
		inds = ind(WIND_FN.eq.fn)
		if ismissing(inds(0)) then
			if .not.isfilepresent(fn) then
;				print("no wind file for hour "+MAP_FN(times(itime)) + \
;				", "+_julian2yyyymmddhhnn(JULIANS(times(itime))))
;				print(""+fn)
;				print("")
				return -1
;				print("**** procedure aborted ****")
;				exit
			end if
			WIND_FN(WIND_FNK) = fn
			k = WIND_FNK
			WIND_FNK = WIND_FNK+1
		else 
			k = inds(0)
		end if
		return k
	end
	
	NO_WIND_DATA = False

	do itime = 0,ntime-1
		WIND_FNAK(itime) = _wind_files(WIND_FNA(itime))
		WIND_FNBK(itime) = _wind_files(WIND_FNB(itime))
		if (WIND_FNAK(itime).lt.0).or.(WIND_FNBK(itime).lt.0) then
			NO_WIND_DATA = True
			MAP_FG(times(itime)) = False
			print("no wind data for hour "+MAP_FN(times(itime)))
;			print("")
		end if
	end do
	
	if NO_WIND_DATA then
		ntime = num(MAP_FG)
		tmp = times(0:ntime-1)
		delete(times)
		times = tmp
		delete(tmp)
		print("")
		print("last hour to process is changed to "+MAP_FN(times(ntime-1)))
		print("")
	end if
	
;do i=0,NTIME
;print(""+i+TAB+MAP_FG(i))
;end do
;exit

	if ntime.eq.0 then
		print("no hours to process")
		print("**** procedure aborted ****")
		exit
	end if

	tmp1 = asciiread(WIND_FNA(0),-1,"string")
	tmp2 = str_split(tmp1(2)," ")
	WIND_LONA = stringtofloat(tmp2(0))
	WIND_LONB = stringtofloat(tmp2(1))
	WIND_LATA = stringtofloat(tmp2(2))
	WIND_LATB = stringtofloat(tmp2(3))
	WIND_LONN = stringtointeger(tmp2(4))
	WIND_LATN = stringtointeger(tmp2(5))
	WIND_LOND = (WIND_LONB-WIND_LONA)/(WIND_LONN-1)
	WIND_LATD = (WIND_LATB-WIND_LATA)/(WIND_LATN-1)
	delete(tmp2)
	tmp2 = str_split(tmp1(3)," ")
	WIND_N = stringtointeger(tmp2(0))
	delete(tmp2)
	tmp2 = _get_nums(tmp1(5:WIND_N+4),(/11,11/))
	WIND_LAT = tmp2(0,:)
	WIND_LON = tmp2(1,:)
	delete(tmp2)
	delete(tmp1)
	
	if VERBOSE then
		print("wind files to process = "+WIND_FNK)
		print("latitude :  min = "+sprintf("%7.4f",WIND_LATA)+"  max = "+sprintf("%7.4f",WIND_LATB)+ \
		"  range = "+sprintf("%6.4f",WIND_LATB-WIND_LATA)+"  step = "+WIND_LATD+"  n = "+WIND_LATN)
		print("longitude:  min = "+sprintf("%7.4f",WIND_LONA)+"  max = "+sprintf("%7.4f",WIND_LONB)+ \
		"  range = "+sprintf("%6.4f",WIND_LONB-WIND_LONA)+"  step = "+WIND_LOND+"  n = "+WIND_LONN)
		print("")
	end if

end if






if False.and.CURRENTS then

	juld = 24.d
	MODEL_FNA = new(NTIME+1,string)
	MODEL_FNB = new(NTIME+1,string)
	MODEL_WEIGHTS = new(NTIME+1,float)
	jula = floor(JULIANS*juld)/juld
	julb = jula+(1.d/24.d)
	jula@units = time_units
	julb@units = time_units
	MODEL_WEIGHTS = doubletofloat((julb-JULIANS)*juld)
;do i=0,NTIME
;print(""+_julian2yymmddhh(JULIANS(i))+TAB+_julian2yymmddhh(jula(i)) \
;+TAB+_julian2yymmddhh(julb(i))+TAB+MODEL_WEIGHTS(i))
;end do
	MODEL_FNA = MODEL_DIR+_julian2yymmddhh(jula)+MODEL_EXT
	MODEL_FNB = MODEL_DIR+_julian2yymmddhh(julb)+MODEL_EXT
	delete(jula)
	delete(julb)
	delete(juld)
		
	do jtime = 0,times-1
		itime = times(jtime)
		if (.not.isfilepresent(MODEL_FNA(itime))).or. \
			 (.not.isfilepresent(MODEL_FNB(itime))) then
			print("no currents file for hour"+MAP_FN(itime))
			print(""+_julian2yyyymmddhhnn(JULIANS(itime)))
			if .not.isfilepresent(MODEL_FNB(itime)) then
				print(""+MODEL_FNA(itime))
			end if
			if .not.isfilepresent(MODEL_FNB(itime)) then
				print(""+MODEL_FNB(itime))
			end if
			print("**** procedure aborted ****")
			_log_finish
			exit
		end if
	end do

end if

;;;;;;;;;;;;;;;; check wind files exist ;;;;;;;;;;;;;;;;;

if False.and.WIND then

	juld = 4.d
	WIND_FNA = new(NTIME+1,string)
	WIND_FNA = new(NTIME+1,string)
	WIND_IA = new(NTIME+1,integer)
	WIND_IB = new(NTIME+1,integer)
	WIND_WEIGHTS = new(NTIME+1,float)
	jula = floor(JULIANS*juld)/juld
	julb = jula+(1.d/4.d)
	jula@units = time_units
	julb@units = time_units
	WIND_WEIGHTS = doubletofloat((julb-JULIANS)*juld)
print("")
do i=0,NTIME
print(""+_julian2yymmddhh(JULIANS(i))+TAB+_julian2yymmddhh(jula(i)) \
+TAB+_julian2yymmddhh(julb(i))+TAB+WIND_WEIGHTS(i))
end do
	MODEL_FNA = MODEL_DIR+_julian2yymmddhh(jula)+MODEL_EXT
	MODEL_FNB = MODEL_DIR+_julian2yymmddhh(julb)+MODEL_EXT
	delete(jula)
	delete(julb)
	delete(juld)	
	do itime = 0,NTIME
		if MAP_FG(itime) then
			if (.not.isfilepresent(MODEL_FNA(itime))).or. \
				 (.not.isfilepresent(MODEL_FNB(itime))) then
				print("no wind file for hour"+MAP_FN(itime))
				print(""+_julian2yyyymmddhhnn(JULIANS(itime)))
				if .not.isfilepresent(MODEL_FNB(itime)) then
					print(""+MODEL_FNA(itime))
				end if
				if .not.isfilepresent(MODEL_FNB(itime)) then
					print(""+MODEL_FNB(itime))
				end if
				print("**** procedure aborted ****")
				_log_finish
				exit
			end if
		end if
	end do

end if

