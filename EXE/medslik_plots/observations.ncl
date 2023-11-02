OBSMAXOBJ   = 100
OBSMAXPOINT = 10000

OBSPOINT = 0
OBSPOLY  = 1

NOBS = 0

OBS_MAX_N = 0

tmp = systemfunc("ls "+DIR+OBS+"??????????.txt 2> /dev/null")
if .not.ismissing(tmp(0)) then
	NOBS = dimsizes(tmp)
	_log("observation files found: "+NOBS)
	OBS_N     = new(NOBS,integer)
	OBS_YMDHN = new(NOBS,integer)
	OBS_JUL   = new(NOBS,double)
	OBS_IND   = new(NOBS,integer)
	OBS_NM    = new((/NOBS,OBSMAXOBJ/),integer)
	OBS_TYPE  = new((/NOBS,OBSMAXOBJ/),integer)
	OBS_LAT   = new((/NOBS,OBSMAXOBJ,OBSMAXPOINT/),float)
	OBS_LON   = new((/NOBS,OBSMAXOBJ,OBSMAXPOINT/),float)
	OBS_N  = 0
	OBS_NM = 0
	do iobs = 0,NOBS-1
		ffn = tmp(iobs)
		n = strlen(ffn)
		fn = str_get_cols(ffn,n-26,n-1)
		str = str_get_cols(ffn,n-14,n-5)
		k = stringtointeger(str)
		OBS_YMDHN(iobs) = k
		OBS_JUL(iobs) = ut_inv_calendar(2000+k/100000000,k/1000000%100,k/10000%100,k/100%100,k%100,0,time_units,0)
		OBS_IND(iobs) = closest_val(OBS_JUL(iobs),JULIANS)
		lines = asciiread(ffn,-1,"string")
		nline = dimsizes(lines)
		do iline = 0,nline-1
			str = str_lower(str_strip(lines(iline)))
			if isStrSubset(str,"lat").or.isStrSubset(str,"/").or.(strlen(str).le.1) then 
				continue
			end if
			if isStrSubset(str,"polygon") then
				i = OBS_N(iobs)
				OBS_N(iobs) = i+1
				OBS_TYPE(iobs,i) = OBSPOLY
				continue
			end if
			if isStrSubset(str,"point") then
				i = OBS_N(iobs)
				OBS_N(iobs) = i+1
				OBS_TYPE(iobs,i) = OBSPOINT
				continue
			end if
			items = str_strip(str_split(str," "))
			_check(dimsizes(items).ge.2,fn+LBR+"line "+(iline+1)+LBR+str+LBR+"must contain lat/lon")
			j = OBS_NM(iobs,i)
			OBS_NM(iobs,i) = j+1
			OBS_LAT(iobs,i,j) = stringtofloat(items(0))
			OBS_LON(iobs,i,j) = stringtofloat(items(1))
			delete(items)
		end do
		delete(lines)
		n = OBS_N(iobs)
		_check(n.gt.0,fn+LBR+"no valid data")
		do i = 0,n-1
			if OBS_TYPE(iobs,i).eq.OBSPOINT then
				_check(OBS_NM(iobs,i).gt.0,fn+LBR+"object Points must have at least one point")
				OBS_MAX_N = max((/OBS_MAX_N,OBS_NM(iobs,i)/))
			end if
			if OBS_TYPE(iobs,i).eq.OBSPOLY then
				_check(OBS_NM(iobs,i).gt.2,fn+LBR+"object Polygon must have at least three points")
				j = OBS_NM(iobs,i)
				OBS_LAT(iobs,i,j) = OBS_LAT(iobs,i,0)
				OBS_LON(iobs,i,j) = OBS_LON(iobs,i,0)
				OBS_NM(iobs,i) = j+1
			end if
		end do
		_log(""+fn)
	end do
	OBS_MAX_LAT = max(OBS_LAT)
	OBS_MAX_LON = max(OBS_LON)
	OBS_MIN_LAT = min(OBS_LAT)
	OBS_MIN_LON = min(OBS_LON)
	OBS_IND = where(OBS_IND.le.0,1,OBS_IND)
	OBS_IND = where(OBS_IND.gt.NTIME,NTIME,OBS_IND)
	print("")
end if
delete(tmp)
