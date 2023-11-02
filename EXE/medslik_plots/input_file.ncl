;;;;;;;;;;;;;;;; parsing medslik_inputfile.txt ;;;;;;;;;;;;;;;;;

SLICKS_MAX_LAT = LAT
SLICKS_MAX_LON = LON
SLICKS_MIN_LAT = LAT
SLICKS_MIN_LON = LON

fn = DIR+TXT
_check_file_exists(fn)
lines = asciiread(fn,-1,"string")
nline = dimsizes(lines)

do iline = 0,nline-1
	if isStrSubset(lines(iline),"ContourSlick=YES") then
		INI_TYPE  = INI_POLY
	end if
	if (INI_TYPE.eq.INI_POLY).and.isStrSubset(lines(iline),"NSlick=") then
		tmp1 = str_split(lines(iline)," ")
		tmp2 = str_split(tmp1(0),"=")
		NSLICK = stringtointeger(tmp2(1))
		delete(tmp1)
		delete(tmp2)
		SLICKS_LAT = new((/NSLICK,nline/),float)
		SLICKS_LON = new((/NSLICK,nline/),float)
		SLICKS_N = new(NSLICK,integer)
	end if
	if (NSLICK.gt.0).and.isStrSubset(lines(iline),"]=") then
		tmp1 = str_sub_str(lines(iline),"[","=")
		tmp2 = str_sub_str(tmp1,"]","=")
		items = str_split(tmp2,"=")
		delete(tmp1)
		delete(tmp2)
		i = stringtointeger(items(1))
		tmp = stringtofloat(items(2))
		chr = stringtocharacter(items(0))
		n = dimsizes(chr)
		j = stringtointeger(charactertostring(chr(1:n-4)))-1
		str = charactertostring(chr(n-4:n-1))
		delete(items)
		delete(chr)
		SLICKS_N(j) = max((/SLICKS_N(j),i/))
		if str.eq."lat" then
			SLICKS_LAT(j,i) = tmp
		else
			SLICKS_LON(j,i) = tmp
		end if
		delete(tmp)
	end if
end do
delete(lines)

if NSLICK.gt.0 then
	do i = 0,NSLICK-1
		n = SLICKS_N(i)
		_check( (.not.any(ismissing(SLICKS_LAT(i,1:n)))).and. \
						(.not.any(ismissing(SLICKS_LON(i,1:n)))).and. \	
						(n.ge.3), TXT+", wrong structure of polygon "+i)
		SLICKS_LAT(i,0) = SLICKS_LAT(i,n)
		SLICKS_LON(i,0) = SLICKS_LON(i,n)
	end do
	SLICKS_MAX_LAT = max(SLICKS_LAT)
	SLICKS_MAX_LON = max(SLICKS_LON)
	SLICKS_MIN_LAT = min(SLICKS_LAT)
	SLICKS_MIN_LON = min(SLICKS_LON)
end if
