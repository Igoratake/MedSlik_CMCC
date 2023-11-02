if NOBS.gt.0 then

obs_map = new((/NOBS,oil_latn,oil_lonn/),float)
obs_map!1 = "lat"
obs_map!2 = "lon"
obs_map&lat = oil_lats
obs_map&lon = oil_lons

obs_oil_n = new(NOBS,integer)
obs_oil_n = 0
obs_oil_val = new((/NOBS,MAX_OIL_N/),float)
obs_oil_lon = new((/NOBS,MAX_OIL_N/),float)
obs_oil_lat = new((/NOBS,MAX_OIL_N/),float)

do iobs = 0,NOBS-1
	i = OBS_IND(iobs)
	k = ind(times.eq.i)
	if .not.ismissing(k(0)) then
		obs_map(iobs,:,:) = oil_map(k,:,:)
		n = oil_n(k)
		obs_oil_n(iobs) = n
		obs_oil_val(iobs,0:n-1) = _oil_val(k,0:n-1)
		obs_oil_lon(iobs,0:n-1) = _oil_lon(k,0:n-1)
		obs_oil_lat(iobs,0:n-1) = _oil_lat(k,0:n-1)
	else
		if (OBS_IND(iobs).gt.0).and.(OBS_IND(iobs).le.NTIME) then
			fn = DIR+PREFIX+MAP_FN(OBS_IND(iobs))+SRF
			lines = asciiread(fn,-1,"string")
			m = dimsizes(lines)
			n = stringtointeger(_get_num(lines(7),0))
			obs_oil_n(iobs) = n
			items = _get_nums(lines(9:n+8),(/9,9,20/))
			obs_oil_val(iobs,0:n-1) = items(2,:)
			obs_oil_lon(iobs,0:n-1) = items(1,:)
			obs_oil_lat(iobs,0:n-1) = items(0,:)
			ii = floattointeger((items(1,:)-oil_lona)/oil_dx)
			jj = floattointeger((items(0,:)-oil_lata)/oil_dy)
			do i = 0,n-1
				obs_map(iobs,jj(i),ii(i)) = (/items(2,i)/)
			end do
			delete(ii)
			delete(jj)
			delete(items)
			delete(lines)
		end if
	end if
end do

end if


