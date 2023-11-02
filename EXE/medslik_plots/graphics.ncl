COLOR_END    = -4
COLOR_START  =  3
PALETTE = "BlAqGrYeOrReVi200"
NCOLOR = 202

NSCALE = dimsizes(oil_scale)
OIL_SCALE = new(NSCALE+2,float)
OIL_SCALE(0) = -99.
OIL_SCALE(NSCALE+1) = 1e10
OIL_SCALE(1:NSCALE) = oil_scale

CENTER_SIZE  = 0.005

OBS_SIZE     = 0.005
OBS_SIZE1    = 0.005
OBS_SIZE10   = 0.002
OBS_SIZE100  = 0.001
OBS_SIZE1000 = 0.0005

if OBS_MAX_N.ge.1 then
	OBS_SIZE = OBS_SIZE1
end if
if OBS_MAX_N.ge.10 then
	OBS_SIZE = OBS_SIZE10
end if
if OBS_MAX_N.ge.100 then
	OBS_SIZE = OBS_SIZE100
end if
if OBS_MAX_N.ge.1000 then
	OBS_SIZE = OBS_SIZE1000
end if

dummy = new((/3,3/),float)
dummy(0,:) = 0.1
dummy(1,:) = 10.
dummy(2,:) = 1000.
dummy(0,:) = 0.001*(vb-va)/2.
dummy(1,:) = (vb-va)/2.
dummy(2,:) = 1000*(vb-va)/2.
dummy!0 = "lat"
dummy!1 = "lon"
dummy&lat = (/0.,1.,2./)
dummy&lon = (/0.,1.,2./)

if COAST_MAP.and.(REAL_CST_N.gt.0) then
	NSCALE_CST = dimsizes(cst_scale)
	CST_SCALE = new(NSCALE_CST+2,float)
	CST_SCALE(0) = -99.
	CST_SCALE(NSCALE_CST+1) = 1e10
	CST_SCALE(1:NSCALE_CST) = cst_scale
	dummy_cst = new((/3,3/),float)
	dummy_cst(0,:) = 0.1
	dummy_cst(1,:) = 10.
	dummy_cst(2,:) = 1000.
	dummy_cst(0,:) = 0.001*(vb_cst-va_cst)/2.
	dummy_cst(1,:) = (vb_cst-va_cst)/2.
	dummy_cst(2,:) = 1000*(vb_cst-va_cst)/2.
	dummy_cst!0 = "lat"
	dummy_cst!1 = "lon"
	dummy_cst&lat = (/0.,1.,2./)
	dummy_cst&lon = (/0.,1.,2./)
end if


function _get_res_base()
local res
begin

	res = True

;  res@gsnMaximize  = True
;  gsnPaperOrientation = "landscape"
	res@gsnFrame = False
	res@gsnDraw  = False
  res@gsnRightString = ""
  res@gsnLeftString  = ""
  res@gsnLeftStringFontHeightF  = FONT_TITLE
  res@gsnRightStringFontHeightF = FONT_TITLE

	res@mpProjection       = "Mercator"
	res@mpOutlineOn        = True
	res@mpFillOn           = True
	res@mpDataBaseVersion  = "HighRes"
	res@mpDataResolution   = "Finest"
	res@mpLandFillColor    = COLOR_LAND
;	res@mpLandFillColor    = -1
  res@mpLimitMode        = "LatLon"
  res@mpFillDrawOrder    = "PreDraw"
  res@mpOutlineDrawOrder = "PreDraw"

  res@mpMinLatF = min_map_lat
  res@mpMaxLatF = max_map_lat
  res@mpMinLonF = min_map_lon
  res@mpMaxLonF = max_map_lon

;print("res@mpMinLatF = "+min_map_lat)
;print("res@mpMaxLatF = "+max_map_lat)
;print("res@mpMinLonF = "+min_map_lon)
;print("res@mpMaxLonF = "+max_map_lon)

  res@mpGridAndLimbOn        = GRID
  res@mpGridLineThicknessF   = GRID_THICKNESS
  res@mpGridLineColor        = COLOR_GRID
  if isvar("GRID_LON_SPACING") then
		res@mpGridLonSpacingF      = GRID_LON_SPACING/60.
	end if
  if isvar("GRID_LAT_SPACING") then
		res@mpGridLatSpacingF      = GRID_LAT_SPACING/60.
	end if

  res@tmXTOn = True
  res@tmYROn = True
  res@pmTickMarkDisplayMode  = "Always"
  res@tmXBMajorLengthF       = 0.002
  res@tmYLMajorLengthF       = 0.002
  res@tmXBMinorLengthF       = 0.001
  res@tmYLMinorLengthF       = 0.001
  res@tmXTMajorLengthF       = 0.002
  res@tmYRMajorLengthF       = 0.002
  res@tmXTMinorLengthF       = 0.001
  res@tmYRMinorLengthF       = 0.001
  res@tmXBLabelFontHeightF   = FONT_AXIS
  res@tmYLLabelFontHeightF   = FONT_AXIS
  res@tmBorderThicknessF     = 1.0

	return res

end

function _get_res_overlay()
local res
begin
	res = True
	res@gsnFrame  = False
	res@gsnDraw   = False
	res@gsnAddCyclic = False
  res@gsnTickMarksOn = False
  res@gsnRightString = ""
  res@gsnLeftString  = ""
  res@tmBorderThicknessF = 0.
	return res
end

function _get_res_vector()
local res
begin
	res = _get_res_overlay
  res@vcRefAnnoPerimOn     = False
  res@vcPositionMode       = "ArrowTail"
  res@vcRefLengthF         = VECTOR_ANNOTATION_LENGTH
  res@vcRefAnnoFontHeightF = FONT_VECTOR_ANNOTATION
	return res
end

function _get_res_model()
local res
begin
	res = _get_res_vector
  res@vcLineArrowColor = COLOR_CURRENTS
  res@vcLineArrowThicknessF = 0.5
  if CURRENTS_CENTER_ONLY then
  	res@vcRefMagnitudeF = model_center_scale
	  res@vcLineArrowThicknessF   = ARROW_CENTER_CURRENTS_WIDTH
	  res@vcLineArrowHeadMaxSizeF   = ARROW_CENTER_CURRENTS_HEAD
  else
  	res@vcRefMagnitudeF         = model_main_scale
	  res@vcLineArrowThicknessF   = ARROW_CURRENTS_WIDTH
	  res@vcLineArrowHeadMaxSizeF = ARROW_CURRENTS_HEAD
	  if isvar("CURLY_CURRENTS").and.CURLY_CURRENTS then
	  	res@vcGlyphStyle  = "CurlyVector"
	  end if
  res@vcMinDistanceF = CURRENTS_MINIMUM_DISTANCE
  end if
  res@vcRefAnnoString1 = ""+res@vcRefMagnitudeF+" m/s"
  res@vcRefAnnoString2 = "Currents"
	return res
end

function _get_res_wind()
local res
begin
	res = _get_res_vector
  res@vcLineArrowColor = COLOR_WIND
  if WIND_CENTER_ONLY then
  	res@vcRefMagnitudeF         = wind_center_scale
	  res@vcLineArrowThicknessF   = ARROW_CENTER_WIND_WIDTH
	  res@vcLineArrowHeadMaxSizeF = ARROW_CENTER_WIND_HEAD
  else
  	res@vcRefMagnitudeF         = wind_main_scale
	  res@vcLineArrowThicknessF   = ARROW_WIND_WIDTH
	  res@vcLineArrowHeadMaxSizeF = ARROW_WIND_HEAD
;;  res@vcMinDistanceF = WIND_MINIMUM_DISTANCE
  end if
  res@vcRefAnnoString1 = ""+res@vcRefMagnitudeF+" m/s"
  res@vcRefAnnoString2 = "Wind"
  res@vcRefAnnoParallelPosF = 0.1
	return res
end

function _get_res_oil()
local res
begin
	res = _get_res_overlay
	res@gsnSpreadColors      = True
	res@gsnSpreadColorEnd    = COLOR_END
	res@gsnSpreadColorStart  = COLOR_START
  res@pmLabelBarSide   = "Right"
  res@pmLabelBarWidthF = 0.07
  res@lbOrientation      = "Vertical"
  res@lbLabelJust        = "CenterLeft"
  res@lbLabelFontHeightF = FONT_COLORBAR
  res@lbLabelStride      = 1
	res@lbLabelStrings     = oil_labels
  res@cnLevelSelectionMode       = "ExplicitLevels"
  res@cnExplicitLabelBarLabelsOn = True
  res@cnLevels                   = oil_scale
	res@cnFillOn              = True
	res@cnLinesOn             = False
	res@cnLineLabelsOn        = False
  res@cnFillMode            = "CellFill"
;  res@cnFillMode            = "AreaFill"
;  res@cnFillMode            = "RasterFill"
;  res@cnRasterSmoothingOn = True
  res@cnMissingValFillColor = -1
  res@cnFillDrawOrder = "PreDraw"
	return res
end

function _get_text_res()
local res
begin
	res = True
	res@txFontHeightF = FONT_VECTOR_TEXT
	res@txJust = "TopLeft"
	return res
end

procedure _plot(fn:string)
local ffn
begin
	ffn = PLOT_DIR+fn
  system("convert -density "+PLOT_DPI+" -trim "+ffn+".ps "+ \
  "+repage -bordercolor White -border 5x5 "+ffn+".png")
  system("rm "+ffn+".ps")
end

procedure _plot_slick(wks:graphic,plot:graphic)
local i,res
begin
	res = True
	if NSLICK.eq.0 then
		res@gsMarkerIndex = 16
		res@gsMarkerSizeF = CENTER_SIZE
		gsn_polymarker(wks,plot,LON,LAT,res)
	else
		do i = 0,NSLICK-1
			n = SLICKS_N(i)
			gsn_polygon(wks,plot,SLICKS_LON(i,0:n),SLICKS_LAT(i,0:n),res)
		end do
	end if
end

procedure _plot_main_map(itime:integer)
local fn,ffn,wks,i,j,k,l,w,c, \
			res_base,plot_base,res_oil,plot_oil,res_model,plot_model,plot_slick
begin
  fn  = "oil_"+sprinti("%0.4i",times(itime))
  fn = "oil_"+MAP_FN(times(itime))
  ffn  = PLOT_DIR+fn
	res_base = _get_res_base
	if times(itime).gt.0 then
		res_base@gsnLeftString = "Surface oil concentration, ton/km~S~2~N~~C~"+TITLE+ \
		" +"+sprinti("%0.4i",times(itime))+" hours"
	else
		res_base@gsnLeftString = "Initial slick position~C~"+TITLE
	end if
;	res_base@gsnRightString = "hour "+MAP_FN(times(itime))
	res_base@gsnRightString = ""
  wks = gsn_open_wks("ps",ffn)
  gsn_define_colormap(wks,PALETTE)
;;;;;;;;;;;;;;;;;;  i = NhlNewColor(wks,COLOR_LAND(0),COLOR_LAND(1),COLOR_LAND(2))
  if GRID then
  	l = NhlNewColor(wks,COLOR_GRID(0),COLOR_GRID(1),COLOR_GRID(2))
  end if
  plot_base = gsn_csm_map(wks,res_base)

  if times(itime).gt.0 then
		res_oil = _get_res_oil
	  plot_oil = gsn_csm_contour(wks,oil_map(itime,:,:),res_oil)
	  overlay(plot_base,plot_oil)
	else
		_plot_slick(wks,plot_base)
	end if

  if MAIN_MAPS_WIND then
	  res_wind = _get_res_wind
  	j = NhlNewColor(wks,COLOR_WIND(0),COLOR_WIND(1),COLOR_WIND(2))
  	w = _get_wind(itime)
	  plot_wind = gsn_csm_vector(wks,w(varU,:,:),w(varV,:,:),res_wind)
	  overlay(plot_base,plot_wind)
	end if
  if MAIN_MAPS_CURRENTS then
	  res_model = _get_res_model
  	k = NhlNewColor(wks,COLOR_CURRENTS(0),COLOR_CURRENTS(1),COLOR_CURRENTS(2))
  	c = _get_model(itime)
;print(""+num(ismissing(c(varU,:,:)))+" "+num(ismissing(c(varV,:,:))))
	  plot_model = gsn_csm_vector(wks,c(varU,:,:),c(varV,:,:),res_model)
	  overlay(plot_base,plot_model)
	end if
  draw(plot_base)
  frame(wks)
	_plot(fn)
  _log(fn)
end

procedure _plot_coast(wks:graphic,plot:graphic,itime:integer)
local i,res
begin
	res = True
	if NSLICK.eq.0 then
		res@gsMarkerIndex = 16
		res@gsMarkerSizeF = CENTER_SIZE
		gsn_polymarker(wks,plot,LON,LAT,res)
	else
		do i = 0,NSLICK-1
			n = SLICKS_N(i)
			gsn_polygon(wks,plot,SLICKS_LON(i,0:n),SLICKS_LAT(i,0:n),res)
		end do
	end if
end





markers = new(1000,graphic)

procedure _plot_obs(wks:graphic,plot:graphic,iobs:integer)
local iobj,i,res,size
begin
	res = True
	res@gsMarkerIndex = 16
	res@gsMarkerSizeF = OBS_SIZE
	do iobj = 0,OBS_N(iobs)-1
		n = OBS_NM(iobs,iobj)
		if OBS_TYPE(iobs,iobj).eq.OBSPOINT then
;			gsn_polymarker(wks,plot,OBS_LON(iobs,iobj,0:n-1),OBS_LAT(iobs,iobj,0:n-1),res)
			markers(iobs) = gsn_add_polymarker(wks,plot,OBS_LON(iobs,iobj,0:n-1),OBS_LAT(iobs,iobj,0:n-1),res)
		else
			gsn_polygon(wks,plot,OBS_LON(iobs,iobj,0:n-1),OBS_LAT(iobs,iobj,0:n-1),res)
		end if
	end do
end

function _get_observation_title(iobs:integer)
local s,k
begin
	k = OBS_YMDHN(iobs)
	s = "Observations and surface oil concentration, ton/km~S~2~N~~C~"+ \
		sprinti("%0.2i",k/10000%100)+"/"+ \
		sprinti("%0.2i",k/1000000%100)+"/"+(2000+k/100000000)+" "+ \
		sprinti("%0.2i",k/100%100)+":"+sprinti("%0.2i",k%100)+" UTC"
	return s
end

procedure _plot_observation_map(iobs:integer)
local fn,ffn,str,wks,i,res_oil,plot_oil,plot_base,str,k, \
			tmp,resmark,rectx,recty,nn,lats,lons,vals,ii,kk,jj,inds
begin
  fn  = "obs_"+sprinti("%0.10i",OBS_YMDHN(iobs))
  fn  = OBS+(iobs+1)
  ffn  = PLOT_DIR+fn
	res_base = _get_res_base
	res_base@gsnLeftString = _get_observation_title(iobs)
	res_base@gsnRightString = "+"+MAP_FN(OBS_IND(iobs))+" hours"

  wks = gsn_open_wks("ps",ffn)
  gsn_define_colormap(wks,PALETTE)
  i = NhlNewColor(wks,COLOR_LAND(0),COLOR_LAND(1),COLOR_LAND(2))
  if GRID then
  	l = NhlNewColor(wks,COLOR_GRID(0),COLOR_GRID(1),COLOR_GRID(2))
  end if
  plot_base = gsn_csm_map(wks,res_base)

  if OBS_IND(iobs).gt.0 then
		res_oil = _get_res_oil
		if False then
		  plot_oil = gsn_csm_contour(wks,obs_map(iobs,:,:),res_oil)
		  overlay(plot_base,plot_oil)
		else
  		plot_oil = gsn_csm_contour(wks,dummy,res_oil)
  		overlay(plot_base,plot_oil)
			tmp = gsn_retrieve_colormap(wks)
			nn = obs_oil_n(iobs)
			if nn.gt.0 then
			  resmark = True
	  		rectx = (/-1,-1,1,1,-1/)*oil_dx*0.5*1
	  		recty = (/-1,1,1,-1,-1/)*oil_dy*0.5*1
				vals = obs_oil_val(iobs,0:nn-1)
				lons = obs_oil_lon(iobs,0:nn-1)
				lats = obs_oil_lat(iobs,0:nn-1)
				do ii = 0,NSCALE
					inds = ind( (OIL_SCALE(ii).lt.vals) .and. (vals.le.(OIL_SCALE(ii+1))) )
					if .not.ismissing(inds(0)) then
						nn = dimsizes(inds)
						kk = COLOR_START+(NCOLOR-COLOR_START+COLOR_END)*(ii+0)/NSCALE
						resmark@gsFillColor = kk
						do jj = 0, nn-1
							gsn_polygon(wks,plot_base,rectx+lons(inds(jj)),recty+lats(inds(jj)),resmark)
						end do
					end if
					delete(inds)
				end do
			end if
		end if
	else
		_plot_slick(wks,plot_base)
	end if
	_plot_obs(wks,plot_base,iobs)

  draw(plot_base)
  frame(wks)
	_plot(fn)
  _log(fn)
end




procedure _test_plot_coast_map(itime:integer)
local fn,ffn,wks,i,j,k,l,w,c,res_base,
plot_base,res_oil,plot_oil,res_model,plot_model \
,res_coast,lona,lonb,lata,latb,poly \
,nn,lats,lons,vals,ii,jj,inds,lon,lat,resmark,colors,tmp,kk \
,ioil
begin
  fn  = "coast_"+sprinti("%0.4i",times(itime))
  ffn  = PLOT_DIR+fn
	res_base = _get_res_base
	if times(itime).gt.-1 then
		res_base@gsnLeftString = "Coast oil concentration, ton/km~C~"+TITLE+ \
		" +"+sprinti("%0.4i",times(itime))+" hours"
	else
		res_base@gsnLeftString = "Initial slick position~C~"+TITLE
	end if
;	res_base@gsnRightString = "hour "+MAP_FN(times(itime))
	res_base@gsnRightString = ""
  wks = gsn_open_wks("ps",ffn)
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  i = NhlNewColor(wks,COLOR_LAND(0),COLOR_LAND(1),COLOR_LAND(2))
  ioil = NhlNewColor(wks,COAST_MAP_OIL_COLOR(0),COAST_MAP_OIL_COLOR(1),COAST_MAP_OIL_COLOR(2))
  if GRID then
  	l = NhlNewColor(wks,COLOR_GRID(0),COLOR_GRID(1),COLOR_GRID(2))
  end if
  plot_base = gsn_csm_map(wks,res_base)

	res_oil = _get_res_oil
	res_oil@cnFillOn              = True
	res_oil@cnInfoLabelOn = False

;  plot_oil = gsn_csm_contour(wks,oil_map(itime,:,:),res_oil)
  plot_oil = gsn_csm_contour(wks,dummy_cst,res_oil)

  overlay(plot_base,plot_oil)


  resmark = True
	resmark@gsFillColor = COAST_MAP_OIL_COLOR
  rectx = (/-1,-1,1,1,-1/)*oil_dx*0.5*1.3
  recty = (/-1,1,1,-1,-1/)*oil_dy*0.5*1.3

if times(itime).gt.0 then
	nn = oil_n(itime)
	vals = oil_val(itime,0:nn-1)
	lons = oil_lon(itime,0:nn-1)
	lats = oil_lat(itime,0:nn-1)
	do ii = 0,NSCALE
		inds = ind( (OIL_SCALE(ii).lt.vals) .and. (vals.le.(OIL_SCALE(ii+1))) )
		if .not.ismissing(inds(0)) then
			nn = dimsizes(inds)
;			do jj = 0, nn-1
;				gsn_polygon(wks,plot_base,rectx+lons(inds(jj)),recty+lats(inds(jj)),resmark)
;			end do
		end if
		delete(inds)
	end do
end if

	if oil_cst_n(itime).gt.0 then
		res_coast = True
		res_coast@gsLineThicknessF = COAST_MAP_THICKNESS
		k = oil_cst_n(itime)
		poly = new(k,graphic)
		do i = 0,k-1
			ii = closest_val(oil_cst(itime,i),CST_SCALE)
			kk = COLOR_START+(NCOLOR-COLOR_START+COLOR_END)*(ii+0)/NSCALE_CST
			res_coast@gsLineColor = kk
			lona = oil_cst_lona(itime,i)
			lonb = oil_cst_lonb(itime,i)
			lata = oil_cst_lata(itime,i)
			latb = oil_cst_latb(itime,i)
			poly(i) = gsn_add_polyline(wks,plot_base,(/lona,lonb/),(/lata,latb/),res_coast)
		end do
	end if

if times(itime).eq.0 then
	_plot_slick(wks,plot_base)
end if

  if MAIN_MAPS_WIND then
	  res_wind = _get_res_wind
  	j = NhlNewColor(wks,COLOR_WIND(0),COLOR_WIND(1),COLOR_WIND(2))
  	w = _get_wind(itime)
	  plot_wind = gsn_csm_vector(wks,w(varU,:,:),w(varV,:,:),res_wind)
	  overlay(plot_base,plot_wind)
	end if
  if MAIN_MAPS_CURRENTS then
	  res_model = _get_res_model
  	k = NhlNewColor(wks,COLOR_CURRENTS(0),COLOR_CURRENTS(1),COLOR_CURRENTS(2))
  	c = _get_model(itime)
;print(""+num(ismissing(c(varU,:,:)))+" "+num(ismissing(c(varV,:,:))))
	  plot_model = gsn_csm_vector(wks,c(varU,:,:),c(varV,:,:),res_model)
	  overlay(plot_base,plot_model)
	end if


  draw(plot_base)
  frame(wks)
	_plot(fn)
  _log(fn)
return

  if True

tmp = gsn_retrieve_colormap(wks)

  resmark = True
  rectx = (/-1,-1,1,1,-1/)*oil_dx*0.5*1
  recty = (/-1,1,1,-1,-1/)*oil_dy*0.5*1

if times(itime).gt.0 then
	nn = oil_n(itime)
	vals = oil_val(itime,0:nn-1)
	lons = oil_lon(itime,0:nn-1)
	lats = oil_lat(itime,0:nn-1)
	do ii = 0,NSCALE
		inds = ind( (OIL_SCALE(ii).lt.vals) .and. (vals.le.(OIL_SCALE(ii+1))) )
		if .not.ismissing(inds(0)) then
			nn = dimsizes(inds)
			kk = COLOR_START+(NCOLOR-COLOR_START+COLOR_END)*(ii+0)/NSCALE
;;;;;print(""+nn+TAB+OIL_SCALE(ii)+TAB+kk+TAB+tmp(kk,0)+TAB+tmp(kk,1)+TAB+tmp(kk,2))
			resmark@gsFillColor = kk
			do jj = 0, nn-1
				gsn_polygon(wks,plot_base,rectx+lons(inds(jj)),recty+lats(inds(jj)),resmark)
			end do
		end if
		delete(inds)
	end do
end if


	end if

if times(itime).eq.0 then
	_plot_slick(wks,plot_base)
end if

	if oil_cst_n(itime).gt.0 then
		res_coast = True
		res_coast@gsLineColor = "black"
		res_coast@gsLineThicknessF = 4.
		k = oil_cst_n(itime)
		poly = new(k,graphic)
		do i = 0,k-1
			lona = oil_cst_lona(itime,i)
			lonb = oil_cst_lonb(itime,i)
			lata = oil_cst_lata(itime,i)
			latb = oil_cst_latb(itime,i)
			poly(i) = gsn_add_polyline(wks,plot_base,(/lona,lonb/),(/lata,latb/),res_coast)
		end do
	end if

  if MAIN_MAPS_WIND then
	  res_wind = _get_res_wind
  	j = NhlNewColor(wks,COLOR_WIND(0),COLOR_WIND(1),COLOR_WIND(2))
  	w = _get_wind(itime)
	  plot_wind = gsn_csm_vector(wks,w(varU,:,:),w(varV,:,:),res_wind)
	  overlay(plot_base,plot_wind)
	end if
  if MAIN_MAPS_CURRENTS then
	  res_model = _get_res_model
  	k = NhlNewColor(wks,COLOR_CURRENTS(0),COLOR_CURRENTS(1),COLOR_CURRENTS(2))
  	c = _get_model(itime)
;print(""+num(ismissing(c(varU,:,:)))+" "+num(ismissing(c(varV,:,:))))
	  plot_model = gsn_csm_vector(wks,c(varU,:,:),c(varV,:,:),res_model)
	  overlay(plot_base,plot_model)
	end if
  draw(plot_base)
  frame(wks)
	_plot(fn)
  _log(fn)

end






procedure _test_plot_main_map(itime:integer)
local fn,ffn,wks,i,j,k,l,w,c,res_base,
plot_base,res_oil,plot_oil,res_model,plot_model \
,res_coast,lona,lonb,lata,latb,poly \
,nn,lats,lons,vals,ii,jj,inds,lon,lat,resmark,colors,tmp,kk
begin
  fn  = "oil_"+sprinti("%0.4i",times(itime))
  ffn  = PLOT_DIR+fn
	res_base = _get_res_base
	if times(itime).gt.0 then
		res_base@gsnLeftString = "Surface oil concentration, ton/km~S~2~N~~C~"+TITLE+ \
		" +"+sprinti("%0.4i",times(itime))+" hours"
	else
		res_base@gsnLeftString = "Initial slick position~C~"+TITLE
	end if
;	res_base@gsnRightString = "hour "+MAP_FN(times(itime))
	res_base@gsnRightString = ""
  wks = gsn_open_wks("ps",ffn)
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  i = NhlNewColor(wks,COLOR_LAND(0),COLOR_LAND(1),COLOR_LAND(2))
  if GRID then
  	l = NhlNewColor(wks,COLOR_GRID(0),COLOR_GRID(1),COLOR_GRID(2))
  end if
  plot_base = gsn_csm_map(wks,res_base)

	res_oil = _get_res_oil
	res_oil@cnFillOn              = True
	res_oil@cnInfoLabelOn = False

;  plot_oil = gsn_csm_contour(wks,oil_map(itime,:,:),res_oil)
  plot_oil = gsn_csm_contour(wks,dummy,res_oil)

  overlay(plot_base,plot_oil)

    if True

tmp = gsn_retrieve_colormap(wks)

  resmark = True
  rectx = (/-1,-1,1,1,-1/)*oil_dx*0.5*1
  recty = (/-1,1,1,-1,-1/)*oil_dy*0.5*1

if times(itime).gt.0 then
	nn = oil_n(itime)
	vals = oil_val(itime,0:nn-1)
	lons = oil_lon(itime,0:nn-1)
	lats = oil_lat(itime,0:nn-1)
	do ii = 0,NSCALE
		inds = ind( (OIL_SCALE(ii).lt.vals) .and. (vals.le.(OIL_SCALE(ii+1))) )
		if .not.ismissing(inds(0)) then
			nn = dimsizes(inds)
			kk = COLOR_START+(NCOLOR-COLOR_START+COLOR_END)*(ii+0)/NSCALE
;;;;;print(""+nn+TAB+OIL_SCALE(ii)+TAB+kk+TAB+tmp(kk,0)+TAB+tmp(kk,1)+TAB+tmp(kk,2))
			resmark@gsFillColor = kk
			do jj = 0, nn-1
				gsn_polygon(wks,plot_base,rectx+lons(inds(jj)),recty+lats(inds(jj)),resmark)
			end do
		end if
		delete(inds)
	end do
end if


	end if

if times(itime).eq.0 then
	_plot_slick(wks,plot_base)
end if

	if oil_cst_n(itime).gt.0 then
		res_coast = True
		res_coast@gsLineColor = "black"
		res_coast@gsLineThicknessF = 4.
		k = oil_cst_n(itime)
		poly = new(k,graphic)
		do i = 0,k-1
			lona = oil_cst_lona(itime,i)
			lonb = oil_cst_lonb(itime,i)
			lata = oil_cst_lata(itime,i)
			latb = oil_cst_latb(itime,i)
			poly(i) = gsn_add_polyline(wks,plot_base,(/lona,lonb/),(/lata,latb/),res_coast)
		end do
	end if

  if MAIN_MAPS_WIND then
	  res_wind = _get_res_wind
  	j = NhlNewColor(wks,COLOR_WIND(0),COLOR_WIND(1),COLOR_WIND(2))
  	w = _get_wind(itime)
	  plot_wind = gsn_csm_vector(wks,w(varU,:,:),w(varV,:,:),res_wind)
	  overlay(plot_base,plot_wind)
	end if
  if MAIN_MAPS_CURRENTS then
	  res_model = _get_res_model
  	k = NhlNewColor(wks,COLOR_CURRENTS(0),COLOR_CURRENTS(1),COLOR_CURRENTS(2))
  	c = _get_model(itime)
;print(""+num(ismissing(c(varU,:,:)))+" "+num(ismissing(c(varV,:,:))))
	  plot_model = gsn_csm_vector(wks,c(varU,:,:),c(varV,:,:),res_model)
	  overlay(plot_base,plot_model)
	end if
  draw(plot_base)
  frame(wks)
	_plot(fn)
  _log(fn)

end


procedure _new_plot_main_map(itime:integer)
local fn,ffn,wks,i,j,k,l,w,c,res_base,plot_base,res_oil,plot_oil,res_model,plot_model,tmp
begin
  fn  = "oil_"+sprinti("%0.4i",times(itime))
  ffn  = PLOT_DIR+fn
	res_base = _get_res_base
	res_base@gsnLeftString = "Surface oil concentration, ton/km~S~2~N~~C~"+TITLE+ \
	" +"+sprinti("%0.4i",times(itime))+" hours"
	res_base@gsnRightString = "hour "+MAP_FN(times(itime))
	res_base@gsnRightString = ""
  wks = gsn_open_wks("ps",ffn)
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  i = NhlNewColor(wks,COLOR_LAND(0),COLOR_LAND(1),COLOR_LAND(2))
  if GRID then
  	l = NhlNewColor(wks,COLOR_GRID(0),COLOR_GRID(1),COLOR_GRID(2))
  end if
  plot_base = gsn_csm_map(wks,res_base)
	res_oil = _get_res_oil
;  plot_oil = gsn_csm_contour(wks,oil_map(itime,:,:),res_oil)
res_oil@cnFillMode            = "AreaFill"

	tmp = _get_oil(itime)
  plot_oil = gsn_csm_contour(wks,tmp,res_oil)
  overlay(plot_base,plot_oil)
  if MAIN_MAPS_WIND then
	  res_wind = _get_res_wind
  	j = NhlNewColor(wks,COLOR_WIND(0),COLOR_WIND(1),COLOR_WIND(2))
  	w = _get_wind(itime)
	  plot_wind = gsn_csm_vector(wks,w(varU,:,:),w(varV,:,:),res_wind)
	  overlay(plot_base,plot_wind)
	end if
  if MAIN_MAPS_CURRENTS then
	  res_model = _get_res_model
  	k = NhlNewColor(wks,COLOR_CURRENTS(0),COLOR_CURRENTS(1),COLOR_CURRENTS(2))
  	c = _get_model(itime)
;print(""+num(ismissing(c(varU,:,:)))+" "+num(ismissing(c(varV,:,:))))
	  plot_model = gsn_csm_vector(wks,c(varU,:,:),c(varV,:,:),res_model)
	  overlay(plot_base,plot_model)
	end if
  draw(plot_base)
  frame(wks)
	_plot(fn)
  _log(fn)
end

procedure _animate_main_maps()
local itime,str
begin
	_log("start producing animation")
	str = "convert -loop 0 -delay "+ANIMATION_DELAY+" +repage"
	do itime = 0,ntime-1
		if times(itime).gt.0 then
			str = str+" "+PLOT_DIR+"oil_"+sprinti("%0.4i",times(itime))+".png"
		end if
	end do
	str = str+" "+PLOT_DIR+"oil_animation.gif"
	system(str)
	_log("oil_animation.gif OK")
	print("")
end

FIRST_IMG = ""

function _option(str:string)
local s,fg
begin
	fg = FIRST_IMG.eq.""
	if fg.and.(NOBS.eq.0) then
		FIRST_IMG = OIL+str
		if str.eq.ANM then
			FIRST_IMG = FIRST_IMG+".gif"
		else
			FIRST_IMG = FIRST_IMG+".png"
		end if
	end if
	if fg.and.(NOBS.gt.0).and..not.ANIMATION then
		FIRST_IMG = "obs_"+str+".png"
	end if
	if fg.and.(NOBS.gt.0).and.ANIMATION then
		FIRST_IMG = OIL+str+".gif"
	end if
	s = TAB+TAB+TAB+"<option"
	if fg then
		s = s+" selected='selected'"
	end if
	s = s+">"+str+"</option>"+LBR
	return s
end

procedure _make_html()
local itime,fn,str
begin
	str = ""
	if SLICK_MAP then
		str = str+_option(SLK)
	end if
	if MAIN_MAPS then
		if ANIMATION then
			str = str+_option(ANM)
		end if
		if NOBS.gt.0 then
			do iobs = 0,NOBS-1
;				str = str+_option(sprinti("%0.10i",OBS_YMDHN(iobs)))
				str = str+_option(OBS+(iobs+1))
			end do
		end if
		do itime = 0,ntime-1
			str = str+_option(MAP_FN(times(itime)))
		end do
	end if
	str = str+TAB+TAB+"</select></td>"+LBR
	str = str+TAB+"<td valign='top'><img src='"+FIRST_IMG+"' name='map' hspace='20' /></td>"+LBR
	str = str+TAB+"</tr>"+LBR+"</table>"+LBR+"</body>"+LBR+"</html>"+LBR
	fn = HTML
	asciiwrite("tmp1.html",str)
	system("cat index_head.html tmp1.html > tmp2.html")
	system("rm tmp1.html")
	system("mv tmp2.html "+PLOT_DIR+HTML)
	_log(HTML+" OK")
end



slick_colors2 = (/"blue4","red4","green4","violetred4","midnightblue","brown4","forestgreen", \
					"purple4","darkslateblue","firebrick4","seagreen4","magenta4","gold4","gray25", \
					"deepskyblue4","coral4","limegreen","deeppink3","red3","darkgreen","darkorchid3"/)

slick_colors = (/ (/0.4,0.4,0.4/), \
(/0.75,0.,0.75/),(/0.,0.75,0./),(/0.,0.,1./), \
(/1.,0.,0./),(/0.,0.4,0.5/),(/0.5,0.,0.5/), \
(/0.7,0.2,0./),(/0.,0.6,0.3/),(/0.3,0.,0.7/) /)

slick_colors = \
(/ (/0.5,0.5,0.5/), \
(/0.8,0.,0.8/), \
(/0.,0.9,0./),\
(/0.,0.,1./), \
(/1.,0.,0./), \
(/0.,0.6,0.6/), \
(/0.6,0.,0.6/), \
(/0.8,0.3,0./), \
(/0.,0.7,0.4/), \
(/0.4,0.,0.8/) /)

procedure _plot_slick_map_old()
local fn,ffn,wks,i,k,l,m,n,res_base,plot_base,color,res,font,x,y,str,tmp1,tmp2,tmp3,resanno, \
lats,lons,lat,lon,j,ns,ks
begin
  fn  = "slicks"
  ffn  = PLOT_DIR+fn
	res_base = _get_res_base
	res_base@gsnLeftString = "Oil slicks (density > 0.05 m~S~3~N~/km~S~2~N~)~C~"+TITLE
  wks = gsn_open_wks("ps",ffn)
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  i = NhlNewColor(wks,COLOR_LAND(0),COLOR_LAND(1),COLOR_LAND(2))
  if GRID then
  	l = NhlNewColor(wks,COLOR_GRID(0),COLOR_GRID(1),COLOR_GRID(2))
  end if
  plot_base = gsn_csm_map(wks,res_base)

  resanno = True
  resanno@amParallelPosF   =  0.52
  resanno@amOrthogonalPosF =  1.0
  resanno@amJust ="TopLeft"

  res = True
	res@txFontHeightF = FONT_SLICKS
  res@txBackgroundFillColor = "white"

  tmp1 = new(ntime,graphic)
  tmp2 = new(ntime,graphic)
  tmp3 = new(ntime,integer)

;  n = dimsizes(slick_colors)
  n = 9
  do i = 0,n-1
	  tmp3(i) = NhlNewColor(wks,slick_colors(i,0),slick_colors(i,1),slick_colors(i,2))
  end do

  font = FONT_SLICKS
	do i = 0,ntime-1
		k = ntime-1-i
		m = k%n
		str = MAP_FN(times(k))
		color = slick_colors(m,:)
		res@txFontColor = color
		res@gsLineColor = color
		res@gsFillColor = color
;		gsn_text_ndc(wks,str,0.75,0.8-k*font*2,res)
    tmp1(i)= gsn_create_text(wks, str, res)
  	resanno@amOrthogonalPosF =  -0.5+k*font*3.
		tmp2(i) = gsn_add_annotation(plot_base, tmp1(i), resanno)
		m = oil_slick_n(k)
	ms = 0
  ns = new(m,integer)
  ks = new((/10,m/),integer)
		lats = oil_slick_lat(k,0:m-1)
		lons = oil_slick_lon(k,0:m-1)
		j0 = 0
		lat = lats(j0)
		lon = lons(j0)
		do j = 1,m-1
			if (j-j0.gt.2).and.(lats(j0).eq.lats(j)).and.(lons(j0).eq.lons(j)) then
				gsn_polygon(wks,plot_base,lons(j0:j),lats(j0:j),res)
print(MAP_FN(times(k))+TAB+j0+TAB+j)
				j0 = j+1
			end if
		end do
		delete(lons)
		delete(lats)
		delete(ns)
		delete(ks)
	end do

  draw(plot_base)
  frame(wks)
	_plot(fn)
  _log(fn)
end

procedure _plot_slick_map()
local fn,ffn,wks,i,k,l,m,n,res_base,plot_base,color,res,font,x,y,str,tmp1,tmp2,tmp3,resanno, \
lats,lons,lat,lon,j,ns,ks,jj,ja,jb,mm,fg
begin
  fn  = "slicks"
  ffn  = PLOT_DIR+fn
	res_base = _get_res_base
	res_base@gsnLeftString = "Oil slicks (density > 0.05 m~S~3~N~/km~S~2~N~)~C~"+TITLE
  wks = gsn_open_wks("ps",ffn)
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  i = NhlNewColor(wks,COLOR_LAND(0),COLOR_LAND(1),COLOR_LAND(2))
  if GRID then
  	l = NhlNewColor(wks,COLOR_GRID(0),COLOR_GRID(1),COLOR_GRID(2))
  end if
  plot_base = gsn_csm_map(wks,res_base)

  resanno = True
  resanno@amParallelPosF   =  0.52
  resanno@amOrthogonalPosF =  1.0
  resanno@amJust ="TopLeft"

  res = True
	res@txFontHeightF = FONT_SLICKS
  res@txBackgroundFillColor = "white"

  tmp1 = new(ntime,graphic)
  tmp2 = new(ntime,graphic)
  tmp3 = new(9,integer)

;  n = dimsizes(slick_colors)
  n = 9
  do i = 0,n-1
	  tmp3(i) = NhlNewColor(wks,slick_colors(i,0),slick_colors(i,1),slick_colors(i,2))
  end do

  font = FONT_SLICKS
	do i = 0,ntime-1
		k = ntime-1-i
		m = k%n
		str = MAP_FN(times(k))
		color = slick_colors(m,:)
		res@txFontColor = color
		res@gsLineColor = color
		res@gsFillColor = color
;		gsn_text_ndc(wks,str,0.75,0.8-k*font*2,res)
    tmp1(i)= gsn_create_text(wks, str, res)
  	resanno@amOrthogonalPosF =  -0.5+k*font*3.
		tmp2(i) = gsn_add_annotation(plot_base, tmp1(i), resanno)
		m = oil_slick_n(k)
	ms = 1
  ns = new(100,integer)
  ns = 0
  ks = new((/100,m/),integer)
		lats = oil_slick_lat(k,0:m-1)
		lons = oil_slick_lon(k,0:m-1)
		ns(0) = 1
		ks = 0

		do j = 1,m-1
			fg = False
			do mm = 0,ms-1
;print(""+k+" "+j+" "+mm+" "+ns(m))
				fg = (abs(oil_slick_lat(k,j)-oil_slick_lat(k,ks(mm,ns(mm)-1))).lt.0.01).and. \
						(abs (oil_slick_lon(k,j)-oil_slick_lon(k,ks(mm,ns(mm)-1))).lt.0.01)
				if fg then
;print(""+k+" "+j+" "+mm+" "+ns(mm)+" "+oil_slick_lon(k,j)+" "+oil_slick_lon(k,ks(mm,ns(mm)-1)) \
;+" "+oil_slick_lat(k,j)+" "+oil_slick_lat(k,ks(mm,ns(mm)-1)))

					ks(mm,ns(mm)) = j
					ns(mm) = ns(mm)+1
;print(""+k+" "+j+" "+mm+" "+ns(mm))
				end if
			end do
			if .not.fg then
				ns(ms) = 1
				ks(ms,0) = j
				ms = ms+1
;print(""+k+" "+j+" "+mm)
			end if
		end do

printVarSummary(ks)
print(ms)

		do mm = 0,ms-1
			jj = ns(mm)
			ja = ks(mm,0)
			jb = ks(mm,jj-1)
print(""+mm+" "+jj+" "+oil_slick_lon(k,ja)+" "+oil_slick_lon(k,jb)+" "+\
oil_slick_lat(k,ja)+" "+oil_slick_lat(k,jb))

;if (jj.lt.4).or.(oil_slick_lon(k,ja).ne.oil_slick_lon(k,jb)).or. \
;(oil_slick_lat(k,ja).ne.oil_slick_lat(k,jb)) then
if (jj.lt.4) then
continue
end if
;;;;print(""+jj)
;;;;ks(mm,jj) = ja
;print(""+oil_slick_lon(k,ks(mm,0:jj-1))+TAB+oil_slick_lat(k,ks(mm,0:jj-1)))
			gsn_polygon(wks,plot_base,oil_slick_lon(k,ks(mm,0:jj)), \
			                          oil_slick_lat(k,ks(mm,0:jj)),res)
		end do

		delete(lons)
		delete(lats)
		delete(ns)
		delete(ks)

;break
	end do

  draw(plot_base)
  frame(wks)
	_plot(fn)
  _log(fn)
end



