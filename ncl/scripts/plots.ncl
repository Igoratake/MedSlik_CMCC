land_color = (/0.85,0.85,0.85/)
wind_color = (/0.4,0.4,0.4/)
model_color = (/0.0,0.0,0.0/)

;log_values = (/0.003, 0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.05, 0.07, 0.1/)
;log_values = (/0.001, 0.0015, 0.002, 0.003, 0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.05/)
;log_values = (/0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2., 5., 10., 20., 50., 100., 200., 500., 1000., 2000./)
;log_values = (/0.5, 1.0, 2., 5., 10., 20., 50., 100., 200., 500., 1000., 2000./)
;log_values = (/0.5, 0.7, 1.0, 1.5, 2., 3., 5., 7.0, 10., 15.,20., 30., 50., 70.,100., \
;150.,200., 300.,500., 700.,1000., 1500.,2000./)
;log_values = (/0.5, 0.7, 1.0, 1.5, 2., 3., 5., 7.0, 10., 15.,20., 30., 50., 70.,100., \
;150.,200., 300./)
log_values = oil_scale
log_labels = ""+log_values

function _get_res_base(itime)
local res
begin
	res = True
	res@gsnFrame   = False
	res@gsnDraw   = False

	res@mpProjection = "Mercator"
	res@mpOutlineOn  = True
	res@mpDataBaseVersion = "HighRes" 
	res@mpDataResolution  = "Finest"
	res@mpLandFillColor      = land_color  ; set land to be gray
  res@mpLimitMode = "LatLon"
  res@mpFillDrawOrder = "PreDraw"
  res@mpOutlineDrawOrder = "PreDraw"

  d = 0.0
  res@mpMinLatF = min_map_lat-d
  res@mpMaxLatF = max_map_lat+d
  res@mpMinLonF = min_map_lon-d
  res@mpMaxLonF = max_map_lon+d
  res@mpGridAndLimbOn =  True
  
  res@pmTickMarkDisplayMode    = "Always"
    
  res@mpGridLineThicknessF   =  0.2   
  res@mpGridLineColor        = "black"
  res@mpGridLonSpacingF      = 1./12.
  res@mpGridLatSpacingF      = 1./12.

  res@tmXTOn       = True
  res@tmYROn       = True
  res@tmXBMajorLengthF       = 0.002
  res@tmYLMajorLengthF       = 0.002
  res@tmXBMinorLengthF       = 0.001
  res@tmYLMinorLengthF       = 0.001
  res@tmXTMajorLengthF       = 0.002
  res@tmYRMajorLengthF       = 0.002
  res@tmXTMinorLengthF       = 0.001
  res@tmYRMinorLengthF       = 0.001
  res@tmXBLabelFontHeightF   = 0.010
  res@tmYLLabelFontHeightF   = 0.010
  
  res@tmBorderThicknessF = 1.0
  res@gsnCenterString = ""
  res@gsnRightString = ""
  res@gsnLeftString = ""
  
  res@gsnCenterStringFontHeightF     = 0.012                     
  res@gsnCenterStringOrthogonalPosF  = 0.014
  
  res@gsnLeftStringFontHeightF       = 0.012                     
  res@gsnLeftStringOrthogonalPosF    = 0.014
  
  res@gsnRightStringFontHeightF      = 0.012                    
  res@gsnRightStringOrthogonalPosF   = 0.014
  
;  res@gsnRightString = "hour "+sprinti("%0.2i",itime)
 ; res@gsnLeftString = _julian2DDMMYYYY(TIMEA)+" "+ \
  ;		sprinti("%0.2i",HOUR)+":"+sprinti("%0.2i",MINUTE)+" UTC"

	return res
end

function _get_res_overlay()
local res
begin
	res = True
	res@gsnAddCyclic = False
  res@gsnTickMarksOn   = False	; no tickmarks
  res@tmBorderThicknessF = 0.
  res@gsnCenterString = ""
  res@gsnRightString = ""
  res@gsnLeftString = ""
  
  res@gsnCenterStringFontHeightF     = 0.014                     
  res@gsnCenterStringOrthogonalPosF  = 0.014
  
  res@gsnLeftStringFontHeightF       = 0.014                     
  res@gsnLeftStringOrthogonalPosF    = 0.014
  
  res@gsnRightStringFontHeightF      = 0.014                    
  res@gsnRightStringOrthogonalPosF   = 0.014
	return res
end

function _get_res_model()
local res
begin
	res = _get_res_overlay
  res@vcRefMagnitudeF           = 0.05             ; define vector ref mag
  res@vcRefLengthF              = 0.05           ; define length of vec ref
  res@vcPositionMode          = "ArrowTail"      
  res@vcLineArrowHeadMaxSizeF = 0.01 
  res@vcLineArrowColor = model_color
  res@vcRefAnnoString2 = "Currents"
;  res@vcRefAnnoPerimSpaceF = 0.5
  res@vcRefAnnoString1 = ""+res@vcRefMagnitudeF+" m/s"
  res@vcRefAnnoPerimOn = False
  res@vcLineArrowThicknessF = 2.
  res@vcRefAnnoFontHeightF = 0.011
	return res
end

function _get_res_wind()
local res
begin
	res = _get_res_overlay
  res@vcRefMagnitudeF           = 5.0             ; define vector ref mag
  res@vcRefLengthF              = 0.05           ; define length of vec ref
  res@vcPositionMode          = "ArrowTail"     
  res@vcLineArrowHeadMaxSizeF = 0.01 
  res@vcRefAnnoParallelPosF = 0.1
  res@vcLineArrowColor = wind_color
  res@vcRefAnnoString2 = "Wind"
;  res@vcRefAnnoPerimSpaceF = 0.5
  res@vcRefAnnoPerimOn = False
  res@vcLineArrowThicknessF = 2.
  res@vcRefAnnoFontHeightF = 0.011
  res@vcRefAnnoString1 = ""+res@vcRefMagnitudeF+" m/s"

	return res
end

function _get_res_oil(itime:integer)
local res
begin
	res = _get_res_overlay
	res@cnFillOn             = True             ; turn on color for contours
	res@cnLinesOn            = False            ; turn off contour lines
	res@cnLineLabelsOn       = False            ; turn off contour line labels
	res@gsnSpreadColors      = True             ; use full color map
	res@gsnSpreadColorEnd    = -4               ; don't use land color
	res@gsnSpreadColorStart  = 3                ; start at color 17
  res@cnFillMode = "CellFill"
  res@cnMissingValFillColor = -1
  res@pmLabelBarSide = "Right"
  res@lbOrientation  = "Vertical"
  res@lbLabelJust    = "CenterLeft"
  res@pmLabelBarWidthF   = 0.07
  res@lbLabelFontHeightF = 0.011
  res@lbLabelStride   = 1
  res@cnLevelSelectionMode = "ExplicitLevels"      
  res@cnLevels = log_values      
  res@cnExplicitLabelBarLabelsOn = True
	res@lbLabelStrings = log_labels
  res@gsnLeftString = "Surface oil concentration, ton/km~S~2~N~~C~"+ \
  _julian2DDMMYYYY(TIMEA)+" "+ \
  		sprinti("%0.2i",HOUR)+":"+sprinti("%0.2i",MINUTE)+" UTC"
  res@gsnRightString = "hour "+sprinti("%0.2i",itime)
;  res@gsnLeftString = _julian2DDMMYYYY(TIMEA)+" "+ \
;  		sprinti("%0.2i",HOUR)+":"+sprinti("%0.2i",MINUTE)+" UTC"

  res@gsnCenterString = ""
  res@gsnCenterStringOrthogonalPosF   = 0.1
;  res@tiMainFontHeightF = 0.014
;  res@tiMainString = "Surface oil concentration, ton/km~S~2"
;  res@gsnCenterStringOrthogonalPosF   = 0.1

	return res
end

procedure _plot_map(itime:integer)
local i,j,k,res_base,plot_base,res_oil,plot_oil,res_model,plot_model,fn,ffn,wks
begin
  fn  = "oil_map_"+sprinti("%0.3i",itime+1)
  ffn  = DIR_PLOTS+fn
	res_base = _get_res_base(itime)
  wks = gsn_open_wks("ps",ffn)  
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200") 
  i = NhlNewColor(wks,land_color(0),land_color(1),land_color(2))
  j = NhlNewColor(wks,wind_color(0),wind_color(1),wind_color(2))
  j = NhlNewColor(wks,model_color(0),model_color(1),model_color(2))
  plot_base = gsn_csm_map(wks,res_base)
	res_oil = _get_res_oil(itime)
	_process_oil(itime)
  plot_oil = gsn_csm_contour(wks,oil_map,res_oil)
  overlay(plot_base,plot_oil)
  res_wind = _get_res_wind
  plot_wind = gsn_csm_vector(wks,wind_model_u(itime,model_ja:model_jb,model_ia:model_ib), \
		wind_model_v(itime,model_ja:model_jb,model_ia:model_ib),res_wind)
  overlay(plot_base,plot_wind)
  res_model = _get_res_model
  plot_model = gsn_csm_vector(wks,model_u(itime,model_ja:model_jb,model_ia:model_ib), \
  	model_v(itime,model_ja:model_jb,model_ia:model_ib),res_model)
  overlay(plot_base,plot_model)
  draw(plot_base)
  frame(wks)
  delete(wks)
  delete(res_base)
  delete(plot_base)
  system("convert -density 120 -trim "+ffn+".ps'[0]' +repage -bordercolor White -border 5x5 "+ffn+"0.png")
  system("convert -density 120 -trim "+ffn+".ps'[1]' +repage -bordercolor White -border 5x5 "+ffn+"1.png")
  system("convert -density 120 -trim "+ffn+".ps[2] +repage -bordercolor White -border 5x5 "+ffn+"2.png")
  system("convert -density 120 -trim "+ffn+".ps[3] +repage -bordercolor White -border 5x5 "+ffn+"3.png")
; system("rm "+ffn+".ps")
  _log(fn)
end














function _get_res_map(itime)
local res
begin
	res = True
	res@gsnFrame   = False

	res@mpProjection = "Mercator"
	res@gsnAddCyclic = False
	res@mpOutlineOn  = True
	res@mpDataBaseVersion = "HighRes" 
	res@mpDataResolution  = "Finest"
	res@mpDataResolution  = "Finest"
	
	res@cnFillOn             = True             ; turn on color for contours
	res@cnLinesOn            = False            ; turn off contour lines
	res@cnLineLabelsOn       = False            ; turn off contour line labels
;	res@gsnScalarContour     = True             ; contours desired
	res@gsnSpreadColors      = True             ; use full color map
	res@gsnSpreadColorEnd    = -3               ; don't use land color
	res@gsnSpreadColorStart  = 3                ; start at color 17
	res@mpLandFillColor      = land_color  ; set land to be gray
  res@mpLimitMode = "LatLon"
  
  res@mpFillDrawOrder = "PreDraw"
  res@mpOutlineDrawOrder = "PreDraw"
  res@cnFillDrawOrder  = "PostDraw"
; res@mpFillDrawOrder = "PostDraw"

  d = 0.0
  res@mpMinLatF = min_map_lat-d
  res@mpMaxLatF = max_map_lat+d
  res@mpMinLonF = min_map_lon-d
  res@mpMaxLonF = max_map_lon+d
  res@mpGridAndLimbOn =  True
  
;      res@cnFillMode = "RasterFill"    
      res@cnFillMode = "CellFill"
      res@cnMissingValFillColor = -1
          

   res@pmTickMarkDisplayMode    = "Always"
 
    res@pmLabelBarSide = "Right"
    res@lbOrientation  = "Vertical"
    res@lbLabelJust    = "CenterLeft"
    res@pmLabelBarWidthF   = 0.07
    res@lbLabelFontHeightF = 0.012
;    res@lbTitleString = "Surface oil concentration, t/km~S~2"
;    res@lbTitlePosition = "Left"
;    res@lbTitleFontHeightF = 0.014
;    res@lbTitleAngleF = 90
    
  res@mpGridLineThicknessF   =  0.2   
  res@mpGridLineColor        = "black"
  res@mpGridLonSpacingF      = 1./12.
  res@mpGridLatSpacingF      = 1./12.
;  res@mpGridLonSpacingF      = 1./2.
;  res@mpGridLatSpacingF      = 1./2.


  res@tmXTOn       = True
  res@tmYROn       = True
  res@tmXBMajorLengthF       = 0.002
  res@tmYLMajorLengthF       = 0.002
  res@tmXBMinorLengthF       = 0.001
  res@tmYLMinorLengthF       = 0.001
  res@tmXTMajorLengthF       = 0.002
  res@tmYRMajorLengthF       = 0.002
  res@tmXTMinorLengthF       = 0.001
  res@tmYRMinorLengthF       = 0.001
  res@tmXBLabelFontHeightF   = 0.010
  res@tmYLLabelFontHeightF   = 0.010

;  res@cnLevelSelectionMode = "ManualLevels"      
;  res@cnMinLevelValF  = oil_map_a
;  res@cnMaxLevelValF  = oil_map_b
 ; res@lbLabelStride   = oil_map_stride
  ;res@cnLevelSpacingF = oil_map_spacing
  
  res@lbLabelStride   = 1
    res@cnLevelSelectionMode = "ExplicitLevels"      
    res@cnLevels = log_values      
  res@cnExplicitLabelBarLabelsOn = True
res@lbLabelStrings = log_labels


	return res
end

procedure _plot_oil_map(itime:integer)
local i,res,wks,plot,tmp,fn,ffn
begin
	_process_oil(itime)
	res = _get_res_map(itime)
  fn  = "oil_map_"+sprinti("%0.3i",itime+1)
  ffn  = DIR_PLOTS+fn
  wks = gsn_open_wks("ps",ffn)  
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200") 
;  i = NhlNewColor(wks,.96,.87,.7)
  i = NhlNewColor(wks,land_color(0),land_color(1),land_color(2))
  plot = gsn_csm_contour_map(wks,oil_map,res)
;  plot = gsn_csm_contour_map_ce(wks,oil_map,res)
  frame(wks)
  delete(wks)
  delete(res)
  delete(plot)
  system("convert -density 120 -trim "+ffn+".ps +repage -bordercolor White -border 5x5 "+ffn+".png")
  system("rm "+ffn+".ps")
  _log(fn)
end
