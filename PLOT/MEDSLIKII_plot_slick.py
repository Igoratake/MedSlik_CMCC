#!/usr/bin/env python

# =================================================================================
#                            PLOT_SLICK.PY
#
# This software read NetCDF Medslik-II outputs and plots:
#
#        1. Oil concentration fields on a 150x150 m regular grid
#
#        2. Ocean currents velocities fields on the original ocean currents grid
#
#        3. Stokes' Drift velocities fields on the original waves grid
#
#-----------------------------------------------------------------------------------

import sys
import numpy as np
import os
import Nio
import Ngl

# -------------- Function nearest --------------------
def nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx
# -----------------------------------------------------

dir  = sys.argv[1]
file = sys.argv[2]
iwave = sys.argv[3]

if iwave <> "000":
   iwave = "Y"
else:
   iwave = "N"

name = dir + file 

f_nc = Nio.open_file(name,"r")
cst_file = dir + "coast.map"

# ----------------------------------------------------------------
#                         TIME variable
# ----------------------------------------------------------------

time = f_nc.variables["time"]

# ----------------------------------------------------------------
#                  OIL CONCENTRATION VARIABLE
# ----------------------------------------------------------------

oil = f_nc.variables["surface"]

if hasattr(oil,"_FillValue"):
   NaN_oil = oil._FillValue
   print "oil FillValue = ", NaN_oil

if hasattr(oil,"scale_factor"):
   sc_f_oil = oil.scale_factor
   print "oil scale_factor = ", sc_f_oil
else:
   sc_f_oil = 1.
if hasattr(oil,"add_offset"):
   ad_o_oil = oil.add_offset
   print "conc add offset = ", ad_o_oil
else:
   ad_o_oil = 0.

C1 = (ad_o_oil + (oil[1:len(time[:]),:,:] * sc_f_oil))

Cmin = np.amin(C1)
print "Cmin: ", Cmin
Cmax = np.amax(C1)
print "Cmax: ", Cmax

C = (ad_o_oil + (oil[:,:,:] * sc_f_oil))

# ----------------------------------------------------------------
#                  OIL CONCENTRATION ON COAST VARIABLE
# ----------------------------------------------------------------

cst_seg  = open(cst_file).readlines()

seg = np.zeros((len(cst_seg),4))

if len(cst_seg) != 0:
   for i in range(0,len(cst_seg)):
       if float(cst_seg[i].split()[0]) !=  999999:
           for j in range(0,4):
               seg[i,j] = float(cst_seg[i].split()[j])

varNames = f_nc.variables.keys()
if varNames[9] == "coast":

   cst = f_nc.variables["coast"]

   if hasattr(cst,"_FillValue"):
      NaN_cst = cst._FillValue
      print "cst FillValue = ", NaN_cst
 
   if hasattr(cst,"scale_factor"):
      sc_f_cst = cst.scale_factor
      print "cst scale_factor = ", sc_f_cst
   else:
      sc_f_cst = 1.
   if hasattr(cst,"add_offset"):
      ad_o_cst = cst.add_offset
      print "conc add offset = ", ad_o_cst
   else:
      ad_o_cst = 0.

   CST1 = (ad_o_cst + (cst[1:len(time[:]),:] * sc_f_cst))

   COASTmin = np.amin(CST1)
   print "COASTmin: ", COASTmin
   COASTmax = np.amax(CST1)
   print "COASTmax: ", COASTmax

   CST = (ad_o_cst + (cst[:,:] * sc_f_cst))
else:
   CST = np.zeros((len(time),1))

# ---------------------------------------------------------------
#                    CURRENTS  U and V variables
# ---------------------------------------------------------------

uu = f_nc.variables["u_cur"]
vv = f_nc.variables["v_cur"]

if hasattr(uu,"_FillValue"):
   NaN_c = uu._FillValue
   print "u FillValue = ", NaN_c

if hasattr(uu,"scale_factor"):
   sc_f_c = uu.scale_factor
   print "conc scale_factor = ", sc_f_c
else:
   sc_f_c = 1.
if hasattr(uu,"add_offset"):
   ad_o_c = uu.add_offset
   print "u add offset = ", ad_o_c
else:
   ad_o_c = 0.

u = (ad_o_c + (uu[:,:,:] * sc_f_c))
v = (ad_o_c + (vv[:,:,:] * sc_f_c))

umin = np.amin(u)      
print "umin: ", umin
umax = np.amax(u)      
print "umax: ", umax

vmin = np.amin(v)
print "vmin: ", vmin
vmax = np.amax(u)
print "vmax: ", vmax

# ---------------------------------------------------------------
#                    WAVE  U and V variables
# ---------------------------------------------------------------

if iwave == "Y":
   wu = f_nc.variables["u_wav"]
   wv = f_nc.variables["v_wav"]

   if hasattr(wu,"_FillValue"):
      NaN_w = wu._FillValue
      print "u FillValue = ", NaN_w

   if hasattr(wu,"scale_factor"):
      sc_f_w = wu.scale_factor
      print "conc scale_factor = ", sc_f_w
   else:
      sc_f_w = 1.
   if hasattr(wu,"add_offset"):
      ad_o_w = wu.add_offset
      print "u add offset = ", ad_o_w
   else:
      ad_o_w = 0.

   uw = (ad_o_w + (wu[:,:,:] * sc_f_w))
   vw = (ad_o_w + (wv[:,:,:] * sc_f_w))

   uwmin = np.amin(uw)
   print "uwmin: ", uwmin
   uwmax = np.amax(uw)
   print "uwmax: ", uwmax

   vwmin = np.amin(vw)
   print "vwmin: ", vwmin
   vwmax = np.amax(vw)
   print "vwmax: ", vwmax

# ------------------------------------------------------------------------------
#                         WIND WINX and WINY variables
# ------------------------------------------------------------------------------

xwin = f_nc.variables["x_wnd"]
ywin = f_nc.variables["y_wnd"]

if hasattr(xwin,"_FillValue"):
   NaN_wd = xwin._FillValue
   print "wind FillValue = ", NaN_wd

if hasattr(xwin,"scale_factor"):
   sc_f_wd = xwin.scale_factor
   print "wind scale_factor = ", sc_f_wd
else:
   sc_f_wd = 1.
if hasattr(xwin,"add_offset"):
   ad_o_wd = xwin.add_offset
   print "wind add offset = ", ad_o_wd
else:
   ad_o_wd = 0.

winx = (ad_o_wd + (xwin[:,:] * sc_f_wd))
winy = (ad_o_wd + (ywin[:,:] * sc_f_wd))

winxmin = np.amin(winx)
print "winx_min: ", winxmin
winxmax = np.amax(winx)
print "winx_max: ", winxmax

winymin = np.amin(winy)
print "winy_min: ", winymin
winymax = np.amax(winy)
print "winy_max: ", winymax

# -------------------------------------------------------------------------------
#                           OIL LON and LAT variables
# -------------------------------------------------------------------------------

lon = f_nc.variables["lon_oil"]
lat = f_nc.variables["lat_oil"]

print "lon shape = ", lon.shape
print "lon dimensions = ", lon.dimensions
print "lat shape = ", lat.shape
print "lat dimensions = ", lat.dimensions

nlon = len(lon[:])
nlat = len(lat[:])

# -------------------------------------------------------------------------------
#                         CURRENTS LON and LAT variables
# -------------------------------------------------------------------------------

lon_cur = f_nc.variables["lon_cur"]
lat_cur = f_nc.variables["lat_cur"]

print "lon_cur shape = ", lon_cur.shape
print "lon_cur dimensions = ", lon_cur.dimensions

print "lat_cur shape = ", lat_cur.shape
print "lat_cur dimensions = ", lat_cur.dimensions

nlon_cur = len(lon_cur[:])
nlat_cur = len(lat_cur[:])

# -------------------------------------------------------------------------------
#                         WAVES LON and LAT variables
# -------------------------------------------------------------------------------

if iwave == "Y":
   lon_wav = f_nc.variables["lon_wav"]
   lat_wav = f_nc.variables["lat_wav"]

   print "lon_wav shape = ", lon_wav.shape
   print "lon_wav dimensions = ", lon_wav.dimensions

   print "lat_wav shape = ", lat_wav.shape
   print "lat_wav dimensions = ", lat_wav.dimensions

   nlon_wav = len(lon_wav[:])
   nlat_wav = len(lat_wav[:])

# --------------------------------------------------------------------------------
#                          CoG LON and LAT variables
# --------------------------------------------------------------------------------

lon_cog = f_nc.variables["lon_cog"]
lat_cog = f_nc.variables["lat_cog"]

print "lon_cog shape = ", lon_cog.shape
print "lon_cog dimensions = ", lon_cog.dimensions
print "lat_cog shape = ", lat_cog.shape
print "lat_cog dimensions = ", lat_cog.dimensions

nlon_cog = lon_cog.shape[1]
nlat_cog = lat_cog.shape[1]

# --------------------------------------------------------------------------------
# Sampling Position During Serious Game:
#
#       SAMPLE ID    Time (UTC)    Latitude              Longitude
#       -------------------------------------------------------------
#          M1         7:08       42 58.01'N           009 59.26'E
#          M2         7:35       42 58.82'N           009 59.34'E
#          M3         8:27       42 57.90'N           009 58.33'E
#          M4         8:49       42 58.21'N           009 59.37'E
# --------------------------------------------------------------------------------

m1_lat = 42. + 58.01 / 60.
m1_lon = 9. + 59.26 / 60.

m2_lat = 42. + 58.82 / 60.
m2_lon = 9. + 59.34 / 60.

m3_lat = 42. + 57.90 / 60.
m3_lon =9. + 58.33 / 60.

m4_lat = 42. + 58.21 / 60.
m4_lon =9. + 59.37 / 60.

m_array = np.array([[m1_lat, m1_lon],
                    [m2_lat, m2_lon],
                    [m3_lat, m3_lon],
                    [m4_lat, m4_lon]])
print "m_array: ", m_array.shape
print m_array

# -------------------------------------------------------------------------------------
# Drifter position:

#      Drifter_Type       Drifter_IMEI     Time(UTC)       Latitude         Longitude
#      --------------------------------------------------------------------------------
#       I-SPHERE        300234061615020     11:38:26     42 58.232' N    9 59.029' E
#       TOSCA           393667769989        11:38:41     42 58.222' N    9 59.016' E
#       CLS(Beacon)     125677              11:39:01     42 58.197' N    9 59.025' E
#       I-SPHERE        300034012659810     11:39:17     42 58.181' N    9 59.030' E
#       TOSCA           393667769966        11:39:35     42 58.160' N    9 59.019' E
#       CLS(Beacon)     125725              11:39:50     42 58.144' N    9 59.010' E
#       CODE Elcon      3356986728          11:08:39     42 58.294' N    9 58.833' E
#       I-SLDMB(ITCG)   300234061640130     11:00:15     42 58.166' N    9 58.863' E
#       CODE Elcon      3351513355          11:04:47     42 58.163' N    9 58.785' E
# --------------------------------------------------------------------------------------

d1_lat = 42 + 58.232 / 60.
d1_lon = 9 + 59.029 / 60.

d2_lat = 42 + 58.222 / 60.
d2_lon = 9 + 59.016 / 60.

d3_lat = 42 + 58.197 / 60.
d3_lon = 9 + 59.025 / 60.

d4_lat = 42 + 58.181 / 60.
d4_lon = 9 + 59.030 / 60.

d5_lat = 42 + 58.160 / 60.
d5_lon = 9 + 59.019 / 60.

d6_lat = 42 + 58.144 / 60.
d6_lon = 9 + 59.010 / 60.

d7_lat = 42 + 58.294 / 60.
d7_lon = 9 + 58.833 / 60.

d8_lat = 42 + 58.166 / 60.
d8_lon = 9 + 58.863 / 60.

d9_lat = 42 + 58.163 / 60.
d9_lon = 9 + 59.785 / 60.

d_array = np.array([[d1_lat, d1_lon],
                    [d2_lat, d2_lon],
                    [d3_lat, d3_lon],
                    [d4_lat, d4_lon],
                    [d5_lat, d5_lon],
                    [d6_lat, d6_lon],
                    [d7_lat, d7_lon],
                    [d8_lat, d8_lon],
                    [d9_lat, d9_lon]])
print d_array

# --------------------------------------------------------------------------------
#                   LOOP on time to plot scalar and vector fields
# --------------------------------------------------------------------------------

for t in range(0, len(time[:])):

    C1 = C[t,:,:]

    CST1 = CST[t,:]

# Changing NaN value in ocean forcing fields ------------------------

    ut = u[t,:,:]
    vt = v[t,:,:] 

    u_log = ut == 9999.
    v_log = vt == 9999.
    ut[u_log] = 0.
    vt[v_log] = 0.    

    if iwave == "Y":
       uwt = uw[t,:,:]
       vwt = vw[t,:,:]

       uw_log = uwt == 9999.
       vw_log = vwt == 9999.
       uwt[uw_log] = 0.
       vwt[vw_log] = 0.
    
# Regridding Wind into a 2D grid -------------------------------------

    winxt = winx[t,:]
    winyt = winy[t,:]

    nxout = 200    # Define 2D output grid 
    nyout = 200
    xmin    = np.amin(lon_cog[t,:])-0.01
    ymin    = np.amin(lat_cog[t,:])-0.01
    xmax    = np.amax(lon_cog[t,:])+0.01
    ymax    = np.amax(lat_cog[t,:])+0.01
    delx    = (xmax-xmin)/(nxout-1)
    dely    = (ymax-ymin)/(nyout-1)
    xo      = xmin + delx*np.arange(0,nxout)
    yo      = ymin + dely*np.arange(0,nyout)

    winx2D = np.zeros((nxout,nyout))
    winy2D = np.zeros((nxout,nyout))

    for k in range(0,nlon_cog):
        j = nearest(xo, lon_cog[t,k])
        print lon_cog[t,k], xo[j]
        i = nearest(yo, lat_cog[t,k])
        print lat_cog[t,k], yo[i]
        winx2D[i,j]  = winxt[k]   # Interpolate.
        winy2D[i,j]  = winyt[k]   # Interpolate.

# Plot ---------------------------------------------------------------

    if t < 10:
       out_name = "oil_000" + str(t)
    elif (t >= 10 and t < 100):
       out_name = "oil_00" + str(t)
    elif t >= 100:
       out_name = "oil_0" + str(t)

# 1. Open the WorKStation
    
    rlist            = Ngl.Resources()
    rlist.wkColorMap = "BlAqGrYeOrReVi200" #cmap
    wks_type = "png"
    wks = Ngl.open_wks(wks_type,out_name,rlist)
    cmap = Ngl.retrieve_colormap(wks)
    Ngl.destroy(wks)

    cmap[0] = [1.,1.,1.] # White
    cmap[1] = [1.,1.,1.] # White
    cmap[2] = [1.,1.,1.] # White
#    cmap[3] = [1.,1.,1.] # White

    rlist.wkColorMap = cmap
    wks = Ngl.open_wks(wks_type,out_name,rlist)

# 2. Set up resource lists will apply to vector, line contour, and
#    filled contour plots.

    mpres   = Ngl.Resources()       # map resources
    c_vcres = Ngl.Resources()       # ocean currents vector resources
    if iwave == "Y":
       w_vcres = Ngl.Resources()    # Stokes' Drift vector resources
    wd_vres = Ngl.Resources()       # Wind vector resources
    cfres   = Ngl.Resources()       # Contour line/fill resources
    m1res   = Ngl.Resources()       # 1 Sampling marker resources
    m2res   = Ngl.Resources()       # 2 Sampling marker resources
    cstres  = Ngl.Resources()       # Coastal segments polyline resources

# 3. Turn off nglDraw and nglFrame because we don't want to draw all
#    these plots until they are all overlaid on the map plot.

    mpres.nglDraw    = False
    mpres.nglFrame   = False

    c_vcres.nglDraw  = False
    c_vcres.nglFrame = False

    if iwave == "Y":
       w_vcres.nglDraw  = False
       w_vcres.nglFrame = False

    wd_vres.nglDraw  = False
    wd_vres.nglFrame = False

    cfres.nglDraw    = False
    cfres.nglFrame   = False

    m1res.nglDraw    = False
    m1res.nglFrame   = False

    m2res.nglDraw    = False
    m2res.nglFrame   = False

    cstres.nglDraw   = False
    cstres.nglFrame  = False

# 4. Set up coordinates of X and Y axes for all plots. This
#    is necessary in order for the Ngl.overlay calls to work
#    later.
 
    cfres.sfXCStartV              = float(lon[0])                   # Define X/Y axes range
    cfres.sfXCEndV                = float(lon[nlon-1])              # for contour plot.
    cfres.sfYCStartV              = float(lat[0])
    cfres.sfYCEndV                = float(lat[nlat-1])

    c_vcres.vfXCStartV              = float(lon_cur[0])             # Define X/Y axes range
    c_vcres.vfXCEndV                = float(lon_cur[nlon_cur-1])    # for vector plot.
    c_vcres.vfYCStartV              = float(lat_cur[0])
    c_vcres.vfYCEndV                = float(lat_cur[nlat_cur-1])

    if iwave == "Y":
       w_vcres.vfXCStartV              = float(lon_wav[0])          # Define X/Y axes range
       w_vcres.vfXCEndV                = float(lon_wav[nlon_wav-1]) # for vector plot.
       w_vcres.vfYCStartV              = float(lat_wav[0])
       w_vcres.vfYCEndV                = float(lat_wav[nlat_wav-1])

# 5. CURRENTS VECTOR FIELD SETUP -----------------------------------------------------------
    c_vcres.vcGlyphStyle              = "LineArrow"
    c_vcres.vcMonoLineArrowColor      = True
    c_vcres.vcLineArrowColor          = "blue" 
    c_vcres.vcRefLengthF              = 0.05
    c_vcres.vcRefMagnitudeF           = 0.2
    c_vcres.vcPositionMode            = "ArrowTail"
    c_vcres.vcLineArrowThicknessF     = 2.              # Double the thickness.

#    c_vcres.vcRefAnnoOrthogonalPosF = -0.179           # Move reference annotation up
    c_vcres.vcRefAnnoParallelPosF   =  0.19             # and over to left.
    c_vcres.vcRefAnnoString1        = "0.2 m/s"
    c_vcres.vcRefAnnoString2        = "Ocean Currents"

# 6. WAVES VECTOR FIELD SETUP --------------------------------------------------------------

    if iwave == "Y":
       w_vcres.vcGlyphStyle              = "LineArrow"
       w_vcres.vcMonoLineArrowColor      = True
       w_vcres.vcLineArrowColor          = "red"
       w_vcres.vcRefLengthF              = 0.05
       w_vcres.vcRefMagnitudeF           = 0.1
       w_vcres.vcPositionMode            = "ArrowTail"
       w_vcres.vcLineArrowThicknessF     = 2.           # Double the thickness.

#    vcres.vcRefAnnoOrthogonalPosF = -0.179             # Move reference annotation up
#    vcres.vcRefAnnoParallelPosF   =  0.19              # and over to left.
       w_vcres.vcRefAnnoString1        = "0.1 m/s"
       w_vcres.vcRefAnnoString2        = "Wave-induced Currents"

# 7. WIND VECTOR FIELD SETUP ---------------------------------------------------------------

    wd_vres.vfXArray                 = xo
    wd_vres.vfYArray                 = yo
    wd_vres.vcGlyphStyle              = "LineArrow"
    wd_vres.vcMonoLineArrowColor      = True
    wd_vres.vcLineArrowColor          = "green"
    wd_vres.vcRefLengthF              = 0.05
    wd_vres.vcRefMagnitudeF           = 2
    wd_vres.vcPositionMode            = "ArrowTail"
    wd_vres.vcLineArrowThicknessF     = 4.                # Quadruple the thickness.

#    wd_vres.vcRefAnnoOrthogonalPosF   =  0.159           # Move reference annotation up
    wd_vres.vcRefAnnoParallelPosF     =  0.59             # and over to left.
    wd_vres.vcRefAnnoString1          = "2 m/s"
    wd_vres.vcRefAnnoString2          = "Wind velocity"


# 8. SCALAR FIELD SETUP ---------------------------------------------------------------------

    cfres.cnFillOn                  = True            # Turn on contour fill.
    cfres.cnLinesOn                 = False           # Turn on contour lines.
    cfres.cnLineLabelsOn            = False           # Turn off contour line labels.
    cfres.lbOrientation             = "Horizontal"    # horizontal labelbar
    cfres.lbLabelFontHeightF        = 0.012           # Decrease font size.
    cfres.pmLabelBarOrthogonalPosF  = -0.05           # Move labelbar up.
    cfres.cnLineDashPatterns        = 3               # dashed contour lines
    cfres.cnLineThicknessF          = 3.              # triple thick contour lines
    cfres.cnInfoLabelOrthogonalPosF = -0.15           # Move info label up.
    cfres.cnFillMode                = "RasterFill"    # or "AreaFill" or "CellFill" or "RasterFill"
    cfres.cnRasterSmoothingOn       = True
    cfres.cnLevelSelectionMode      = "ExplicitLevels"    # Define your own contour levels
    cfres.cnLevels                  = np.arange(0.01,Cmax,(Cmax-0.01)/199)
#    cfres.cnLevels                  = [ 0.001,0.002,0.005,0.1,0.2,0.5,1.,2.,5.,10.,20.,
#                                        50.,100.,200.,500.,1000.,1100.,1200.,1500.,
#                                        2000.,2100.]

    cfres.lbOrientation              = "Horizontal"
    cfres.lbTitleString              = "Surface oil concentration, Kg/m^2"
    cfres.lbTitleFontHeightF        = 0.012
    cfres.lbLabelFontHeightF        = 0.008
    cfres.lbTitleOffsetF            = -0.27
    cfres.lbBoxMinorExtentF         = 0.15
    cfres.pmLabelBarOrthogonalPosF  = -0.01


# 9. SAMPLE MARKER SETUP --------------------------------------------------------------------

    m1res.gsMarkerSizeF             = 8.0
    m1res.gsMarkerColor             = "red"
    m1res.gsMarkerIndex             = 16

    m2res.gsMarkerSizeF             = 8.0
    m2res.gsMarkerColor             = "red"
    m2res.gsMarkerIndex             = 16

# 10. COASTAL POLYLINE SETUP ----------------------------------------------------------------

    cstres.gsLineThicknessF         = 8.0
    cstres.gsLineColor              = "black"
    cstres.gsFillColor              = True

# 11. Draw VECTOR and SCALAR FIELD ----------------------------------------------------------
    
    cf    = Ngl.contour(wks,C1,cfres)
    vccur = Ngl.vector(wks,ut,vt,c_vcres)
    if iwave == "Y":
       vcwav = Ngl.vector(wks,uwt,vwt,w_vcres)
    vcwnd = Ngl.vector(wks,winx2D,winy2D,wd_vres)

# 12. MAP SETUP -----------------------------------------------------------------------------

    xs = Ngl.get_float(cf.sffield,"sfXCActualStartF") - 0.2
    xe = Ngl.get_float(cf.sffield,"sfXCActualEndF") + 0.2
    ys = Ngl.get_float(cf.sffield,"sfYCActualStartF") - 0.2
    ye = Ngl.get_float(cf.sffield,"sfYCActualEndF") + 0.2

#    xs  = 9.7336 
#    xe  = 10.1997
#    ys  = 42.7996 
#    ye  = 43.1503

    mpres.mpProjection                = "CylindricalEquidistant"
    mpres.mpDataBaseVersion           = "HighRes"
    mpres.mpLimitMode                 = "LatLon"
    mpres.mpMinLonF                   = xs
    mpres.mpMaxLonF                   = xe
    mpres.mpMinLatF                   = ys
    mpres.mpMaxLatF                   = ye
    mpres.mpPerimOn                   = True    # Turn on map perimeter.
    mpres.mpGridAndLimbOn             = True
    mpres.mpPerimDrawOrder            = "PostDraw"
    mpres.mpFillDrawOrder             = "PostDraw"
    mpres.mpOutlineBoundarySets       = "GeophysicalAndUSStates"
    mpres.mpGeophysicalLineThicknessF = 5.0          # thickness of outlines
    mpres.mpFillOn                    = True
#    mpres.mpFillColors               = ["background","transparent","LightGray","LightGray"]   # Fill land and inland water 
                                                                                               # and leave ocean transparent
    mpres.mpFillColors               = ["background","transparent","transparent","transparent"]

    mpres.pmTitleDisplayMode         = "Always"       # Turn on map title.
    mpres.pmTickMarkDisplayMode      = "Always"
    mpres.pmLabelBarOrthogonalPosF   = -0.05

# 13. TITLE resources ------------------------------------------------------------------------

    mpres.tiMainString               = "Surface OIL Concentration"
    mpres.tiMainFontHeightF          = 0.020
    mpres.tiMainFont                 = "Helvetica-bold"
    mpres.tiMainOffsetYF             = 0.025

# 14. PLOT THE MAP ---------------------------------------------------------------------------

    mp = Ngl.map(wks,mpres)                           # Draw a plot of map

# 15. OVERLAY EVERYTHING ON THE MAP PLOT ------------------------------------------------------

    Ngl.overlay(mp,cf)
    
    Ngl.overlay(mp,vccur)
    if iwave == "Y":
       Ngl.overlay(mp,vcwav)
    Ngl.overlay(mp,vcwnd)

    Ngl.maximize_plot(wks,mp) # Maximize size of plot in frame.

# 16. IF PRESENT, PLOT SAMPLE OBSERVATIONS ----------------------------------------------------

#    if t == 3:
#       for i in range(0,4):
#           Ngl.add_polymarker(wks,mp,m_array[i,1],m_array[i,0],m1res)

#    if t == 6:
#       for i in range(0,9):
#           Ngl.add_polymarker(wks,mp,d_array[i,1],d_array[i,0],m2res)


# 17. ADD THE COASTAL SEGMENTS IF PRESENT ------------------------------------------------------

    for i in range(0,len(CST1)):
        if CST1[i] != 0.:
           cstres.gsFillIndex = CST1[i]
           x1 = seg[i,0]
           x2 = seg[i,2]
           y1 = seg[i,1]
           y2 = seg[i,3]
           Ngl.add_polyline(wks,mp,[x1,x2],[y1,y2],cstres)

# 18. MAXIMIZE SIZE OF PLOT IN FRAME -------------------------------------------------------------

    Ngl.draw(mp)

# 19. Draw text strings after the plot is drawn to make sure plot gets maximized properly --------

    txres               = Ngl.Resources()
    txres.txFontHeightF = 0.015
    
    count_date = t + 6
    if (count_date <= 24):
       string1 = "17-05-2014"
    else:
       string1 = "18-05-2014"
    if (t < 10) :
        string2 = "+ 000" + str(t)
    else:
        string2 = "+ 00" + str(t)
    string = string1 + " 05:38 UTC " + string2 + " hours "
    Ngl.text_ndc(wks,string, 0.335, 0.920, txres)
#    Ngl.text_ndc(wks,string2, 0.85, 0.915, txres)
#    Ngl.text_ndc(wks,"Max Elevation: 4322", 0.85, 0.775, txres)

    Ngl.frame(wks)
    Ngl.destroy(wks)

Ngl.end()
