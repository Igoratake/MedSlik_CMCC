#!/usr/bin/env python

# ============================================================================
#                           plot_traj_BEACON.py
#------------------------------------------------------------------------------

import numpy as np
import os
import Nio
import Ngl

drift_dir = "../RUN/drifters_observations/"

d1_file = drift_dir + "a300234061615020.txt"
d2_file = drift_dir + "d300034012659810.txt"

struct1 = np.dtype([('lat',np.float32),('lon',np.float32),('date', 'S13'),('hour','S8')])

data1 = open(d1_file).readlines()
drift1 = np.loadtxt(d1_file,dtype=struct1)

data2 = open(d2_file).readlines()
drift2 = np.loadtxt(d2_file,dtype=struct1)

for i in xrange(2):

    if i == 0:
       out_dir = "../OUTPUT/EULER_SG_drifter_test_case/"
       out_name         = "eul_ISPHERE_trajectory"
       SIM = [1,2,3]
       LEG = [0,1,2,3]
    if i == 1:
       out_dir = "../OUTPUT/RUNGE_SG_drifter_test_case/"
       out_name         = "run_ISPHERE_trajectory"
       SIM = [1,3]
       LEG = [0,1,3]
# Plot -------------------------

# 1. Open the WorKStation

    rlist            = Ngl.Resources()
    rlist.wkColorMap = "rainbow+gray" #cmap
    wks_type = "png"
    wks = Ngl.open_wks(wks_type,out_name,rlist)
   
# 2. Set up resource lists will apply to vector, line contour, and
#    filled contour plots.

    mpres   = Ngl.Resources()     # map resources
    plres   = Ngl.Resources()     # polyline resources

# 3. Turn off nglDraw and nglFrame because we don't want to draw all
#    these plots until they are all overlaid on the map plot.

    mpres.nglDraw    = False
    mpres.nglFrame   = False

    plres.nglDraw    = False
    plres.nglFrame   = False

# 4. POLYLINE SETUP

    plres.gsLineThicknessF     = 15.0

# 5. MAP SETUP AND DRAW

    xs = float(9.78)
    xe = float(10.20)
    ys = float(42.92)
    ye = float(43.61)

    mpres.tiMainString               = "ISPHERE"
    mpres.tiMainFontHeightF          = 0.020
    mpres.tiMainFont                 = "Helvetica-bold"
    mpres.tiMainOffsetYF             = 0.025

    mpres.mpProjection               = "CylindricalEquidistant"
    mpres.mpDataBaseVersion          = "HighRes"
    mpres.mpLimitMode                = "LatLon"
    mpres.mpMinLonF                  = xs
    mpres.mpMaxLonF                  = xe
    mpres.mpMinLatF                  = ys
    mpres.mpMaxLatF                  = ye
    mpres.mpLabelsOn                 = True
    mpres.mpOutlineOn                = True
    mpres.mpPerimOn                  = True    # Turn on map perimeter.
    mpres.mpGridAndLimbOn            = True
    mpres.mpPerimDrawOrder           = "PostDraw"
    mpres.mpGridAndLimbDrawOrder     = "PostDraw"
    mpres.mpFillDrawOrder            = "PostDraw"
    mpres.mpOutlineBoundarySets      = "GeophysicalAndUSStates"
    mpres.mpFillOn                   = True
    mpres.mpFillColors               = ["background","transparent","LightGray","LightGray"] # Fill land and inland water 
                                                                                            # and leave oceans transparent.
    mpres.tmYROn                     = True
    mpres.tmXTOn                     = True
    mpres.lbOrientation              = "Horizontal"
    mpres.pmTitleDisplayMode         = "Always"       # Turn on map title.
#    mpres.pmTickMarkDisplayMode      = "Never"
    mpres.pmLabelBarOrthogonalPosF   = -0.05

    mp = Ngl.map(wks,mpres)                           # Draw a plot of map

# 6. TRAJECTORY PLOT

    plres.gsLineColor           = "black"
    for j in xrange(len(data1)-1):
        Ngl.add_polyline(wks,mp,[drift1[j][1],drift1[j+1][1]],[drift1[j][0],drift1[j+1][0]],plres)

    plres.gsLineColor           = "black"
    for j in xrange(len(data2)-1):
        Ngl.add_polyline(wks,mp,[drift2[j][1],drift2[j+1][1]],[drift2[j][0],drift2[j+1][0]],plres)
# -------------------------------------------------------------------------------------------------------
    
    for s in SIM:
 
        if s == 1:
           exp_dir = "cur"
        if s == 2:
           exp_dir = "cur_stok"
        if s == 3:
           exp_dir = "cur_wav" 

        exp= [1,4]
        for drft in exp:
            if drft == 1:
               start = "1138"
            if drft == 4:
               start = "1139"

            m_name = out_dir+"MEDSLIKII_2014_05_17_"+start+"_test_case_drft_"+str(drft)+"_v1.02_"+exp_dir+"/dr"+str(drft)+"_140517_"+start+"_F.trj"
            print m_name
            if s == 1: 
               plres.gsLineColor           = 224 #red
               plres.gsLineThicknessF      = 8
               plres.gsLineDashPattern     = 0
            if s == 2: # blu
               plres.gsLineColor           = 48
               plres.gsLineThicknessF      = 8
               plres.gsLineDashPattern     = 0 
            if s == 3: # green
               plres.gsLineColor           = 128
               plres.gsLineThicknessF      = 8
               plres.gsLineDashPattern     = 0

            traj = np.zeros((24,2))        
            data = open(m_name).readlines()

            traj[0][1] = 42.96907
            traj[0][0] = 9.98350
       
            for t in range(1, 24):
                model = data[t+1]
#                print model 
                lat_m = float(model.split('    ')[1])
#                print lat_m
                lon_m = float(model.split('    ')[2])
#                print lon_m
                traj[t][0] = lon_m
                traj[t][1] = lat_m
            for t in xrange(len(traj)-1):
                Ngl.add_polyline(wks,mp,[traj[t][0],traj[t+1][0]],[traj[t][1],traj[t+1][1]],plres)
 
#--------------------------------------------------------------------------------------------------------
# LEGEND

    x_l1 = 10.02
    x_l2 = 10.08
    y_l1 = 43.55
    y_l2 = 43.55

    txres = Ngl.Resources()
    txres.txFontHeightF = 0.015
    delx = 0.05
    txres.txFontColor = "black"
    x_l  = 10.08
    y_l  = 43.55

    for s in LEG:

         if s == 0:
            exp_dir = "I-SPHERE"
            plres.gsLineColor           = 2
            plres.gsLineThicknessF      = 15
            plres.gsLineDashPattern     = 0
         if s == 1:
            exp_dir = "cur"
            plres.gsLineColor           = 224 #red
            plres.gsLineThicknessF      = 8
            plres.gsLineDashPattern     = 0
         if s == 2:
            exp_dir = "cur_stok"
            plres.gsLineColor           = 48 #blue
            plres.gsLineThicknessF      = 8
            plres.gsLineDashPattern     = 0
         if s == 3:
            exp_dir = "cur_wav"
            plres.gsLineColor           = 128 #green
            plres.gsLineThicknessF      = 8
            plres.gsLineDashPattern     = 0
   
         Ngl.add_polyline(wks,mp,[x_l1,x_l2],[y_l1,y_l2],plres)
         y_l1 = y_l1 - 0.0150
         y_l2 = y_l2 - 0.0150
         Ngl.add_text(wks, mp, exp_dir, x_l + delx, y_l, txres)
         y_l = y_l - 0.0150
#--------------------------------------------------------------------------------------------------------
    Ngl.maximize_plot(wks,mp)                         # Maximize size of plot in frame.
    Ngl.draw(mp)
    Ngl.frame(wks)
    Ngl.destroy(wks)

Ngl.end()

