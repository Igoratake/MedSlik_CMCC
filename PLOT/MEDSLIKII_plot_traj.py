#!/usr/bin/env python

# ============================================================================
#                            MEDESS_PLOT_TRAJ.PY
#------------------------------------------------------------------------------

import sys
import numpy as np
import os
import Nio
import Ngl

model_dir = sys.argv[1]
exp_name  = sys.argv[2]
max_time = int(sys.argv[3]) + 1

# Plot -------------------------

# 1. Open the WorKStation

out_name = "trajectory"
rlist            = Ngl.Resources()
rlist.wkColorMap = "WhBlGrYeRe" #cmap
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

plres.gsLineThicknessF        = 8.0

# 5. MAP SETUP AND DRAW

limit_file = model_dir + "medslik.tmp"
coord = open(limit_file).readlines()

xs = float(coord[1].split('   ')[1])
xe = float(coord[1].split('   ')[4])
ys = float(coord[2].split('   ')[1])
ye = float(coord[2].split('   ')[4])

mpres.tiMainString               = "MEDSLIK-II trajectories simulation"
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
mpres.mpFillColors               = ["background","transparent","LightGray","LightGray"]   # Fill land and inland water 
                                                                                          # and leave oceans transparent.
mpres.tmYROn                     = True
mpres.tmXTOn                     = True 
mpres.lbOrientation              = "Horizontal"
mpres.pmTitleDisplayMode         = "Always"       # Turn on map title.
#mpres.pmTickMarkDisplayMode      = "Never"
mpres.pmLabelBarOrthogonalPosF   = -0.05

mp = Ngl.map(wks,mpres)                           # Draw a plot of map

# -------------------------------------------------------------------------------------------------------

plres.gsLineColor           = "red"


m_name = model_dir + exp_name
traj = np.zeros((max_time,2))        
data = open(m_name).readlines()

traj[0][1] = 42.97053
traj[0][0] = 9.98382
    
for t in range(1, max_time):
    model = data[t+1]
    print model
    lat_m = float(model.split('    ')[1])
    print lat_m
    lon_m = float(model.split('    ')[2])
    print lon_m
    traj[t][0] = lon_m
    traj[t][1] = lat_m

for t in xrange(len(traj)-1):
    Ngl.add_polyline(wks,mp,[traj[t][0],traj[t+1][0]],[traj[t][1],traj[t+1][1]],plres)


#--------------------------------------------------------------------------------------------------------
Ngl.maximize_plot(wks,mp)                         # Maximize size of plot in frame.
Ngl.draw(mp)
Ngl.frame(wks)
Ngl.destroy(wks)

Ngl.end()

