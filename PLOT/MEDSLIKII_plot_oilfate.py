#!/usr/bin/env python

# ============================================================================
#                            PLOT_OILFATE.PY
#
# This software read NetCDF Medslik-II outputs and plots:
#
#        1. Oil concentration fields on a 150x150 m regular grid
#
#        2. Ocean currents velocities fields on the original ocean currents grid
#
#        3. Stokes' Drift velocities fields on the original waves grid
#
#------------------------------------------------------------------------------

import sys
import numpy as np
import Ngl

dir  = sys.argv[1]
file = dir + "medslik.fte"

fte_file  = open(file).readlines()

x = np.zeros((len(fte_file)-7))
y = np.zeros((4,len(fte_file)-6))
#y = np.zeros((len(fte_file)-7))

for i in range(7,len(fte_file)):
    x[i-7] = float(fte_file[i].split()[0])
#    y[i-7] = float(fte_file[i].split()[2])
    y[0][i-7] = float(fte_file[i].split()[2])
    y[1][i-7] = float(fte_file[i].split()[3])
    y[2][i-7] = float(fte_file[i].split()[5])
    y[3][i-7] = float(fte_file[i].split()[7])

out_name = "oilfate"
wks_type = "png"
wks = Ngl.open_wks(wks_type,out_name)

res = Ngl.Resources()
res.xyLineColors        = ["green","red","black","blue"]  # Define line colors.
res.xyLineThicknesses   = [5.,5.,5.,5.]    # Define line thicknesses
                                              # (1.0 is the default).

res.tiMainString    = "Oil Fate Parameters"  # Title for the XY plot
res.tiXAxisString   = "Hours after the simulation start"    # Label for the X axis
res.tiYAxisString   = "Percentage"    # Label for the Y axis
res.tiMainFont      = "Helvetica" # Font for title
res.tiXAxisFont     = "Helvetica" # Font for X axis label
res.tiYAxisFont     = "Helvetica" # Font for Y axis label
res.tiXAxisFontHeightF     = 0.02        # Change the font size.
res.tiYAxisFontHeightF     = 0.02

res.pmLegendDisplayMode    = "Always"     # Turn on the drawing
res.pmLegendZone           = 0            # Change the location
res.pmLegendOrthogonalPosF = 0.31         # of the legend
res.lgJustification        = "BottomRight"
res.pmLegendWidthF         = 0.3          # Change width and
res.pmLegendHeightF        = 0.10         # height of legend.
res.pmLegendSide           = "Top"        # Change location of
res.lgPerimOn              = True        # legend and turn off
                                                  # the perimeter.

res.xyExplicitLegendLabels = ["% of OIL evaporated","% of OIL on the Sea Surface","% of OIL dispersed in the WATER Column","% of OIL on the Coast"]


plot = Ngl.xy(wks,x,y,res)         
