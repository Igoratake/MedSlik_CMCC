; The main NCL script for MEDSLIK results visualisation
;
;-----------------------------------------------------------------------------------
;  Copyright (C) <2012>
;  This program was originally written
;  by INGV, Slava Lyubartsev, 2011.02.22 
;----------------------------------------------------------------------------------
;  The development of the MEDSLIK-II model is supported by a formal agreement
;  Memorandum of Agreement for the Operation and Continued Development of MEDSLIK-II
;  signed by the following institutions:
;  INGV - Istituto Nazionale di Geofisica e Vulcanologia
;  OC-UCY - Oceanography Center at the University of Cyprus
;  CNR-IAMC - Consiglio Nazionale delle Ricerche – Istituto per 
;  lo Studio dell’Ambiente Marino Costiero
;  CMCC - Centro Euro-Mediterraneo sui Cambiamenti Climatici
; 
;  This program is free software: you can redistribute it and/or modify
;  it under the terms of the GNU General Public License as published by
;  the Free Software Foundation, either version 3 of the License, or
;  any later version.
;
;  This program is distributed in the hope that it will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;  You should have received a copy of the GNU General Public License
;  along with this program.  If not, see <http://www.gnu.org/licenses/>.
;-----------------------------------------------------------------------------------
; Must be launched from medslik_plots.sh
; Environmental variable MEDSLIK_OUTPUT_DIR is necessary

loadscript("lib.ncl") ; common library, time conversion, benchmarking, etc

INP  = "medslik5.inp"
FTE  = "medslik.fte"
OBS  = "observation_"
INI  = "initial.txt"
TXT  = "medslik_inputfile.txt"
HTML = "index.html"
HEAD = "index_head.html"
ANM = "animation"
OIL = "oil_"
SLK = "slicks"
SRF = ".srf"
CST = ".cst"
DSP = ".dsp"
ECM = "ecm"

DIR = getenv("MEDSLIK_OUTPUT_DIR")

_check_file_exists(DIR+TXT)

PLOT_DIR = getenv("PLOT_DIR")

if ismissing(PLOT_DIR) then
	PLOT_DIR = DIR+"plots/"
else
	PLOT_DIR = DIR+"plots/"+PLOT_DIR
end if

if .not.isfilepresent(DIR+"plots") then
	system("mkdir "+DIR+"plots")
end if

if .not.isfilepresent(PLOT_DIR) then
	system("mkdir "+PLOT_DIR)
end if

TASK_NCL = getenv("TASK_NCL")

if ismissing(TASK_NCL) then
	TASK_NCL = DIR+"medslik_plots.ncl"
	if .not.isfilepresent(TASK_NCL) then
		system("cp default_medslik_plots.ncl "+TASK_NCL)
	end if
end if

COAST_MAP = True
COAST_MAP_OIL_COLOR = (/0.3,0.3,0.3/)
COAST_MAP_THICKNESS = 4.0

loadscript(TASK_NCL)  ; visualisation parameters

if (.not.MAIN_MAPS).and.(.not.CURRENTS_MAPS).and.(.not.SLICK_MAP) then
	print("")
	print("No maps assigned to plot (main, currents, slicks)")
	print("See file "+TASK_NCL)
	print("")
	exit
end if

WIND     = (MAIN_MAPS.and.MAIN_MAPS_WIND)
CURRENTS = (MAIN_MAPS.and.MAIN_MAPS_CURRENTS).or.CURRENTS_MAPS

if WIND then
	if isfilepresent(DIR+"E25/") then
		WIND_DIR = DIR+"E25/"
	else
		if isfilepresent(DIR+"ECM/") then
			WIND_DIR = DIR+"ECM/"
		else
			print("No wind data directory found")
			print("**** procedure aborted ****")
			exit
		end if
	end if
end if

if CURRENTS then
	modelMFS = 0
	modelAFS = 1
        modelTIRR = 2
        modelWME = 3
        modelSCRM = 4
        modelMFS24 = 5
	MODEL_NAMES = (/"MFS","AFS","TYRRHENIAN","WEST_MED","SICILY","MFS24"/)
	
	MODEL_DIR = ""
	
	if isfilepresent(DIR+"A1h") then
		MODEL = modelAFS
		MODEL_DIR = "A1h/adri"
		MODEL_EXT = ".adr"
	end if
	if isfilepresent(DIR+"O1h") then
		MODEL = modelMFS
		MODEL_DIR = "O1h/medf"
		MODEL_EXT = ".opa"
	end if
        if isfilepresent(DIR+"T1h") then
		MODEL = modelTIRR
		MODEL_DIR = "T1h/tyrr"
		MODEL_EXT = ".tyr"
	end if
        if isfilepresent(DIR+"WME") then
		MODEL = modelWME
		MODEL_DIR = "WME/wmed"
		MODEL_EXT = ".wme"
	end if
        if isfilepresent(DIR+"S1h") then
		MODEL = modelSCRM
		MODEL_DIR = "S1h/sici"
		MODEL_EXT = ".sic"
	end if
        if isfilepresent(DIR+"OPA") then
		MODEL = modelMFS24
		MODEL_DIR = "OPA/medf"
		MODEL_EXT = ".opa"
	end if
	if isvar("MODEL") then 
		print("")
		print("model = "+MODEL_NAMES(MODEL))
		print("")
	else
		print("")
		print("No currents directory found")
		print("**** procedure aborted ****")
		exit
	end if

MODEL_DIR = DIR+MODEL_DIR

varU = 0
varV = 1
varW = 2

end if

str = "hours to process = from "+HOUR_BEGIN+" to "
if HOUR_END.le.0 then
	str = str+"the end"
else
	str = str+HOUR_END
end if
print(str+", step "+HOUR_STEP)
print("")

; initial position type
INI_POINT = 0
INI_SAT   = 1
INI_POLY  = 2
INI_TYPE  = INI_POINT
NSLICK = 0

loadscript("files.ncl")        ; read MEDSLIK output metadata, and check all files exist

loadscript("input_file.ncl")   ; read medslik_input.txt

loadscript("observations.ncl") ; read obseravation_YYYYMMDDHHNN.txt files if exist

loadscript("oil.ncl")          ; read oil files

loadscript("obs_oil.ncl")

if CURRENTS.or.WIND then
	loadscript("grid.ncl")
end if

if CURRENTS then
	loadscript("model.ncl")
end if

if WIND then
	loadscript("wind.ncl")
end if

loadscript("graphics.ncl")

if SLICK_MAP then
	_plot_slick_map
end if

if NOBS.gt.0 then
	do iobs = 0,NOBS-1
		_plot_observation_map(iobs)
	end do
end if

if MAIN_MAPS then
	do itime = 0,ntime-1
		_test_plot_main_map(itime)
		if COAST_MAP.and.(REAL_CST_N.gt.0) then
			_test_plot_coast_map(itime)
		end if
;		_plot_main_map(itime)
	end do
	if ANIMATION then
		_animate_main_maps
	end if
	_make_html
end if

print("")
_log_finish
