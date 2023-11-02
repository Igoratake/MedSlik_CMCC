loadscript("lib.ncl")

modelMFS = 0
modelAFS = 1

MAXN = 40000

PLOTS = "medslik.plots"
DIR = getenv("MEDSLIK_OUTPUT_DIR")
DIR_E25 = DIR+"E25/"
DIR_PLOTS = DIR+"plots/"
if .not.isfilepresent(DIR_PLOTS) then
	system("mkdir "+DIR_PLOTS)
end if
if .not.isfilepresent(DIR+PLOTS) then
	system("cp "+PLOTS+" "+DIR)
end if

DPI = 120

loadscript(DIR+PLOTS)

if isfilepresent(DIR+"A1h") then
	MODEL = modelAFS
	model_dir = "A1h"
	MODEL_EXT = "adr"
end if
if isfilepresent(DIR+"O1h") then
	MODEL = modelMFS
	model_dir = "O1h"
	MODEL_EXT = "opa"
end if
MODEL_DIR = DIR+model_dir+"/"

loadscript("oil.ncl")

loadscript("model.ncl")

loadscript("wind.ncl")

loadscript("plots.ncl")

_plot_map(22)
_plot_map(40)
exit

;_plot_oil_map(0)
;_plot_oil_map(6)
;_plot_oil_map(33)
;_plot_oil_map(66)
;_plot_oil_map(70)

do i = 0,NTIME-1
	_plot_oil_map(i)
end do

_log_finish

exit

