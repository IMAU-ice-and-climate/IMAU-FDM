; ***********************************************
; Written by Stefan Ligtenberg
; Adapted by Peter Kuipers Munneke, June 2012
; Adapted by Max Brils, July 2021
; Adapted by Elizabeth Case, 2024
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin
;************************************************
; read in data
;************************************************

variab = "zs" ;"FirnAir"  ;"zs"   ;"FirnAir"
dt = 0      ;(0=no detrending, 1=detrending)
atzero = 0  ;(0=don't start at 0, 1=force variab start at 0)

project_name="run-1957-2023-FGRN055-era055"

dir_pointlist = "$SCRATCH/FGRN055_era055/reference/"
dir_mask =  "$SCRATCH/FGRN055_era055/mask/"
dir_out = "$SCRATCH/run-1957-2023-FGRN055-era055/"
fname_out = "_FGRN055_1957-2023_GrIS_GIC"

; variables concerning files run
name_pointlist1 = dir_pointlist + "IN_ll_FGRN055_GrIS_GIC_implicit.txt"
dir_in1 = "$SCRATCH/"+project_name+"/output/"
name_input1 = "FGRN055_era055"
ncols1 = numAsciiCol(name_pointlist1)
Grid1 = readAsciiTable(name_pointlist1,ncols1,"float",0)
dim_grid1 = dimsizes(Grid1)
numpoints1 = dim_grid1(0)

;****************************************
;****************************************


;****************************************
;****************************************

; load one file to extract the dimensions from
fname = dir_in1 + name_input1 + "_1D_1.nc"

do nnn = 0, numpoints1-1
 numb = nnn + 1
 fname = dir_in1 + name_input1 + "_1D_" + numb + ".nc"

 if (isfilepresent(fname)) then
  ec_out = addfile(fname,"r")


 else

  print(numb)
 
 end if
end do


;***************************************





end