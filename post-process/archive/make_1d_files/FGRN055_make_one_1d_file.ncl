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

;variab = "zs" ;"FirnAir"  ;"zs"   ;"FirnAir"
;dt = 0      ;(0=no detrending, 1=detrending)
;atzero = 1  ;(0=don't start at 0, 1=force variab start at 0)

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

printVarSummary(Grid1)
print(numpoints1+0)

;****************************************
;****************************************

fname_lsm = dir_mask + "FGRN055_Masks.nc"
model_lsm = addfile(fname_lsm,"r")
LSM = model_lsm->IceMask
Lat = model_lsm->lat
Lon = model_lsm->lon
dimsize_lsm = dimsizes(LSM)
printVarSummary(LSM)
n_lat = dimsize_lsm(0)
n_lon = dimsize_lsm(1)

;****************************************
;****************************************

; load one file to extract the dimensions from
fname = dir_in1 + name_input1 + "_1D_1.nc"
forT = addfile(fname,"r")
Time_zs = forT->$variab$

ds_time = dimsizes(Time_zs)
ntime = ds_time(0)

print(ntime)
print(" ")

;****************************************
;****************************************

Matrix_055 = new((/ntime,n_lat,n_lon/),float)

Matrix_055 = -999.
Matrix_055@_FillValue = -999.

; calculate time of each time step
year = 1957.
month = 10.
day = 1.
TimeY = new((/ntime/),float)
do ttt = 0, ntime-1
 ;TimeY(ttt) = 1957. + ((month-1.)/12.) + ((day-1.)/365.25) + (ttt/ntime*(62.+3./12.))
 TimeY(ttt) = 1957. + ((month-1.)/12.) + ((day-1.)/365.25) + (ttt/ntime*(66.+3./12.))
end do 
print(TimeY(ntime-1)+0)
print(" ")

;*****************************************


do nnn = 0, numpoints1-1
 numb = nnn + 1
 fname = dir_in1 + name_input1 + "_1D_" + numb + ".nc"

 if (isfilepresent(fname)) then
  ec_out = addfile(fname,"r")
  Tussen = ec_out->$variab$(0:ntime-1)

  xx = floattointeger(Grid1(nnn,5))
  yy = floattointeger(Grid1(nnn,6))
  
  print(numb + " " + xx + " " + yy)
  
  if (atzero .eq. 1) then
    Tussen(:) = Tussen(:) - Tussen(0) 
  end if

  if (dt .eq. 0) then
   Matrix_055(:,xx,yy) = Tussen(:) 
  end if
 
  if (dt .eq. 1) then
   rc0 = (Tussen(ntime-1)-Tussen(0))/ntime
   do sss = 0, ntime-1
    Matrix_055(sss,xx,yy) = Tussen(sss) - rc0*sss
   end do
  end if   


 else

  print("------------------------------> Deze is er niet: " + numb)
  print("------------------------------> Versie: " + name_input1)
 
 end if
end do


;***************************************


fout = dir_out + "FDM_" + variab + fname_out + ".nc"
system("rm -rf " + fout)
ncdf_out = addfile(fout,"c")

dimNames = (/"time","rlat","rlon"/)  
dimSizes = (/ -1   , n_lat, n_lon/) 
dimUnlim = (/ True , False, False/)   
filedimdef(ncdf_out,dimNames,dimSizes,dimUnlim)

filevardef(ncdf_out,variab,typeof(Matrix_055),(/"time","rlat","rlon"/))
filevardef(ncdf_out,"lsm",typeof(LSM),(/"rlat","rlon"/))
filevardef(ncdf_out,"lat",typeof(Lat),(/"rlat","rlon"/))
filevardef(ncdf_out,"lon",typeof(Lon),(/"rlat","rlon"/))
filevardef(ncdf_out,"time",typeof(TimeY),(/"time"/))

if (variab.eq."FirnAir") then
  Matrix_055@long_name = "Firn air content"             
  Matrix_055@units     = "m"
else if (variab.eq."zs") then
  Matrix_055@long_name = "Elevation change"             
  Matrix_055@units     = "m"
else if (variab.eq."vice") then
  Matrix_055@long_name = "Ice velocity"             
  Matrix_055@units     = "m y-1"
else if (variab.eq."vacc") then
  Matrix_055@long_name = "Accumulation velocity"             
  Matrix_055@units     = "m y-1"
else if (variab.eq."vfc") then
  Matrix_055@long_name = "Firn compaction velocity"             
  Matrix_055@units     = "m y-1"
else if (variab.eq."vmelt") then
  Matrix_055@long_name = "Melt velocity"             
  Matrix_055@units     = "m y-1"
else if (variab.eq."vbuoy") then
  Matrix_055@long_name = "Buoyancy velocity"             
  Matrix_055@units     = "m y-1"
else if (variab.eq."vtotal") then
  Matrix_055@long_name = "Total velocity"             
  Matrix_055@units     = "m y-1"
else if (variab.eq."Runoff") then
  Matrix_055@long_name = "Runoff"             
  Matrix_055@units     = "mm w.e. y-1"
else if (variab.eq."vsub") then
  Matrix_055@long_name = "Sublimation velocity"             
  Matrix_055@units     = "m y-1"
else if (variab.eq."vsnd") then
  Matrix_055@long_name = "Snowdrift velocity"             
  Matrix_055@units     = "m y-1"
else if (variab.eq."totlwc") then
  Matrix_055@long_name = "Vertically integrated liquid water content"             
  Matrix_055@units     = "mm"
  
end if
end if
end if
end if
end if
end if
end if
end if
end if
end if
end if
end if

TimeY@long_name = "Time in years, starting on 01-10-1957"           ; assign attributes, DD-MM-YYYY 
TimeY@units     = "years" 

printVarSummary(Lat)
printVarSummary(Lon)
printVarSummary(TimeY)
printVarSummary(Matrix_055)

ncdf_out->time = TimeY
ncdf_out->lat = Lat
ncdf_out->lon = Lon
ncdf_out->lsm = LSM
ncdf_out->$variab$ = (/Matrix_055(:,:,:)/)



end