#!/bin/ksh -x

raw_dir="$SCRATCH/data/input/ssp126/rawfiles/"
ave_dir="$SCRATCH/data/input/ssp126/averages/"

nco_dir="/usr/local/apps/nco/4.3.7/bin/"
ncra=${nco_dir}"ncra"
ncrcat=${nco_dir}"ncrcat"

netcdf_dir="/opt/cray/netcdf/4.3.1/bin/"
nccopy=${netcdf_dir}"nccopy"

varlist=$1 #"precip snowfall evap tskin sndiv snowmelt"

temp1=${raw_dir}"temp1.nc"
temp2=${raw_dir}"temp2.nc"

for varname in $varlist; do
  (( year = 1950 ))
  while [ $year -lt 1970 ]; do
    echo $year

    fname=${raw_dir}${varname}"_FGRN11_forFDM_Year"${year}".nc"
    fname_yave=${raw_dir}${varname}"_Yave_"${year}".nc"
    ${ncra} ${fname} ${fname_yave}
   
    (( year = year + 1 ))
  done

  fname_out=${ave_dir}${varname}"_FGRN11_50-70_ave.nc"
  
  ${ncrcat} ${raw_dir}${varname}"_Yave"*".nc" ${temp1}
  ${ncra} ${temp1} ${temp2}
  ${nccopy} -k classic ${temp2} ${fname_out}
  
  rm ${raw_dir}${varname}"_Yave"*".nc" ${temp1} ${temp2}

done
