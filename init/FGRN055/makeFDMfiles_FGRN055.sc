#!/bin/ksh -x

raw_dir="$SCRATCH/data/input/era055_files/rawfiles/"
ave_dir="$SCRATCH/data/input/era055_files/averages/"
files_dir="$SCRATCH/data/input/era055_files/files/"

nco_dir="/usr/local/apps/nco/4.3.7/bin/"
ncks=${nco_dir}"ncks"
ncrcat=${nco_dir}"ncrcat"

netcdf_dir="/opt/cray/netcdf/4.3.1/bin/"
nccopy=${netcdf_dir}"nccopy"

varlist=$1 #"precip snowfall evap tskin sndiv snowmelt"

for varname in $varlist; do
  temp1=${raw_dir}${varname}"_temp1.nc"

  (( year = 1957 ))
  while [ $year -lt 2021 ]; do
    echo $varname" "$year
    fname_in=${raw_dir}${varname}"_FGRN055_forFDM_Year"${year}".nc"

    (( part = 1 ))
    (( start_t = 0 ))
    (( end_t = -1 ))
    while [ $part -lt 74 ]; do
      (( start_t = end_t + 1 ))
      (( end_t = start_t + 5 ))    
      fname_out=${raw_dir}${varname}"_"${year}"_part"${part}".nc"
      
      if [[ -f ${fname_out} ]]; then
        echo "File is already present"
      else      
        ${ncks} -d rlon,${start_t},${end_t} ${fname_in} ${fname_out}
      fi
      
      (( part = part + 1 ))
    done 
    (( year = year + 1 ))
  done
  
  (( part2 = 1 ))
  while [ $part2 -lt 74 ]; do
    fname_final=${files_dir}${varname}"_FGRN055_57-20_p"${part2}".nc"

      if [[ -f ${fname_final} ]]; then
        echo "File is already present"
      else  
        ${ncrcat} ${raw_dir}${varname}"_"*"_part"${part2}".nc" ${temp1}
        ${nccopy} -k classic ${temp1} ${fname_final} 
        rm ${temp1}
      fi	
    
    (( part2 = part2 + 1 ))  
  done
  
  rm ${raw_dir}${varname}"_"*"part"*".nc"
done
