#!/bin/ksh -x

raw_dir="$SCRATCH/data/input/ssp585/rawfiles/"
ave_dir="$SCRATCH/data/input/ssp585/averages/"
files_dir="$SCRATCH/data/input/ssp585/files/"

varlist=$1 #"precip snowfall evap tskin sndiv snowmelt"

for varname in $varlist; do
  temp1=${raw_dir}${varname}"_temp1.nc"

  (( year = 1950 ))
  while [ $year -lt 2100 ]; do
    echo $varname" "$year
    fname_in=${raw_dir}${varname}"_FGRN11_forFDM_Year"${year}".nc"

    (( part = 1 ))
    (( start_t = 0 ))
    (( end_t = -1 ))
    while [ $part -le 65 ]; do
      (( start_t = end_t + 1 ))
      (( end_t = start_t + 4 ))    
      fname_out=${raw_dir}${varname}"_"${year}"_part"${part}".nc"
      
      ncks -d rlon,${start_t},${end_t} ${fname_in} ${fname_out}
      
      (( part = part + 1 ))
    done 
    (( year = year + 1 ))
  done
  
  (( part2 = 1 ))
  while [ $part2 -le 65 ]; do
    fname_final=${files_dir}${varname}"_FGRN11_1950-2099_p"${part2}".nc"
  
    ncrcat ${raw_dir}${varname}"_"*"_part"${part2}".nc" ${temp1}
    nccopy -k classic ${temp1} ${fname_final} 
    rm ${temp1}
    
    (( part2 = part2 + 1 ))  
  done
  
  rm ${raw_dir}${varname}"_"*"part"*".nc"
done

