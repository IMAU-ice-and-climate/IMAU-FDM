#!/bin/ksh -x

raw_dir="$SCRATCH/FGRN055_era055/raw/historical/years/"
parts_dir="$SCRATCH/FGRN055_era055/raw/historical/parts-temp/"
files_dir="$SCRATCH/FGRN055_era055/input/timeseries/"
project_name="FGRN055-era055"
start_year=1957
end_year=2023
num_rlon_strips=73

varlist=$1 #"precip snowfall evap tskin sndiv snowmelt"
# varlist="precip"

for varname in $varlist; do
  

  # (( year = ${start_year} ))
  # while [ $year -le ${end_year} ]; do
  #   echo $varname" "$year
  #   fname_in="${raw_dir}${varname}_${project_name}_forFDM_Year${year}.nc"

  #   (( part = 1 ))
  #   (( start_lon = 0 ))
  #   (( end_lon = -1 ))

  #   while [ $part -le ${num_rlon_strips} ]; do
  #     (( start_lon = end_lon + 1 ))
  #     (( end_lon = start_lon + 5 ))

  #     fname_out="${parts_dir}${varname}_parts/${varname}_${year}_part${part}.nc"
  #     if [[ -f ${fname_out} ]]; then
  #       echo "Part ${part} already exists"
  #     else
  #       ncks -O -d rlon,${start_lon},${end_lon} ${fname_in} ${fname_out}
  #     fi
      
  #     (( part = part + 1 ))
  #   done 
  #   (( year = year + 1 ))
  # done
  
  temp1="${files_dir}${varname}_temp1.nc"
  (( part2 = 1 ))
  while [ ${part2} -le ${num_rlon_strips} ]; do
    fname_part="${parts_dir}${varname}_parts/${varname}_*_part${part2}.nc"
    fname_final="${files_dir}${varname}_${project_name}_${start_year}-${end_year}_p${part2}.nc"
    
    if [[ -f ${fname_final} ]]; then
        echo "Part ${part2} is already present"
    else
      ncrcat ${fname_part} ${temp1}
      nccopy -k classic ${temp1} ${fname_final} 
      rm ${temp1}
    fi
    (( part2 = part2 + 1 ))  
  done
  
  #rm ${parts_dir}${varname}"_"*"part"*".nc"
done

# check to see that all parts from all years have been made
# ii=1
# while [ $ii -le ${num_rlon_strips} ]; do
#   echo "Part $ii"
#   ls -U1q evap*part${ii}.nc | wc -l 
#   (( ii = ii + 1 ))
# done