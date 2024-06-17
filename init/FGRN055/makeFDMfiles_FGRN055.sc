#!/bin/ksh -x

## slices yearly var into longitude strips and                  ##
## sews together longitude files into time series               ##
## ------------------------------------------------------------ ##

## set project name & filepaths                                 ##
## ------------------------------------------------------------ ##
#1-script 2-var 3-project_name 4-RACMO_3H 5-years_dir 6-ts_dir 7-ave_dir
project_name=$2
base_dir=$3
years_dir="${base_dir}/raw/historical/years/"
files_dir="${base_dir}/input/timeseries/"

## set variable in `all_*`                                      ##
## will be one of: precip snowfall evap tskin sndiv snowmelt    ##
## ------------------------------------------------------------ ##
varlist=$1

## set start and end years
## TK: update to 2023?
## ------------------------------------------------------------ ##
start_year= $4 #1957
end_year= $5 #2020

## sets number of longitude bands and cell-width of strips      ##
## TK: cell_width must match openNetCDF pref; make auto
## ------------------------------------------------------------ ##
num_long_bands=$6 #74
cell_width=$7 #5

## create timeseries for each variable
## ------------------------------------------------------------ ##
for varname in $varlist; do
  temp1=${years_dir}${varname}"_temp1.nc"

  (( year = ${start_year} ))
  while [ $year -le ${end_year} ]; do
    # echo $varname" "$year
    fname_in="${years_dir}${varname}_${project_name}_forFDM_Year${year}.nc"

    ## loop through all years
    ## ------------------------------------------------------------ ##
    (( part = 1 ))
    (( start_t = 0 ))
    (( end_t = -1 ))
    
    ## loop through all longitude bands
    ## ------------------------------------------------------------ ##
    while [ $part -lt ${num_long_bands} ]; do
      (( start_t = end_t + 1 ))
      (( end_t = start_t + ${cell_width} )) 
      fname_out="${years_dir}${varname}_${year}_part${part}.nc"
      
      ## check if file exists, otherwise slices yearlies into bands   ##
      ## ------------------------------------------------------------ ##
      if [[ -f ${fname_out} ]]; then
        echo "File is already present"
      else    
        ncks -d rlon,${start_t},${end_t} ${fname_in} ${fname_out}
      fi
      
      (( part = part + 1 ))
    done 
    (( year = year + 1 ))
  done
  
  (( part2 = 1 ))

  ## sew yearly files, sliced by longitude, into single timeseries
  ## ------------------------------------------------------------ ##
  while [ $part2 -lt ${num_long_bands} ]; do
    fname_final="${files_dir}${varname}_${project_name}_${start_year}-${end_year}_p${part2}.nc"
    touch $fname_final

      if [[ -f ${fname_final} ]]; then
        echo "File is already present"
      else  
        ncrcat ${years_dir}${varname}"_"*"_part"${part2}".nc" ${temp1}
        nccopy -k classic ${temp1} ${fname_final} 
        rm ${temp1}
      fi  
    
    (( part2 = part2 + 1 ))  
  done
  
  rm ${years_dir}${varname}"_"*"part"*".nc"
done
