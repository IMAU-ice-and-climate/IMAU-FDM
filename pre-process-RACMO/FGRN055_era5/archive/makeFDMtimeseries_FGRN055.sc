#!/bin/ksh -x

## slices yearly var into longitude strips and                  ##
## sews together longitude files into time series               ##

## 1: variable 2: project name 3: base directory                ##
## 4: years_dir 5: files_dir 6: ave_dir                         ##
## 7: number of longitudinal bands 8: cell width                ##
## 9: ts start year 10: ts end year                             ##
## 11: avg start year 12: avg end year                          ##

## ------------------------------------------------------------ ##

## set start and end years
## TK: update to 2023?
## ------------------------------------------------------------ ##
start_year=$9
end_year=${10}

## set project name & filepaths                                 ##
## ------------------------------------------------------------ ##
project_name=$2
base_dir=$3
years_dir=$4
files_dir=$5

## set variable in `all_*`                                      ##
## will be one of: precip snowfall evap tskin sndiv snowmelt    ##
## ------------------------------------------------------------ ##
varlist=$1

## sets number of longitude bands and cell-width of strips      ##
## TK: cell_width must match openNetCDF pref; make auto
## ------------------------------------------------------------ ##
num_long_bands=$7
cell_width=$8

do_parts=1 # set to 0 if parts already created; set to 1 to recreate parts

## create timeseries for each variable
## ------------------------------------------------------------ ##
for varname in $varlist; do
  temp1="${years_dir}/parts/${varname}_temp1.nc"

  if (( do_parts -eq 1 )); then #only split years into parts if requested 
   
    (( year = ${start_year} ))
    
    while [ $year -le ${end_year} ]; do
      fname_in="${years_dir}/${varname}_${project_name}_forFDM_Year${year}.nc"

      ## loop through all years
      ## ------------------------------------------------------------ ##
      (( part = 1 ))
      (( start_t = 0 ))
      (( end_t = -1 ))
    
      ## loop through all longitude bands
      ## ------------------------------------------------------------ ##
      while [ $part -le ${num_long_bands} ]; do
        
        (( start_t = end_t + 1 ))
        (( end_t = start_t + ${cell_width} - 1)) 
        fname_out="${years_dir}parts/${varname}_${year}_part${part}.nc"
      
        ## check if file exists, then slice year files into bands       ##
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
  
  fi

  (( part2 = 1 ))

  ## sew slices from each year into single timeseries             ##
  ## ------------------------------------------------------------ ##
  
  while [ $part2 -le ${num_long_bands} ]; do
    fname_final="${files_dir}/${varname}_${project_name}_${start_year}-${end_year}_p${part2}.nc"

      if [[ -f ${fname_final} ]]; then
        echo "File is already present"
      else  
        ncrcat ${years_dir}"/parts/"${varname}"_"*"_part"${part2}".nc" ${temp1}
        nccopy -k classic ${temp1} ${fname_final} 
        rm ${temp1}
      fi  
    
    (( part2 = part2 + 1 ))  
  done
  
  #rm ${years_dir}"parts/"${varname}"_"*"part"*".nc"
done
