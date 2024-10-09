#!/bin/ksh -x

## make average inputs for FDM                                ##
## requirements: `makeFDMyears_FGRN055.sc` must be run first  ##
## ---------------------------------------------------------- ##

## one var per job, set in `all_*`                            ##
## will be one of: precip snowfall evap tskin sndiv snowmelt  ##
## ---------------------------------------------------------- ##
varname=$1 

## set project name & filepaths                               ##
## ---------------------------------------------------------- ##
project_name=$2
base_dir=$3
years_dir=$4
ave_dir=$5

## set start and end years for averaging                      ##
## ---------------------------------------------------------- ##
start_year=${10} #1960
end_year=${11} #1981

## takes averages across year files                           ##
## ---------------------------------------------------------- ##

temp1="${years_dir}temp1.nc"
temp2="${years_dir}temp2.nc"

(( year = $start_year ))
while [ $year -lt $end_year ]; do
  echo $year

  fname="${years_dir}${varname}_${project_name}_forFDM_Year${year}.nc"
  fname_yave="${years_dir}${varname}_Yave_${year}.nc"

  ncra ${fname} ${fname_yave}
 
  (( year = year + 1 ))
done

fname_out="${ave_dir}${varname}_${project_name}_${start_year}-${end_year}_ave.nc"

ncrcat ${years_dir}${varname}"_Yave"*".nc" ${temp1}
ncra ${temp1} ${temp2}
nccopy -k classic ${temp2} ${fname_out}

rm ${years_dir}${varname}"_Yave"*".nc" ${temp1} ${temp2}