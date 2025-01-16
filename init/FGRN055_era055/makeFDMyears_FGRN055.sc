#!/bin/ksh -x

## makeFDMyears slices decadal RACMO output into year files     ##
## for each variable                                            ##
## ------------------------------------------------------------ ##

## 
varname=$1
project_name=$2
base_dir=$3
years_dir=$4
start_year=$9

RACMO_3H_dir="${base_dir}/raw/historical-${start_year}"

#for era5 runs
fname_extra_1=".FGRN055.BN_RACMO2.3p2_ERA5_3h_1940_FGRN055.3H.nc"
fname_extra_2=".FGRN055.BN_RACMO2.3p2_ERA5_3h_1940_FGRN055.3H.nc"
fname_extra_3=".FGRN055.BN_RACMO2.3p2_ERA5_3h_1940_FGRN055.3H.nc"

# for interim runs
#fname_extra_1=".FGRN055.BN_RACMO2.3p2_FGRN055.3H.nc"
#fname_extra_2=".FGRN055.BN_RACMO2.3p2_ERA5_3h_FGRN055.3H.nc"
#fname_extra_3=".FGRN055.BN_RACMO2.3p2_ERA5_3h_1940_FGRN055.3H.nc"

## sets number of data points in each year (#d/yr * #pts/d)     ##
## TK: update to 2023?                                          ##
## in korn, -A assigns parametr to Y01list                      ##
## https://www.mkssoftware.com/docs/man1/set.1.asp              ##
## ------------------------------------------------------------ ##

if [ "${start_year}" = 1957 ]; then

  set -A yearlist 1957 1961 1971 1981 1991 2001 2011 2021 2023
  set -A nyears   4    10   10   10   10   10   10   2    1

  set -A Y57list 736 2920 2920 2928
  set -A Y61list 2920 2920 2920 2928 2920 2920 2920 2928 2920 2920
  set -A Y71list 2920 2928 2920 2920 2920 2928 2920 2920 2920 2928
  set -A Y81list 2920 2920 2920 2928 2920 2920 2920 2928 2920 2920
  set -A Y91list 2920 2928 2920 2920 2920 2928 2920 2920 2920 2928 
  set -A Y01list 2920 2920 2920 2928 2920 2920 2920 2928 2920 2920
  set -A Y11list 2920 2928 2920 2920 2920 2928 2920 2920 2920 2928
  set -A Y21list 2920 2920 
  set -A Y23list 2920

  set -A data_in_years_list Y57list Y61list Y71list Y81list Y91list Y01list Y11list Y21list Y23list

elif [ "${start_year}" = 1939 ]; then

  set -A yearlist 1939 1941 1951 1961 1971 1981 1991 2001 2011 2021 2023
  set -A nyears   2    10   10   10   10   10   10   10   10   2    1

  set -A Y39list 976 2928
  set -A Y41list 2920 2920 2920 2928 2920 2920 2920 2928 2920 2920
  set -A Y51list 2920 2928 2920 2920 2920 2928 2920 2920 2920 2928
  set -A Y61list 2920 2920 2920 2928 2920 2920 2920 2928 2920 2920
  set -A Y71list 2920 2928 2920 2920 2920 2928 2920 2920 2920 2928
  set -A Y81list 2920 2920 2920 2928 2920 2920 2920 2928 2920 2920
  set -A Y91list 2920 2928 2920 2920 2920 2928 2920 2920 2920 2928 
  set -A Y01list 2920 2920 2920 2928 2920 2920 2920 2928 2920 2920
  set -A Y11list 2920 2928 2920 2920 2920 2928 2920 2920 2920 2928
  set -A Y21list 2920 2920 
  set -A Y23list 2920

  set -A data_in_years_list Y39list Y41list Y51list Y61list Y71list Y81list Y91list Y01list Y11list Y21list Y23list

else
  
  echo "Start year -- ${start_year} -- does not match either 1939 or 1957"

fi

## for each variable and each year, create year file            ##
## ------------------------------------------------------------ ##

(( years_i = ${#yearlist[@]} - 1 ))

for ii in {0..$years_i}; do

   (( fileyear = yearlist[ii] ))
  
  if [[ ${fileyear} -lt 1990 ]]; then
    infile="${RACMO_3H_dir}/${varname}.KNMI-${fileyear}${fname_extra_1}"
  elif [[ ${fileyear} -gt 1990 ]] && [[ ${fileyear} -le 2020 ]]; then
    infile="${RACMO_3H_dir}/${varname}.KNMI-${fileyear}${fname_extra_2}"
  elif [[ ${fileyear} -gt 2020 ]]; then 
    infile="${RACMO_3H_dir}/${varname}.KNMI-${fileyear}${fname_extra_3}"
  fi
  
   start_t=0
   end_t=-1

   (( nyears_jj = ${nyears[$ii]} - 1))

   for jj in {0..$nyears_jj}; do
     
     (( year = fileyear + jj ))
    
     outfile=${years_dir}/${varname}"_${project_name}_forFDM_Year"${year}".nc"
    
     list_of_years=${data_in_years_list[$ii]}

     num_data_in_years=$(eval echo \${${list_of_years}[$jj]})

     #echo "Year: " ${year}
     #echo "Num data points: " $num_data_in_years
     
     (( start_t = end_t + 1 ))
     (( end_t = start_t - 1 + num_data_in_years)) 

     echo ${varname} ${fileyear}" "${start_t}" "${end_t}" ..... "${year}
     if [[ -f ${outfile} ]]; then
       echo "File is already present"
     else
      #echo "ncks for "${varname}", year "${year}", yearlist " ${list_of_years} ", with " $num_data_in_years
      ncks -d time,${start_t},${end_t} ${infile} ${outfile}
     fi  
     
   done

 done
