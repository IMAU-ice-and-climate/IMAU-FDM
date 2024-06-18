#!/bin/ksh -x

## slices decadal racmo output into year files                  ##
## ------------------------------------------------------------ ##

project_name=$2 #"FGRN055-era055"
base_dir=$3
RACMO_3H_dir="${base_dir}/raw/historical/" #"${HPCPERM}/${project_name}/raw/historical/"
years_dir="${base_dir}/raw/historical/years/"
fname_extra_1=".FGRN055.BN_RACMO2.3p2_FGRN055.3H.nc"
fname_extra_2=".FGRN055.BN_RACMO2.3p2_ERA5_3h_FGRN055.3H.nc"


## will be one of: precip snowfall evap tskin sndiv snowmelt    ##
## ------------------------------------------------------------ ##
varlist=$1

## sets number of data points in each year (#d/yr * #pts/d)     ##
## TK: update to 2023?                                          ##
## in korn, -A assigns parametr to Y01list                      ##
## https://www.mkssoftware.com/docs/man1/set.1.asp              ##
## ------------------------------------------------------------ ##
#set -A yearlist 1957 1961 1971 1981 1991 2001 2011
#set -A nyears   4    10   10   10   10   10   10
#set -A yearlist 1939 1941 1951 1961 1971 1981 1991 2001 2011 2021
#set -A nyears   2    10   10   10   10   10   10   10   10   3
set -A yearlist 2021
set -A nyears   2

#set -A Y57list 736 2920 2920 2928
#set -A Y39list 2920 2928
#set -A Y41list 2920 2920 2920 2928 2920 2920 2920 2928 2920 2920
#set -A Y51list 2920 2928 2920 2920 2920 2928 2920 2920 2920 2928
#set -A Y61list 2920 2920 2920 2928 2920 2920 2920 2928 2920 2920
#set -A Y71list 2920 2928 2920 2920 2920 2928 2920 2920 2920 2928
#set -A Y81list 2920 2920 2920 2928 2920 2920 2920 2928 2920 2920
#set -A Y91list 2920 2928 2920 2920 2920 2928 2920 2920 2920 2928 
#set -A Y01list 2920 2920 2920 2928 2920 2920 2920 2928 2920 2920
#set -A Y11list 2920 2928 2920 2920 2920 2928 2920 2920 2920 2928
set -A Y21list 2920 2920 

set -A data_in_years_list Y21list

## for each variable and each year, create year file            ##
## TK: update to 2023?
## ------------------------------------------------------------ ##

(( years_i = ${#yearlist[@]} - 1 ))

for varname in ${varlist}; do

  for ii in {0..$years_i}; do
  
     (( fileyear = yearlist[ii] ))
    
    if [[ ${fileyear} -lt 1990 ]]; then
      infile="${RACMO_3H_dir}${varname}.KNMI-${fileyear}${fname_extra_1}"
    elif [[ ${fileyear} -gt 1990 ]]; then
      infile="${RACMO_3H_dir}${varname}.KNMI-${fileyear}${fname_extra_2}"
    fi
    
     start_t=0
     end_t=-1

     for jj in {1..${nyears[$ii]}}; do
       (( year = fileyear + jj - 1 ))
      
       outfile=${years_dir}${varname}"_${project_name}_forFDM_Year"${year}".nc"
      
       list_of_years=\$${data_in_years_list[$ii]}

       num_data_in_years=${list_of_years[jj-1]}

       echo $list_of_years
       echo $num_data_in_years
       echo ${Y21list[jj-1]}
#       (( start_t = end_t + 1 ))
#       (( end_t = start_t - 1 + num_data_in_years)) 

#       #echo ${fileyear}" "${start_t}" "${end_t}" ..... "${year}
#       if [[ -f ${outfile} ]]; then
#         echo "File is already present"
#       else
#       echo "ncks for "${varname}", year "${year}", yearlist" ${list_of_years},"with "$num_data_in_years
#       #ncks -d time,${start_t},${end_t} ${infile} ${outfile}
#       fi  
       
     done

   done
done





### ARCHIVE ####
# for varname in ${varlist}; do

#     ## TK: update so that 1) new era runs are used 2) number of years arrays counted automatically 3) individual yearlists don't have to be specified
  
#   for ii in {0..6}; do
#     (( fileyear = yearlist[ii] ))
    
#     if [[ ${fileyear} -lt 1990 ]]; then
#       infile="${RACMO_3H_dir}${varname}.KNMI-${fileyear}${fname_extra_1}"
#     elif [[ ${fileyear} -gt 1990 ]]; then
#       infile="${RACMO_3H_dir}${varname}.KNMI-${fileyear}${fname_extra_2}"
#     fi
    
#     start_t=0
#     end_t=-1
    


#     for jj in {1..${nyears[$ii]}}; do
#       (( year = fileyear + jj - 1 ))
      
#       outfile=${years_dir}${varname}"_${project_name}_forFDM_Year"${year}".nc"
#       (( start_t = end_t + 1 ))
#       if [[ ${ii} -eq 0 ]]; then
#         (( end_t = start_t - 1 + Y57list[jj - 1] )) 
#       elif [[ ${ii} -eq 1 ]]; then
#         (( end_t = start_t - 1 + Y61list[jj - 1] ))
#       elif [[ ${ii} -eq 2 ]]; then
#         (( end_t = start_t - 1 + Y71list[jj - 1] )) 
#       elif [[ ${ii} -eq 3 ]]; then
#         (( end_t = start_t - 1 + Y81list[jj - 1] ))
#       elif [[ ${ii} -eq 4 ]]; then
#         (( end_t = start_t - 1 + Y91list[jj - 1] ))     
#       elif [[ ${ii} -eq 5 ]]; then
#         (( end_t = start_t - 1 + Y01list[jj - 1] ))
#       elif [[ ${ii} -eq 6 ]]; then
#         (( end_t = start_t - 1 + Y11list[jj - 1] ))

#       fi 

#       echo ${fileyear}" "${start_t}" "${end_t}" ..... "${year}
#       if [[ -f ${outfile} ]]; then
#         echo "File is already present"
#       else
#       echo "ncks for "${varname}", year "${year}
#       ncks -d time,${start_t},${end_t} ${infile} ${outfile}
#       fi  
       
#     done

#   done

# done