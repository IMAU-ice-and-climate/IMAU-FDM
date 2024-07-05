#!/bin/ksh -x

raw_dir="$SCRATCH/data/input/era055_files/rawfiles/"
ave_dir="$SCRATCH/data/input/era055_files/averages/"
files_dir="$SCRATCH/data/input/era055_files/files/"
#ectape_dir="ec:/nmg/FM_Data/INPUT_PEN055/"

ecp="/usr/local/apps/ecfs/2.0.13rc2/bin/ecp"
nco_dir="/usr/local/apps/nco/4.3.7/bin/"
ncks=${nco_dir}"ncks"

fname_extra=".FGRN055.BN_RACMO2.3p2_FGRN055.3H.nc"

#############################################

varlist=$1 #"snowmelt evap precip tskin sndiv snowfall"
#varlist="snowmelt evap precip tskin sndiv snowfall ff10m"

set -A yearlist 1957 1961 1971 1981 1991 2001 2011
set -A nyears 4 10 10 10 10 10 10

set -A Y57list 736 2920 2920 2928
set -A Y61list 2920 2920 2920 2928 2920 2920 2920 2928 2920 2920
set -A Y71list 2920 2928 2920 2920 2920 2928 2920 2920 2920 2928
set -A Y81list 2920 2920 2920 2928 2920 2920 2920 2928 2920 2920
set -A Y91list 2920 2928 2920 2920 2920 2928 2920 2920 2920 2928 
set -A Y01list 2920 2920 2920 2928 2920 2920 2920 2928 2920 2920
set -A Y11list 2920 2928 2920 2920 2920 2928 2920 2920 2920 2928

for varname in ${varlist}; do
  
  for ii in {0..6}; do
    (( fileyear = yearlist[ii] ))
#    tapefile=${ectape_dir}${varname}".KNMI-"${fileyear}${fname_extra}
    infile=${raw_dir}${varname}".KNMI-"${fileyear}${fname_extra}   
#    ${ecp} ${tapefile} ${infile}
    
    start_t=0
    end_t=-1
    for jj in {1..${nyears[$ii]}}; do
#      echo $fileyear 
      (( year = fileyear + jj - 1 ))
      
      outfile=${raw_dir}${varname}"_FGRN055_forFDM_Year"${year}".nc"
      (( start_t = end_t + 1 ))
      if [[ ${ii} -eq 0 ]]; then
        (( end_t = start_t - 1 + Y57list[jj - 1] )) 
      elif [[ ${ii} -eq 1 ]]; then
        (( end_t = start_t - 1 + Y61list[jj - 1] ))
      elif [[ ${ii} -eq 2 ]]; then
        (( end_t = start_t - 1 + Y71list[jj - 1] )) 
      elif [[ ${ii} -eq 3 ]]; then
        (( end_t = start_t - 1 + Y81list[jj - 1] ))
      elif [[ ${ii} -eq 4 ]]; then
        (( end_t = start_t - 1 + Y91list[jj - 1] ))     
      elif [[ ${ii} -eq 5 ]]; then
        (( end_t = start_t - 1 + Y01list[jj - 1] ))
      elif [[ ${ii} -eq 6 ]]; then
        (( end_t = start_t - 1 + Y11list[jj - 1] ))
      fi 

#      echo ${fileyear}" "${start_t}" "${end_t}"         "${year}
      if [[ -f ${outfile} ]]; then
        echo "File is already present"
      else
      	echo "ncks for "${varname}", year "${year}
	${ncks} -d time,${start_t},${end_t} ${infile} ${outfile}
      fi	
       
    done

  done
done

