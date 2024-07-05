#!/bin/ksh -x

raw_dir="$SCRATCH/data/input/ssp126/rawfiles/"

fname_extra=".FGRN11.BN_RACMO2.3p2_CESM2_FGRN11.3H.nc"

#############################################

varlist=$1 #"snowmelt evap precip tskin sndiv snowfall"
#varlist="snowmelt evap precip tskin sndiv snowfall ff10m"

set -A yearlist 1950 1951 1961 1971 1981 1991 2001 2011 2015 2021 2031 2041 2051 2061 2071 2081 2091
set -A nyears   1    10   10   10   10   10   10   4    6    10   10   10   10   10   10   10   9

set -A Y1950list 2920
set -A Y1951list 2920 2920 2920 2920 2920 2920 2920 2920 2920 2920
set -A Y1961list 2920 2920 2920 2920 2920 2920 2920 2920 2920 2920
set -A Y1971list 2920 2920 2920 2920 2920 2920 2920 2920 2920 2920
set -A Y1981list 2920 2920 2920 2920 2920 2920 2920 2920 2920 2920
set -A Y1991list 2920 2920 2920 2920 2920 2920 2920 2920 2920 2920
set -A Y2001list 2920 2920 2920 2920 2920 2920 2920 2920 2920 2920
set -A Y2011list 2920 2920 2920 2920
set -A Y2015list 2920 2920 2920 2920 2920 2920
set -A Y2021list 2920 2920 2920 2920 2920 2920 2920 2920 2920 2920
set -A Y2031list 2920 2920 2920 2920 2920 2920 2920 2920 2920 2920
set -A Y2041list 2920 2920 2920 2920 2920 2920 2920 2920 2920 2920
set -A Y2051list 2920 2920 2920 2920 2920 2920 2920 2920 2920 2920
set -A Y2061list 2920 2920 2920 2920 2920 2920 2920 2920 2920 2920
set -A Y2071list 2920 2920 2920 2920 2920 2920 2920 2920 2920 2920
set -A Y2081list 2920 2920 2920 2920 2920 2920 2920 2920 2920 2920
set -A Y2091list 2920 2920 2920 2920 2920 2920 2920 2920 2920


for varname in ${varlist}; do
  
  for ii in {0..16}; do			# loop over yearlist
    (( fileyear = yearlist[ii] ))

    infile=${raw_dir}${varname}".KNMI-"${fileyear}${fname_extra}   

    start_t=0
    end_t=-1
    if [[ ${ii} -eq 8 ]]; then
        start_1=2920*4
        end_t=2920*4-1
    fi

    for jj in {1..${nyears[$ii]}}; do	# loop over years in yearlist

      (( year = fileyear + jj - 1 ))

      outfile=${raw_dir}${varname}"_FGRN11_forFDM_Year"${year}".nc"

      (( start_t = end_t + 1 ))		# start at the index of the end of previous year + 1
      if [[ ${ii} -eq 0 ]]; then
        (( end_t = start_t - 1 + Y1950list[jj - 1] )) 	# end at the last index of that year
      elif [[ ${ii} -eq 1 ]]; then
        (( end_t = start_t - 1 + Y1951list[jj - 1] ))
      elif [[ ${ii} -eq 2 ]]; then
        (( end_t = start_t - 1 + Y1961list[jj - 1] )) 
      elif [[ ${ii} -eq 3 ]]; then
        (( end_t = start_t - 1 + Y1971list[jj - 1] ))
      elif [[ ${ii} -eq 4 ]]; then
        (( end_t = start_t - 1 + Y1981list[jj - 1] ))     
      elif [[ ${ii} -eq 5 ]]; then
        (( end_t = start_t - 1 + Y1991list[jj - 1] ))
      elif [[ ${ii} -eq 6 ]]; then
        (( end_t = start_t - 1 + Y2001list[jj - 1] ))
      elif [[ ${ii} -eq 7 ]]; then
        (( end_t = start_t - 1 + Y2011list[jj - 1] ))
      elif [[ ${ii} -eq 8 ]]; then
        (( end_t = start_t - 1 + Y2015list[jj - 1] ))
      elif [[ ${ii} -eq 9 ]]; then
        (( end_t = start_t - 1 + Y2021list[jj - 1] ))
      elif [[ ${ii} -eq 10 ]]; then
        (( end_t = start_t - 1 + Y2031list[jj - 1] ))
      elif [[ ${ii} -eq 11 ]]; then
        (( end_t = start_t - 1 + Y2041list[jj - 1] ))
      elif [[ ${ii} -eq 12 ]]; then
        (( end_t = start_t - 1 + Y2051list[jj - 1] ))
      elif [[ ${ii} -eq 13 ]]; then
        (( end_t = start_t - 1 + Y2061list[jj - 1] ))
      elif [[ ${ii} -eq 14 ]]; then
        (( end_t = start_t - 1 + Y2071list[jj - 1] ))
      elif [[ ${ii} -eq 15 ]]; then
        (( end_t = start_t - 1 + Y2081list[jj - 1] ))
      elif [[ ${ii} -eq 16 ]]; then
        (( end_t = start_t - 1 + Y2091list[jj - 1] ))
     fi 

     echo $year
     echo $start_t
     echo $end_t

      if [[ -f ${outfile} ]]; then
        echo "File is already present"
      else
      	echo "ncks for "${varname}", year "${year}
	    ncks -d time,${start_t},${end_t} ${infile} ${outfile}
      fi	
       
    done

  done
done

