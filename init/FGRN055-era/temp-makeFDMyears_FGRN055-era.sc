#!/bin/ksh -x

raw_dir="$HPCPERM/FGRN055_era055/raw/historical/"
years_dir="${raw_dir}years/"
fname_extra=".FGRN055.BN_RACMO2.3p2_ERA5_3h_1940_FGRN055.3H.nc"
project_name="FGRN055-era055"
#############################################

#varlist="snowmelt evap precip tskin sndiv snowfall"
varlist="ff10m" # snowmelt evap precip tskin sndiv snowfall 

set -A yearlist 2021
set -A nyears 2

set -A Y21list 2920 2920


for varname in ${varlist}; do
  
  for ii in {0..0}; do			# loop over yearlist
    (( fileyear = yearlist[ii] ))

    infile=${raw_dir}${varname}".KNMI-"${fileyear}${fname_extra}   

    start_t=0
    end_t=-1

    for jj in {1..${nyears[$ii]}}; do	# loop over years in yearlist

      (( year = fileyear + jj - 1 ))

      outfile=${years_dir}${varname}"_${project_name}_forFDM_Year"${year}".nc"

      (( start_t = end_t + 1 ))		# start at the index of the end of previous year + 1
      (( end_t = start_t - 1 + Y21list[jj - 1] ))
   
      echo "year: $year"
      echo "start_t: $start_t"
      echo "end_t: $end_t"

      if [[ -f ${outfile} ]]; then
        echo "File is already present"
      else
      	echo "ncks for "${varname}", year "${year}
	    ncks -d time,${start_t},${end_t} ${infile} ${outfile}
      fi	
       
    done

  done
done

