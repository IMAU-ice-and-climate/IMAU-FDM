#!/bin/ksh -x

## runs makeFDMfiles_FGRN055 for each variable				 ##
## structure: ./run_sc_job <FILE_PREFIX> <CODE_SOURCE> <VAR> ##
## everything else is set in run_sc_job						 ##
## --------------------------------------------------------- ##
script_filename="makeFDMfiles"

# 1:job name 2:script name 3: variable

for var in "evap" "snowfall" "snowmelt" "precip" "sndiv" "tskin" "ff10m"; do
	./run_sc_job.sc ts_${var} ${script_filename} ${var}
done

## ./run_sc_job.sc ts_evap makeFDMfiles_FGRN055.sc evap
## ./run_sc_job.sc ts_snowfall makeFDMfiles_FGRN055.sc snowfall
## ./run_sc_job.sc ts_snowmelt makeFDMfiles_FGRN055.sc snowmelt
## ./run_sc_job.sc ts_precip makeFDMfiles_FGRN055.sc precip
## ./run_sc_job.sc ts_sndiv makeFDMfiles_FGRN055.sc sndiv
## ./run_sc_job.sc ts_tskin makeFDMfiles_FGRN11.sc tskin
## ./run_sc_job.sc ts_ff10m makeFDMfiles_FGRN11.sc ff10m
