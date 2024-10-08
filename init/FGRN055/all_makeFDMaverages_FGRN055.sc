#!/bin/ksh -x

## runs makeFDMaverages_FGRN055 for each variable 			 ##
## structure: ./run_sc_job <FILE_PREFIX> <CODE_SOURCE> <VAR> ##
## --------------------------------------------------------- ##

script_filename="makeFDMaverages_${domain}.sc"

## 1:job name 2:script name 3: variable

for var in "evap" "snowfall" "snowmelt" "precip" "sndiv" "tskin" "ff10m"; do
	./run_sc_job.sc ave_${var} ${script_filename} ${var}
done

## ./run_sc_job.sc ave_evap makeFDMaverages_FGRN055.sc evap
## ./run_sc_job.sc ave_snowfall makeFDMaverages_FGRN055.sc snowfall
## ./run_sc_job.sc ave_snowmelt makeFDMaverages_FGRN055.sc snowmelt
## ./run_sc_job.sc ave_precip makeFDMaverages_FGRN055.sc precip
## ./run_sc_job.sc ave_sndiv makeFDMaverages_FGRN055.sc sndiv
## ./run_sc_job.sc ave_tskin makeFDMaverages_FGRN055.sc tskin
## ./run_sc_job.sc ave_ff10m makeFDMaverages_FGRN055.sc ff10m
