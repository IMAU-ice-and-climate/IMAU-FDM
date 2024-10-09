#!/bin/ksh -x

## runs makeFDMyears_FGRN055 for each variable				 ##
## structure: ./run_sc_job <FILE_PREFIX> <CODE_SOURCE> <VAR> ##
## --------------------------------------------------------- ##

script_filename="makeFDMyears_${domain}.sc"

## 1:job name 2:script name 3: variable 4:project name 5:base directory 6: script directory 7: start year 8: end year 9: number of longitudinal bands 10: cell width

for var in "evap" "snowfall" "snowmelt" "precip" "sndiv" "tskin" "ff10m"; do
	./run_sc_job.sc yrs_${var} ${script_filename} ${var}
done