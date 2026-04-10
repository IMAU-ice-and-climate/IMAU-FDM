#!/bin/ksh -x

## runs makeFDMyears_FGRN055 for each variable				 ##
## structure: ./run_sc_job <FILE_PREFIX> <CODE_SOURCE> <VAR> ##
## --------------------------------------------------------- ##

script_filename="makeFDMyears"

## 1:job name 2:script name 3: variable

for var in "snowfall" "snowmelt" "sndiv" "precip" "tskin" "ff10m" "evap"; do
	./run_sc_job.sc yrs_${var} ${script_filename} ${var}
done