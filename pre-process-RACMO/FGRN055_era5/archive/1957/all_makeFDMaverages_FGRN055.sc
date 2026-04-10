#!/bin/ksh -x

## runs makeFDMaverages_FGRN055 for each variable 			 ##
## structure: ./run_sc_job <FILE_PREFIX> <CODE_SOURCE> <VAR> ##
## --------------------------------------------------------- ##

script_filename="makeFDMaverages"

## 1:job name 2:script name 3: variable

for var in "evap" "snowfall" "snowmelt" "precip" "sndiv" "tskin" "ff10m"; do
	./run_sc_job.sc ave_${var} ${script_filename} ${var}
done
