#!/bin/ksh -x

## runs makeFDMfiles_FGRN055 for each variable				 ##
## structure: ./run_sc_job <FILE_PREFIX> <CODE_SOURCE> <VAR> ##
## --------------------------------------------------------- ##
script_filename="makeFDMfiles"

# 1:job name 2:script name 3: variable 4:project name 5:base directory 6: script directory 7: start year 8: end year 9: number of longitudinal bands 10: cell width

for var in "evap"; do #"snowfall" "snowmelt" "precip" "sndiv" "tskin" "ff10m"; do
	./run_sc_job.sc ts_${var} ${script_filename} ${var} ${project_name} ${base_dir} ${script_dir} ${ts_start_year} ${ts_end_year} ${num_long_bands} ${cell_width}
done

## ./run_sc_job.sc ts_evap makeFDMfiles_FGRN055.sc evap
## ./run_sc_job.sc ts_snowfall makeFDMfiles_FGRN055.sc snowfall
## ./run_sc_job.sc ts_snowmelt makeFDMfiles_FGRN055.sc snowmelt
## ./run_sc_job.sc ts_precip makeFDMfiles_FGRN055.sc precip
## ./run_sc_job.sc ts_sndiv makeFDMfiles_FGRN055.sc sndiv
## ./run_sc_job.sc ts_tskin makeFDMfiles_FGRN11.sc tskin
## ./run_sc_job.sc ts_ff10m makeFDMfiles_FGRN11.sc ff10m
