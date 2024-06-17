#!/bin/ksh -x

## runs makeFDMyears_FGRN055 for each variable				 ##
## structure: ./run_sc_job <FILE_PREFIX> <CODE_SOURCE> <VAR> ##
## --------------------------------------------------------- ##

project_name="FGRN055-era055"
base_dir="${HPCPERM}/${project_name}/"
script_dir="$PERM/code/init-scripts/FGRN055/"

make_years="makeFDMyears_FGRN055.sc"
make_files="makeFDMfiles_FGRN055.sc"
make_avgs="makeFDMaverages_FGRN055.sc"

ts_start_year=1957
ts_end_year=2020
avg_start_year=1960
avg_end_year=1981
num_long_bands=74
cell_width=5

## 1:job name 2:script name 3: variable 4:project name 5:base directory 6: script directory 7: start year 8: end year 9: number of longitudinal bands 10: cell width

for var in "evap" "snowfall" "snowmelt" "precip" "sndiv" "tskin" "ff10m"; do
	./run_sc_job.sc ave_${var} makeFDMyears_FGRN055.sc ${var} ${project_name} ${base_dir} ${script_dir}
done