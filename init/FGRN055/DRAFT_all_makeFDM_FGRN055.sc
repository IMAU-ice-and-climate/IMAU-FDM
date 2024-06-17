#!/bin/ksh -x

## runs makeFDM* for each variable				 			 ##
## structure: ./run_sc_job <FILE_PREFIX> <CODE_SOURCE> <VAR> ##
## --------------------------------------------------------- ##

project_name="FGRN055-era055"
base_dir="${HPCPERM}/${project_name}/"
script_dir="$PERM/code/init-scripts/FGRN055/"

# defined by base_dir in each makeFDM* file as:
# RACMO_3H_dir="${base_dir}/raw/historical/"
# years_dir="${base_dir}/raw/historical/years/"
# ts_dir="${base_dir}/input/timeseries/"
# ave_dir="${base_dir}/input/averages/"

make_years="makeFDMyears_FGRN055.sc"
make_files="makeFDMfiles_FGRN055.sc"
make_avgs="makeFDMaverages_FGRN055.sc"

ts_start_year=1957
ts_end_year=2020
avg_start_year=1960
avg_end_year=1981
num_long_bands=74
cell_width=5

do_years=false
do_files=false
do_averages=true
varlist=("evap" "snowfall" "snowmelt" "precip" "sndiv" "tskin" "ff10m")
jobids=


## 1:job name 2:script name 3: variable 4:project name 5:base directory 6: script directory 7: start year 8: end year 9: number of longitudinal bands 10: cell width

if [ do_years ]; then

	## for var in "evap" "snowfall" "snowmelt" "precip" "sndiv" "tskin" "ff10m"; do
	for var in ${varlist[@]}; do
		initRacmoFile=${./run_sc_job.sc years_${var} ${make_years} ${var} ${project_name} ${base_dir} ${script_dir}}
		job_ids=$jobids":"$(sbatch --parsable ${initRacmoFile})
	done

	sbatch_command_files="sbatch --dependency=aftercorr${job_ids} "
	sbatch_command_avgs="sbatch --dependency=aftercorr${job_ids} "

else

	sbatch_command_files="sbatch "
	sbatch_command_avgs="sbatch "

fi

for var in ${varlist[@]}; do

	if [ do_files ]; then
		initRacmoFile_files=${./run_sc_job.sc ts_${var} ${make_files} ${var} ${project_name} ${base_dir} ${script_dir} ${ts_start_year} ${ts_end_year} ${num_long_bands} ${cell_width}}
		echo ${sbatch_command_files}${initRacmoFile_files}
	fi

	if [ do_averages ]; then
		initRacmoFile_avgs=${./run_sc_job.sc ave_${var} ${make_avgs} ${var} ${project_name} ${base_dir} ${script_dir} ${avg_start_year} ${avg_end_year}}
		echo ${sbatch_command_avgs}${initRacmoFile_avgs}
	fi
	
done

