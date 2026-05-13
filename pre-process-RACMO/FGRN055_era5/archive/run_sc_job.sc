#!/bin/ksh -x

## submits job to turn decadal variable outputs into sliced timeseries	##
##																		##
## 1:job name 2:script name 3: variable				 					##
##																		##
## --------------------------------------------------------------------	##

domain="FGRN055"
forcing="era055"

vars=("evap" "ff10m" "precip" "sndiv" "snowfall" "snowmelt" "tskin")

jobname=$1

if [ "$jobname" = "years" ]; then
	scriptname="makeFDMyears_${domain}.sc"
elif [ "$jobname" = "timeseries" ]; then
	scriptname="makeFDMtimeseries_${domain}.sc"
elif [ "$jobname" = "averages" ]; then
	scriptname="makeFDMaverages_${domain}.sc"
else
	echo "Invalid script name: ${jobname}_${domain}.sc"
	exit 1
fi

data_base_dir="$SCRATCH" #where the data should be stored
scripts_base_dir="$PERM/code" #where IMAU-FDM directory is

project_name="${domain}_${forcing}"
base_dir="${data_base_dir}/${project_name}"
script_dir="${scripts_base_dir}/IMAU-FDM/init/${domain}_${forcing}"

ts_start_year=1939
ts_end_year=2023
ave_start_year=1940
ave_end_year=1970
num_long_bands=73
cell_width=6

years_dir="${base_dir}/process-RACMO/years-${ts_start_year}"
files_dir="${base_dir}/input/timeseries-${ts_start_year}"
ave_dir="${base_dir}/input/averages-${ts_start_year}_${ave_start_year}-${ave_end_year}"
jobfile_dir="${base_dir}/process-RACMO/jobs"
logfile_dir="${base_dir}/process-RACMO/logs"

mkdir -p ${years_dir}
mkdir -p ${files_dir}
mkdir -p ${ave_dir}
mkdir -p ${logfile_dir}
mkdir -p ${jobfile_dir}

## creates and writes initRacmoFile, then submits with sbatch			##
##																		##
## 1: variable 2: project name 3: base directory 						##
## 4: years_dir 5: files_dir 6: ave_dir									##
## 7: number of longitudinal bands 8: cell width 						##
## 9: ts start year 10: ts end year 									##
## 11: avg start year 12: avg end year 									##
##																		##
## -------------------------------------------------------------------	##

for var in ${vars[@]}; do

initRacmoFile="${jobfile_dir}/preprocess-RACMO_job_${jobname}_${var}"

cat <<EOS1 > ${initRacmoFile} 
#!/bin/ksh -f

#SBATCH -q nf
#SBATCH -J ${jobname}_${var}_preprocess
#SBATCH --time=48:00:00
#SBATCH -o ${logfile_dir}/${jobname}_${var}_preprocess.log
#SBATCH --mem-per-cpu=4000mb

module load nco

echo "Start year: " ${ts_start_year}
echo "End year: " ${ts_end_year}
echo "Spinup start year: " ${ave_start_year}
echo "Spinup end year: " ${ave_end_year}
echo "Number of lon bands: " ${num_long_bands}
echo "Cell width: " ${cell_width}

cd ${script_dir}

./${scriptname} ${var} ${project_name} ${base_dir} ${years_dir} ${files_dir} ${ave_dir} ${num_long_bands} ${cell_width} ${ts_start_year} ${ts_end_year} ${ave_start_year} ${ave_end_year}

EOS1

sbatch ${initRacmoFile}

done
