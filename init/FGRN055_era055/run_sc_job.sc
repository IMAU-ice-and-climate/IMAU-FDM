#!/bin/ksh -x

## submits job to turn decadal variable outputs into sliced timeseries	##
##																		##
## 1:job name 2:script name 3: variable				 					##
##																		##
## --------------------------------------------------------------------	##

jobname=$1
scriptname=$2
vars=$3

domain="FGRN055"
forcing="era055"
project_name="${domain}_${forcing}"
base_dir="${SCRATCH}/${project_name}"
script_dir="${PERM}/IMAU-FDM/init/${domain}"

ts_start_year=1957
ts_end_year=2023
avg_start_year=1960
avg_end_year=1981
num_long_bands=74
cell_width=5

years_dir="${base_dir}/process-RACMO/years/"
ave_dir="${base_dir}/input/averages/"
jobfile_dir="${base_dir}/process-RACMO/jobs/"
logfile_dir="${base_dir}/logfiles/process-RACMO/"

mkdir -p years_dir
mkdir -p avg_dir
mkdir -p logfile_dir
mkdir -p jobfile_dir

## creates and writes initRacmoFile, then submits with sbatch			##
##																		##
## 1: variable 2: project name 3: base directory 						##
## 4: years_dir 5: ave_dir												##
## 5: number of longitudinal bands 6: cell width 						##
## 8: ts start year 9: ts end year 										##
## 10: avg start year 11: avg end year 									##
##																		##
## -------------------------------------------------------------------	##

initRacmoFile="${jobfile_dir}init_job_${jobname}"

cat <<EOS1 > ${initRacmoFile} 
#!/bin/ksh -f

#SBATCH -q nf
#SBATCH -J ${jobname}_preprocess
#SBATCH --time=48:00:00
#SBATCH -o ${logfile_dir}${jobname}_preprocess.log
#SBATCH --mem-per-cpu=4000mb

module load nco

cd ${script_dir}

./${scriptname} ${vars} ${project_name} ${base_dir} ${years_dir} ${ave_dir} ${num_long_bands} ${cell_width} ${ts_start_year} ${ts_end_year} ${avg_start_year} ${avg_end_year}


EOS1

sbatch ${initRacmoFile}
