#!/bin/ksh -x

## submits job to turn decadal variable outputs into sliced timeseries ##
## 1:job name 2:script name 3: variable 4:project name 5:base directory 6: script directory 7: start year 8: end year 9: number of longitudinal bands 10: cell width
## ------------------------------------------------------------------- ##

jobname=$1
scriptname=$2
vars=$3

project_name=$4 # "FGRN055-era055"
base_dir=$5 # "${HPCPERM}/${project_name}/"
script_dir=$6 # "$PERM/code/init-scripts/FGRN055/"

jobfile_dir="${base_dir}raw/jobs/"
logfile_dir="${base_dir}logfiles/init/"

start_year=$7
end_year=$8

num_long_bands=$9
cell_width=${10}

## prints job command to initRacmoFile, then submits with sbatch	   ##
## ------------------------------------------------------------------- ##

initRacmoFile="${jobfile_dir}init_job_${jobname}"

cat <<EOS1 > ${initRacmoFile} 
#!/bin/ksh -f

#SBATCH -q nf
#SBATCH -J SCjob_${jobname}
#SBATCH --time=48:00:00
#SBATCH -o ${logfile_dir}"log_SCjob_"${jobname}".log"
#SBATCH --mem-per-cpu=4000mb

module load nco

cd ${script_dir}

./${scriptname} ${vars} ${project_name} ${base_dir} ${start_year} ${end_year} ${num_long_bands} ${cell_width}


EOS1

sbatch ${initRacmoFile}
