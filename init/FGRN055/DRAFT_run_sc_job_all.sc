#!/bin/ksh -x

## submits job to turn decadal variable outputs into sliced timeseries ##
## ------------------------------------------------------------------- ##

jobname=$1
scriptname=$2
vars=$3

project_name=$4 # "FGRN055-era055"
base_dir=$5 # "${HPCPERM}/${project_name}/"
script_dir=$6 # "$PERM/code/init-scripts/FGRN055/"

jobfile_dir="${base_dir}raw/jobs/"
logfile_dir="${base_dir}logfiles/init/"

RACMO_3H_dir="${base_dir}/raw/historical/"
years_dir="${base_dir}/raw/historical/years/"
ts_dir="${base_dir}/input/timeseries/"
ave_dir="${base_dir}/input/averages/"

start_year=$8
end_year=$9
num_long_bands=${10}
cell_width=${11}


#echo ${jobname}"   "${scriptname}"    "${vars}

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
## 1-var 2-project_name 3-base_dir 4-start-year 5-end-year 6-longitude-bands 7-cell-width

date

EOS1

echo ${initRacmoFile}
#sbatch ${initRacmoFile}
#job_id=$(sbatch --parsable ${initRacmoFile})
#echo ${job_id}