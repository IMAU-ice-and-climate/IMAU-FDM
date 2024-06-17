#!/bin/ksh -x

## submits job to turn decadal variable outputs into sliced timeseries ##
## ------------------------------------------------------------------- ##

project_name="FGRN055-era055"
base_dir="${HPCPERM}/${project_name}/"
jobfile_dir="${base_dir}raw/jobs/"
logfile_dir="${base_dir}logfiles/init/"
script_dir="$PERM/code/init-scripts/FGRN055/"

timestamp=$(date +%H%M)

jobname=$1
scriptname=$2
vars=$3

echo ${jobname}"   "${scriptname}"    "${vars}

## prints job command to initRacmoFile, then submits with sbatch	   ##
## ------------------------------------------------------------------- ##

initRacmoFile="${jobfile_dir}init_job_${jobname}"

cat <<EOS1 > ${initRacmoFile} 
#!/bin/ksh -f

#SBATCH -q nf
#SBATCH -J SCjob_${jobname}
#SBATCH --time=48:00:00
#SBATCH -o ${logfile_dir}${timestamp}"log_SCjob_"${jobname}".log"
#SBATCH --mem-per-cpu=4000mb

module load nco

cd ${script_dir}

./${scriptname} "${vars}"

date

EOS1

sbatch ${initRacmoFile}
