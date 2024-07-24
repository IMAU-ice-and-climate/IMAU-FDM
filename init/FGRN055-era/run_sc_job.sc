#!/bin/ksh -x

log_dir="$HPCPERM/FGRN055_era055/logfiles/init/"
job_script_dir="$HPCPERM/FGRN055_era055/raw/jobs/"
script_dir="$PERM/code/IMAU-FDM/init/FGRN055-era/"

jobname=$1
scriptname=$2
var=$3

echo ${jobname}"   "${scriptname}"    "${var}
loadscript="${job_script_dir}init_job_ts_${vars}"

cat <<EOS1 > ${loadscript}
#!/bin/ksh -f

#SBATCH -q nf
#SBATCH -J SCjob_${jobname}
#SBATCH --time=48:00:00
#SBATCH -o ${log_dir}log_SCjob_${jobname}.log
#SBATCH --mem=15G

module load nco

cd ${script_dir}

./${scriptname} "${var}"

date

EOS1

sbatch ${loadscript}

