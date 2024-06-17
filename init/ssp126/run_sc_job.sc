#!/bin/ksh -x

perm_dir="$PERM/pbs/"
script_dir=$PERM"/scripts/input/ssp126/"

jobname=$1
scriptname=$2
vars=$3

echo ${jobname}"   "${scriptname}"    "${vars}
loadscript=${perm_dir}"loadscripts/input/loadscript_"${jobname}

cat <<EOS1 > ${loadscript}
#!/bin/ksh -f

#SBATCH -q nf
#SBATCH -J SCjob_${jobname}
#SBATCH --time=48:00:00
#SBATCH -o ${perm_dir}"/logfiles/input/log_SCjob_"${jobname}".log"
#SBATCH --mem-per-cpu=4000mb

module load nco

cd ${script_dir}

./${scriptname} "${vars}"

date

EOS1

sbatch ${loadscript}

