#!/bin/ksh -x

perm_dir="$PERM/pbs/"

jobname=$1
scriptname=$2
vars=$3

echo ${jobname}"   "${scriptname}"    "${vars}
loadscript=${perm_dir}"loadscripts/loadscript_"${jobname}

cat <<EOS1 > ${loadscript}
#!/bin/ksh -f

# #PBS -S /usr/bin/ksh
# #PBS -N SCjob_${jobname}
# #PBS -q ns
# #PBS -j oe
# #PBS -V
# #PBS -l EC_billing_account=spnlberg
# #PBS -o ${perm_dir}"logfiles/log_SCjob_"${jobname}".out"
# #PBS -l EC_memory_per_task=2000mb
# #PBS -l EC_threads_per_task=1
# #PBS -m ae
# #PBS -M brils.hpc@gmail.com
# #PBS -l walltime=48:00:00

#SBATCH -q nf
#SBATCH -J SCjob_${jobname}
#SBATCH --time=48:00:00
#SBATCH -o ${perm_dir}"/logfiles/input/log_SCjob_"${jobname}".log"
#SBATCH --mem-per-cpu=2000mb

module load nco

cd $PERM/scripts/input/FGRN11_ssp585/

./${scriptname} "${vars}"

date

EOS1

sbatch ${loadscript}

