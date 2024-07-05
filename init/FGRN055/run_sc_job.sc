#!/bin/ksh -x

perm_dir="/perm/ms/nl/rumb/pbs/"

jobname=$1
scriptname=$2
vars=$3

echo ${jobname}"   "${scriptname}"    "${vars}
loadscript=${perm_dir}"loadscripts/loadscript_"${jobname}

cat <<EOS1 > ${loadscript}
#!/bin/ksh -f
#PBS -S /usr/bin/ksh
#PBS -N SCjob_${jobname}
#PBS -q ns
#PBS -j oe
#PBS -V
#PBS -l EC_billing_account=spnlberg
#PBS -o ${perm_dir}"logfiles/log_SCjob_"${jobname}".out"
#PBS -l EC_memory_per_task=2000mb
#PBS -l EC_threads_per_task=1
#PBS -m ae
#PBS -M brils.hpc@gmail.com
#PBS -l walltime=48:00:00

date
echo 'hostname'

cd /perm/ms/nl/rumb/scripts/input

./${scriptname} "${vars}"

date

EOS1

qsub ${loadscript}
