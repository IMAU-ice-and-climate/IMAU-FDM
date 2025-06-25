#!/bin/ksh -x

#to submit, e.g. ./submit_1d_job.sc NCL_1d_h_surf h_surf 0 0
#can also be submitted through run_submit_1d_jobs to submit multiple


dir="/perm/nld3562/code/IMAU-FDM/post-process/"

jobname=$1
variab=$2
detrend=$3
startzero=$4

loadscript="${dir}loadscript/loadscript_${jobname}.sc"
homedir=`pwd`

cat <<EOFp > ${loadscript}
#!/bin/ksh 

#SBATCH -J ${jobname}
#SBATCH -q nf
#SBATCH -e ${dir}/logs/log_${jobname}.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=e.h.case@uu.nl
#SBATCH --export=ALL
#SBATCH -o ${dir}/logs/log_${jobname}.out
#SBATCH --mem-per-cpu=16000mb
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

ncl 'variab="${variab}"' dt=${detrend} atzero=${startzero} ${homedir}/FGRN055_nld3561_1d.ncl

exit 0
EOFp

chmod u+x ${loadscript}
sbatch ${loadscript}
