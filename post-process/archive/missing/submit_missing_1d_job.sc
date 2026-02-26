#!/bin/ksh -x

#find missing points and write to a file
#to submit, e.g. ./submit_1d_job.sc miss_1D


dir="/perm/nld4814/code/IMAU-FDM/post-process/"

jobname=$1

loadscript="${dir}loadscript/loadscript_${jobname}.sc"
homedir=`pwd`

cat <<EOFp > ${loadscript}
#!/bin/ksh 

#SBATCH -J ${jobname}
#SBATCH -q nf
#SBATCH -e ${dir}logs/log_${jobname}.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=e.h.case@uu.nl
#SBATCH --export=ALL
#SBATCH -o ${dir}logs/log_${jobname}.out
#SBATCH --mem-per-cpu=500mb
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00

ncl ${homedir}/write_out_missing_1D.ncl

exit 0
EOFp

chmod u+x ${loadscript}
sbatch ${loadscript}
