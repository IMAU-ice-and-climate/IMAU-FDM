#!/bin/ksh 

#SBATCH -J 1d_miss2
#SBATCH -q nf
#SBATCH -e /perm/nld4814/code/IMAU-FDM/post-process//logs/log_1d_miss2.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=e.h.case@uu.nl
#SBATCH --export=ALL
#SBATCH -o /perm/nld4814/code/IMAU-FDM/post-process//logs/log_1d_miss2.out
#SBATCH --mem-per-cpu=500mb
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

ncl /perm/nld4814/code/IMAU-FDM/post-process/make_1d_files/write_out_missing_1D.ncl

exit 0
