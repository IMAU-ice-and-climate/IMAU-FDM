#!/bin/ksh 

#SBATCH -J NCLjob_NCL_1d_zs
#SBATCH -q nf
#SBATCH -e /perm/nld4814/code/IMAU-FDM/post-process//logs/log_NCL_1d_zs.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=e.h.case@uu.nl
#SBATCH --export=ALL
#SBATCH -o /perm/nld4814/code/IMAU-FDM/post-process//logs/log_NCL_1d_zs.out
#SBATCH --mem-per-cpu=8000mb
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

ncl 'variab="zs"' dt=0 atzero=0 /perm/nld4814/code/IMAU-FDM/post-process/make_1d_files/FGRN055_make_one_1d_file.ncl

exit 0
