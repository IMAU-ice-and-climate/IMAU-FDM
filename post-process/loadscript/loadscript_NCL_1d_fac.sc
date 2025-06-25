#!/bin/ksh 

#SBATCH -J NCL_1d_fac
#SBATCH -q nf
#SBATCH -e /perm/nld3562/code/IMAU-FDM/post-process//logs/log_NCL_1d_fac.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=e.h.case@uu.nl
#SBATCH --export=ALL
#SBATCH -o /perm/nld3562/code/IMAU-FDM/post-process//logs/log_NCL_1d_fac.out
#SBATCH --mem-per-cpu=16000mb
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

ncl 'variab="FirnAir"' dt=0 atzero=0 /perm/nld3562/code/IMAU-FDM/post-process/make_1d_files/FGRN055_nld3561_1d.ncl

exit 0
