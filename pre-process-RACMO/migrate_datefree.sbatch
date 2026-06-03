#!/bin/bash
# Submit the date-free forcing migration as a SLURM job with parallel dtobs stamping.
#
#   sbatch migrate_datefree.sbatch /path/to/<DOMAIN>_<FORCING>/input
#
# Renames are instant; stamping dtobs on each multi-GB classic-netCDF timeseries
# file is a full rewrite, so it is parallelised across the allocated CPUs.

#SBATCH -q nf
#SBATCH -J fdm_datefree_migrate
#SBATCH --cpus-per-task=32
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=1000mb
#SBATCH -o /home/nld4814/scratch/datefree_migrate_%j.out

set -euo pipefail
module load nco

INPUT_DIR="${1:?usage: sbatch migrate_datefree.sbatch <DOMAIN_FORCING>/input}"
SCRIPT="/home/nld4814/perm/code/IMAU-FDM/pre-process-RACMO/migrate_existing_to_datefree.py"

echo "Input dir : ${INPUT_DIR}"
echo "Workers   : ${SLURM_CPUS_PER_TASK:-32}"
echo "Started   : $(date)"

python3 "${SCRIPT}" "${INPUT_DIR}" --apply --workers "${SLURM_CPUS_PER_TASK:-32}"

echo "Finished  : $(date)"
