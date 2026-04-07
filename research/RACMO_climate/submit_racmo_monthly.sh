#!/bin/bash
#SBATCH --job-name=racmo_monthly
#SBATCH --output=logs/racmo_monthly_%a_%j.out
#SBATCH --error=logs/racmo_monthly_%a_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-6
#SBATCH --qos=nf

# Usage
# -----
#   sbatch submit_racmo_monthly.sh
#
# Runs 7 parallel jobs (one per variable), each reading all longitude-band
# timeseries files, resampling to monthly means, and writing:
#   /scratch/FGRN055_era055/input/monthly/{var}_FGRN055_era055_1939-2023_monthly.nc
#
# Output shape: (time=1014, rlat=566, rlon=438)

set -euo pipefail

VARIABLES=(ff10m precip evap sndiv snowfall snowmelt tskin)
VAR="${VARIABLES[$SLURM_ARRAY_TASK_ID]}"

echo "=================================================="
echo "RACMO monthly means — variable: ${VAR}"
echo "Array task : ${SLURM_ARRAY_TASK_ID}"
echo "Job ID     : ${SLURM_JOB_ID}"
echo "Node       : ${SLURM_NODELIST}"
echo "Start      : $(date)"
echo "=================================================="

module load python3

mkdir -p logs

OUTPUT_DIR="/home/nld4814/scratch/FGRN055_era055/monthlies"
SCRIPT_DIR="/home/nld4814/perm/code/IMAU-FDM/research/RACMO_climate"

python3 "${SCRIPT_DIR}/make_racmo_monthly.py" --var "${VAR}"

echo ""
echo "=================================================="
echo "Finished: $(date)"
echo "=================================================="
