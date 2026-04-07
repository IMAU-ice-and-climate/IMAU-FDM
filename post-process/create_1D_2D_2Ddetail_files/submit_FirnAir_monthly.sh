#!/bin/bash
#SBATCH --job-name=fdm_FirnAir_monthly
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --qos=nf

SCRIPT_DIR="/home/nld4814/perm/code/IMAU-FDM/post-process/create_1D_2D_2Ddetail_files"
WORKERS=${SLURM_CPUS_PER_TASK:-16}

echo "Job ID:     ${SLURM_JOB_ID}"
echo "Node:       ${SLURM_NODELIST}"
echo "Start time: $(date)"

module load python3 2>/dev/null || true

cd "${SCRIPT_DIR}"
mkdir -p logs

python3 make_1d_files.py \
    --var FirnAir \
    --timestep monthly \
    --spinup-start 1940 \
    --spinup-end 1970 \
    --workers ${WORKERS}

EXIT_CODE=$?
echo "Finished at $(date), exit code: ${EXIT_CODE}"
exit ${EXIT_CODE}
