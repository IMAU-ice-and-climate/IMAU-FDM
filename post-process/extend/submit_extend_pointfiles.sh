#!/bin/bash
#SBATCH --job-name=fdm_extend_pointfiles
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --qos=nf

SCRIPT_DIR="/home/nld4814/perm/code/IMAU-FDM/post-process/extend"

echo "Job ID:     ${SLURM_JOB_ID}"
echo "Node:       ${SLURM_NODELIST}"
echo "Start time: $(date)"

module load python3 2>/dev/null || true

cd "${SCRIPT_DIR}"
mkdir -p logs

python3 extend_pointfiles.py --workers ${SLURM_CPUS_PER_TASK}

EXIT_CODE=$?
echo "Finished at $(date), exit code: ${EXIT_CODE}"
exit ${EXIT_CODE}
