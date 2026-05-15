#!/bin/bash
#SBATCH --job-name=fdm_extend_variable_2d
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --qos=nf

SCRIPT_DIR="/home/nld4814/perm/code/IMAU-FDM/post-process/extend"

echo "Job ID:     ${SLURM_JOB_ID}"
echo "Node:       ${SLURM_NODELIST}"
echo "Start time: $(date)"

module load python3 2>/dev/null || true

cd "${SCRIPT_DIR}"
mkdir -p logs

W=${SLURM_CPUS_PER_TASK}

echo "--- z830 ---"
python3 extend_variable_2d.py --output-var z830 --var dens --threshold 830 --workers ${W}

echo "--- z550 ---"
python3 extend_variable_2d.py --output-var z550 --var dens --threshold 550 --workers ${W}

EXIT_CODE=$?
echo "Finished at $(date), exit code: ${EXIT_CODE}"
exit ${EXIT_CODE}
