#!/bin/bash
#SBATCH --job-name=fdm_extend_variable_2ddetail
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

echo "--- T10m ---"
python3 extend_variable_2ddetail.py --output-var T10m --var temp --depth 10 --workers ${W}

# Uncomment additional variables as needed:
# echo "--- SSN ---"
# python3 extend_variable_2ddetail.py --output-var SSN --var dens --depth-begin 0 --depth-end 0.5 --workers ${W}
#
# echo "--- LWC_surf ---"
# python3 extend_variable_2ddetail.py --output-var LWC_surf --var lwc --depth-begin 0 --depth-end 0.2 --workers ${W}
#
# echo "--- rho_2m ---"
# python3 extend_variable_2ddetail.py --output-var rho_2m --var dens --depth-begin 0 --depth-end 2 --workers ${W}

EXIT_CODE=$?
echo "Finished at $(date), exit code: ${EXIT_CODE}"
exit ${EXIT_CODE}
