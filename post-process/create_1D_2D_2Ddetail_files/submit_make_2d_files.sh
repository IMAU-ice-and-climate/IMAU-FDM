#!/bin/bash
#SBATCH --job-name=fdm_2d
#SBATCH --output=logs/slurm-2d-%j.out
#SBATCH --error=logs/slurm-2d-%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --qos=nf

# =============================================================================
# IMAU-FDM Post-Processing: Create gridded files from 2D profile output
# =============================================================================
#
# This script processes 2D profile files (time × depth) to create gridded files
# for variables like critical density depths (z550, z830, z917).
#
# Usage:
#   sbatch submit_make_2d_files.sh                     # Uses defaults below
#   sbatch submit_make_2d_files.sh --var dens --threshold 830 --output-var z830
#
# =============================================================================

set -e

# -----------------------------------------------------------------------------
# Default configuration - modify these or override via command line
# -----------------------------------------------------------------------------
SCRIPT_DIR="/home/nld4814/perm/code/IMAU-FDM/post-process/python"
OUTPUT_DIR="/home/nld4814/scratch/run_FGRN055-era055_1939-2023/output"
WORKERS=${SLURM_CPUS_PER_TASK:-8}

# Default: find depth where density = 830 kg/m³
DEFAULT_VAR="dens"
DEFAULT_THRESHOLD="550" # "830"
DEFAULT_OUTPUT_VAR="z550" #"z830"

# -----------------------------------------------------------------------------
# Job info
# -----------------------------------------------------------------------------
echo "=================================================="
echo "IMAU-FDM 2D Profile Processing Job"
echo "=================================================="
echo "Job ID:       ${SLURM_JOB_ID}"
echo "Node:         ${SLURM_NODELIST}"
echo "CPUs:         ${SLURM_CPUS_PER_TASK}"
echo "Start time:   $(date)"
echo "=================================================="

# Load Python environment
if command -v module &> /dev/null; then
    module purge 2>/dev/null || true
    module load python3 2>/dev/null || module load python 2>/dev/null || true
fi

if [ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
    source "${HOME}/miniconda3/etc/profile.d/conda.sh"
    conda activate base 2>/dev/null || true
fi

if ! command -v python3 &> /dev/null; then
    echo "ERROR: python3 not found"
    exit 1
fi

echo "Python:       $(which python3)"
echo "=================================================="

# -----------------------------------------------------------------------------
# Run processing
# -----------------------------------------------------------------------------
cd "${SCRIPT_DIR}"

if [ $# -eq 0 ]; then
    # No arguments: use defaults
    echo "Using defaults: --var ${DEFAULT_VAR} --threshold ${DEFAULT_THRESHOLD} --output-var ${DEFAULT_OUTPUT_VAR}"
    python3 make_2d_files.py \
        -o "${OUTPUT_DIR}" \
        -v "${DEFAULT_VAR}" \
        -t "${DEFAULT_THRESHOLD}" \
        --output-var "${DEFAULT_OUTPUT_VAR}" \
        -n "${WORKERS}"
else
    # Pass all arguments through
    echo "Arguments: $@"
    python3 make_2d_files.py -o "${OUTPUT_DIR}" -n "${WORKERS}" "$@"
fi

EXIT_CODE=$?

echo ""
echo "=================================================="
echo "Job completed at $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "=================================================="

exit ${EXIT_CODE}
