#!/bin/bash
#SBATCH --job-name=fdm_2ddet
#SBATCH --output=logs/slurm-2ddet-%j.out
#SBATCH --error=logs/slurm-2ddet-%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --qos=nf

# =============================================================================
# IMAU-FDM Post-Processing: Create gridded files from 2Ddetail profile output
# =============================================================================
#
# This script processes 2Ddetail files (high-resolution depth profiles) to
# create gridded output like surface snow density (SSN) or near-surface temp.
#
# Usage:
#   sbatch submit_make_2ddetail_files.sh                    # Uses defaults (SSN)
#   sbatch submit_make_2ddetail_files.sh -v dens --depth-begin 0 --depth-end 0.5 --output-var SSN
#   sbatch submit_make_2ddetail_files.sh -v temp -d 10 --output-var T10m
#
# =============================================================================

set -e

# -----------------------------------------------------------------------------
# Default configuration
# -----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="/home/nld4814/scratch/run_FGRN055-era055_1939-2023/output"
WORKERS=${SLURM_CPUS_PER_TASK:-8}

# Default: Surface snow density (upper 50cm)
# DEFAULT_VAR="dens"
# DEFAULT_DEPTH_BEGIN="0"
# DEFAULT_DEPTH_END="0.5"
# DEFAULT_OUTPUT_VAR="SSN"

DEFAULT_VAR="temp"
DEFAULT_DEPTH_BEGIN="9.95"
DEFAULT_DEPTH_END="10.05"
DEFAULT_OUTPUT_VAR="T10m"

# -----------------------------------------------------------------------------
# Job info
# -----------------------------------------------------------------------------
echo "=================================================="
echo "IMAU-FDM 2Ddetail Profile Processing Job"
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
    # No arguments: use defaults (SSN)
    echo "Using defaults: SSN (density 0-0.5m)"
    python3 make_2Ddetail_files.py \
        -o "${OUTPUT_DIR}" \
        -v "${DEFAULT_VAR}" \
        --depth-begin "${DEFAULT_DEPTH_BEGIN}" \
        --depth-end "${DEFAULT_DEPTH_END}" \
        --output-var "${DEFAULT_OUTPUT_VAR}" \
        -n "${WORKERS}"
else
    # Pass all arguments through
    echo "Arguments: $@"
    python3 make_2Ddetail_files.py -o "${OUTPUT_DIR}" -n "${WORKERS}" "$@"
fi

EXIT_CODE=$?

echo ""
echo "=================================================="
echo "Job completed at $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "=================================================="

exit ${EXIT_CODE}
