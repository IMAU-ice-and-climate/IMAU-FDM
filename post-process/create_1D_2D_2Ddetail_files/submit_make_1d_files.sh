#!/bin/bash
#SBATCH --job-name=fdm_1d_files
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --qos=nf

# =============================================================================
# IMAU-FDM Post-Processing: Convert 1D output files to gridded NetCDF files
# =============================================================================
#
# This script processes all point files and creates gridded files for
# each variable. Output files are written to the post-process directory.
#
# Usage:
#   sbatch submit_make_1d_files.sh                    # Process all variables
#   sbatch submit_make_1d_files.sh h_surf FirnAir    # Process specific variables
#
# Output location: /home/nld4814/scratch/run_FGRN055-era055_1939-2023/post-process/
# =============================================================================

# Exit on error
set -e

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKERS=${SLURM_CPUS_PER_TASK:-16}
TIMESTEP="10day"
SPINUP_START=1940
SPINUP_END=1970

# -----------------------------------------------------------------------------
# Load required modules (adjust for your ECMWF environment)
# -----------------------------------------------------------------------------
echo "=================================================="
echo "IMAU-FDM Post-Processing Job"
echo "=================================================="
echo "Job ID:       ${SLURM_JOB_ID}"
echo "Node:         ${SLURM_NODELIST}"
echo "CPUs:         ${SLURM_CPUS_PER_TASK}"
echo "Memory:       ${SLURM_MEM_PER_NODE}MB"
echo "Start time:   $(date)"
echo "=================================================="

# Load Python environment
# Option 1: If using ECMWF module system
if command -v module &> /dev/null; then
    module purge 2>/dev/null || true
    module load python3 2>/dev/null || module load python 2>/dev/null || true
fi

# Option 2: If using conda/mamba
if [ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
    source "${HOME}/miniconda3/etc/profile.d/conda.sh"
    conda activate base 2>/dev/null || true
fi

# Verify Python is available
if ! command -v python3 &> /dev/null; then
    echo "ERROR: python3 not found. Please load appropriate modules."
    exit 1
fi

echo "Python:       $(which python3)"
echo "Version:      $(python3 --version)"
echo "=================================================="

# -----------------------------------------------------------------------------
# Determine which variables to process
# -----------------------------------------------------------------------------
cd "${SCRIPT_DIR}"

if [ $# -eq 0 ]; then
    # No arguments: process all variables
    VARS="all"
    echo "Processing: ALL variables"
else
    # Process specified variables
    VARS="$@"
    echo "Processing: ${VARS}"
fi

echo "Timestep:     ${TIMESTEP}"
echo "Spinup:       ${SPINUP_START}-${SPINUP_END}"
echo "Workers:      ${WORKERS}"
echo "=================================================="

# -----------------------------------------------------------------------------
# Run the processing script
# -----------------------------------------------------------------------------
echo ""
echo "Starting processing at $(date)"
echo ""

python3 make_1d_files.py \
    --var ${VARS} \
    --timestep ${TIMESTEP} \
    --spinup-start ${SPINUP_START} \
    --spinup-end ${SPINUP_END} \
    --workers ${WORKERS}

EXIT_CODE=$?

# -----------------------------------------------------------------------------
# Report completion
# -----------------------------------------------------------------------------
echo ""
echo "=================================================="
echo "Job completed at $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "=================================================="

# List output files
OUTPUT_DIR="/home/nld4814/scratch/run_FGRN055-era055_1939-2023/post-process"
if [ -d "${OUTPUT_DIR}" ]; then
    echo ""
    echo "Output files in ${OUTPUT_DIR}:"
    ls -lh "${OUTPUT_DIR}"/*.nc 2>/dev/null || echo "  (no NetCDF files found)"
fi

exit ${EXIT_CODE}
