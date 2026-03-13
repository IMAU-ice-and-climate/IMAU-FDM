#!/bin/bash
#SBATCH --job-name=firn_ice_content
#SBATCH --chdir=/home/nld4814/perm/code/IMAU-FDM/research
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6

# Usage
# -----
# Full run (all ~58 265 columns, 6 parallel workers):
#   sbatch submit_firn_ice_content.sh
#
# Test run (first 500 columns, single worker):
#   TEST=1 sbatch submit_firn_ice_content.sh
#
# The output file lands in config.PROCESSED_DIR
# (scratch/run_FGRN055-era055_1939-2023/post-process/) unless overridden.

set -euo pipefail

module load python3

SCRIPT_DIR="/home/nld4814/perm/code/IMAU-FDM/research"

if [[ "${TEST:-0}" == "1" ]]; then
    echo "=== TEST RUN: first 500 columns ==="
    python3 "$SCRIPT_DIR/firn_ice_content_analysis.py" \
        --output     "FDM_firn_ice_content_FGRN055_TEST.nc" \
        --file-limit 500
else
    echo "=== FULL RUN: all columns, 6 workers ==="
    python3 "$SCRIPT_DIR/firn_ice_content_analysis.py" \
        --output    "FDM_firn_ice_content_FGRN055_1939-2023_2D.nc" \
        --n-workers 6
fi

echo "=== Done ==="
