#!/bin/bash
#SBATCH --job-name=ice_lens
#SBATCH --chdir=/home/nld4814/perm/code/IMAU-FDM/research
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6

# Usage
# -----
# Full run (all ~58 265 columns, 6 parallel workers):
#   sbatch submit_ice_lens.sh
#
# Test run (first 500 columns, single worker):
#   sbatch submit_ice_lens.sh --test
#
# The output file lands in config.PROCESSED_DIR
# (scratch/run_FGRN055-era055_1939-2023/post-process/) unless overridden.
#
# Note: 24 GB requested because the output arrays for 2Ddetail files
# are ~7 GB (3080 time steps × 566 × 438 grid, 3 variables).

set -euo pipefail

module load python3

SCRIPT_DIR="/home/nld4814/perm/code/IMAU-FDM/research"

if [[ "${1:-}" == "--test" ]]; then
    echo "=== TEST RUN: first 500 columns ==="
    python3 "$SCRIPT_DIR/ice_lens_detection.py" \
        --output     "FDM_ice_lens_FGRN055_TEST.nc" \
        --file-limit 500
else
    echo "=== FULL RUN: all columns, 6 workers ==="
    python3 "$SCRIPT_DIR/ice_lens_detection.py" \
        --output    "FDM_ice_lens_FGRN055_1939-2023_2Ddetail.nc" \
        --n-workers 6
fi

echo "=== Done ==="
