#!/bin/bash
#SBATCH --job-name=fdm_batch
#SBATCH --chdir=/home/nld4814/perm/code/IMAU-FDM/research
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6

# Usage
# -----
# Full run:
#   JOB=firn_memory       sbatch jobs/submit_batch.sh
#   JOB=firn_ice_content  sbatch jobs/submit_batch.sh
#   JOB=ice_lens          sbatch jobs/submit_batch.sh
#
# Test run (first 500 columns):
#   TEST=1 JOB=firn_memory  sbatch jobs/submit_batch.sh
#
# Memory: 24 GB (ice_lens needs ~24 GB; firn jobs use ~16 GB)
# Workers: 6 (matches --cpus-per-task)

set -euo pipefail

module load python3

SCRIPT_DIR="/home/nld4814/perm/code/IMAU-FDM/research"

: "${JOB:?ERROR: set JOB before calling sbatch (firn_memory | firn_ice_content | ice_lens)}"

if [[ "${TEST:-0}" == "1" ]]; then
    EXTRA="--file-limit 500"
    echo "=== TEST RUN: first 500 columns ==="
else
    EXTRA=""
    echo "=== FULL RUN: all columns, 6 workers ==="
fi

case "$JOB" in
    firn_memory)
        OUTPUT="${TEST:+FDM_firn_memory_FGRN055_TEST.nc}"
        OUTPUT="${OUTPUT:-FDM_firn_memory_FGRN055_1939-2023_2D.nc}"
        python3 "$SCRIPT_DIR/firn_analysis.py" memory \
            --output "$OUTPUT" --n-workers 6 $EXTRA
        ;;
    firn_ice_content)
        OUTPUT="${TEST:+FDM_firn_ice_content_FGRN055_TEST.nc}"
        OUTPUT="${OUTPUT:-FDM_firn_ice_content_FGRN055_1939-2023_2D.nc}"
        python3 "$SCRIPT_DIR/firn_analysis.py" ice_content \
            --output "$OUTPUT" --n-workers 6 $EXTRA
        ;;
    ice_lens)
        OUTPUT="${TEST:+FDM_ice_lens_FGRN055_TEST.nc}"
        OUTPUT="${OUTPUT:-FDM_ice_lens_FGRN055_1939-2023_2Ddetail.nc}"
        python3 "$SCRIPT_DIR/ice_lens_detection.py" \
            --output "$OUTPUT" --n-workers 6 $EXTRA
        ;;
    *)
        echo "ERROR: unknown JOB='$JOB'. Choose: firn_memory | firn_ice_content | ice_lens"
        exit 1
        ;;
esac

echo "=== Done ==="
