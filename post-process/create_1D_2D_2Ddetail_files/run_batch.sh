#!/bin/bash
#SBATCH --job-name=fdm_1d_maps
#SBATCH --qos=nf
#SBATCH --time=06:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/fdm_1d_maps_%j.out
#SBATCH --error=logs/fdm_1d_maps_%j.err

# =============================================================================
# IMAU-FDM 1D to Gridded Maps - ECMWF Batch Job
# =============================================================================
#
# Submit this script with: sbatch run_batch.sh
#
# Or with custom options:
#   sbatch run_batch.sh h_surf           # Single variable
#   sbatch run_batch.sh all              # All variables
#   sbatch run_batch.sh "h_surf FirnAir" # Multiple variables
#
# Monitor job: squeue -u $USER
# View output: tail -f logs/fdm_1d_maps_<jobid>.out
#
# =============================================================================

# Set up environment
echo "=============================================="
echo "IMAU-FDM 1D to Gridded Maps"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=============================================="

# Change to script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Create logs directory if it doesn't exist
mkdir -p logs

# Load required modules (adjust for your ECMWF environment)
module load python3 2>/dev/null || true

# Activate virtual environment if available
if [ -f "$HOME/.venv/fdm/bin/activate" ]; then
    source "$HOME/.venv/fdm/bin/activate"
    echo "Activated virtual environment: $HOME/.venv/fdm"
fi

# Parse variables from command line or use default
VARIABLES="${1:-all}"

echo ""
echo "Processing variables: $VARIABLES"
echo ""

# Run the processing script
if [ "$VARIABLES" = "all" ]; then
    python make_1d_maps.py --var all --workers "$SLURM_CPUS_PER_TASK" --timestep 10day
else
    python make_1d_maps.py --var $VARIABLES --workers "$SLURM_CPUS_PER_TASK" --timestep 10day
fi

# Report completion
EXIT_CODE=$?
echo ""
echo "=============================================="
echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"
echo "=============================================="

exit $EXIT_CODE
