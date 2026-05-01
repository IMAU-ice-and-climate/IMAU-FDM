"""
Configuration for concatenating original FDM per-column output with extended runs.

Each extended file covers the full model period but has NaN for all timesteps
already present in the original run. Concatenation takes the original data as-is
and appends only the new timesteps (tail of the extended file) to produce a single
combined file spanning the full original + extension period.
"""

import os
from pathlib import Path

# =============================================================================
# MODEL METADATA
# =============================================================================

DOMAIN = 'FGRN055'
FORCING = 'era055'

# =============================================================================
# PATHS
# =============================================================================

SCRATCH_DIR = Path('/home/nld4814/scratch')

# Original run — provides data up to end of its simulation period
ORIG_RUN_NAME = 'run_FGRN055-era055_1939-2023'
ORIG_OUTPUT_DIR = SCRATCH_DIR / ORIG_RUN_NAME / 'output'

# Extended run output directories, checked in order — first directory containing
# a file for a given point wins. Put the preferred (most recent/fixed) run first.
EXT_OUTPUT_DIRS = [
    SCRATCH_DIR / 'extend-FGRN055-era055_1939-2023_to_2025-2' / 'output',
    SCRATCH_DIR / 'extend-FGRN055-era055_1939-2023_to_2025' / 'output',
]

# Combined output directory
OUTPUT_DIR = SCRATCH_DIR / 'FDM_FGRN055_output/output/points'

# =============================================================================
# PROCESSING
# =============================================================================

# Total number of grid points (1-based point IDs: 1 to N_POINTS)
N_POINTS = 58265

# File types to concatenate
FILE_TYPES = ['1D', '2D', '2Ddetail']

# Parallel workers — reads SLURM_CPUS_PER_TASK if available, otherwise defaults to 4
NUM_WORKERS = int(os.environ['SLURM_CPUS_PER_TASK']) if 'SLURM_CPUS_PER_TASK' in os.environ else 4
