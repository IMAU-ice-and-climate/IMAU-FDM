"""
Configuration for extend_variable_2d.py and extend_variable_2ddetail.py.

Points at the merged per-column files produced by extend_pointfiles.py,
the existing 1939-2023 gridded files, and the output directory for the
extended 1939-2025 gridded files.
"""

import os
from pathlib import Path
from datetime import datetime

SCRATCH_DIR = Path('/home/nld4814/scratch')
BASE_DIR    = Path('/home/nld4814/perm/code/IMAU-FDM')

DOMAIN  = 'FGRN055'
FORCING = 'era055'

# =============================================================================
# INPUT PATHS
# =============================================================================

# Merged per-column files (output of extend_pointfiles.py), covering 1939-2025
EXT_INPUT_DIR = SCRATCH_DIR / 'FDM_FGRN055_output/output/points'

# Filename patterns for merged per-column files
EXT_2D_FILENAME_PATTERN       = 'FGRN055_era055_2D_{point_id}.nc'
EXT_2DDETAIL_FILENAME_PATTERN = 'FGRN055_era055_2Ddetail_{point_id}.nc'

# Existing 1939-2023 gridded post-process files
ORIG_OUTPUT_DIR = SCRATCH_DIR / 'run_FGRN055-era055_1939-2023/post-process'

# Static reference files (same as create_1D_2D_2Ddetail_files/config.py)
POINTLIST_FILE = BASE_DIR / 'reference' / DOMAIN / f'IN_ll_{DOMAIN}.txt'
MASK_FILE      = BASE_DIR / 'reference' / DOMAIN / f'{DOMAIN}_Masks.nc'
GRID_FILE      = BASE_DIR / 'reference' / DOMAIN / f'{DOMAIN}_grid.nc'

# =============================================================================
# OUTPUT PATH
# =============================================================================

EXT_OUTPUT_DIR = SCRATCH_DIR / 'run_FGRN055-era055_1939-2025/post-process'

# =============================================================================
# TIME SETTINGS
# =============================================================================

MODEL_START    = datetime(1939, 9,  1)
ORIG_MODEL_END = datetime(2023, 12, 31)
EXT_MODEL_END  = datetime(2025, 12, 31)

ORIG_DATE_TAG  = '1939-2023'
EXT_DATE_TAG   = '1939-2025'

# =============================================================================
# PROCESSING
# =============================================================================

# Parallel workers — reads SLURM_CPUS_PER_TASK if available
NUM_WORKERS = int(os.environ['SLURM_CPUS_PER_TASK']) if 'SLURM_CPUS_PER_TASK' in os.environ else 4
