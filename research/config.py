"""
Configuration for percolation zone research analyses.

Points to the processed gridded files produced by the
create_1D_2D_2Ddetail_files pipeline, and sets default parameters
for the three event-frequency analyses in transient_percolation_analysis.py.
"""

from pathlib import Path

# =============================================================================
# PATHS
# =============================================================================

SCRATCH_DIR = Path('/home/nld4814/scratch')
PROJECT_NAME = 'run_FGRN055-era055_1939-2023'
PROCESSED_DIR = SCRATCH_DIR / PROJECT_NAME / 'post-process'

# Static grid metadata: topography, masks, drainage basins.
MASKS_FILE    = Path('/home/nld4814/perm/code/IMAU-FDM/reference/FGRN055/FGRN055_Masks.nc')

# Point list: maps file number N → (rlat_idx, rlon_idx) on the FGRN055 grid.
# Row N-1 (0-indexed) → file FGRN055_era055_2Ddetail_{N}.nc
POINTLIST_FILE = Path('/home/nld4814/perm/code/IMAU-FDM/reference/FGRN055/IN_ll_FGRN055.txt')

# Subsurface gridded diagnostics (produced by post-processing pipeline).
Z830_FILE         = PROCESSED_DIR / 'FDM_z830_FGRN055_1939-2023_2D.nc'
T10M_FILE         = PROCESSED_DIR / 'FDM_T10m_FGRN055_1939-2023_2Ddetail.nc'
FIRN_MEMORY_FILE  = PROCESSED_DIR / 'FDM_firn_memory_FGRN055_1939-2023_2D.nc'
ICE_LENS_FILE     = PROCESSED_DIR / 'FDM_ice_lens_FGRN055_1939-2023_2Ddetail.nc'

# Gridded files for each variable (produced by the post-processing pipeline).
# Add new variables here to make them available to the analysis functions.
VARIABLE_FILES = {
    'surfmelt': PROCESSED_DIR / 'FDM_surfmelt_FGRN055_1939-2023_10day.nc',
    'Runoff':   PROCESSED_DIR / 'FDM_Runoff_FGRN055_1939-2023_10day.nc',
    'vacc':     PROCESSED_DIR / 'FDM_vacc_FGRN055_1939-2023_10day.nc',
    'refreeze': PROCESSED_DIR / 'FDM_refreeze_FGRN055_1939-2023_10day.nc',
    'vsub':     PROCESSED_DIR / 'FDM_vsub_FGRN055_1939-2023_10day.nc',
}

# =============================================================================
# EVENT THRESHOLDS
# =============================================================================

# Annual total (mm w.e.) that must be exceeded for a year to count as an
# 'event' for each variable.  Adjust these to suit your definition.
EVENT_THRESHOLDS = {
    'surfmelt': 10.0,  # mm w.e.
    'Runoff':    10.0,  # mm w.e.
}

# =============================================================================
# ANALYSIS 2 - HISTOGRAM: comparison periods
# =============================================================================

# Each tuple is (start_year, end_year), both inclusive.
PERIOD_1 = (1940, 1970)
PERIOD_2 = (1993, 2023)

# =============================================================================
# ANALYSIS 3 - PERCOLATION MIGRATION: early / late windows and thresholds
# =============================================================================

# Comparison windows (both inclusive, typically ~10 years each)
EARLY_WINDOW = (1940, 1963)
LATE_WINDOW  = (2000, 2023)

# Dry facies: pixel had at most this many events exceeding the threshold in the window
DRY_EVENTS_MAX = 0

# Percolation zone: pixel had at least this many events exceeding the threshold in the window
PERCOLATION_EVENTS_MIN = 1
