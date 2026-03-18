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

# Gridded files for each variable (produced by the post-processing pipeline).
# Add new variables here to make them available to the analysis functions.
VARIABLE_FILES = {
    'surfmelt': PROCESSED_DIR / 'FDM_surfmelt_FGRN055_1939-2023_10day.nc',
    'Runoff':   PROCESSED_DIR / 'FDM_Runoff_FGRN055_1939-2023_10day.nc',
}

# =============================================================================
# EVENT THRESHOLDS
# =============================================================================

# Annual total (mm w.e.) that must be exceeded for a year to count as an
# 'event' for each variable.  Adjust these to suit your definition.
EVENT_THRESHOLDS = {
    'surfmelt': 10.0,  # mm w.e.
    'Runoff':    1.0,  # mm w.e.
}

# =============================================================================
# ANALYSIS 2 - HISTOGRAM: comparison periods
# =============================================================================

# Each tuple is (start_year, end_year), both inclusive.
PERIOD_1 = (1940, 1970)
PERIOD_2 = (1990, 2020)

# =============================================================================
# ANALYSIS 3 - PERCOLATION MIGRATION: early / late windows and thresholds
# =============================================================================

# Comparison windows (both inclusive, typically ~10 years each)
EARLY_WINDOW = (1940, 1963)
LATE_WINDOW  = (2000, 2023)

# Dry facies: pixel had at most this many melt events in the window
DRY_EVENTS_MAX = 2

# Percolation zone: pixel had at least this many melt events in the window
PERCOLATION_EVENTS_MIN = 3
