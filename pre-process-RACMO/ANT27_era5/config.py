"""
Configuration for RACMO pre-processing (ANT27 / ERA5).
All paths and settings are defined here; the processing scripts import this.

=== Normal mode (full pipeline from scratch) ===

    python3 submit_jobs.py all

This runs three steps for all variables (one SLURM job per variable per step):
  1. years      — slice decade files into yearly files
  2. timeseries — slice yearly files into lon-band parts, concatenate into timeseries
  3. averages   — compute spinup-period averages

Run a subset of steps or variables:
    python3 submit_jobs.py years timeseries --skip-averages
    python3 submit_jobs.py all --vars precip snowfall tskin

Check what would be submitted without actually submitting:
    python3 submit_jobs.py all --dry-run

=== Extension mode (adding new years to an existing timeseries) ===

1. Update TS_END_YEAR in this file to the new target year (e.g. 2025).

2. Create yearly files for the new years:
    python3 submit_jobs.py years --skip-averages

3a. Safe first pass — create lon-band parts for new years only, in parts_extend/:
    python3 submit_jobs.py timeseries --extend-start 2024 --extend-end 2025 \\
        --extend-mode parts-only --skip-averages

3b. Once satisfied, append new parts to the existing timeseries files:
    python3 submit_jobs.py timeseries --extend-start 2024 --extend-end 2025 \\
        --extend-mode append --skip-averages

    The old timeseries files are renamed to .bak before being replaced.
    The new files are named with the updated end year, e.g.:
      precip_FGRN055_era055_1939-2025_p1.nc

Notes:
  - Steps 1 and 2 skip files that already exist, so it is safe to rerun.
  - averages are based on the spinup period (AVE_START_YEAR–AVE_END_YEAR) and
    do not need to be recomputed when extending the timeseries.
"""

import os
from pathlib import Path

DOMAIN = "ANT27"
FORCING = "ERA5"
PROJECT_NAME = f"{DOMAIN}_{FORCING}"

VARS = ["evap", "ff10m", "precip", "sndiv", "snowfall", "snowmelt", "tskin"]

SCRATCH = os.environ.get("SCRATCH", "")
PERM   = os.environ.get("PERM",   "")

BASE_DIR    = Path(SCRATCH) / PROJECT_NAME
SCRIPTS_DIR = Path(PERM) / "code" / "IMAU-FDM" / "init" / f"{DOMAIN}_{FORCING}"

# Timeseries extent.
# Update TS_END_YEAR here when extending the timeseries to new years,
# then use --extend-start/--extend-end in submit_jobs.py to specify
# which new years to process in this run.
TS_START_YEAR = 1979
TS_END_YEAR   = 2025

# Spinup averaging period
AVE_START_YEAR = 1979
AVE_END_YEAR   = 2020

# Longitude band parameters
NUM_LONG_BANDS = 17
CELL_WIDTH     = 6   # rlon indices per band

# Raw RACMO file suffix 

FNAME_SUFFIX = f"{DOMAIN}.ERA5-3H_RACMO2.3p2.3H.nc"

# Derived paths
RAW_DIR     = BASE_DIR / "raw"
YEARS_DIR   = BASE_DIR / "process-RACMO" / f"years-{TS_START_YEAR}"
FILES_DIR   = BASE_DIR / "input" / "timeseries"
AVE_DIR     = BASE_DIR / "input" / "averages"
JOBFILE_DIR = BASE_DIR / "process-RACMO" / "jobs"
LOGFILE_DIR = BASE_DIR / "process-RACMO" / "logs"

# SLURM settings
SLURM_QUEUE = "nf"
SLURM_TIME  = "48:00:00"
SLURM_MEM   = "4000mb"
