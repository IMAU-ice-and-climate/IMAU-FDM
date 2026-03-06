"""
Configuration for IMAU-FDM 1D to gridded maps post-processing.

This module contains all configurable paths, parameters, and variable metadata.
Users can modify the settings below or override them via command-line arguments.
"""

from pathlib import Path
from datetime import datetime

# =============================================================================
# USER-CONFIGURABLE SETTINGS
# =============================================================================

# Time aggregation: 'daily', '10day', or 'monthly'
TIME_AGGREGATION = '10day'

# Spinup period for detrending (h_surf and FirnAir only)
# The model runs repeated forcing during this period
SPINUP_START = datetime(1940, 1, 1)
SPINUP_END = datetime(1970, 1, 1)

# Variables to process (None = all variables)
# Example: ['h_surf', 'FirnAir', 'Runoff']
VARIABLES_TO_PROCESS = None

# Variables that require detrending
DETREND_VARIABLES = ['h_surf', 'FirnAir']

# Number of parallel workers (None = use all available CPUs)
NUM_WORKERS = None

# =============================================================================
# PATHS - Adjust these for your environment
# =============================================================================

# Base directories
BASE_DIR = Path('/home/nld4814/perm/code/IMAU-FDM')
SCRATCH_DIR = Path('/home/nld4814/scratch')

# =============================================================================
# MODEL METADATA
# =============================================================================

# Domain name — also controls the expected mask/pointlist file locations:
#   <BASE_DIR>/reference/<DOMAIN>/<DOMAIN>_Masks.nc
#   <BASE_DIR>/reference/<DOMAIN>/IN_ll_<DOMAIN>.txt
DOMAIN = 'FGRN055'

# Input files
PROJECT_NAME = 'run_FGRN055-era055_1939-2023'
INPUT_DIR = SCRATCH_DIR / PROJECT_NAME / 'output'
POINTLIST_FILE = BASE_DIR / 'reference' / DOMAIN / f'IN_ll_{DOMAIN}.txt'
MASK_FILE = BASE_DIR / 'reference' / DOMAIN / f'{DOMAIN}_Masks.nc'

# Output directory
OUTPUT_DIR = SCRATCH_DIR / PROJECT_NAME / 'post-process'

# File naming pattern for input 1D files
INPUT_FILENAME_PATTERN = f'{DOMAIN}_era055_1D_{{point_id}}.nc'

# Model run period
MODEL_START = datetime(1939, 9, 1)
MODEL_END = datetime(2023, 12, 31)

# Timestep in the input files (seconds)
INPUT_TIMESTEP_SECONDS = 86400  # daily

# Grid dimensions are read from the mask file at runtime
# (see utils.get_grid_dimensions())

# =============================================================================
# VARIABLE METADATA
# =============================================================================

# All available variables in 1D output files with their metadata
# Note: 'vbouy' is the actual name in the files (typo in original model)
#
# Units note: In the raw 1D output files:
#   - Velocities are in m/yr (mean rate over output period)
#   - Flux variables (Runoff, refreeze, rain, surfmelt) are TOTALS per output
#     timestep (mm w.e. per day), NOT rates
#   - State variables (h_surf, FirnAir, TotLwc, icemass, Rho0) are instantaneous
#
# Aggregation: 'sum' for flux totals, 'mean' for rates and state variables
VARIABLES = {
    'h_surf': {
        'long_name': 'Surface height',
        'units': 'm',
        'needs_detrend': True,
        'aggregation': 'mean',  # State variable - instantaneous
    },
    'FirnAir': {
        'long_name': 'Firn air content',
        'units': 'm',
        'needs_detrend': True,
        'aggregation': 'mean',  # State variable - instantaneous
    },
    'vice': {
        'long_name': 'Ice velocity at model base',
        'units': 'm/yr',
        'needs_detrend': False,
        'aggregation': 'mean',  # Rate - already annualized in model
    },
    'vacc': {
        'long_name': 'Accumulation velocity',
        'units': 'm/yr',
        'needs_detrend': False,
        'aggregation': 'mean',  # Rate - already annualized in model
    },
    'vfc': {
        'long_name': 'Firn compaction velocity',
        'units': 'm/yr',
        'needs_detrend': False,
        'aggregation': 'mean',  # Rate - already annualized in model
    },
    'vmelt': {
        'long_name': 'Melt velocity',
        'units': 'm/yr',
        'needs_detrend': False,
        'aggregation': 'mean',  # Rate - already annualized in model
    },
    'vbouy': {  # Note: typo in original model output
        'long_name': 'Buoyancy velocity',
        'units': 'm/yr',
        'needs_detrend': False,
        'aggregation': 'mean',  # Rate - already annualized in model
    },
    'vsub': {
        'long_name': 'Sublimation velocity',
        'units': 'm/yr',
        'needs_detrend': False,
        'aggregation': 'mean',  # Rate - already annualized in model
    },
    'vsnd': {
        'long_name': 'Snowdrift velocity',
        'units': 'm/yr',
        'needs_detrend': False,
        'aggregation': 'mean',  # Rate - already annualized in model
    },
    'vtotal': {
        'long_name': 'Total velocity',
        'units': 'm/yr',
        'needs_detrend': False,
        'aggregation': 'mean',  # Rate - already annualized in model
    },
    'Runoff': {
        'long_name': 'Runoff',
        'units': 'mm w.e.',  # Total per output period
        'needs_detrend': False,
        'aggregation': 'sum',  # Flux - sum daily totals
    },
    'TotLwc': {
        'long_name': 'Total liquid water content',
        'units': 'mm',
        'needs_detrend': False,
        'aggregation': 'mean',  # State variable - instantaneous
    },
    'refreeze': {
        'long_name': 'Refreezing',
        'units': 'mm w.e.',  # Total per output period
        'needs_detrend': False,
        'aggregation': 'sum',  # Flux - sum daily totals
    },
    'rain': {
        'long_name': 'Rainfall',
        'units': 'mm w.e.',  # Total per output period
        'needs_detrend': False,
        'aggregation': 'sum',  # Flux - sum daily totals
    },
    'surfmelt': {
        'long_name': 'Surface melt',
        'units': 'mm w.e.',  # Total per output period
        'needs_detrend': False,
        'aggregation': 'sum',  # Flux - sum daily totals
    },
    'solin': {
        'long_name': 'Solar insolation',
        'units': 'W/m2',
        'needs_detrend': False,
        'aggregation': 'mean',  # Intensive variable - mean is appropriate
    },
    'icemass': {
        'long_name': 'Ice mass',
        'units': 'kg/m2',
        'needs_detrend': False,
        'aggregation': 'mean',  # State variable - instantaneous
    },
    'Rho0': {
        'long_name': 'Surface density',
        'units': 'kg/m3',
        'needs_detrend': False,
        'aggregation': 'mean',  # State variable - instantaneous
    },
}


def get_variable_names():
    """Return list of all available variable names."""
    return list(VARIABLES.keys())


def get_detrend_variables():
    """Return list of variables that need detrending."""
    return [var for var, meta in VARIABLES.items() if meta.get('needs_detrend', False)]


def get_aggregation_method(var_name):
    """
    Return the aggregation method for a variable.

    Parameters
    ----------
    var_name : str
        Variable name

    Returns
    -------
    str
        'sum' for flux variables, 'mean' for rates and state variables
    """
    if var_name in VARIABLES:
        return VARIABLES[var_name].get('aggregation', 'mean')
    return 'mean'  # Default to mean for unknown variables


def get_output_filename(var_name, timestep='10day', detrended=False):
    """
    Generate output filename for a gridded variable.

    Parameters
    ----------
    var_name : str
        Variable name
    timestep : str
        Time aggregation ('daily', '10day', 'monthly')
    detrended : bool
        Whether detrending was applied

    Returns
    -------
    str
        Output filename
    """
    detrend_suffix = '_detrended' if detrended else ''
    return f'FDM_{var_name}_{DOMAIN}_1939-2023_{timestep}{detrend_suffix}.nc'
