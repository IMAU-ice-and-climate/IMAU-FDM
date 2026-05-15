# IMAU-FDM Post-Processing

This directory contains two pipelines for converting IMAU-FDM per-column output into
gridded NetCDF files and for extending existing gridded files to a later end date.

```
post-process/
├── create_1D_2D_2Ddetail_files/   # Initial gridded file creation (1939-2023)
├── extend/                         # Extend existing gridded files (1939-2025)
├── visualization/                  # Plotting and data-loading utilities
└── README.md                       # This file
```

---

## Part 1 — Creating gridded files from scratch

Scripts live in `create_1D_2D_2Ddetail_files/`.
All three scripts read per-column output from the FDM run directory and write one
gridded NetCDF per variable to the post-process output directory.

Output filename convention: `FDM_{variable}_{domain}_{date_range}_{timestep}.nc`
e.g. `FDM_Runoff_FGRN055_1939-2023_10day.nc`

### 1a. 1D variables (daily point output → gridded)

These are scalar time series per point (runoff, surface melt, firn air content, etc.).
The script can aggregate to 10-day or monthly means/sums and optionally detrend
selected variables (h_surf, FirnAir) to remove the spinup drift.

**Script:** `make_1d_files.py`

```
# All variables, 10-day aggregation
python3 make_1d_files.py --var all --timestep 10day \
    --spinup-start 1940 --spinup-end 1970 --workers 16

# Specific variables
python3 make_1d_files.py --var Runoff surfmelt --timestep 10day --workers 8

# List available variables
python3 make_1d_files.py --list-vars
```

**SLURM:** `sbatch submit_make_1d_files.sh`
(48 h wall time, 64 GB RAM, 16 CPUs)

Key arguments:

| Flag | Default | Description |
|------|---------|-------------|
| `--var` | required | Variable name(s), or `all` |
| `--timestep` | `10day` | `daily`, `10day`, or `monthly` |
| `--spinup-start` | 1940 | Year spinup detrend window starts |
| `--spinup-end` | 1970 | Year spinup detrend window ends |
| `--workers` | 16 | Parallel worker processes |
| `--max-points` | — | Cap for testing (e.g. 500) |

**Output period slicing:** set `OUTPUT_START` and `OUTPUT_END` in `config.py` to save
only a sub-period of the full model run (e.g. model covers 1939-2025, output only 1958-2023).
Set either to `None` to disable. Applies to all three file types; for 2D and 2Ddetail scripts
the `--start-year`/`--end-year` CLI flags override the config defaults.

### 1b. 2D variables (layered monthly profiles → scalar gridded)

These files contain full depth profiles (3000 layers, ~122 m) at 30-day intervals.
`make_2d_files.py` derives one scalar per timestep per point — typically the depth
at which a density threshold is crossed (e.g. the firn–ice transition depth z830).

**Script:** `make_2d_files.py`

```
# Depth where density first reaches 830 kg/m³ (firn–ice transition)
python3 make_2d_files.py --var dens --threshold 830 --output-var z830 \
    -o /path/to/output --workers 8

# Depth where density first reaches 550 kg/m³
python3 make_2d_files.py --var dens --threshold 550 --output-var z550 \
    -o /path/to/output --workers 8
```

**SLURM:** `sbatch submit_make_2d_files.sh`
(24 h, 32 GB, 8 CPUs)

Key arguments:

| Flag | Description |
|------|-------------|
| `--var` | Variable in the 2D per-column file (`dens`, `temp`, etc.) |
| `--threshold` | Density threshold (kg/m³) to find depth crossing |
| `--output-var` | Name of the output variable (e.g. `z830`) |
| `-o` | Output directory |
| `-n` | Parallel workers |

### 1c. 2Ddetail variables (high-res near-surface profiles → scalar gridded)

These files contain the top 20 m at 4 cm resolution (500 layers), at 10-day intervals.
`make_2Ddetail_files.py` extracts a value at a target depth or averages over a depth
range — used for near-surface temperature (T10m) and surface snow density (SSN).

**Script:** `make_2Ddetail_files.py`

```
# 10 m temperature
python3 make_2Ddetail_files.py --var temp --depth 10 --output-var T10m \
    -o /path/to/output --workers 8

# Surface snow density (0–0.5 m mean)
python3 make_2Ddetail_files.py --var dens --depth-begin 0 --depth-end 0.5 \
    --output-var SSN -o /path/to/output --workers 8

# Near-surface liquid water content (0–0.2 m)
python3 make_2Ddetail_files.py --var lwc --depth-begin 0 --depth-end 0.2 \
    --output-var LWC_surf -o /path/to/output --workers 8
```

**SLURM:** `sbatch submit_make_2ddetail_files.sh`
(24 h, 32 GB, 8 CPUs)

Key arguments:

| Flag | Description |
|------|-------------|
| `--var` | Variable in the 2Ddetail file (`temp`, `dens`, `lwc`, etc.) |
| `--depth` | Single target depth in metres (e.g. `10` for T10m) |
| `--depth-begin`, `--depth-end` | Depth range for averaging (m) |
| `--output-var` | Name of the output variable |
| `-o` | Output directory |
| `-n` | Parallel workers |

---

## Part 2 — Extending existing gridded files to 2025

Scripts live in `extend/`.  There are two distinct steps:

1. **Extend pointfiles** — merge per-column output from the original run (1939-2023)
   with the extension run (2024-2025) into combined per-point files. This is a
   prerequisite for step 2.
2. **Extend variable files** — read the merged per-column files and append the
   2024-2025 timesteps to existing 1939-2023 gridded files, producing
   `1939-2025` domain-wide maps.

> **Note on 1D variables:** 1D gridded files (Runoff, surfmelt, h_surf, etc.)
> cannot be extended with a simple slice — their resampling windows are anchored
> to `MODEL_START` and detrending requires the full post-spinup series. These must
> be remade from scratch using `create_1D_2D_2Ddetail_files/make_1d_files.py` once
> the merged pointfiles are available.

### 2a. Extend pointfiles (prerequisite)

Merges per-column 1D, 2D, and 2Ddetail files from the original and extension
runs into a single file per point covering the full 1939-2025 period.

Configuration: `extend_pointfiles_config.py`

```
# Run once (~1 h with 16 workers)
sbatch submit_extend_pointfiles.sh

# Or directly, e.g. to reprocess a single point:
python3 extend_pointfiles.py --start 56928 --end 56928 --workers 1
```

Key settings in `extend_pointfiles_config.py`:

| Setting | Value |
|---------|-------|
| `ORIG_OUTPUT_DIR` | `scratch/run_FGRN055-era055_1939-2023/output` |
| `EXT_OUTPUT_DIRS` | extension run output directories (checked in order) |
| `OUTPUT_DIR` | `scratch/FDM_FGRN055_output/output/points` |

### 2b. Extend variable files — 2D (z830, z550, etc.)

Reads the merged 2D per-column files and appends new timesteps to the existing
gridded file. Each timestep is independent (no resampling), so only the
2024-2025 slice of each per-column file is actually used.

Configuration: `extend_variable_config.py`

**Script:** `extend_variable_2d.py`

```
python3 extend_variable_2d.py --output-var z830 --var dens --threshold 830 --workers 16
python3 extend_variable_2d.py --output-var z550 --var dens --threshold 550 --workers 16
```

**SLURM:** `sbatch submit_extend_variable_2d.sh`
(12 h, 32 GB, 16 CPUs)

Key arguments:

| Flag | Description |
|------|-------------|
| `--output-var` | Output variable name (e.g. `z830`) |
| `--var` | Variable in the 2D per-column file |
| `--threshold` | Density threshold (kg/m³) |
| `--workers` | Parallel workers |

### 2c. Extend variable files — 2Ddetail (T10m, SSN, etc.)

Same approach as 2b — each 10-day timestep is independent.

**Script:** `extend_variable_2ddetail.py`

```
python3 extend_variable_2ddetail.py --output-var T10m --var temp --depth 10 --workers 16
python3 extend_variable_2ddetail.py --output-var SSN --var dens \
    --depth-begin 0 --depth-end 0.5 --workers 16
python3 extend_variable_2ddetail.py --output-var LWC_surf --var lwc \
    --depth-begin 0 --depth-end 0.2 --workers 16
```

**SLURM:** `sbatch submit_extend_variable_2ddetail.sh`
(12 h, 32 GB, 16 CPUs)

Key settings in `extend_variable_config.py`:

| Setting | Value |
|---------|-------|
| `EXT_INPUT_DIR` | `scratch/FDM_FGRN055_output/output/points` |
| `ORIG_OUTPUT_DIR` | `scratch/run_FGRN055-era055_1939-2023/post-process` |
| `EXT_OUTPUT_DIR` | `scratch/run_FGRN055-era055_1939-2025/post-process` |
| `ORIG_MODEL_END` | 2023-12-31 |
| `EXT_MODEL_END` | 2025-12-31 |

---

## Part 3 — Loading and visualising gridded files

Utilities live in `visualization/`.  Import from a notebook or script:

```python
from visualization.visualization import load_gridded_data, plot_variable_map

# Load 1939-2023 Runoff (default)
ds = load_gridded_data('Runoff', timestep='10day')

# Load extended 1939-2025 file from a different directory
ds = load_gridded_data('Runoff', timestep='10day',
                       output_dir=Path('/path/to/1939-2025/post-process'),
                       date_tag='1939-2025')

# Load detrended h_surf
ds = load_gridded_data('h_surf', timestep='10day', detrended=True)
```

`load_gridded_data` parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `var_name` | required | Variable name (must match filename) |
| `timestep` | `'10day'` | `'daily'`, `'10day'`, or `'monthly'` |
| `detrended` | `None` | Auto-detect from `DETREND_VARIABLES`; or `True`/`False` |
| `output_dir` | `config.OUTPUT_DIR` | Directory containing the `.nc` file |
| `date_tag` | `None` | Date string in filename; defaults to `config.DATE_TAG` (`'1939-2023'`) |

### Filename convention

All gridded files follow:
```
FDM_{variable}_{domain}_{date_range}_{timestep}.nc
```

Examples:
- `FDM_Runoff_FGRN055_1939-2023_10day.nc`
- `FDM_h_surf_FGRN055_1939-2025_10day_detrended.nc`
- `FDM_z830_FGRN055_1939-2023_30day.nc`
- `FDM_T10m_FGRN055_1939-2023_10day.nc`
- `FDM_firn_memory_FGRN055_1939-2023_30day.nc`
- `FDM_ice_lens_FGRN055_1939-2023_10day.nc`

---

## Typical workflow for extending the model in time

```
# 1. Create gridded files for the base run (1939-2023)
cd create_1D_2D_2Ddetail_files/
sbatch submit_make_2d_files.sh
sbatch submit_make_2ddetail_files.sh

# 2. Merge per-column output with extension years into full-period pointfilchkes
cd ../extend/
sbatch submit_extend_pointfiles.sh

# 3a. Remake 1D gridded files for 1939-2025 (full reprocessing required)
cd ../create_1D_2D_2Ddetail_files/
sbatch submit_make_1d_files.sh   # point at merged per-column dir

# 3b. Extend 2D/2Ddetail gridded variable files (append 2024-2025 slice only)
cd ../extend/
sbatch submit_extend_variable_2d.sh
sbatch submit_extend_variable_2ddetail.sh
```

Output directories:
- `scratch/run_FGRN055-era055_1939-2023/post-process/`  — base run
- `scratch/FGRN055-era055_1939-2025/post-process/`  — extended run
