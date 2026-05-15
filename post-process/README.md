# IMAU-FDM Post-Processing

This directory contains two pipelines for converting IMAU-FDM per-column output into
gridded NetCDF files and for extending existing gridded files to a later end date.

```
post-process/
├── create_1D_2D_2Ddetail_files/   # Initial gridded file creation
│   ├── config.py
│   ├── run_config.py
│   ├── utils.py
│   ├── make_1d_files.py
│   ├── make_2d_files.py
│   ├── make_2Ddetail_files.py
│   ├── submit_make_1d_files.sh
│   ├── submit_make_2d_files.sh
│   ├── submit_make_2ddetail_files.sh
│   └── logs/
├── extend/                         # Extend existing gridded files to a later end date
├── visualization/                  # Plotting and data-loading utilities
└── README.md                       # This file
```

---

## Output file structure

All gridded files follow the naming convention:

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

All output files share a common internal structure:

```
Dimensions:  time, rlat, rlon
Coordinates: time (datetime), rlat (degrees), rlon (degrees)
Data vars:   <variable>(time, rlat, rlon), lat(rlat,rlon), lon(rlat,rlon),
             x(rlat,rlon), y(rlat,rlon), x_FDM(rlon), y_FDM(rlat),
             rotated_pole, crs
```

- `rlat`/`rlon`: rotated-pole degree coordinates from the RACMO grid file
- `x`/`y`: EPSG:3413 projection coordinates in metres
- `x_FDM`/`y_FDM`: 0-based integer grid indices for cross-referencing model files
- `time`: CF-compliant datetime (xarray encodes as `days since ...`)
- CRS: EPSG:3413 (WGS 84 / NSIDC Sea Ice Polar Stereographic North)

---

## Part 1 — Creating gridded files from scratch

Scripts live in `create_1D_2D_2Ddetail_files/`.  The model produces one file per
grid column (~58,000 for FGRN055); this pipeline assembles them into 3D gridded
files `(time, rlat, rlon)` in CF-1.8 compliant NetCDF4.

Three types of per-column output are handled by separate scripts:

| Input type | What it contains | Script |
|---|---|---|
| **1D** | Daily surface/integrated variables | `make_1d_files.py` |
| **2D** | Monthly depth profiles (layer × time, ~122 m) | `make_2d_files.py` |
| **2Ddetail** | 10-day high-res near-surface profiles (4 cm layers, upper 20 m) | `make_2Ddetail_files.py` |

### Configuration (`config.py`)

Edit before running. Key settings:

```python
BASE_DIR     = Path('/home/nld4814/perm/code/IMAU-FDM')
SCRATCH_DIR  = Path('/home/nld4814/scratch')
DOMAIN       = 'FGRN055'
PROJECT_NAME = 'run_FGRN055-era055_1939-2023'
MODEL_START  = datetime(1939, 9, 1)
MODEL_END    = datetime(2023, 12, 31)
```

Expected reference files:

```
<BASE_DIR>/reference/<DOMAIN>/<DOMAIN>_Masks.nc
<BASE_DIR>/reference/<DOMAIN>/<DOMAIN>_grid.nc
<BASE_DIR>/reference/<DOMAIN>/IN_ll_<DOMAIN>.txt
```

Full `config.py` settings reference:

| Setting | Description |
|---|---|
| `BASE_DIR` | Root of the IMAU-FDM code/reference directory |
| `SCRATCH_DIR` | Root of the model output directory |
| `DOMAIN` | Domain name; controls mask/grid/pointlist file paths |
| `PROJECT_NAME` | Run name; combined with `SCRATCH_DIR` to find output |
| `MODEL_START` / `MODEL_END` | Simulation period |
| `INPUT_TIMESTEP_SECONDS` | Timestep of 1D input files in seconds (default: 86400) |
| `TIME_AGGREGATION_1D` | Output time aggregation for 1D resampling: `daily`, `10day`, `monthly` |
| `TIME_AGGREGATION_2D` | Expected aggregation of 2D files (must match model output) |
| `TIME_AGGREGATION_2Ddetail` | Expected aggregation of 2Ddetail files (must match model output) |
| `SPINUP_START` / `SPINUP_END` | Period used for spinup detrending of `h_surf` and `FirnAir` |
| `NUM_WORKERS` | Default number of parallel workers (`None` = all CPUs) |
| `OUTPUT_START` / `OUTPUT_END` | Sub-period to save (set either to `None` to disable) |

> **Note on time aggregation:** `TIME_AGGREGATION_1D` controls resampling and can be freely
> changed. `TIME_AGGREGATION_2D` and `TIME_AGGREGATION_2Ddetail` are fixed by the model run —
> changing them without re-running the model will raise an error at startup.

> **Output period slicing:** `OUTPUT_START`/`OUTPUT_END` let you save only a sub-period of the
> full model run (e.g. model covers 1939-2025, output only 1958-2023). For 2D and 2Ddetail
> scripts the `--start-year`/`--end-year` CLI flags override config defaults.

### 1a. 1D variables

Daily surface/integrated scalar time series per point (runoff, surface melt, firn air
content, etc.). Aggregates to 10-day or monthly means/sums; optionally detrends
`h_surf` and `FirnAir` to remove spinup drift.

**Script:** `make_1d_files.py` &nbsp; **SLURM:** `sbatch submit_make_1d_files.sh` (48 h, 64 GB, 16 CPUs)

```bash
# All variables, 10-day aggregation
python3 make_1d_files.py --var all --timestep 10day \
    --spinup-start 1940 --spinup-end 1970 --workers 16

# Specific variables
python3 make_1d_files.py --var Runoff surfmelt --timestep 10day --workers 8

# List available variables
python3 make_1d_files.py --list-vars
```

Arguments:

```
--var, -v         Variable(s) to process, or "all"
--timestep, -t    Output time aggregation: daily | 10day | monthly
--spinup-start    Spinup start year for detrending (default: 1940)
--spinup-end      Spinup end year for detrending (default: 1970)
--workers, -w     Number of parallel workers
--output-dir, -o  Output directory (default: from config.py)
--max-points      Limit to N columns (useful for testing)
--list-vars       Print available variables and exit
--quiet, -q       Suppress progress output
```

Available variables:

| Variable | Long name | Units | Aggregation |
|---|---|---|---|
| `h_surf` | Surface height | m | mean |
| `FirnAir` | Firn air content | m | mean |
| `Runoff` | Runoff | mm w.e. | sum |
| `refreeze` | Refreezing | mm w.e. | sum |
| `rain` | Rainfall | mm w.e. | sum |
| `surfmelt` | Surface melt | mm w.e. | sum |
| `TotLwc` | Total liquid water content | mm | mean |
| `icemass` | Ice mass | kg/m² | mean |
| `Rho0` | Surface density | kg/m³ | mean |
| `solin` | Solar insolation | W/m² | mean |
| `vice` | Ice velocity at base | m/yr | mean |
| `vacc` | Accumulation velocity | m/yr | mean |
| `vfc` | Firn compaction velocity | m/yr | mean |
| `vmelt` | Melt velocity | m/yr | mean |
| `vsub` | Sublimation velocity | m/yr | mean |
| `vsnd` | Snowdrift velocity | m/yr | mean |
| `vtotal` | Total velocity | m/yr | mean |
| `vbouy` | Buoyancy velocity | m/yr | mean |

### 1b. 2D variables

Monthly depth profiles (3000 layers, ~122 m). `make_2d_files.py` derives one scalar per
timestep per point — typically the depth at which a density threshold is crossed (e.g.
the firn–ice transition depth z830).

**Script:** `make_2d_files.py` &nbsp; **SLURM:** `sbatch submit_make_2d_files.sh` (24 h, 32 GB, 8 CPUs)

```bash
# Depth where density first reaches 830 kg/m³ (firn–ice transition)
python3 make_2d_files.py -o OUTPUT_DIR -v dens -t 830 --output-var z830

# Depth where density first reaches 550 kg/m³
python3 make_2d_files.py -o OUTPUT_DIR -v dens -t 550 --output-var z550

# Test on a single year
python3 make_2d_files.py -o OUTPUT_DIR -v dens -t 830 --output-var z830 \
    --start-year 2020 --end-year 2020
```

Arguments:

```
-o, --output-dir    Directory containing model output files (required)
-v, --var           Input variable: dens | temp | lwc | ... (default: dens)
-t, --threshold     Find depth where variable crosses this value
-d, --depth         Extract value at this depth in metres
--output-var        Name for the output variable (required, e.g. z830)
--reference-dir     Directory with mask and pointlist files
--processed-dir     Output directory for gridded files
--start-year        Only output timesteps >= this year
--end-year          Only output timesteps <= this year
-n, --num-workers   Number of parallel workers
--dry-run           Print detected configuration and exit
```

### 1c. 2Ddetail variables

High-resolution near-surface profiles (4 cm layers, upper ~20 m) at 10-day intervals.
`make_2Ddetail_files.py` extracts a value at a target depth or averages over a depth
range — used for near-surface temperature (T10m) and surface snow density (SSN).

**Script:** `make_2Ddetail_files.py` &nbsp; **SLURM:** `sbatch submit_make_2ddetail_files.sh` (24 h, 32 GB, 8 CPUs)

```bash
# Surface snow density (upper 50 cm)
python3 make_2Ddetail_files.py -o OUTPUT_DIR -v dens \
    --depth-begin 0 --depth-end 0.5 --output-var SSN

# 10 m firn temperature
python3 make_2Ddetail_files.py -o OUTPUT_DIR -v temp -d 10 --output-var T10m

# Liquid water content in upper 20 cm
python3 make_2Ddetail_files.py -o OUTPUT_DIR -v lwc \
    --depth-begin 0 --depth-end 0.2 --output-var LWC_surf
```

Arguments:

```
-o, --output-dir    Directory containing model output files (required)
-v, --var           Input variable: dens | temp | lwc | refreeze (required)
--output-var        Name for the output variable (required)

Depth selection (choose one):
  --z-begin / --z-end         Layer indices (0-based from surface)
  --depth-begin / --depth-end Depth range in metres
  -d, --depth                 Single depth in metres

--operation         Operation over depth range: average | sum (default: average)
--layer-thickness   Layer thickness in metres (default: auto-detect from files)
--reference-dir     Directory with mask and pointlist files
--processed-dir     Output directory for gridded files
--start-year        Only output timesteps >= this year
--end-year          Only output timesteps <= this year
-n, --num-workers   Number of parallel workers
--dry-run           Print detected configuration and exit
```

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

```bash
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

**Script:** `extend_variable_2d.py` &nbsp; **SLURM:** `sbatch submit_extend_variable_2d.sh` (12 h, 32 GB, 16 CPUs)

```bash
python3 extend_variable_2d.py --output-var z830 --var dens --threshold 830 --workers 16
python3 extend_variable_2d.py --output-var z550 --var dens --threshold 550 --workers 16
```

Arguments:

| Flag | Description |
|------|-------------|
| `--output-var` | Output variable name (e.g. `z830`) |
| `--var` | Variable in the 2D per-column file |
| `--threshold` | Density threshold (kg/m³) |
| `--workers` | Parallel workers |

### 2c. Extend variable files — 2Ddetail (T10m, SSN, etc.)

Same approach as 2b — each 10-day timestep is independent.

**Script:** `extend_variable_2ddetail.py` &nbsp; **SLURM:** `sbatch submit_extend_variable_2ddetail.sh` (12 h, 32 GB, 16 CPUs)

```bash
python3 extend_variable_2ddetail.py --output-var T10m --var temp --depth 10 --workers 16
python3 extend_variable_2ddetail.py --output-var SSN --var dens \
    --depth-begin 0 --depth-end 0.5 --workers 16
python3 extend_variable_2ddetail.py --output-var LWC_surf --var lwc \
    --depth-begin 0 --depth-end 0.2 --workers 16
```

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

### Inspect output

```python
import xarray as xr

ds = xr.open_dataset("path/to/output.nc")
print(ds)

# Select by year (time is datetime)
ds_2022 = ds.sel(time="2022")

# Select nearest timestep
da = ds['FirnAir'].sel(time="2022-06", method="nearest")
```

---

## Typical workflow for extending the model in time

```bash
# 1. Create gridded files for the base run (1939-2023)
cd create_1D_2D_2Ddetail_files/
sbatch submit_make_2d_files.sh
sbatch submit_make_2ddetail_files.sh

# 2. Merge per-column output with extension years into full-period pointfiles
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
