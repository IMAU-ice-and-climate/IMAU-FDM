# IMAU-FDM Post-Processing: 1D → Gridded NetCDF

Converts IMAU-FDM model output into gridded NetCDF files on the RACMO rotated-pole
grid (rlat × rlon), ready for analysis and visualisation.

---

## Overview

The model produces one file per grid column (~58,000 files for FGRN055). This
pipeline assembles those into 3D gridded files `(time, rlat, rlon)` in CF-1.8
compliant NetCDF4.

Three types of input file are handled separately, each with its own script:

| Input type | What it contains | Script |
|---|---|---|
| **1D** | Daily surface/integrated variables (one per column) | `make_1d_files.py` |
| **2D** | Monthly depth profiles (layer × time) | `make_2d_files.py` |
| **2Ddetail** | 10-day high-res near-surface profiles (4 cm layers, upper 20 m) | `make_2Ddetail_files.py` |

---

## Output file structure

All output files share a common structure:

```
Dimensions:  time, rlat, rlon
Coordinates: time (datetime), rlat (degrees), rlon (degrees)
Data vars:   <variable>(time, rlat, rlon), lat(rlat,rlon), lon(rlat,rlon),
             x(rlat,rlon), y(rlat,rlon), x_FDM(rlon), y_FDM(rlat),
             rotated_pole, crs
```

- `rlat`/`rlon` are rotated-pole degree coordinates from the RACMO grid file
- `x`/`y` are EPSG:3413 projection coordinates in metres
- `x_FDM`/`y_FDM` are 0-based integer grid indices for cross-referencing model files
- `time` is CF-compliant datetime (xarray encodes as `days since ...`)
- CRS is EPSG:3413 (WGS 84 / NSIDC Sea Ice Polar Stereographic North)

---

## Quick start

### 1. Configure paths

Edit `config.py` before running. The key settings are:

```python
BASE_DIR    = Path('/home/nld4814/perm/code/IMAU-FDM')
SCRATCH_DIR = Path('/home/nld4814/scratch')
DOMAIN      = 'FGRN055'
PROJECT_NAME = 'run_FGRN055-era055_1939-2023'
MODEL_START = datetime(1939, 9, 1)
MODEL_END   = datetime(2023, 12, 31)
```

The script expects the following reference files to exist:
```
<BASE_DIR>/reference/<DOMAIN>/<DOMAIN>_Masks.nc
<BASE_DIR>/reference/<DOMAIN>/<DOMAIN>_grid.nc
<BASE_DIR>/reference/<DOMAIN>/IN_ll_<DOMAIN>.txt
```

### 2. Run interactively (small test)

```bash
cd /home/nld4814/perm/code/IMAU-FDM/post-process/create_1D_2D_2Ddetail_files

# 1D: process Runoff for 100 columns only (fast test)
python3 make_1d_files.py --var Runoff --max-points 100 --workers 4

# 2D: z830 for one year only
python3 make_2d_files.py \
  -o /home/nld4814/scratch/run_FGRN055-era055_1939-2023/output \
  -v dens --threshold 830 --output-var z830 \
  --start-year 2020 --end-year 2020

# 2Ddetail: surface snow density for one year
python3 make_2Ddetail_files.py \
  -o /home/nld4814/scratch/run_FGRN055-era055_1939-2023/output \
  -v dens --depth-begin 0 --depth-end 0.5 --output-var SSN \
  --start-year 2020 --end-year 2020
```

### 3. Submit as SLURM batch jobs (full run)

```bash
# 1D: all variables
sbatch submit_make_1d_files.sh

# 1D: specific variables
sbatch submit_make_1d_files.sh h_surf FirnAir Runoff

# 2D: density at 830 kg/m³
sbatch submit_make_2d_files.sh --var dens --threshold 830 --output-var z830

# 2D: temperature at 10 m depth
sbatch submit_make_2d_files.sh --var temp --depth 10 --output-var T10m

# 2Ddetail: surface snow density
sbatch submit_make_2ddetail_files.sh -v dens --depth-begin 0 --depth-end 0.5 --output-var SSN

# 2Ddetail: 10 m temperature
sbatch submit_make_2ddetail_files.sh -v temp --depth 10 --output-var T10m
```

Logs are written to `logs/slurm-<jobid>.out`.

---

## Script reference

### `make_1d_files.py` — 1D surface/integrated variables

Reads daily 1D files, optionally detrends spinup artefacts, resamples to the
configured time aggregation, and assembles a gridded output file.

```
--var, -v         Variable(s) to process, or "all"
--timestep, -t    Output time aggregation: daily | 10day | monthly
                  (default: TIME_AGGREGATION_1D from config.py)
--spinup-start    Spinup start year for detrending (default: 1940)
--spinup-end      Spinup end year for detrending (default: 1970)
--workers, -w     Number of parallel workers (default: all CPUs)
--output-dir, -o  Output directory (default: from config.py)
--max-points      Limit to N columns (useful for testing)
--list-vars       Print available variables and exit
--quiet, -q       Suppress progress output
```

Variables that are detrended (spinup-aware): `h_surf`, `FirnAir`

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

---

### `make_2d_files.py` — 2D depth profiles (monthly)

Extracts a scalar per column per timestep from the monthly depth-profile files.
Two extraction modes:

```
-o, --output-dir    Directory containing model output files (required)
-v, --var           Input variable: dens | temp | lwc | ... (default: dens)
-t, --threshold     Find depth where variable crosses this value
-d, --depth         Extract value at this depth in metres
--output-var        Name for the output variable (required, e.g. z830, T10m)
--reference-dir     Directory with mask and pointlist files
--processed-dir     Output directory for gridded files
--start-year        Only output timesteps >= this year
--end-year          Only output timesteps <= this year
-n, --num-workers   Number of parallel workers
--dry-run           Print detected configuration and exit
```

Common examples:

```bash
# Critical density horizons
python3 make_2d_files.py -o OUTPUT_DIR -v dens -t 550 --output-var z550
python3 make_2d_files.py -o OUTPUT_DIR -v dens -t 830 --output-var z830
python3 make_2d_files.py -o OUTPUT_DIR -v dens -t 917 --output-var z917

# Temperature at depth
python3 make_2d_files.py -o OUTPUT_DIR -v temp -d 10 --output-var T10m
```

---

### `make_2Ddetail_files.py` — 2Ddetail near-surface profiles (10-day)

Extracts a scalar per column per timestep from the high-resolution near-surface
profiles (4 cm layers, upper ~20 m). Three depth selection modes:

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

Common examples:

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

---

## Configuration reference (`config.py`)

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

> **Note on time aggregation:** `TIME_AGGREGATION_1D` controls resampling and
> can be freely changed. `TIME_AGGREGATION_2D` and `TIME_AGGREGATION_2Ddetail`
> are fixed by the model run — changing them without re-running the model will
> raise an error at startup.

---

## File structure

```
create_1D_2D_2Ddetail_files/
├── config.py                    # Paths, variables, time aggregation settings
├── run_config.py                # Auto-detecting RunConfig (for 2D/2Ddetail)
├── utils.py                     # Shared utilities: resampling, dataset creation
├── make_1d_files.py             # 1D → gridded surface/integrated variables
├── make_2d_files.py             # 2D profiles → gridded depth diagnostics
├── make_2Ddetail_files.py       # 2Ddetail profiles → gridded near-surface fields
├── submit_make_1d_files.sh      # SLURM batch script for 1D processing
├── submit_make_2d_files.sh      # SLURM batch script for 2D processing
├── submit_make_2ddetail_files.sh # SLURM batch script for 2Ddetail processing
└── logs/                        # SLURM job logs (created automatically)
```

---

## Inspect output

```python
import xarray as xr

ds = xr.open_dataset("path/to/output.nc")
print(ds)

# Select by year (time is datetime)
ds_2022 = ds.sel(time="2022")

# Select nearest timestep
da = ds['FirnAir'].sel(time="2022-06", method="nearest")
```
