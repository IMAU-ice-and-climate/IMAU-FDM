# Post-Processing

Converts IMAU-FDM per-column output into gridded NetCDF files.

Source: `post-process/`

## Output filename convention

```
FDM_{variable}_{domain}_{date_range}_{timestep}.nc
```

Examples:
- `FDM_Runoff_FGRN055_1939-2023_10day.nc`
- `FDM_h_surf_FGRN055_1939-2025_10day_detrended.nc`
- `FDM_z830_FGRN055_1939-2023_30day.nc`
- `FDM_T10m_FGRN055_1939-2023_10day.nc`

## Two pipelines

**Part 1 — Create from scratch** (`create_1D_2D_2Ddetail_files/`)

Reads raw per-column output and writes one gridded file per variable.

**Part 2 — Extend existing files** (`extend/`)

Merges per-column output from an extension run with the original, then appends
new timesteps to existing gridded files.

## Output period slicing

Set `OUTPUT_START` and `OUTPUT_END` in
`post-process/create_1D_2D_2Ddetail_files/config.py` to save only a
sub-period (e.g. model covers 1939–2025, output only 1958–2023). Set either to
`None` to output the full period. For 2D and 2Ddetail scripts the
`--start-year`/`--end-year` CLI flags override the config defaults.

## Loading gridded files

```python
from post-process.visualization.visualization import load_gridded_data

ds = load_gridded_data('Runoff', timestep='10day')
ds = load_gridded_data('h_surf', timestep='10day', detrended=True)
ds = load_gridded_data('Runoff', timestep='10day',
                       output_dir=Path('{OUTPUT_PATH}/post-process'),
                       date_tag='1939-2025')
```

## Typical workflow

```bash
# 1. Create gridded files for base run (1939-2023)
cd create_1D_2D_2Ddetail_files/
sbatch submit_make_1d_files.sh
sbatch submit_make_2d_files.sh
sbatch submit_make_2ddetail_files.sh
```

## Typical extension workflow

```bash

#1 Merge per-column output into full-period pointfiles (1939-2025)
cd ../extend/
sbatch submit_extend_pointfiles.sh

# 2a. Remake 1D files for 1939-2025 (full reprocess required)
cd ../create_1D_2D_2Ddetail_files/
sbatch submit_make_1d_files.sh   # point INPUT_DIR at merged files

# 2b. Extend 2D/2Ddetail files (append 2024-2025 only)
cd ../extend/
sbatch submit_extend_variable_2d.sh
sbatch submit_extend_variable_2ddetail.sh
```
