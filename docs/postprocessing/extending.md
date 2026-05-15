# Extending Gridded Files

When the model run is extended to a later end date (e.g. 1939–2023 → 1939–2025),
existing gridded files can be appended rather than reprocessed from scratch.

Source: `post-process/extend/`

## Two-step process

### Step 1 — Extend pointfiles (prerequisite)

Merges per-column 1D, 2D, and 2Ddetail files from the original run and the
extension run into a single file per point covering the full period.

```bash
# ~1 h with 16 workers
sbatch submit_extend_pointfiles.sh

# Reprocess a single point (for testing)
python3 extend_pointfiles.py --start 56928 --end 56928 --workers 1
```

Configure paths in `extend_pointfiles_config.py`:

| Setting | Value |
|---------|-------|
| `ORIG_OUTPUT_DIR` | `scratch/run_FGRN055-era055_1939-2023/output` |
| `EXT_OUTPUT_DIRS` | Extension run output directories (checked in order) |
| `OUTPUT_DIR` | `scratch/FDM_FGRN055_output/output/points` |

### Step 2a — Extend 2D variable files (z830, z550, etc.)

Each timestep in the 2D files is independent, so only the new 2024–2025 slice
of each per-column file is needed.

```bash
python3 extend_variable_2d.py --output-var z830 --var dens --threshold 830 --workers 16
python3 extend_variable_2d.py --output-var z550 --var dens --threshold 550 --workers 16
sbatch submit_extend_variable_2d.sh
```

### Step 2b — Extend 2Ddetail variable files (T10m, SSN, etc.)

```bash
python3 extend_variable_2ddetail.py --output-var T10m --var temp --depth 10 --workers 16
python3 extend_variable_2ddetail.py --output-var SSN --var dens \
    --depth-begin 0 --depth-end 0.5 --workers 16
sbatch submit_extend_variable_2ddetail.sh
```

Configure in `extend_variable_config.py`:

| Setting | Value |
|---------|-------|
| `EXT_INPUT_DIR` | Merged per-column pointfiles |
| `ORIG_OUTPUT_DIR` | Existing 1939–2023 gridded files |
| `EXT_OUTPUT_DIR` | Output directory for 1939–2025 gridded files |
| `ORIG_MODEL_END` | `2023-12-31` |
| `EXT_MODEL_END` | `2025-12-31` |

## Important note on 1D variables

1D gridded files (`Runoff`, `surfmelt`, `h_surf`, etc.) **cannot** be extended
with a simple slice — their resampling windows are anchored to `MODEL_START` and
detrending requires the full post-spinup series. These must be remade from
scratch using `make_1d_files.py` once the merged pointfiles are available.
