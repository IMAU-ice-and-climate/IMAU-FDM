# RACMO Pre-processing for IMAU-FDM

Converts raw RACMO2.3p2 output (decadal NetCDF files) into the timeseries and
spinup-average files expected by IMAU-FDM. Jobs are submitted to SLURM; one job
is submitted per variable per step.

---

## Pipeline overview

```
Raw RACMO decade files
        │
        ▼  make_fdm_years.py
Yearly files  (one .nc per variable per year)
        │
        ▼  make_fdm_timeseries.py
Longitude-band timeseries  (one .nc per variable per lon band)
        │
        ▼  make_fdm_averages.py
Spinup-period averages  (one .nc per variable)
```

All three steps are orchestrated by `submit_jobs.py`, which chains them as
SLURM `afterok` dependencies so they run in order without manual intervention.

---

## Directory structure

```
pre-process-RACMO/
├── submit_jobs.py          # master job submission script (start here)
├── make_fdm_years.py       # step 1: slice decades into yearly files
├── make_fdm_timeseries.py  # step 2: build lon-band timeseries
├── make_fdm_averages.py    # step 3: compute spinup-period averages
│
├── FGRN055_era5/           # domain-specific config (Greenland, ERA5)
│   └── config.py
├── ANT27_era5/             # domain-specific config (Antarctica, ERA5)
│   └── config.py
│
├── QAQC_pre-processing/    # manual quality-check notebook (not part of pipeline)
│   └── check_pre-processing.ipynb
│
└── FGRN055_era5/archive/   # legacy ksh scripts — superseded, kept for reference
```

---

## Prerequisites

- Python 3 with `xarray` and `numpy`
- NCO tools (`ncks`, `ncrcat`, `ncra`) — load via `module load nco`
- SLURM (`sbatch`)
- `$SCRATCH` and `$PERM` environment variables set (used in config.py for paths)

---

## Running the pipeline

Always run from within the domain subdirectory (so `config.py` is importable):

```bash
cd FGRN055_era5          # or ANT27_era5
python3 ../submit_jobs.py all
```

This submits one SLURM job per variable per step (7 variables × 3 steps = 21 jobs),
chained so each step waits for all jobs in the previous step to finish.

### Common usage patterns

```bash
# Full pipeline for all variables
python3 ../submit_jobs.py all

# Preview what would be submitted, without actually submitting
python3 ../submit_jobs.py all --dry-run

# Only the first step (yearly files)
python3 ../submit_jobs.py years

# Steps 1 and 2 only, skip averages
python3 ../submit_jobs.py years timeseries --skip-averages

# Subset of variables
python3 ../submit_jobs.py all --vars precip snowfall tskin
```

---

## Extending an existing timeseries (adding new years)

When new RACMO data arrives and you need to append years to an existing run:

**1. Update `TS_END_YEAR` in `config.py`** to the new target year (e.g. 2025).

**2. Create yearly files for the new years:**
```bash
python3 ../submit_jobs.py years --skip-averages
```

**3a. Safe first pass** — write lon-band parts for the new years only into
`parts_extend/` (non-destructive, existing timeseries files untouched):
```bash
python3 ../submit_jobs.py timeseries \
    --extend-start 2024 --extend-end 2025 \
    --extend-mode parts-only --skip-averages
```

**3b. Append** — once satisfied with the parts, append them to the existing
timeseries files (old files are backed up as `.bak`):
```bash
python3 ../submit_jobs.py timeseries \
    --extend-start 2024 --extend-end 2025 \
    --extend-mode append --skip-averages
```

The averages step does not need to be rerun when extending a timeseries.

---

## Adding a new domain

1. Copy an existing domain directory (e.g. `cp -r FGRN055_era5 MYNEWDOMAIN/`).
2. Edit `config.py` in the new directory:
   - Set `DOMAIN`, `FORCING`, and `PROJECT_NAME`
   - Update `TS_START_YEAR`, `TS_END_YEAR`, `AVE_START_YEAR`, `AVE_END_YEAR`
   - Set `NUM_LONG_BANDS` and `CELL_WIDTH` to match the grid
   - Update `FNAME_SUFFIX` to match your raw RACMO filenames
3. Run from within the new directory as above.

---

## Output file locations

All output is written under `$SCRATCH/<PROJECT_NAME>/`:

| Step | Output path |
|------|-------------|
| Yearly files | `process-RACMO/years-<START>/` |
| Timeseries | `input/timeseries/` |
| Averages | `input/averages/` |
| SLURM job files | `process-RACMO/jobs/` |
| SLURM log files | `process-RACMO/logs/` |

---

## Quality checking

`QAQC_pre-processing/check_pre-processing.ipynb` is a Jupyter notebook for
manually inspecting the output. It is not part of the automated pipeline.
Paths inside the notebook are hardcoded — update them before use.
