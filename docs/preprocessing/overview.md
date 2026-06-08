# Pre-Processing

Converts RACMO output (decadal NetCDF files) into the
timeseries and spinup-average files expected by IMAU-FDM.

Source: `pre-process-RACMO/`

## Pipeline

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

## Running the pipeline

```bash
cd pre-process-RACMO/

# Submit all three steps for Greenland ERA5
python submit_jobs.py --domain FGRN055_era5 --start-year 1939 --end-year 2023

# Extend an existing run (step 1 and 2 only, for adding new years if new data becomes available)
python submit_jobs.py --domain FGRN055_era5 --start-year 2024 --end-year 2025 \
    --skip-averages
```

See `pre-process-RACMO/README.md` for full argument reference.

## Domain configuration

Each domain has a `config.py` in its own subdirectory:

```
pre-process-RACMO/
├── FGRN055_era5/config.py   # Greenland, ERA5
└── ANT27_era5/config.py     # Antarctica, ERA5
```

Key settings in `config.py`:

| Setting | Description |
|---------|-------------|
| `DOMAIN` | Domain name (e.g. `FGRN055`) |
| `FORCING` | Forcing name (e.g. `era5`) |
| `INPUT_DIR` | Raw RACMO/ERA5 files |
| `OUTPUT_DIR` | Pre-processed output |
| `SPINUP_START` / `SPINUP_END` | Years used for spinup averages |
| `NUM_LONG_BANDS` | Number of longitude bands to split into |

## Output structure

```
PATH_TO_OUTPUT/{domain}_{forcing}/input/
├── timeseries/           # lon-band timeseries (used during model run)
│   ├── {domain}_{forcing}_{start_year}-{end_year}_p1.nc
│   └── ...
└── averages/             # spinup-period averages (used during spinup)
    └── {domain}_{forcing}_{spinup_start_year}-{spinup_end_year}_ave.nc
```

## QA notebook

`pre-process-RACMO/QAQC_pre-processing/check_pre-processing.ipynb` checks:

1. **File completeness** — all expected yearly, timeseries, and average files exist
2. **Time axis integrity** — no gaps or duplicates in timeseries
3. **NaN check** — no unexpected NaNs in output variables
4. **Spatial coverage** — all land-ice grid points have valid forcing data

````{admonition} Walkthrough notebook
:class: tip
To embed the QA notebook here, copy it to `docs/notebooks/preprocessing_qaqc.ipynb`
and add `- file: notebooks/preprocessing_qaqc` to `_toc.yml`.
````
