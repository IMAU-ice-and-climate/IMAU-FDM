# Repository Architecture

## Directory layout

```
IMAU-FDM/
├── source/             # Fortran source code (compiled to imau-fdm.x)
├── settings/           # TOML configuration files
├── rundir/             # SLURM submission scripts and run setup
├── pre-process-RACMO/  # Convert RACMO/ERA5 NetCDF → model input
├── post-process/       # Convert per-column output → gridded NetCDF
├── MO_fit/             # Calibrate densification coefficients
├── QAQC/               # Run monitoring and log inspection
├── reference/          # Static grid files and point lists (per domain)
├── distributor/        # Job distribution scripts (HPC)
├── example/            # Minimal worked example
└── docs/               # This documentation (Jupyter Book)
```

## End-to-end data flow

```
ERA5 / RACMO2.3p2 (raw decadal NetCDF)
        │
        ▼  pre-process-RACMO/
Yearly files → lon-band timeseries → spinup averages
        │
        ▼  rundir/ + source/
Per-column output: 1D (daily), 2D (monthly profiles), 2Ddetail (10-day near-surface)
        │
        ▼  post-process/
Gridded NetCDF maps  (FDM_{var}_{domain}_{years}_{timestep}.nc)
```

## Key configuration files

| File | Purpose |
|------|---------|
| `settings/model_settings.toml` | Domain, restart type, output layer counts |
| `settings/run_settings.toml` | Submission iteration counter |
| `settings/paths.toml` | All input/output directory paths |
| `settings/ecmwf_settings.toml` | HPC-specific settings |
| `settings/model_variables.toml` | Variables to include in output |

## Reference files (per domain)

Stored in `reference/{DOMAIN}/`:

| File | Contents |
|------|---------|
| `{DOMAIN}_Masks.nc` | Grid mask, lat/lon, x/y in EPSG:3413 |
| `IN_ll_{DOMAIN}.txt` | Point list: lon, lat, ..., rlat\_idx, rlon\_idx |
| `{DOMAIN}_grid.nc` | Rotated-pole grid coordinates |
