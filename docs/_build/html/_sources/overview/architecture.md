# Repository Architecture

## Directory layout

```
IMAU-FDM/
├── source/             # Fortran source code (compiled to imau-fdm.x)
├── rundir/             # SLURM submission scripts and run setup
├── pre-process-RACMO/  # Convert RACMO/ERA5 NetCDF → model input
├── post-process/       # Convert per-column output → gridded NetCDF
├── MO_fit/             # Calibrate densification coefficients
├── QAQC/               # Run monitoring and log inspection
├── reference/          # Static grid files and point lists (per domain)
├── example/            # Minimal worked example
└── docs/               # This documentation (Jupyter Book)
```

## End-to-end data flow

```
RACMO: raw decadal NetCDF files of 7 variables (10m wind speed, precipitation, snowfall, snowmelt, snowdrift, evaporation, skin temperature)
        │
        ▼  pre-process-RACMO/
Yearly files → lon-band timeseries → spinup averages
        │
        ▼  rundir/ + source/
Model output: 1D (daily), 2D (monthly profiles), 2Ddetail (10-day near-surface) for each cell
        │
        ▼  post-process/
Gridded NetCDF maps: specified variables across whole domain/timeseries (e.g., firn air content, integrated liquid water content, etc)
```

## Key configuration files

| File | Contents |
|------|---------|
|`rundir/launch_new_job.sc`| Set project_name, domain, pointlist, restart_type, SLURM options|
|`rundir/start_model_ccab.sc`| Sets various model parameters, including timestep & dimensions of input and output|
|`source/model_settings.sc`| Constants, model metadata, and paths are set in this file|

_Note that once the new [distributor](development/distributor) comes online, this structure will change*

## Reference files (per domain)

Stored in `reference/{DOMAIN}/`:

| File | Contents |
|------|---------|
| `{DOMAIN}_Masks.nc` | Grid mask, lat/lon, x/y in EPSG:3413 |
| `IN_ll_{DOMAIN}.txt` | Point list: lon, lat, ..., rlat\_idx, rlon\_idx |