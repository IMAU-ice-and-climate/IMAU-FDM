# Repository Architecture

## Directory layout

```
IMAU-FDM/
├── source/             # Fortran source (compiled to imau-fdm via fpm)
├── launch_job/         # launch_job.sh + submit_job.sh + MPI distributor launch
│   └── pointlists/     # point lists (which grid cells to run)
├── settings/           # per-domain TOML configuration (settings/<DOMAIN>/)
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
RACMO/ERA5: raw decadal NetCDF of 7 variables
        (10 m wind, precip, snowfall, snowmelt, snowdrift, evaporation, skin temperature)
        │
        ▼  pre-process-RACMO/
Lon-band timeseries files + spin-up averages (date-free names; dims & dates in NetCDF metadata)
        │
        ▼  launch_job/ + source/
Per-column model output: 1D (daily), 2D (monthly profiles), 2Ddetail (10-day near-surface)
        │
        ▼  post-process/
Gridded NetCDF maps of selected variables across the whole domain/timeseries
```

## Key configuration files

All run configuration lives in `settings/<DOMAIN>/` as TOML (read at startup —
no recompilation needed to change a value). See [Settings](../running/settings).

| File | Contents |
|------|---------|
| `settings/<DOMAIN>/run.toml` | Per-run: project name, domain/forcing, restart type, directories, SLURM/offline options |
| `settings/<DOMAIN>/model.toml` | Physics, numerics, and output dimensions (shared across runs with the same setup) |
| `settings/<DOMAIN>/constants.toml` | Physical constants |
| `launch_job/launch_job.sh` | Entry point — reads `run.toml`, builds the working directory, compiles, and submits |

Forcing dimensions, start/end years, and spin-up averaging years are **not** set
here — they are read from the forcing NetCDF metadata at runtime.

## Reference files (per domain)

Stored in `reference/{DOMAIN}/`:

| File | Contents |
|------|---------|
| `{DOMAIN}_Masks.nc` | Grid mask, lat/lon, x/y in EPSG:3413 |
| `IN_ll_{DOMAIN}.txt` | Point list: lon, lat, ..., rlat\_idx, rlon\_idx |
