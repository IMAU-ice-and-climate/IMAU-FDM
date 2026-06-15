# IMAU-FDM

**IMAU Firn Densification Model** is a one-dimensional column model for firn
compaction and melt–refreeze processes over ice sheets and glaciers.

This documentation covers the complete workflow:

| Section | What it covers |
|---------|----------------|
| [Architecture](overview/architecture) | Repository layout and data flow |
| [Physics](overview/physics) | Firn densification equations and key parameters |
| [Pre-Processing](preprocessing/overview) | Converting RACMO/ERA5 output into model input |
| [Running the Model](running/overview) | Settings, launching, and the MPI distributor |
| [Post-Processing](postprocessing/overview) | Gridded NetCDF output from per-column files |
| [MO Fitting](mo_fitting/overview) | Calibrating densification rate coefficients |
| [Source Code](source/overview) | Fortran module descriptions and compiling |
| [Design Notes](development/overview) | The distributor architecture and run structure |

## Quick start on ECMWF

```bash
# 1. Pre-process forcing (Greenland, ERA5)
cd pre-process-RACMO/
python submit_jobs.py --domain FGRN055 --forcing era5 \
    --start-year 1939 --end-year 2025

# 2. Configure the run
#    Edit settings/FGRN055/run.toml   (project_name, domain, forcing, restart_type, paths, SLURM)
#    Edit settings/FGRN055/model.toml (physics / numerics / output) if needed
#    Put your point list in launch_job/pointlists/

# 3. Launch — compiles first if recompile = true in run.toml
cd ../launch_job/
./launch_job.sh

# 4. Post-process per-column output into gridded NetCDF
cd ../post-process/create_1D_2D_2Ddetail_files/
sbatch submit_make_1d_files.sh
```

The model reads **all** of its configuration from TOML files in
`settings/<DOMAIN>/` and reads its forcing dimensions and dates directly from
the NetCDF forcing metadata — nothing about the run is hard-coded in the Fortran.
See [Running the Model](running/overview).

## Supported domains

| Domain | Grid | Forcing |
|--------|------|---------|
| FGRN055 | Greenland 5.5 km | ERA5 / RACMO2.3p2 |
| ANT27 | Antarctica 27 km | ERA5 / RACMO2.3p2 |
