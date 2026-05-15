# IMAU-FDM

**IMAU Firn Densification Model** is a one-dimensional column model for firn
compaction and melt–refreeze processes over ice sheets and glaciers.

This documentation covers the complete workflow:

| Section | What it covers |
|---------|----------------|
| [Architecture](overview/architecture) | Repository layout and data flow |
| [Physics](overview/physics) | Firn densification equations and key parameters |
| [Pre-Processing](preprocessing/overview) | Converting RACMO/ERA5 output into model input |
| [Running the Model](running/overview) | Domains, settings, SLURM submission |
| [Post-Processing](postprocessing/overview) | Gridded NetCDF output from per-column files |
| [MO Fitting](mo_fitting/overview) | Calibrating densification rate coefficients |
| [Source Code](source/overview) | Fortran module descriptions |
| [Development](development/distributor) | New distributor architecture |

## Quick start on the ECMWF

```bash
# 1. Pre-process forcing (Greenland, ERA5)
cd pre-process-RACMO/
python submit_jobs.py --domain FGRN055 --forcing era5 \
    --start-year 1939 --end-year 2025

# 2. Compile the model
cd ../
./compile_hpc.sh

# 3 Update model settings in rundir/launch_new_job, rundir/start_model_ccab, and source/model_settings

# 4. Launch model run
cd ../rundir/
sbatch submit_job.sc

# 5. Update post-processing config files

# 6. Post-process output 
cd ../post-process/create_1D_2D_2Ddetail_files/
sbatch submit_make_1d_files.sh
sbatch submit_make_2d_files.sh
sbatch submit_make_2ddetail_files.sh
```

## Supported domains

| Domain | Grid | Forcing |
|--------|------|---------|
| FGRN055 | Greenland 5.5 km | ERA5 / RACMO2.3p2 |
| ANT27 | Antarctica 27 km | ERA5 / RACMO2.3p2 |
