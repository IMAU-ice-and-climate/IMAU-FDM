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

## Quick start on the ECMWF

```bash
# 1. Pre-process forcing (Greenland, ERA5)
#    Edit FGRN055_era5/config.py first, then run from inside the domain dir:
cd pre-process-RACMO/FGRN055_era5
python3 ../submit_jobs.py all

# 2. Compile the model
cd ../..
./compile_hpc.sh

# 3. Edit run settings in rundir/launch_new_job.sc, rundir/start_model_ccab.sc,
#    and source/model_settings.f90  (see "Running the Model")

# 4. Launch the run (interactive: echoes all paths and prompts before submitting)
cd rundir
./launch_new_job.sc

# 5. Edit post-processing config (create_1D_2D_2Ddetail_files/config.py)

# 6. Post-process output into gridded NetCDF
cd ../post-process/create_1D_2D_2Ddetail_files
sbatch submit_make_1d_files.sh
sbatch submit_make_2d_files.sh
sbatch submit_make_2ddetail_files.sh
```

## Supported domains

| Domain | Grid | Forcing |
|--------|------|---------|
| FGRN055 | Greenland 5.5 km | ERA5 / RACMO2.3p2 |
| ANT27 | Antarctica 27 km | ERA5 / RACMO2.3p2 |
