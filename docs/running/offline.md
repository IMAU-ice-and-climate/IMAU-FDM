# Running Offline

You can run IMAU-FDM on a workstation (or any machine with MPI + NetCDF) instead
of the ECMWF HPC. Set `run_type = "offline"` in `settings/<DOMAIN>/run.toml`.

## Prerequisites

- `mpifort` / `mpirun` (MPICH or OpenMPI) on your `PATH`
- `gfortran`
- `netcdf4-fortran` (see [Compiling](../source/compiling))
- `fpm`

The provided conda environment installs these:

```bash
conda env create -f environment.yml
conda activate imau-fdm
```

## Configure

```toml
# settings/<DOMAIN>/run.toml
run_type = "offline"

[offline]
n_procs = 2      # 1 distributor + (n_procs - 1) workers
```

Point `forcing_root` and `output_root` at local paths.

## Launch

```bash
cd launch_job/
./launch_job.sh
```

With `run_type = "offline"`, `launch_job.sh` recompiles with `fpm` (when
`recompile = true`) and `submit_job.sh` runs the job directly with
`mpirun -n <n_procs>` instead of `sbatch`. Everything else — the working
directory layout, the [distributor](../development/distributor), restart
handling, and resubmission — behaves exactly as on the HPC.
