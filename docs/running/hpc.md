# HPC Setup (ECMWF)

The model normally runs on ECMWF's HPC (Atos/Bologna); see [offline](offline) to
run on your own machine.

## Environment

| Variable | Path |
|----------|------|
| `$PERM` | `/home/$USER/perm` — persistent storage (code, reference files) |
| `$SCRATCH` | `/home/$USER/scratch` — fast scratch (forcing input, model output) |

`compile_hpc.sh` loads the required modules itself (`prgenv/gnu`,
`hpcx-openmpi`, `netcdf4`, `nco`, `cdo`, `python3`, …), so you normally don't
need to load anything by hand before compiling.

## Compiling

```bash
cd /perm/nld4814/code/IMAU-FDM/
./compile_hpc.sh             # release build (optimised)
./compile_hpc.sh debug       # debug build
```

This runs `fpm install --prefix <repo>` and places the `imau-fdm` executable in
the repo's `bin/`. With `recompile = true` in `run.toml`, `launch_job.sh` runs
this for you on each launch.

## Submitting a run

```bash
cd launch_job/
./launch_job.sh
```

`launch_job.sh` reads `run.toml`, sets up the working directory, and calls
`submit_job.sh`, which `sbatch`es an MPI job. The job runs the
[distributor](../development/distributor): one rank hands out points, the rest
run the model.

### SLURM sizing (`[ecmwf]` in `run.toml`)

| Key | Meaning |
|-----|---------|
| `max_nodes` | Node budget; workers = `min(n_points, max_nodes × 128)` |
| `walltime` | Job wall-clock limit |
| `account_no` | Charge account |
| `memory_per_task` | `--mem-per-cpu` for each rank |

The number of nodes is derived automatically from the worker count (128
cores/node, +1 rank for the distributor).

## Checking job status

```bash
squeue -u $USER
sacct -j <JOBID> --format=JobID,State,ExitCode,Elapsed
```

The distributor's SLURM log is `<project>/launch_<project>_<iteration>.log`;
per-point and distributor logs are under `<project>/logfiles/`.

## Tape archive

Raw RACMO/ERA5 forcing is on ECMWF tape (MARS). See
`misc/copy_files_from_tape_FGRN055.sc` for the retrieval script.
