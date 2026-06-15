# HPC Setup (ECMWF)

## Environment

The model normally runs on ECMWF's HPC (Atos/Bologna), though you can also run it offline. Key environment variables:

| Variable | Path |
|----------|------|
| `$PERM` | `/perm/$USER` — persistent storage (code, reference files) |
| `$SCRATCH` | `/ec/res4/scratch/$USER` — fast scratch (model output) |

You don't normally need to load modules by hand — `compile_hpc.sh` loads what
it needs (`prgenv/gnu`, `hpcx-openmpi`, `netcdf4`, `nco`, `cdo`, …) and installs
`fpm` from conda if it isn't already on your path.

## Compiling

```bash
cd $PERM/code/IMAU-FDM/
./compile_hpc.sh             # release build (optimised)
./compile_hpc.sh debug       # debug build
```

This runs `fpm install` and places the binary at `~/.local/bin/imau-fdm`.
The `rundir/` scripts expect `imau-fdm.x` in the repo root — the compile
script copies it there.

## Submitting a run

```bash
cd rundir/
./launch_new_job.sc
```

`launch_new_job.sc` is where you set the run parameters (project name, domain,
forcing, restart type, point list, SLURM options). It is interactive: it echoes
all resolved paths and prompts for confirmation, then calls `submit_job.sc`,
which constructs the point-list range and submits an array job via
`npnf_outer_script.sc`.

## Typical resource requirements

| Job | Nodes | CPUs | RAM | Wall time |
|-----|-------|------|-----|-----------|
| Full FGRN055 run | 8 | 128/node | — | 48 h |
| Pre-processing | 1 | 8 | 32 GB | 6 h |
| Post-processing 1D | 1 | 16 | 64 GB | 48 h |
| Post-processing 2D | 1 | 8 | 32 GB | 24 h |

## Checking job status

```bash
squeue -u $USER               # all running jobs
sacct -j <JOBID> --format=JobID,State,ExitCode,Elapsed
```

Run logs are written under the run's scratch directory:
`$SCRATCH/<project_name>/nplogs/` (SLURM/distributor logs) and
`$SCRATCH/<project_name>/logfiles/` (per-point model logs).

## Tape archive

Raw RACMO/ERA5 forcing data is on ECMWF tape (MARS).
See `misc/copy_files_from_tape_FGRN055.sc` for the retrieval script.
