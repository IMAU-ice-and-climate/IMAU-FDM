# HPC Setup (ECMWF)

## Environment

The model runs on ECMWF's HPC (Atos/Bologna). Key environment variables:

| Variable | Path |
|----------|------|
| `$PERM` | `/home/nld4814/perm` — persistent storage (code, reference files) |
| `$SCRATCH` | `/home/nld4814/scratch` — fast scratch (model output) |

Load required modules before compiling or running:

```bash
module load netcdf4-fortran
module load fpm              # if not installed locally
```

## Compiling

```bash
cd /perm/nld4814/code/IMAU-FDM/
./compile_hpc.sh             # release build (optimised)
./compile_hpc.sh debug       # debug build
```

This runs `fpm install` and places the binary at `~/.local/bin/imau-fdm`.
The `rundir/` scripts expect `imau-fdm.x` in the repo root — the compile
script copies it there.

## Submitting a run

```bash
cd rundir/
sbatch submit_job.sc
```

`submit_job.sc` reads `settings/paths.toml` and `settings/model_settings.toml`,
constructs the point-list range, and submits an array job via `npnf_outer_script.sc`.

## Typical resource requirements

| Job | Nodes | CPUs | RAM | Wall time |
|-----|-------|------|-----|-----------|
| Full FGRN055 run | ~10 | 128/node | 256 GB | 48 h |
| Pre-processing | 1 | 8 | 32 GB | 6 h |
| Post-processing 1D | 1 | 16 | 64 GB | 48 h |
| Post-processing 2D | 1 | 8 | 32 GB | 24 h |

## Checking job status

```bash
squeue -u $USER               # all running jobs
sacct -j <JOBID> --format=JobID,State,ExitCode,Elapsed
```

Log files are written to `rundir/logs/ecsbatch.log.*`.

## Tape archive

Raw RACMO/ERA5 forcing data is on ECMWF tape (MARS).
See `misc/copy_files_from_tape_FGRN055.sc` for the retrieval script.
