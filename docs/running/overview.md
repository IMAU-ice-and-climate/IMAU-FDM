# Running the Model

Runs are launched from `launch_job/`. A single script, `launch_job.sh`, reads
`settings/<DOMAIN>/run.toml`, optionally recompiles, builds the working
directory, and submits the job via the MPI
[distributor](../development/distributor).

## Quick start

```bash
cd launch_job/
# 1. edit settings/<DOMAIN>/run.toml          (see Settings)
# 2. set SETTINGS_FOLDER at the top of launch_job.sh to your domain folder
# 3. put your point list in launch_job/pointlists/
./launch_job.sh
```

`launch_job.sh` will:

1. recompile if `recompile = true` (`compile_hpc.sh` for ECMWF, `fpm` offline);
2. create `output_root/<project_name>/` and its subdirectories;
3. copy the source, settings, executable, and launch scripts into `localcode/`
   — a frozen snapshot for reproducibility;
4. hand the point list and settings to `submit_job.sh`, which submits the job.

## Run types

Set `run_type` in `run.toml`:

| Value | Launch | Parallelism |
|-------|--------|-------------|
| `ECMWF` | `sbatch` an MPI job | `[ecmwf]` settings; workers = min(points, nodes × 128) |
| `offline` | `mpirun` locally | `[offline] n_procs` |

`restart_type` controls spin-up vs. continuation — see [Settings](settings).

## Point lists

The model runs one grid cell ("point") at a time. Reference lists live in
`reference/{DOMAIN}/`; the run-specific list (point indices) lives in
`launch_job/pointlists/` and is named by `pointlist_name` in `run.toml`.

## Output directory structure

```
<project_name>/
├── output/                    # per-point output (1D, 2D, 2Ddetail)
├── post-process/              # gridded output (created during post-processing)
├── restart/
│   ├── spinup/                # restart after spin-up
│   └── run/                   # restart from the end of this run
├── logfiles/
│   ├── model_logfiles/        # per-point logs
│   └── distributor_logfiles/  # distributor logs
├── localcode/                 # frozen snapshot: source/, settings/, executable, launch scripts
└── pointlist_<N>.txt          # point list for submission iteration N
```

## Resuming an unfinished run

If a job doesn't complete, the distributor writes the remaining points to
`pointlist_<N+1>.txt`. With `relaunch = "yes"` it resubmits automatically. To
resume manually, set `submission_iteration = N` in `run.toml` and re-run
`launch_job.sh` — it reuses the existing `pointlist_N.txt` and keeps the
original `localcode` snapshot.

## Monitoring a run

```bash
squeue -u $USER                                 # SLURM queue

cd QAQC/
jupyter notebook plot_ongoing_run.ipynb         # inspect FDM log output
jupyter notebook check_run_is_completed.ipynb   # check completion
```

Per-point logs are in `logfiles/model_logfiles/`; the distributor log is in
`logfiles/distributor_logfiles/`. Check that `remove_*` counts are reasonable,
file paths are correct, and the spin-up converges.
