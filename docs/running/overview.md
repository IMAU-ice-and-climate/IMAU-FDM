# Running the Model

Source: `rundir/`

## Quick start

```bash
cd rundir/

# Edit settings first (see Settings page)
# Then submit:
sbatch submit_job.sc
```

## Run types

Set `restart_type` in `settings/model_settings.toml`:

| Value | Behaviour |
|-------|-----------|
| `none` | Start from scratch (spinup) |
| `spinup` | Restart from end of spinup (begin main run) |
| `run` | Continue from end of last main run |

## Point lists

The model runs one column per MPI rank. A point list specifies which grid
points to run. Reference point lists live in `reference/{DOMAIN}/`:

- `IN_ll_FGRN055.txt` — all 58 265 ice-covered Greenland points

Custom subsets live in `rundir/pointlists/`. To generate a new subset:

```bash
sbatch rundir/make_pointlist.sc
```

Point list format (space-separated, 7 columns):

```
lon  lat  ??  ??  ??  rlat_idx  rlon_idx
```

Rows are 1-based; row N corresponds to output file `{prefix}_1D_N.nc`.

## SLURM scripts

| Script | Purpose |
|--------|---------|
| `submit_job.sc` | Main submission script |
| `launch_new_job.sc` | Launch a fresh run |
| `launch_example_job.sc` | Run the minimal example |
| `npnf_outer_script.sc` | Outer SLURM array wrapper |
| `npnf_inner_script.sc` | Inner per-point launcher |
| `ns_script.sc` | Node-sharing variant |
| `offline.sc` | Run without HPC scheduler |
| `CancelMyJob.sc` | Cancel all running jobs |

## Monitoring a run

```bash
# Check SLURM queue
squeue -u $USER

# Inspect FDM log output
cd QAQC/
jupyter notebook plot_ongoing_run.ipynb

# Check completion
jupyter notebook check_run_is_completed.ipynb
```

Log files are written to `rundir/logs/`.
