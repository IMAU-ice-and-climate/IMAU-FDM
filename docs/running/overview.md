# Running the Model

Source: `rundir/`

## Quick start

```bash
cd rundir/

# Edit settings first (see Settings page)
# Then submit:
sbatch launch_new_job.sc
```

## Run types

Set `restart_type` in `settings/model_settings.toml`:

| Value | Behaviour |
|-------|-----------|
| `none` | Start from scratch (spinup) |
| `spinup` | Restart from end of spinup (begin main run) |
| `run` | Continue from end of last main run - ensure restart directory is set correctly|

## Point lists

The model runs one column per MPI rank. Reference point lists live in `reference/{DOMAIN}/`. A customized pointlist (list of indices) specifies which grid
points to run. These live in `rundir/pointlists/`.

## SLURM scripts

| Script | Purpose |
|--------|---------|
| `launch_new_job.sc` | Main script; launch a fresh run |
| `submit_job.sc` | Starts job |
| `npnf_outer_script.sc` | Outer SLURM array wrapper |
| `npnf_inner_script.sc` | Inner per-point launcher |
| `CancelMyJob.sc` | Cancel all running jobs - updated when new job is launched|

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
