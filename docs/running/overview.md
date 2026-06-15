# Running the Model

Always re-compile the model first. See [compiling](../source/compiling).

Source: `rundir/`

## Quick start on ECMWF

```bash
cd rundir/

# Edit settings first (see the Settings Files page)
# Then launch. launch_new_job.sc is interactive: it echoes all paths and
# prompts for confirmation, then calls submit_job.sc to submit the array job.
./launch_new_job.sc
```

## Run types

Set `restart_type` in `rundir/launch_new_job.sc`:

| Value | Behaviour |
|-------|-----------|
| `none` | Start from scratch (spinup) |
| `spinup` | Restart from end of spinup (begin main run) |
| `run` | Continue from end of last main run - ensure restart directory is set correctly|

## Point lists

The model runs one column/point/cell at a time. Reference point lists live in `reference/{DOMAIN}/`. A customized pointlist (list of indices) specifies which grid points to run. These live in `rundir/pointlists/`.

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
chk
#or 
squeue -u $USER

# see info about ongoing run or whether it actually succeeded
cd QAQC/

#Inspect FDM log output
jupyter notebook plot_ongoing_run.ipynb

# Check completion
jupyter notebook check_run_is_completed.ipynb
```

You can also look at logs in the output directory. There are multiple types of logs:

1. `output/logfiles` - individual logfiles for each point; check remove_ values are resasonable (~<100), filepaths are correct, spinup converges, etc
2. `nplogs/*_np.log` - log about the SLURM job itself
3. `nplogs/*_runlog.log` - log about the distributor
4. `nplogs/Threads_iter_#` - log about each core