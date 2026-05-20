# Proposal: Job resubmission and submission_iteration

## How the old system works

The old system has three layers:

1. **`launch_new_job.sc`** — sets all paths, calls `submit_job.sc` which writes a
   SLURM script and submits it with `sbatch`.

2. **`npnf_outer_script.sc`** (runs inside SLURM) — launches the external
   `distribute_points.x` distributor as a background process, then uses `srun`
   to run worker tasks in parallel. After all workers finish, it polls for a
   signal file at `$requestdir/DP`. If the signal is `"continue"` and
   `relaunch != "no"`, it increments `submission_iteration` and calls
   `submit_job.sc` again, which submits a new SLURM job with the remaining
   points.

3. **`start_model_ccab.sc`** (runs inside each worker task) — runs the model,
   then writes a completion file to `$readydir/$cpoint`.

The distributor (`distribute_points.x`) monitors elapsed time relative to
walltime. When approaching the end, it stops dispatching new points, waits for
in-flight workers to finish, writes a new pointlist of the remaining points,
then signals `"continue"` or `"stop"` via `$requestdir/DP`.

**Key mechanisms:**
- Each completed point writes a file to `readydir/`
- Distributor writes remaining pointlist and signals outer script
- Outer script resubmits from within the running SLURM job


## Challenges in the new architecture

The new system collapses all three layers into one MPI program
(`distributor.f90`). This creates two problems:

1. **No outer wrapper** to detect the "continue" signal and resubmit — the MPI
   job itself is the whole job. When SLURM kills it at walltime, there is
   nothing left running to resubmit.

2. **SLURM walltime kill is unclean** — if the job runs out of walltime, all
   MPI ranks are killed immediately with no graceful shutdown. The distributor
   cannot write a remaining pointlist at that moment.

The solution is to separate "tracking completion" from "deciding to resubmit":
track completion continuously during the run (so the state survives a kill),
and resubmit via a SLURM job dependency that runs after the main job ends.


## Proposed architecture

### 1. Track completed points with done-files (in `distributor.f90`)

Each worker rank writes a small file to `READY_DIR` immediately after
`Run_Model` returns, one file per completed point:

```fortran
! After Run_Model returns, worker writes done-file
character(len=512) :: done_fname
write(done_fname, '(2A,I0)') trim(READY_DIR), "/", pointlist(cur_i_point)
open(newunit=done_unit, file=trim(done_fname), action='write', status='replace')
write(done_unit, '(A,I0,A)') "done"
close(done_unit)
```

This mirrors the old `echo ... > $readydir/$cpoint` in `start_model_ccab.sc`.
Since it is written immediately on completion, it survives a walltime kill.

### 2. Relaunch via SLURM job dependency (in `launch_job.sh`)

Instead of resubmitting from within the running job, submit a lightweight
relaunch script at the same time as the main job, using
`--dependency=afterany:$MAIN_JOBID`. SLURM will run it after the main job
ends for any reason (completion, walltime, failure).

```bash
# Submit main MPI job
MAIN_JOBID=$(sbatch --parsable main_job.sc)

# Submit relaunch script, runs after main job finishes
sbatch --dependency=afterany:${MAIN_JOBID} relaunch.sc
```

`relaunch.sh` (a new script in `launch_job/`):
1. Reads `relaunch` from `run.toml`
2. Reads the original pointlist for this iteration
3. Checks which points have done-files in `READY_DIR`
4. If remaining points > 0 and `relaunch = true`:
   - Writes a new pointlist `pointlist_<next_iteration>.txt` with remaining points
   - Increments `submission_iteration` (stored as a file in `WORK_DIR`, not
     in `run.toml` — avoids modifying shared config mid-run)
   - Resubmits the main job + a new relaunch dependency

```bash
# relaunch.sh sketch
ITER=$(cat "${WORK_DIR}/submission_iteration.txt")
RELAUNCH=$(toml_get "$SETTINGS_PATH/run.toml" job relaunch)

# Build remaining pointlist by diffing original vs done-files
python3 -c "
import os, sys
orig = [int(l.strip()) for l in open(sys.argv[1])]
done = set(int(f) for f in os.listdir(sys.argv[2]) if f.isdigit())
remaining = [p for p in orig if p not in done]
print(f'{len(remaining)} points remaining')
for p in remaining:
    print(p, file=open(sys.argv[3], 'a'))
" "${WORK_POINTLIST_PATH}" "${READY_DIR}" "${NEXT_POINTLIST}"

N_REMAINING=$(wc -l < "${NEXT_POINTLIST}")

if [[ $N_REMAINING -gt 0 && "$RELAUNCH" == "True" ]]; then
    NEXT_ITER=$((ITER + 1))
    echo $NEXT_ITER > "${WORK_DIR}/submission_iteration.txt"
    # resubmit main job and new relaunch dependency
    NEW_JOBID=$(sbatch --parsable main_job.sc)
    sbatch --dependency=afterany:${NEW_JOBID} relaunch.sc
else
    echo "All points done, or relaunch=false. Stopping."
fi
```

### 3. Add `relaunch` to `run.toml`

```toml
[job]
relaunch           = true    # automatically resubmit if job hits walltime
```

### 4. Track `submission_iteration` in `WORK_DIR`

`submission_iteration` in `run.toml` is the *starting* value (always 1 for a
new run). The *current* value lives in `WORK_DIR/submission_iteration.txt`,
written and incremented by `relaunch.sh`. This avoids modifying the shared
TOML config between iterations.

`launch_job.sh` initialises it:
```bash
echo "1" > "${WORK_DIR}/submission_iteration.txt"
```

### 5. Per-iteration pointlist naming

Pointlists in `WORK_DIR`:
```
pointlist_1.txt   <- copy of original on first launch
pointlist_2.txt   <- remaining points after iteration 1, written by relaunch.sh
pointlist_3.txt   <- etc.
```

The active pointlist is always `pointlist_<current_iter>.txt`.  
`distributor.f90` reads the path from its command-line arguments (set by
`launch_job.sh` to the correct iteration's file).


## Summary of files changed or added

| File | Change |
|---|---|
| `source/distributor.f90` | Worker writes done-file after each `Run_Model` call |
| `launch_job/launch_job.sh` | Initialise `submission_iteration.txt`; submit main job + relaunch dependency |
| `launch_job/relaunch.sh` | New script: diffs pointlist vs done-files, writes remaining list, resubmits |
| `settings/default/run.toml` | Add `relaunch = true` field |

The `distributor.f90` itself needs no walltime-awareness — all resubmission
logic lives in bash, which is easier to modify without recompiling.
