# The Distributor

Source: `source/distributor.f90` (the model's **main program**) + `launch_job/`.

IMAU-FDM runs many independent firn columns. The **distributor** spreads those
columns across MPI ranks: one rank hands out work, the rest run the model.

## How it works

```
launch_job.sh → submit_job.sh → srun / mpirun -n (workers + 1) imau-fdm
                                          │
        ┌─────────────────────────────────┴──────────────────────────────┐
        │                                                                  │
   rank 0: distributor                                   ranks 1..N: workers
   hands out point numbers               each pulls the next point, runs the
   from the point list                   full single-column model, writes
                                         output, then asks for another
```

- **Rank 0 (distributor)** reads the point list and hands out point indices on
  request.
- **Workers (ranks 1…N)** each run the complete single-column model
  (`model_main`, formerly `main.f90`) for one point, then request the next.
  Because workers pull work on demand, fast columns don't wait for slow ones —
  load balancing is automatic.
- The run finishes when every point has been handed out and completed.

`submit_job.sh` launches this with `srun` on ECMWF (`run_type = "ECMWF"`) or
`mpirun` offline. The worker count comes from `[ecmwf] max_nodes` or
`[offline] n_procs` in `run.toml`.

## Resubmission

When the job ends, `launch_job/make_resubmit_pointlist.py` scans the `output/`
directory, determines which points completed, and writes the remainder to
`pointlist_<N+1>.txt`. If `relaunch = "yes"` and progress was made,
`submit_job.sh` resubmits itself for the next iteration. This replaces the old
manual "find the failed points and resubmit" loop.

## Reproducibility

On launch, the source, settings, executable, and launch scripts are copied into
`<project>/localcode/`, and the submitted job runs from that frozen snapshot. A
run's exact code and configuration are therefore always recoverable, even if the
main repository changes afterwards.
