# New Workflow (In Development)

Source: `distributor/` (development branch)

<!-- | File | Purpose |
|------|---------|
| `settings/model_settings.toml` | Domain, restart type, output layer counts |
| `settings/run_settings.toml` | Submission iteration counter |
| `settings/paths.toml` | All input/output directory paths |
| `settings/ecmwf_settings.toml` | HPC-specific settings |
| `settings/model_variables.toml` | Variables to include in output | -->

The **distributor** is a new job management layer that replaces the current
SLURM array approach with a more flexible work-queue pattern.

````{admonition} Status
:class: warning
The distributor is under active development on the `development` branch.
This page will be expanded once the design stabilises.
````

## Motivation

The current SLURM array approach (`npnf_outer_script.sc` / `npnf_inner_script.sc`)
submits one job per point and relies on SLURM's array job mechanism for
parallelism. Limitations:

- Failed points must be re-identified and resubmitted manually.
- No dynamic load balancing between fast and slow columns.

## Design

The distributor separates the work queue from the compute workers:

```
launch_job.sh
     │
     ├─ distributor process (1×)   ← hands out point numbers
     │
     └─ worker processes (N×)      ← each pulls next point, runs model
```

Workers request points from the distributor over a socket. If a worker fails,
the distributor re-queues that point. The run finishes when all points are
complete.

## Current state

`distributor/launch_job_archive.sh` — archived version of the old launcher,
kept for reference.

The new distributor implementation is being developed on the `development`
branch. See that branch for current code and tests.
