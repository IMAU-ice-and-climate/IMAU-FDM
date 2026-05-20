# Proposal: Per-point and distributor-level logfiles

## Background

The old system ran each grid point as a separate process launched by a bash script, which redirected that process's stdout to a per-point file:
```
log_IMAU-FDM_<jobname>_<point_numb>.out
```

The new system uses MPI: all worker ranks share stdout, so per-point output must be written explicitly to per-point file units from within Fortran.

## Overview

Two log streams:

| Stream | Content | Written by |
|---|---|---|
| Per-point log | Full model output for a single grid point | Worker rank running that point |
| Distributor log | Job management: points dispatched, ranks, timings, errors | Rank 0 only |

## Implementation

### Step 1: Add `log_unit` to `model_settings.f90`

Add one public variable to the module-level declarations:

```fortran
integer, public :: log_unit = 6   ! 6 = stdout (Fortran default)
```

The default of `6` means the model still works correctly when run outside MPI (e.g. offline/testing) without any log file being opened. All source files already `use model_settings`, so `log_unit` is immediately available everywhere.

### Step 2: Replace `print *` with `write(log_unit, *)` in source files

Files to update (~140 `print` statements total):
- `model_main.f90`
- `model_settings.f90`
- `initialise_model.f90`
- `initialise_variables.f90`
- `openNetCDF.f90`
- `output.f90`
- `time_loop.f90`
- `firn_physics.f90`

Run from `source/`:
```bash
sed -i 's/print \*,/write(log_unit, *)/g; s/print \*$/write(log_unit, *)/g' \
    model_main.f90 model_settings.f90 initialise_model.f90 \
    initialise_variables.f90 openNetCDF.f90 output.f90 \
    time_loop.f90 firn_physics.f90
```

**Exceptions — leave unchanged:**
- `write(stderr, ...)` calls in `model_settings.f90` — error messages should always go to stderr
- `print *` in `distributor.f90` — not per-point model output

### Step 3: Per-point log in `distributor.f90` (worker block)

Before calling `Run_Model`, the worker opens a point-specific log file and sets `log_unit`. After `Run_Model` returns, it closes the file and resets `log_unit` to stdout.

```fortran
! Inside worker do-loop, before calling Run_Model
character(len=512) :: log_fname
integer :: log_unit_local

write(log_fname, '(3A,I0,A)') &
    trim(POINT_LOG_DIR), "/log_IMAU-FDM_", &
    trim(PROJECT_NAME), "_", pointlist(cur_i_point), ".out"

open(newunit=log_unit_local, file=trim(log_fname), &
     action='write', status='replace', iostat=rc)
if (rc /= 0) then
    write(dist_log_unit, *) "ERROR: could not open log file: ", trim(log_fname)
end if

log_unit = log_unit_local   ! set module-level variable

call Run_Model(cur_lat, cur_lon, cur_int_lat, cur_int_lon, &
               point_numb, SETTINGS_PATH)

close(log_unit_local)
log_unit = 6                ! reset to stdout
```

`POINT_LOG_DIR` is passed to the distributor as command-line argument 5 (see below).

### Step 4: Distributor-level log in `distributor.f90` (rank 0)

Rank 0 opens a single log at startup and writes one line per dispatch event.

```fortran
! Rank 0 initialisation
integer :: dist_log_unit
character(len=512) :: dist_log_fname

write(dist_log_fname, '(3A)') &
    trim(DIST_LOG_DIR), "/distributor_", trim(PROJECT_NAME), ".log"

open(newunit=dist_log_unit, file=trim(dist_log_fname), &
     action='write', status='replace')

write(dist_log_unit, '(2A)')   "Distributor started: ", trim(PROJECT_NAME)
write(dist_log_unit, '(A,I0)') "Total points: ", i_point
write(dist_log_unit, '(A,I0)') "MPI size: ", size

! Inside distributor dispatch loop
write(dist_log_unit, '(A,I0,A,I0,A,F8.4,A,F8.4)') &
    "Rank ", status(MPI_Source), " <- point ", pointlist(cur_i_point), &
    "  lat=", new_lat, "  lon=", new_lon

! When all done
write(dist_log_unit, '(A)') "All points dispatched. Distributor exiting."
close(dist_log_unit)
```

## File naming

| File | Path |
|---|---|
| Per-point log | `<WORK_DIR>/model_logfiles/log_IMAU-FDM_<PROJECT_NAME>_<point_numb>.out` |
| Distributor log | `<WORK_DIR>/distributor_logfiles/distributor_<PROJECT_NAME>.log` |

Both directories are already created by `launch_job.sh` (`mkdir -p`).

## Command-line arguments to distributor

| # | Current | Proposed |
|---|---|---|
| 1 | `PROJECT_NAME` | `PROJECT_NAME` |
| 2 | `DOMAIN` | `DOMAIN` |
| 3 | `WORK_POINTLIST_PATH` | `WORK_POINTLIST_PATH` |
| 4 | `SETTINGS_PATH` | `SETTINGS_PATH` |
| 5 | — | `POINT_LOG_DIR` |
| 6 | — | `DIST_LOG_DIR` |

In `launch_job.sh`, the model call becomes:
```bash
"${WORKEXE_DIR}/${FDM_EXECUTABLE}" \
    "$PROJECT_NAME" "$DOMAIN" "$WORK_POINTLIST_PATH" "$SETTINGS_PATH" \
    "$POINT_LOG_DIR" "$DISTRIBUTOR_LOG_DIR"
```
