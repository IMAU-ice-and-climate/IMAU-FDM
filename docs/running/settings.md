# Settings Files

All model configuration lives in `settings/<DOMAIN>/` as TOML files, read at
startup. Changing a value requires **no recompilation**.

For each substantially new run, copy an existing domain folder (e.g.
`settings/FGRN055/`) and edit the copy, so each run's choices are preserved.
There are three files:

| File | Scope | Edit per run? |
|------|-------|---------------|
| `run.toml` | This specific run (paths, project name, SLURM) | **Yes** |
| `model.toml` | Physics, numerics, output dimensions | Only when changing the setup |
| `constants.toml` | Physical constants | Rarely |

## `run.toml`

```toml
[job]
project_name         = "my_run"   # unique; output goes to output_root/project_name
pointlist_name       = "pointlist_1"   # file in launch_job/pointlists/
domain               = "FGRN055"  # FGRN055 | ANT27
forcing              = "era5"     # used in input/output file naming
recompile            = true       # recompile on each launch?
run_type             = "ECMWF"    # ECMWF | offline
restart_type         = "none"     # none | spinup | run
restart_project      = ""         # empty = this project; else read restart from another project
relaunch             = "yes"      # resubmit remaining points if the job doesn't finish
submission_iteration = 1          # 1 for a fresh launch; N to resume from pointlist_N.txt

[metadata]
model_version = "1.4"

[directories]
code_dir     = "/perm/nld4814/code/IMAU-FDM/"
forcing_root = "/home/nld4814/scratch/"   # inputs: <forcing_root>/<domain>_<forcing>/input/{averages,timeseries}
output_root  = "/home/nld4814/scratch/"   # outputs: <output_root>/<project_name>/

[code-subdirectories]
fdm_executable = "imau-fdm"

[ecmwf]                 # used when run_type = "ECMWF"
max_nodes       = 1     # workers = min(n_points, max_nodes * 128)
walltime        = "48:00:00"
account_no      = "spnlberg"
memory_per_task = "999Mb"

[offline]               # used when run_type = "offline"
n_procs = 2             # 1 distributor + (n_procs - 1) workers
```

### Run types (`restart_type`)

| Value | Behaviour |
|-------|-----------|
| `none` | Start from scratch (run the spin-up) |
| `spinup` | Restart from the end of spin-up (begin the main run) |
| `run` | Continue from the end of a previous main run |

When `restart_project` is set, the **initial** restart files are read from that
project; output (including new restart files) is always written to *this* run's
`project_name`.

## `model.toml`

Physics, numerics, and output cadence/dimensions — shared across runs with the
same setup.

```toml
[model_physics]
do_mo_fit = false                 # if true, sets MO factors to 1 for refitting
LWC_avail = "Coleou1998_1p2"      # Coleou1998_corr | Coleou1998_1p2

[initialization]
startasice = 1     # 1=linear density profile, 2=ice
beginT     = 3     # 1=winter, 2=summer, other=isothermal (T=tsav)
initdepth  = 1.0   # initial firn profile depth [m]

[model_choices]
ImpExp     = 1        # 1=implicit (fast), 2=explicit
dtmodelImp = 900      # implicit timestep [s]
dtmodelExp = 180      # explicit timestep [s]
dtSnow     = 31557600 # snow-parameterisation averaging window [s]
dzmax      = 0.15     # vertical resolution [m]
th         = 0.5      # theta (0.5 = Crank-Nicolson)

[output_dimensions]
writeinspeed  = 86400    # speed-component output frequency [s]
writeinprof   = 2592000  # firn-profile output frequency [s] (~30 days)
writeindetail = 864000   # detailed-profile output frequency [s] (~10 days)
proflayers    = 3000     # layers in 2D profile output
detlayers     = 500      # layers in 2Ddetail output
detthick      = 0.04     # 2Ddetail layer thickness [m]
```

## `constants.toml`

Physical constants (gas constant, gravity, ice density, melting point,
activation energies, latent heat, `days_per_year`, fill value). Edit only with
good reason.

## Forcing dimensions are *not* set here

The number of timesteps, forcing years, grid dimensions (`Nlon`/`Nlat`), the
timestep `dtobs`, and the spin-up averaging years are all read from the **forcing
NetCDF metadata** at runtime (`Set_Forcing_Dimensions` in `model_settings.f90`).
Keeping these out of the TOML means the model adapts automatically when the
forcing record is extended.
