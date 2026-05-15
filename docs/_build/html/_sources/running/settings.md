# Settings Files

All configuration is in `settings/` as TOML files. The model reads them at
startup; changing a value requires no recompilation.

## `model_settings.toml`

Core model configuration.

```toml
run_location = "ecmwf"     # "ecmwf" | "offline"
restart_type = "none"      # "none" | "spinup" | "run"

domain  = "FGRN055"        # domain name (used in file naming)
forcing = "era055"         # forcing name (used in file naming)

[output_dimensions]
proflayers    = 3000       # layers in 2D profile output (~122 m)
writeindetail = 864000     # seconds between 2Ddetail writes (10 days)
detlayers     = 500        # layers in 2Ddetail output (top 20 m)
detthick      = 0.04       # 2Ddetail layer thickness (m)

[minimum_values]
ts_minimum    = 1.0e-04    # floor for 1D timeseries values
det2d_minimum = 1.0e-05    # floor for 2Ddetail refreezing
```

## `paths.toml`

All file and directory paths. Edit this when setting up a new run.

```toml
[directories]
project_name = ""          # unique run name; output goes to output_dir/project_name
code_dir     = "/perm/nld4814/code/IMAU-FDM/"
input_dir    = "/scratch/FGRN055_era055/input/"
output_dir   = "/scratch/"

[input-subdirectories]
averages   = "averages/"
timeseries = "timeseries/"

[output-subdirectories]
output  = "output/"
settings = "ms_files/"
restart = "restart/"

[reference-files]
mask          = "FGRN055_Masks.nc"
pointlist_ref = "IN_ll_FGRN055.txt"
```

## `ecmwf_settings.toml`

ECMWF HPC-specific settings (module names, scratch paths). See
[HPC Setup](hpc) for details.

## `model_variables.toml`

Selects which variables appear in the 1D output files. Set a variable to
`false` to suppress it (reduces output size).

## `run_settings.toml`

Tracks the current submission iteration (`submission_iteration`). The SLURM
scripts increment this automatically; rarely need to edit by hand.
