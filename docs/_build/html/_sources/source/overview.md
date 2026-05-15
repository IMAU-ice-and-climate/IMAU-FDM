# Source Code

The model is written in Fortran 90 and compiled with `gfortran` / `mpifort`.
Source files live in `source/`; compiled objects in `objects/` and module
interfaces in `modules/`.

## Module overview

| Module | File | Description |
|--------|------|-------------|
| `model_settings` | `model_settings.f90` | Reads and stores all configuration (TOML files, namelist) |
| `initialise_variables` | `initialise_variables.f90` | Allocates and zeroes all model arrays |
| `initialise_model` | `initialise_model.f90` | Reads forcing files, sets up grid, applies spinup |
| `grid_routines` | `grid_routines.f90` | Layer insertion/removal, mass redistribution |
| `firn_physics` | `firn_physics.f90` | Densification, temperature diffusion, compaction |
| `water_physics` | `water_physics.f90` | Melt, percolation, refreezing, liquid water transport |
| `time_loop` | `time_loop.f90` | Outer time integration loop |
| `output` | `output.f90` | Writes 1D, 2D, 2Ddetail output NetCDF files |
| `openNetCDF` | `openNetCDF.f90` | NetCDF read/write wrappers |
| `main` | `main.f90` | Entry point; MPI initialisation |

## Build system

The model is built with [`fpm`](https://fpm.fortran-lang.org/):

```bash
fpm build --profile release    # optimised
fpm build --profile debug      # with bounds-checking
fpm install                    # installs binary to ~/.local/bin/imau-fdm
```

On ECMWF HPC, use the provided wrapper:

```bash
./compile_hpc.sh               # release
./compile_hpc.sh debug         # debug
```

## Key arrays

| Array | Shape | Description |
|-------|-------|-------------|
| `dens` | `(proflayers)` | Density profile (kg m⁻³) |
| `temp` | `(proflayers)` | Temperature profile (K) |
| `depth` | `(proflayers)` | Layer mid-point depth (m) |
| `dz` | `(proflayers)` | Layer thickness (m) |
| `lwc` | `(proflayers)` | Liquid water content (m w.e.) |
| `dRho` | `(proflayers)` | Density change per timestep |

Layer 0 is the bottom of the column (~55–122 m); the surface layer has the
highest index. New snow is inserted at the top each timestep.

## Adding a new output variable

1. Declare the variable in `initialise_variables.f90`.
2. Compute it inside `time_loop.f90` (or the relevant physics module).
3. Register it for output in `output.f90` and add it to `settings/model_variables.toml`.
