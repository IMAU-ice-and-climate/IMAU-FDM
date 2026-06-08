# Source Code

The model is written in Fortran 90 and compiled with `gfortran` / `mpifort`.
Source files live in `source/`; compiled objects in `objects/` and module
interfaces in `modules/`.

## Module overview

| Module/file  | Description |
|--------|-------------|
| `model_settings`  | Reads and stores all configurations (pathnames, public constants, etc) |
| `initialise_variables` | Allocates and zeroes all model arrays |
| `initialise_model` |  Reads forcing files, sets up grid |
| `grid_routines` |  Layer insertion/removal, mass redistribution |
| `firn_physics` |  Densification, temperature diffusion |
| `water_physics` |  Melt, percolation, refreezing, liquid water transport |
| `time_loop` |  Spinup and model run |
| `output` |  Writes 1D, 2D, 2Ddetail output NetCDF files |
| `openNetCDF` | NetCDF read/write wrappers |
| `model_main` | Per-point single-column driver (formerly `main.f90`) |
| `distributor` | MPI main program — hands points out to worker ranks |

## Compile
(see [compiling](compiling.md) for more details)

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

Layer 0 is the bottom of the column; the surface layer has the
highest index.

## Adding a new output variable

1. Declare the variable in `initialise_variables.f90` and across any other modules.
2. Compute it inside `time_loop.f90` (or the relevant physics module).
3. Register it for output in `output.f90`.
