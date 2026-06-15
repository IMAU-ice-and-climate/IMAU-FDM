# IMAU-FDM

**IMAU Firn Densification Model** — a one-dimensional firn compaction and
melt–refreeze model for ice sheets and glaciers. The `main` branch always holds
a working version.

📖 **Documentation:** <https://imau-ice-and-climate.github.io/IMAU-FDM/> (or build
the [`docs/`](./docs) Jupyter Book locally). The end-to-end workflow is
pre-processing → running → post-processing.

# Repository layout

| Directory | Contents |
|-----------|----------|
| `source/` | Fortran model source (compiled to `imau-fdm.x`) |
| `rundir/` | Run setup and SLURM submission scripts |
| `pre-process-RACMO/` | Build model forcing from RACMO/ERA5 input |
| `post-process/` | Gridded NetCDF output from per-column files |
| `MO_fit/` | Calibrate densification coefficients against firn cores |
| `QAQC/` | Run monitoring and log inspection |
| `reference/` | Static grid files and point lists (per domain) |
| `docs/` | Documentation (Jupyter Book) |
 
# Development

`mpifort` and `gfortran` are used for compilation.

## Prerequisites

- [`mpifort`](https://www.mpich.org/downloads/) must be available on your path
- `gfortran` must be available on your path
- `netcdf4-fortran` must be installed
  - `libnetcdf` is a dependency, install it first
  -  to check if `netcdf` is installed properly, try: `nf-config --all`
  - `netcdf4-fortran` can then be built from source

## Compilation

### Using `fpm`

The easiest way to compile this project is using [fpm](https://github.com/fortran-lang/fpm):

1. Make sure `fpm` is installed and on your path
2. From this repository's root directory, run: `fpm install`
3. This will compile and move the `imau-fdm` binary to `~/.local/bin/`
4. To run, execute `~/.local/bin/imau-fdm` (alternatively, add `~/.local/bin/` to your `PATH` and just run `imau-fdm`)

See the [fpm documentation](https://fpm.fortran-lang.org/) for more info, e.g. on cleaning the build.

#### On HPC

For running on ECMWF, a convenience script `compile_hpc.sh` is provided. Simply run it with `./compile_hpc.sh`. This will:

1. Install `fpm` from `conda`, if it's not already found in `~/.local/bin/`
2. Activate the right modules required for compilation.
3. Clean all existing build artifacts.
4. Set some extra flags (e.g. where to find `netcdf4`).
5. Run `fpm install`

By default, `compile_hpc.sh` will build the project using `fpm`'s `--profile release` option, which optimizes for production. Run `./compile_hpc.sh debug` to use the debug profile instead.

#### Locally

You can run `fpm build` or `fpm install`. `fpm` will automatically download certain dependencies, and attempt to find certain others (e.g. `netcdf`). To ensure that you have all required dependencies, you can use the provided `conda` `environment.yml`:

1. `conda env create -f environment.yml`
2. `conda activate imau-fdm` (this will also install `fpm` automatically)
3. `fpm build`

