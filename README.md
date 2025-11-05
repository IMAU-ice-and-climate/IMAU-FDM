####This repo is under construction####

# IMAU-FDM
Working repository for the IMAU-FDM code. The main branch will always have a working version. 

# Code structure
- main: was @mbrils' cleanup-code-working; now main
- v1.2.0: archived old code structure; also perserved as release 1.2.0G
- ehc: elizabeth case current working branch
- grain-size: based on main, currently being tested
- documentation: documentation being developed on this branch
 
## Main branch (2024-07-24)

![rundir](https://github.com/user-attachments/assets/7733a84f-9a20-484d-8332-ab31489426bb)

![sourcedir](https://github.com/user-attachments/assets/4c4e2eb3-fdf9-4282-9565-c88cebc1050d)


# Development

`mpifort` and `gfortran` are used for compilation.

## Prerequisites

- [`mpifort`](https://www.mpich.org/downloads/) must be available on your path
- `gfortran` must be available on your path
- `netcdf4-fortran` must be installed
  - `libnetcdf` is a dependency, install it first
  -  to check if `netcdf` is installed properly, try: `nf-config --all`
  - `netcdf4-fortran` can then be built from source

See the example [Containerfile](./Containerfile) for how to install `libnetcdf` and `netcdf4-fortran` on Ubuntu.

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

### Manual

The project can be compiled using the [`Makefile`](./Makefile). We need to set at least two environment variables:

- `NETCDF4_LIB`: `netcdf4-fortran` library path
- `NETCDF4_INCLUDE`: `netcdf4-fortran` include path

There is also an optional variable :

- `DEBUG`: if set to `true`, will configure the compiler to enable debugging options (backtrace etc.) while optimizing for compilation speed. This is recommended for the edit-compile-debug cycle during development. *For production, do net set this flag---that will ensure the compiler optimizes for performance.*

These paths can be obtained using the `nc-config` tool. So to compile, the following command suffices:

- `NETCDF4_LIB=$(nc-config --flibs) NETCDF4_INCLUDE=$(nc-config --fflags) make`

After compilation succeeds, you can run `./imau-fdm.x`.
