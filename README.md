# IMAU-FDM

**IMAU Firn Densification Model** — a one-dimensional firn compaction and
melt–refreeze model for ice sheets and glaciers. `main` always holds a working
version; day-to-day development happens on the `development` branch.

📖 **Documentation:** <https://imau-ice-and-climate.github.io/IMAU-FDM/> (or build
the [`docs/`](./docs) Jupyter Book locally). The end-to-end workflow is
pre-processing → running → post-processing.

## Repository layout

| Directory | Contents |
|-----------|----------|
| `source/` | Fortran model source (built to `bin/imau-fdm`) |
| `launch_job/` | Run launcher + SLURM/MPI submission scripts |
| `settings/` | Per-domain run config: `<domain>/{run,model,constants}.toml` |
| `pre-process-RACMO/` | Build model forcing from RACMO/ERA5 input |
| `post-process/` | Gridded NetCDF output from per-column files |
| `MO_fit/` | Calibrate densification coefficients against firn cores |
| `QAQC/` | Run monitoring and log inspection |
| `reference/` | Static grid files and point lists (per domain) |
| `docs/` | Documentation (Jupyter Book) |

## How a run works

A model run is an **MPI job**: one rank is the **distributor**
(`source/distributor.f90`), which hands grid points out to the remaining **worker**
ranks; each worker runs the one-dimensional column model
(`source/model_main.f90`) point by point. Everything about a run is configured
through TOML files in `settings/<domain>/`.

## Compiling

Built with [`fpm`](https://fpm.fortran-lang.org/); the executable is installed to
`bin/imau-fdm`.

```bash
./compile_hpc.sh            # release build (ECMWF: loads modules, installs fpm via conda if needed)
./compile_hpc.sh debug      # debug build
```

`launch_job.sh` recompiles for you on launch when `recompile = true` in `run.toml`.

## Running

Runs are configured entirely through TOML and launched from `launch_job/`:

1. **Edit the settings** for your domain in `settings/<DOMAIN>/`:
   - `run.toml` — project name, point list, `domain`/`forcing`, `run_type`
     (`ecmwf` or `offline`), `restart_type`, and the `code_dir` / `forcing_root` /
     `output_root` directories.
   - `model.toml` — numerical & physics options. `constants.toml` — physical constants.
2. **Set `SETTINGS_FOLDER`** in `launch_job/launch_job.sh` to your domain folder.
3. **Launch:**
   ```bash
   cd launch_job/
   ./launch_job.sh
   ```
   It reads `run.toml`, (re)compiles if requested, and submits the run — an MPI
   `sbatch` job on ECMWF, or a local MPI run using `[offline] n_procs` when
   `run_type = "offline"`.

## Prerequisites

- `gfortran` and an MPI library (`hpcx-openmpi` on ECMWF; `mpifort` locally)
- `netcdf4-fortran` (install `libnetcdf` first; verify with `nf-config --all`)
- `fpm` — installed automatically by `compile_hpc.sh`, or via the conda
  [`environment.yml`](./environment.yml)

## Supported domains

| Domain | Grid | Forcing |
|--------|------|---------|
| FGRN055 | Greenland 5.5 km | ERA5 / RACMO2.3p2 |
| ANT27 | Antarctica 27 km | ERA5 / RACMO2.3p2 |

Each domain has its own `settings/<domain>/` folder (currently provided: `FGRN055`).
