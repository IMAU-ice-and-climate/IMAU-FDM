# Gridded File Types

All scripts live in `post-process/create_1D_2D_2Ddetail_files/`.

## 1D variables — daily scalars → gridded

**Script:** `make_1d_files.py` | **SLURM:** `submit_make_1d_files.sh`

These are scalar time series per point (runoff, surface melt, firn air content,
etc.). The script aggregates to 10-day or monthly means/sums and optionally
detrends selected variables (`h_surf`, `FirnAir`) to remove spinup drift.

```bash
# All variables, 10-day aggregation
python3 make_1d_files.py --var all --timestep 10day \
    --spinup-start 1940 --spinup-end 1970 --workers 16

# Specific variables
python3 make_1d_files.py --var Runoff surfmelt --timestep 10day --workers 8

# List available variables
python3 make_1d_files.py --list-vars
```

Key arguments:

| Flag | Default | Description |
|------|---------|-------------|
| `--var` | required | Variable name(s), or `all` |
| `--timestep` | `10day` | `daily`, `10day`, or `monthly` |
| `--spinup-start` | 1940 | Detrend window start year |
| `--spinup-end` | 1970 | Detrend window end year |
| `--workers` | 16 | Parallel worker processes |
| `--max-points` | — | Cap for testing (e.g. 500) |

**Available 1D variables:**

| Variable | Units | Description |
|----------|-------|-------------|
| `h_surf` | m | Surface height (detrended) |
| `FirnAir` | m | Firn air content (detrended) |
| `Runoff` | mm w.e. | Runoff (sum per period) |
| `surfmelt` | mm w.e. | Surface melt (sum per period) |
| `refreeze` | mm w.e. | Refreezing (sum per period) |
| `rain` | mm w.e. | Rainfall (sum per period) |
| `TotLwc` | mm | Total liquid water content |
| `icemass` | kg m⁻² | Ice mass |
| `Rho0` | kg m⁻³ | Surface snow density |
| `vice` | m yr⁻¹ | Ice velocity at base |
| `vacc` | m yr⁻¹ | Accumulation velocity |
| `vfc` | m yr⁻¹ | Firn compaction velocity |
| `vmelt` | m yr⁻¹ | Melt velocity |
| `vtotal` | m yr⁻¹ | Total velocity |

---

## 2D variables — monthly depth profiles → scalar gridded

**Script:** `make_2d_files.py` | **SLURM:** `submit_make_2d_files.sh`

Full depth profiles (3000 layers, ~122 m) at 30-day intervals. The script
extracts the depth at which a density threshold is crossed (e.g. firn–ice
transition `z830`).

```bash
# Depth where density reaches 830 kg/m³ (firn–ice transition)
python3 make_2d_files.py --var dens --threshold 830 --output-var z830 \
    -o /path/to/output --workers 8

# Depth where density reaches 550 kg/m³
python3 make_2d_files.py --var dens --threshold 550 --output-var z550 \
    -o /path/to/output --workers 8
```

| Flag | Description |
|------|-------------|
| `--var` | Variable in 2D file (`dens`, `temp`, etc.) |
| `--threshold` | Density threshold (kg m⁻³) |
| `--output-var` | Name for output variable (e.g. `z830`) |
| `--start-year`, `--end-year` | Override output period from config |

---

## 2Ddetail variables — high-res near-surface → scalar gridded

**Script:** `make_2Ddetail_files.py` | **SLURM:** `submit_make_2ddetail_files.sh`

Top 20 m at 4 cm resolution (500 layers), at 10-day intervals. Extracts a value
at a target depth or averages over a depth range — used for 10 m temperature
and surface snow density.

```bash
# 10 m temperature
python3 make_2Ddetail_files.py --var temp --depth 10 --output-var T10m \
    -o /path/to/output --workers 8

# Surface snow density (0–0.5 m mean)
python3 make_2Ddetail_files.py --var dens --depth-begin 0 --depth-end 0.5 \
    --output-var SSN -o /path/to/output --workers 8

# Near-surface liquid water content (0–0.2 m)
python3 make_2Ddetail_files.py --var lwc --depth-begin 0 --depth-end 0.2 \
    --output-var LWC_surf -o /path/to/output --workers 8
```

| Flag | Description |
|------|-------------|
| `--var` | Variable in 2Ddetail file (`temp`, `dens`, `lwc`, etc.) |
| `--depth` | Single target depth in metres |
| `--depth-begin`, `--depth-end` | Depth range for averaging |
| `--output-var` | Name for output variable |
| `--start-year`, `--end-year` | Override output period from config |
