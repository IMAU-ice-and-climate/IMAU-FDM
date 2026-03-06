# IMAU-FDM Percolation Zone Research

**Run:** `FGRN055-era055_1939-2023`
**Goal:** Quantify how the runoff limit and percolation zone have changed over 1939–2023, and diagnose why IMAU-FDM's runoff limit sits lower than observations.

---

## Scientific context

- **Machguth et al. 2026** (*The Cryosphere*) — Comparison of IMAU-FDM, MAR, and MODIS runoff limits over Greenland.
  Key numbers:
  | Model/Obs | Mean runoff limit (m a.s.l.) |
  |---|---|
  | IMAU-FDM | 1545 ± 31 |
  | MODIS    | 1613 ± 99 |
  | MAR      | 1816 ± 94 |

- **Machguth et al. 2016** (*Nature Climate Change*) — Near-surface ice formation limits firn water storage. Key finding: meltwater retention is far lower than expected because ice lenses/slabs impede deep percolation and force lateral runoff.

### Why IMAU-FDM's runoff limit is too low

1. **Bucket runoff condition** — the bucket scheme only generates runoff at the base of the firn column; ice lenses do not trigger lateral runoff.
2. **Cold firn** — IMAU-FDM firn is colder than observations at the K-transect, so less meltwater refreezes near-surface.
3. **Deep firn column** (up to 105 m) — large cold-content reservoir absorbs meltwater that would run off laterally in reality.
4. **RACMO albedo discontinuity** at the ELA — forcing has an abrupt step that pins the runoff limit near the ELA.

### Appropriate runoff fractions at ice structures (MAR approach)
| Firn structure | Runoff fraction |
|---|---|
| Thick ice slab (>0.5 m) | ~75% (field observations) |
| Thin ice lens (<0.1 m) | ~33% (MAR default) |
| Porous firn above the limit | ~0% |

---

## Repository layout

```
research/
├── README.md                        ← this file
├── config.py                        ← all paths and default parameters
├── transient_percolation_analysis.py ← Analyses 1–3 (event frequency, histogram, zone migration)
├── ice_lens_detection.py            ← Ice lens / slab detection from 2Ddetail profiles
├── percolation_zone_research.ipynb  ← Main analysis notebook
└── [runoff_limit_analysis.py]       ← NOT YET WRITTEN — see Pending Work below
```

---

## Data paths

### Static grid (never changes)
```
/home/nld4814/perm/code/IMAU-FDM/reference/FGRN055/
  FGRN055_Masks.nc        # Topography (m), IceMask, Mouginot_basins (1-7), lat, lon
  IN_ll_FGRN055.txt       # Point list: 58 265 rows, cols = lon, lat, ??, ??, ??, rlat_idx, rlon_idx
                           # File number N = row index N-1 (1-indexed files)
```

### Per-column raw output (one file per grid column)
```
/home/nld4814/scratch/run_FGRN055-era055_1939-2023/output/
  FGRN055_era055_1D_{N}.nc        # 1D surface/integrated vars, daily timestep (~30 803 steps)
  FGRN055_era055_2D_{N}.nc        # Full firn column (3000 layers, Lagrangian, monthly)
  FGRN055_era055_2Ddetail_{N}.nc  # Top 20 m, 500 layers × 4 cm, 10-day timestep (3080 steps)
                                   # Variables: dens, temp, lwc, depth, dz, refreeze
```

### Gridded post-processed files (2D space × time)
```
/home/nld4814/scratch/run_FGRN055-era055_1939-2023/post-process/
  FDM_Runoff_FGRN055_1939-2023_10day.nc    # (3081, 566, 438) — mm w.e., 10-day
  FDM_surfmelt_FGRN055_1939-2023_10day.nc  # (3081, 566, 438) — mm w.e., 10-day
  FDM_T10m_FGRN055_1939-2023_2Ddetail.nc  # (3080, 566, 438) — K, 10-day (reference time axis)
  FDM_z830_FGRN055_1939-2023_2D.nc        # (1026, 566, 438) — m, monthly (depth to rho=830)
  FDM_z550_FGRN055_1939-2023_2D.nc        # (1026, 566, 438) — m, monthly (depth to rho=550)
  FDM_refreeze_FGRN055_1939-2023_10day.nc # (3081, 566, 438) — mm w.e., 10-day
  ... [other 10-day gridded variables: vfc, vbouy, rain, solin, TotLwc, vacc, vice,
       vmelt, vsnd, vsub, vtotal, icemass, Rho0, h_surf_detrended, FirnAir_detrended]
```

**Grid dimensions:** rlat 0–565 (566 pts), rlon 0–437 (438 pts)
**Time coordinate:** fractional year float (e.g., 1939.666 = ~1 Sep 1939)

---

## Modules

### `config.py`

Central configuration. All analysis functions read defaults from here.

| Parameter | Value | Meaning |
|---|---|---|
| `SCRATCH_DIR` | `/home/nld4814/scratch` | |
| `PROJECT_NAME` | `run_FGRN055-era055_1939-2023` | |
| `MASKS_FILE` | `.../FGRN055_Masks.nc` | Static grid |
| `POINTLIST_FILE` | `.../IN_ll_FGRN055.txt` | Column → grid mapping |
| `Z830_FILE` | `.../FDM_z830_...2D.nc` | Monthly depth-to-ice |
| `T10M_FILE` | `.../FDM_T10m_...2Ddetail.nc` | 10-day 10m temp (also reference time axis) |
| `EVENT_THRESHOLDS['surfmelt']` | 10.0 mm w.e. | Annual melt threshold |
| `EVENT_THRESHOLDS['Runoff']` | 10.0 mm w.e. | Annual runoff threshold |
| `PERIOD_1` | (1940, 1970) | Histogram early period |
| `PERIOD_2` | (1993, 2023) | Histogram late period |
| `EARLY_WINDOW` | (1940, 1963) | Migration map early window |
| `LATE_WINDOW` | (2000, 2023) | Migration map late window |
| `DRY_EVENTS_MAX` | 1 | Max event years → low-event category |
| `PERCOLATION_EVENTS_MIN` | 2 | Min event years → high-event category |

---

### `transient_percolation_analysis.py`

Three analyses on any gridded variable in `config.VARIABLE_FILES`.

**`_load_annual_totals(variable='surfmelt')`**
Internal helper. Converts fractional-year time → integer year, groups and sums within each year. Returns `(year, rlat, rlon)` DataArray. Preserves off-ice NaN.

**`map_melt_frequency(variable, threshold, years_back, ...)`**
Map: how many years each pixel exceeded the annual threshold. Optional `years_back` to restrict to the last N years.

**`plot_melt_histogram(variable, period1, period2, threshold, ...)`**
Side-by-side bar chart: per-pixel event count distribution for two periods.

**`plot_percolation_migration(variable, early_window, late_window, ...)`**
4-category map: stayed dry / stayed percolation / expanded (dry→percolation) / contracted.
Categories: 1=stayed low, 2=stayed high, 3=expanded, 4=contracted.

---

### `ice_lens_detection.py`

Detects ice lenses/slabs in 2Ddetail subsurface density profiles.

**Definition:** a layer qualifies as an ice lens if `dens[i] >= 900 kg/m³` AND at least one layer below it has `dens < 900 kg/m³` (ice sitting on top of permeable firn — the configuration that blocks downward percolation).

**`_compute_ice_lens_diagnostics(dens, depth, dz, min_ice_dens=900.0, require_firn_above=False)`**
Core vectorised detection. Inputs are raw numpy arrays `(n_layers, n_time)`.
- `require_firn_above=False`: includes superimposed ice at the surface
- `require_firn_above=True`: only embedded lenses (firn → ice → firn geometry)
Returns dict: `has_ice_lens`, `depth_top_lens`, `ice_lens_thickness` — all `(n_time,)`.

**`detect_ice_lenses_in_column(nc_file, min_ice_dens=900.0, require_firn_above=False, time_coord=None)`**
Loads one `2Ddetail_N.nc` file, returns `xr.Dataset` with fractional-year time coordinate. Time axis sourced from `config.T10M_FILE`.

**`create_gridded_ice_lens_file(output_file, output_dir=None, min_ice_dens=900.0, require_firn_above=False, verbose=True)`**
Batch loop over all 58 265 per-column 2Ddetail files. Assembles `(time, rlat, rlon)` NetCDF. **Run as a job, not interactively.**

**Test results** (K-transect, ~67°N, ~47–48°W):
| Column | Elevation | Ice lens fraction | Notes |
|---|---|---|---|
| 8840 | 1494 m | 3.8% (embedded) | Occasional superimposed ice |
| 8841 | 1542 m | 18.2% (embedded) | Active percolation zone |
| 8842 | 1590 m | 28.8% (embedded) | More firn above base ice |
| 30000 | Accumulation | 0% | No melt |

The rlat/rlon grid indices for these test columns: 8840→(186,114), 8841→(186,115), 8842→(186,116).

---

## Pending work

### 1. `runoff_limit_analysis.py` — NOT YET WRITTEN

Two functions were fully planned but not implemented (ice lens detection was prioritised):

**`plot_runoff_limit_timeseries(threshold, basin, ax, show_benchmarks)`**
- Load annual Runoff totals via `_load_annual_totals('Runoff')`
- Load `Topography` from `MASKS_FILE`
- For each year: `topo.where((annual_runoff > threshold) & ice_mask).max(dim=['rlat','rlon'])`
- Plot with MODIS (1613±99 m) and MAR (1816±94 m) benchmark bands
- Optional Mouginot basin filter (column 7 of Masks.nc, values 0–7)

**`plot_near_surface_ice_at_runoff_limit(threshold, elev_band, ax)`**
- Compute annual runoff limit elevation per year
- For each year: pixels in `[max_elev, max_elev + elev_band]` — the transition zone
- Get annual-mean z830 at those pixels (from `Z830_FILE`, aggregate monthly → annual mean)
- Plot z830 vs year (y-axis inverted: low z830 = ice near surface = percolation zone expanding)

### 2. Run `create_gridded_ice_lens_file()` as a batch job

```python
from ice_lens_detection import create_gridded_ice_lens_file
create_gridded_ice_lens_file('FDM_ice_lens_FGRN055_1939-2023_2Ddetail.nc')
```

Output will be ~3 GB at `(3080, 566, 438)` int8/float32. Add to `config.VARIABLE_FILES` or a new `config.ICE_LENS_FILE` path once created.

### 3. Combine ice lens + runoff limit

Once both are available:
- At pixels just above the annual runoff limit: what fraction of 10-day timesteps have an ice lens?
- Time series of that fraction → shows whether the transition zone is progressively building ice
- Compare with refreeze totals (`FDM_refreeze_FGRN055_1939-2023_10day.nc`)

---

## Quick-start (notebook)

```python
import sys
sys.path.insert(0, '/home/nld4814/perm/code/IMAU-FDM/research')
import config
from transient_percolation_analysis import map_melt_frequency, plot_melt_histogram, plot_percolation_migration
from ice_lens_detection import detect_ice_lenses_in_column

# --- Existing analyses ---
count, ax = map_melt_frequency(variable='Runoff')
ax2 = plot_melt_histogram(variable='Runoff')
classif, ax3 = plot_percolation_migration(variable='Runoff')

# --- Ice lens on a single column ---
result = detect_ice_lenses_in_column(
    '/home/nld4814/scratch/run_FGRN055-era055_1939-2023/output/FGRN055_era055_2Ddetail_8841.nc',
    require_firn_above=True,   # embedded lenses only
)
result['ice_lens_thickness'].plot()
```
