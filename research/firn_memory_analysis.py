"""
Firn memory analysis from IMAU-FDM 2D Lagrangian profile output.

Firn memory is defined as the time (in years) for a surface firn parcel
to reach a target depth (default 10 m).  It is read directly from the
``year`` variable in the per-column 2D monthly files, which stores the
fractional simulation year at which each Lagrangian layer was deposited.

Physical interpretation
-----------------------
IMAU-FDM tracks each firn layer as a Lagrangian parcel.  The ``year``
variable records when that parcel was first deposited at the surface (as
a fraction of years from the model start, 0 = 1939-09-01).  At any
timestep ``t``, the layer currently sitting at depth ``d`` has

    firn_memory = time_sim[t] - year_at_d[t]

which is the age of the oldest firn at that depth — equivalently, the
time it took the original surface snowfall to be buried to depth ``d``.

Regions where the column does not reach the target depth (ablation zone,
thin firn) are returned as NaN.

Public interface
----------------
compute_firn_memory_column(nc_file, time_sim, target_depth=10.0)
    Process one 2D per-column file.  Returns a (n_time,) float32 array.

create_gridded_firn_memory(output_file, target_depth=10.0, ...)
    Loop over all per-column 2D files, assemble onto the FGRN055 grid,
    and save a NetCDF.  Batch operation (~58 265 files) — run as a job.
"""

import sys
from pathlib import Path
from datetime import datetime

import numpy as np
import xarray as xr
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import config


# =============================================================================
# Parallel worker (module-level so it is picklable)
# =============================================================================

# Shared state initialised once per worker process via ProcessPoolExecutor
# initializer — avoids re-serialising large arrays with every task.
_W = {}   # _W['time_sim'], _W['n_time'], _W['target_depth'], _W['input_dir']


def _fm_worker_init(time_sim, n_time, target_depth, input_dir_str):
    _W['time_sim']     = time_sim
    _W['n_time']       = n_time
    _W['target_depth'] = target_depth
    _W['input_dir']    = Path(input_dir_str)


def _fm_worker(nn_ri_rj):
    """Process one per-column 2D file. Returns (ri, rj, nt, fm_10m, fm_z830) or None."""
    nn, ri, rj        = nn_ri_rj
    time_sim          = _W['time_sim']
    n_time            = _W['n_time']
    target_depth      = _W['target_depth']
    fname             = _W['input_dir'] / f'FGRN055_era055_2D_{nn + 1}.nc'
    if not fname.exists():
        return None
    try:
        ds    = xr.open_dataset(fname, decode_times=False)
        depth = ds['depth'].values
        year  = ds['year'].values
        dens  = ds['dens'].values
        ds.close()
        nt     = min(depth.shape[1], n_time)
        fm_10m = _compute_firn_memory(depth[:, :nt], year[:, :nt], time_sim[:nt], target_depth)
        fm_z830, z830_depth = _compute_firn_memory_z830(
            dens[:, :nt], year[:, :nt], depth[:, :nt], time_sim[:nt])
        fm_10m = np.where(z830_depth >= target_depth, fm_10m, np.nan).astype(np.float32)
        return (ri, rj, nt, fm_10m, fm_z830)
    except Exception as exc:
        print(f'  WARNING col {nn + 1}: {exc}', flush=True)
        return None


# =============================================================================
# Core per-column computation (fully vectorised over layers × time)
# =============================================================================

def _compute_firn_memory(depth, year, time_sim, target_depth=10.0):
    """
    Compute the firn memory time series for one firn column.

    Parameters
    ----------
    depth    : (n_layers, n_time) float array
        Layer centre depths (m).  Layer 0 is the deepest (oldest); the
        highest-index active layer is the most recently deposited surface
        material.
    year     : (n_layers, n_time) float array
        Fractional simulation year of deposition for each layer
        (0 = model start 1939-09-01).  Increases monotonically from
        deep (old) to shallow (new) within a valid column.
    time_sim : (n_time,) float array
        Fractional simulation year at each timestep (0 at t=0).
    target_depth : float
        Depth (m) at which to evaluate firn memory (default 10 m).

    Returns
    -------
    (n_time,) float32 array
        Firn memory in years.  NaN where the column does not span
        target_depth, or where data are invalid.
    """
    FILL = 9.9e35

    depth = np.where(np.abs(depth) > FILL, np.nan, depth).astype(np.float64)
    year  = np.where(np.abs(year)  > FILL, np.nan, year).astype(np.float64)

    n_layers, n_time = depth.shape

    # Distance from target_depth for each layer × time; inf where invalid
    dist = np.abs(depth - target_depth)                           # (n_layers, n_time)
    dist = np.where(np.isfinite(depth) & np.isfinite(year), dist, np.inf)

    # Layer index closest to target depth, per timestep
    idx      = np.argmin(dist, axis=0)                            # (n_time,)
    min_dist = dist.min(axis=0)                                   # (n_time,)

    # Year value of the nearest layer
    yr_at_target = year[idx, np.arange(n_time)]                   # (n_time,)

    # Firn memory = elapsed time since that layer was deposited
    firn_memory = (time_sim - yr_at_target).astype(np.float32)

    # Mask where:
    #   • nearest valid layer is >1 m from target (column doesn't span it)
    #   • result is non-positive (fill-value contamination or model artefact)
    valid = (min_dist <= 1.0) & (firn_memory > 0)
    return np.where(valid, firn_memory, np.nan).astype(np.float32)


def _compute_firn_memory_z830(dens, year, depth, time_sim):
    """
    Compute the time (years) for a surface firn parcel to reach the
    firn-ice transition (830 kg/m³ density horizon, z830), and the depth
    of that horizon.

    The z830 horizon is defined as the top of the *continuous* ice body
    that extends uninterrupted from the column base (layer 0) upward.
    Isolated ice lenses or slabs embedded in the firn column are ignored —
    a layer only qualifies if every layer below it is also ≥ 830 kg/m³.

    This is implemented via a cumulative product from the base upward:
    a single non-ice layer anywhere below zeroes out all layers above it,
    so only the unbroken basal ice body contributes.

    Parameters
    ----------
    dens     : (n_layers, n_time) float array
        Layer density (kg/m³).  Layer 0 = deepest; layer N-1 = surface.
    year     : (n_layers, n_time) float array
        Fractional simulation year of deposition for each layer.
    depth    : (n_layers, n_time) float array
        Layer centre depths (m).
    time_sim : (n_time,) float array
        Fractional simulation year at each timestep.

    Returns
    -------
    fm_z830 : (n_time,) float32
        Firn memory at z830 in years.  NaN where the column base is not
        ice or where data are invalid.
    z830_depth : (n_time,) float32
        Depth (m) of the z830 horizon.  NaN where no continuous basal
        ice exists.
    """
    FILL = 9.9e35

    dens  = np.where(np.abs(dens)  > FILL, np.nan, dens ).astype(np.float64)
    year  = np.where(np.abs(year)  > FILL, np.nan, year ).astype(np.float64)
    depth = np.where(np.abs(depth) > FILL, np.nan, depth).astype(np.float64)

    n_layers, n_time = dens.shape

    # Ice mask: layer is part of continuous basal ice only if dens ≥ 830
    # and data are valid.
    ice_mask = (dens >= 830.0) & np.isfinite(dens) & np.isfinite(year)

    # Cumulative AND from the base (layer 0) upward.
    # cum_ice[i, t] = 1 iff ALL layers 0 … i are ≥ 830 kg/m³.
    # A single non-ice layer zeroes everything above it.
    cum_ice = np.cumprod(ice_mask.astype(np.int8), axis=0)   # (n_layers, n_time)

    # Number of consecutive ice layers from the base = sum of 1s.
    # Top of continuous ice = that count − 1 (0-indexed layer).
    # A value of -1 means no ice at the base → NaN.
    z830_layer = cum_ice.sum(axis=0).astype(int) - 1          # (n_time,)

    has_continuous_ice = z830_layer >= 0

    # Clip to valid index range before fancy-indexing (avoids −1 wraparound)
    idx         = np.clip(z830_layer, 0, n_layers - 1)
    yr_at_z830  = year [idx, np.arange(n_time)]
    dep_at_z830 = depth[idx, np.arange(n_time)]
    fm_z830     = (time_sim - yr_at_z830).astype(np.float32)

    # Three output categories:
    #
    #   > 0        Physically meaningful: z830 formed during the main simulation
    #              (year >= 0 ↔ deposited after 1939-09-01).
    #
    #   SPINUP_FILL (-9999)  z830 exists in the column but was already present
    #              before the simulation started.  The layer was created during
    #              the spin-up (42 × 30-yr repeated ERA5 cycles) so its age is a
    #              model artefact, not a historically meaningful deposition date.
    #              Distinct from NaN so downstream plots can colour it separately.
    #
    #   NaN        No continuous ice body at the column base (no z830), or data
    #              quality issues.  Represents truly ice-free / firn-only columns.
    SPINUP_FILL = np.float32(-9999.0)

    valid_hist  = (has_continuous_ice & np.isfinite(yr_at_z830)
                   & (fm_z830 > 0) & (yr_at_z830 >= 0.0))
    spinup_era  = (has_continuous_ice & np.isfinite(yr_at_z830)
                   & (yr_at_z830 < 0.0))

    fm_out = np.where(valid_hist, fm_z830,    np.nan).astype(np.float32)
    fm_out = np.where(spinup_era, SPINUP_FILL, fm_out)           # override

    # z830_depth: real depth where z830 exists (either era), NaN otherwise
    depth_out = np.where(has_continuous_ice & np.isfinite(dep_at_z830),
                         dep_at_z830, np.nan).astype(np.float32)

    return fm_out, depth_out


# =============================================================================
# Single-column convenience wrapper
# =============================================================================

def compute_firn_memory_column(nc_file, time_sim, target_depth=10.0):
    """
    Compute firn memory for one per-column 2D output file.

    Parameters
    ----------
    nc_file : Path or str
        Path to a ``FGRN055_era055_2D_N.nc`` per-column output file.
    time_sim : (n_time,) array
        Fractional simulation year at each of the 1026 monthly timesteps
        (0 = model start 1939-09-01).
    target_depth : float
        Depth (m) at which to measure firn memory (default 10 m).

    Returns
    -------
    (n_time,) float32 array
        Firn memory in years.  NaN where the column does not reach
        target_depth.
    """
    nc_file = Path(nc_file)
    ds      = xr.open_dataset(nc_file, decode_times=False)
    depth   = ds['depth'].values   # (n_layers, n_time)
    year    = ds['year'].values    # (n_layers, n_time)
    dens    = ds['dens'].values    # (n_layers, n_time)
    ds.close()

    nt = min(depth.shape[1], len(time_sim))
    fm_10m = _compute_firn_memory(
        depth[:, :nt], year[:, :nt], time_sim[:nt], target_depth
    )
    # Mask z10m where the firn column is shallower than target_depth
    _, z830_depth = _compute_firn_memory_z830(
        dens[:, :nt], year[:, :nt], depth[:, :nt], time_sim[:nt]
    )
    fm_10m = np.where(z830_depth >= target_depth, fm_10m, np.nan).astype(np.float32)
    return fm_10m


# =============================================================================
# Gridded post-processing (batch — run as a job)
# =============================================================================

def create_gridded_firn_memory(
    output_file,
    output_dir=None,
    target_depth=10.0,
    file_limit=None,
    n_workers=1,
    verbose=True,
):
    """
    Create a gridded (time, rlat, rlon) NetCDF of firn memory.

    Loops over per-column 2D files listed in the FGRN055 pointlist,
    computes firn memory at two depths for each column, and assembles
    the results onto the 2-D model grid.

    This is a batch operation (~58 265 files).  Run it as a job script,
    not interactively.  For a quick test run pass file_limit=500 and
    a test output filename.

    Parameters
    ----------
    output_file : str
        Filename of the output NetCDF file (placed in output_dir).
    output_dir : Path or str, optional
        Directory for the output file.  Defaults to config.PROCESSED_DIR.
    target_depth : float
        Fixed depth (m) at which to measure firn memory (default 10 m).
    file_limit : int, optional
        If given, stop after processing this many columns.  Useful for
        testing.  Default None processes all columns.
    verbose : bool
        Print progress every 1 000 columns (default True).

    Output variables
    ----------------
    firn_memory      (time, rlat, rlon)  — firn memory at target_depth (years)
    firn_memory_mean (rlat, rlon)        — time-mean
    firn_memory_std  (rlat, rlon)        — temporal std
    firn_memory_z830      (time, rlat, rlon)  — firn memory at z830 (years)
    firn_memory_z830_mean (rlat, rlon)        — time-mean
    firn_memory_z830_std  (rlat, rlon)        — temporal std
    lat / lon        (rlat, rlon)        — geographic coordinates
    """
    if output_dir is None:
        output_dir = config.PROCESSED_DIR
    output_path = Path(output_dir) / output_file

    # --- Grid metadata ---
    masks  = xr.open_dataset(config.MASKS_FILE)
    n_rlat = len(masks['rlat'])
    n_rlon = len(masks['rlon'])
    lat_2d = masks['lat'].values
    lon_2d = masks['lon'].values
    masks.close()

    # --- Time axis from the gridded z830 file (same monthly steps as 2D files) ---
    ref      = xr.open_dataset(config.Z830_FILE, decode_times=False)
    time_arr = ref['time'].values   # fractional calendar year, shape (n_time,)
    n_time   = len(time_arr)
    ref.close()
    time_sim = time_arr - time_arr[0]   # simulation years from model start (0 at t=0)

    # --- Point list: maps file number → (rlat_idx, rlon_idx) ---
    grid     = pd.read_csv(config.POINTLIST_FILE, header=None, sep=',')
    n_cols   = len(grid)
    rlat_col = grid.iloc[:, 5].astype(int).values
    rlon_col = grid.iloc[:, 6].astype(int).values

    if file_limit is not None:
        n_cols = min(n_cols, file_limit)
        if verbose:
            print(f'  file_limit={file_limit}: processing first {n_cols} columns only')

    # --- Allocate output arrays ---
    fm_out      = np.full((n_time, n_rlat, n_rlon), np.nan, dtype=np.float32)
    fm_z830_out = np.full((n_time, n_rlat, n_rlon), np.nan, dtype=np.float32)

    input_dir = Path(config.SCRATCH_DIR) / config.PROJECT_NAME / 'output'

    # --- Main loop ---
    n_processed = 0

    if n_workers > 1:
        from concurrent.futures import ProcessPoolExecutor
        tasks = [(nn, int(rlat_col[nn]), int(rlon_col[nn])) for nn in range(n_cols)]
        if verbose:
            print(f'  parallel mode: {n_workers} workers, {n_cols} columns', flush=True)
        with ProcessPoolExecutor(
            max_workers=n_workers,
            initializer=_fm_worker_init,
            initargs=(time_sim, n_time, target_depth, str(input_dir)),
        ) as pool:
            for result in pool.map(_fm_worker, tasks, chunksize=100):
                if result is not None:
                    ri, rj, nt, fm_10m, fm_z830 = result
                    fm_out     [:nt, ri, rj] = fm_10m
                    fm_z830_out[:nt, ri, rj] = fm_z830
                    n_processed += 1
                    if verbose and n_processed % 1000 == 0:
                        print(f'  processed {n_processed} / {n_cols}', flush=True)
    else:
        for nn in range(n_cols):
            file_num = nn + 1
            fname    = input_dir / f'FGRN055_era055_2D_{file_num}.nc'
            if not fname.exists():
                continue
            ri = rlat_col[nn]
            rj = rlon_col[nn]
            try:
                ds    = xr.open_dataset(fname, decode_times=False)
                depth = ds['depth'].values
                year  = ds['year'].values
                dens  = ds['dens'].values
                ds.close()
                nt = min(depth.shape[1], n_time)
                fm_10m = _compute_firn_memory(
                    depth[:, :nt], year[:, :nt], time_sim[:nt], target_depth)
                fm_z830, z830_depth = _compute_firn_memory_z830(
                    dens[:, :nt], year[:, :nt], depth[:, :nt], time_sim[:nt])
                fm_10m = np.where(z830_depth >= target_depth, fm_10m, np.nan).astype(np.float32)
                fm_out     [:nt, ri, rj] = fm_10m
                fm_z830_out[:nt, ri, rj] = fm_z830
            except Exception as exc:
                if verbose:
                    print(f'  WARNING: skipping file {file_num}: {exc}')
                continue
            n_processed += 1
            if verbose and n_processed % 1000 == 0:
                print(f'  processed {n_processed} / {n_cols}')

    if verbose:
        print(f'  done — {n_processed} columns written', flush=True)

    # --- Summary statistics ---
    fm_mean      = np.nanmean(fm_out,      axis=0).astype(np.float32)
    fm_std       = np.nanstd( fm_out,      axis=0).astype(np.float32)
    fm_z830_mean = np.nanmean(fm_z830_out, axis=0).astype(np.float32)
    fm_z830_std  = np.nanstd( fm_z830_out, axis=0).astype(np.float32)

    # --- Build and save dataset ---
    out = xr.Dataset(
        {
            'firn_memory': xr.DataArray(
                fm_out, dims=('time', 'rlat', 'rlon'),
                attrs={
                    'long_name'  : f'Firn memory at {target_depth} m depth',
                    'units'      : 'years',
                    'description': (
                        f'Time elapsed since the firn layer now at {target_depth} m '
                        'was deposited at the surface.  NaN in the ablation zone '
                        'or where the firn column does not reach this depth.'
                    ),
                },
            ),
            'firn_memory_mean': xr.DataArray(
                fm_mean, dims=('rlat', 'rlon'),
                attrs={
                    'long_name': f'Time-mean firn memory at {target_depth} m',
                    'units'    : 'years',
                },
            ),
            'firn_memory_std': xr.DataArray(
                fm_std, dims=('rlat', 'rlon'),
                attrs={
                    'long_name': f'Temporal std of firn memory at {target_depth} m',
                    'units'    : 'years',
                },
            ),
            'firn_memory_z830': xr.DataArray(
                fm_z830_out, dims=('time', 'rlat', 'rlon'),
                attrs={
                    'long_name'  : 'Firn memory at z830 (firn–ice transition)',
                    'units'      : 'years',
                    'description': (
                        'Time elapsed since the firn layer now at the 830 kg/m³ '
                        'density horizon was deposited at the surface.  '
                        'Positive values: z830 formed during the 1939–2023 '
                        'simulation (physically meaningful).  '
                        '-9999: z830 exists but was already present before 1939 '
                        '(spin-up artefact; historical run does not reach z830).  '
                        'NaN: no continuous basal ice in the column (truly no '
                        'firn–ice transition detected).'
                    ),
                    'spinup_fill_value': -9999.0,
                },
            ),
            'firn_memory_z830_mean': xr.DataArray(
                fm_z830_mean, dims=('rlat', 'rlon'),
                attrs={
                    'long_name': 'Time-mean firn memory at z830',
                    'units'    : 'years',
                },
            ),
            'firn_memory_z830_std': xr.DataArray(
                fm_z830_std, dims=('rlat', 'rlon'),
                attrs={
                    'long_name': 'Temporal std of firn memory at z830',
                    'units'    : 'years',
                },
            ),
            'lat': xr.DataArray(
                lat_2d, dims=('rlat', 'rlon'),
                attrs={'units': 'degrees_north', 'standard_name': 'latitude'},
            ),
            'lon': xr.DataArray(
                lon_2d, dims=('rlat', 'rlon'),
                attrs={'units': 'degrees_east', 'standard_name': 'longitude'},
            ),
        },
        coords={
            'time': xr.DataArray(
                time_arr,
                attrs={'long_name': 'time', 'units': 'fractional calendar year'},
            ),
            'rlat': np.arange(n_rlat, dtype=np.float32),
            'rlon': np.arange(n_rlon, dtype=np.float32),
        },
        attrs={
            'title'       : f'IMAU-FDM firn memory at {target_depth} m and z830',
            'source'      : 'IMAU-FDM version 1.2+',
            'domain'      : 'FGRN055',
            'institution' : 'IMAU, Utrecht University',
            'history'     : f'Created on {datetime.utcnow().isoformat()}',
            'Conventions' : 'CF-1.8',
            'target_depth': f'{target_depth} m',
        },
    )

    out.to_netcdf(output_path)
    if verbose:
        print(f'Saved: {output_path}')
    return output_path


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Create gridded firn memory NetCDF')
    parser.add_argument('--output',     default='FDM_firn_memory_FGRN055_1939-2023_2D.nc',
                        help='Output filename (default: FDM_firn_memory_FGRN055_1939-2023_2D.nc)')
    parser.add_argument('--output-dir', default=None,
                        help='Output directory (default: config.PROCESSED_DIR)')
    parser.add_argument('--target-depth', type=float, default=10.0,
                        help='Fixed depth in metres (default: 10.0)')
    parser.add_argument('--file-limit', type=int, default=None,
                        help='Stop after N columns — useful for test runs')
    parser.add_argument('--n-workers', type=int, default=1,
                        help='Number of parallel worker processes (default: 1)')
    args = parser.parse_args()

    create_gridded_firn_memory(
        output_file  = args.output,
        output_dir   = args.output_dir,
        target_depth = args.target_depth,
        file_limit   = args.file_limit,
        n_workers    = args.n_workers,
        verbose      = True,
    )
