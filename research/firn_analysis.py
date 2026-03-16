"""
Firn column analyses from IMAU-FDM 2D monthly Lagrangian output.

Two analyses, both looping over the ~58 265 per-column 2D files:

  firn_memory
    Time (in years) for a surface firn parcel to descend to a target depth
    (default 10 m).  Also computed at the firn-ice transition (z830).

    Public interface
    ----------------
    compute_firn_memory_column(nc_file, time_sim, target_depth=10.0)
        Single column; returns (n_time,) float32 array.

    create_gridded_firn_memory(output_file, target_depth=10.0, ...)
        Batch gridded output; saves firn_memory (z10m) and firn_memory_z830.

  firn_ice_content
    Integrated ice content (kg m⁻²) within the firn column above the
    firn-ice transition (top of the continuous basal ice body).

    Public interface
    ----------------
    compute_firn_ice_content_column(nc_file)
        Single column; returns (n_time,) float32 array.

    create_gridded_firn_ice_content(output_file, ...)
        Batch gridded output; saves firn_ice_content and summary statistics.

Run as a script
---------------
    python3 firn_analysis.py memory      [options]
    python3 firn_analysis.py ice_content [options]
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
# Firn memory — parallel workers
# =============================================================================

_FM_W = {}   # shared state initialised once per worker process


def _fm_worker_init(time_sim, n_time, target_depth, input_dir_str):
    _FM_W['time_sim']     = time_sim
    _FM_W['n_time']       = n_time
    _FM_W['target_depth'] = target_depth
    _FM_W['input_dir']    = Path(input_dir_str)


def _fm_worker(nn_ri_rj):
    """Process one per-column 2D file. Returns (ri, rj, nt, fm_10m, fm_z830) or None."""
    nn, ri, rj        = nn_ri_rj
    time_sim          = _FM_W['time_sim']
    n_time            = _FM_W['n_time']
    target_depth      = _FM_W['target_depth']
    fname             = _FM_W['input_dir'] / f'FGRN055_era055_2D_{nn + 1}.nc'
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
# Firn memory — core computation
# =============================================================================

def _compute_firn_memory(depth, year, time_sim, target_depth=10.0):
    """
    Compute the firn memory time series for one firn column.

    Parameters
    ----------
    depth    : (n_layers, n_time) float array
    year     : (n_layers, n_time) float array — fractional simulation year of deposition
    time_sim : (n_time,) float array — fractional simulation year at each timestep
    target_depth : float

    Returns
    -------
    (n_time,) float32 — firn memory in years; NaN where column < target_depth
    """
    FILL = 9.9e35
    depth = np.where(np.abs(depth) > FILL, np.nan, depth).astype(np.float64)
    year  = np.where(np.abs(year)  > FILL, np.nan, year ).astype(np.float64)

    n_layers, n_time = depth.shape
    dist     = np.abs(depth - target_depth)
    dist     = np.where(np.isfinite(depth) & np.isfinite(year), dist, np.inf)
    idx      = np.argmin(dist, axis=0)
    min_dist = dist.min(axis=0)

    yr_at_target = year[idx, np.arange(n_time)]
    firn_memory  = (time_sim - yr_at_target).astype(np.float32)
    valid        = (min_dist <= 1.0) & (firn_memory > 0)
    return np.where(valid, firn_memory, np.nan).astype(np.float32)


def _compute_firn_memory_z830(dens, year, depth, time_sim):
    """
    Compute firn memory at the z830 (830 kg/m³ density) horizon.

    The z830 horizon is the top of the *continuous* ice body from the column
    base upward.  Isolated ice lenses are ignored.

    Returns
    -------
    fm_z830 : (n_time,) float32
        Positive = valid (year >= 0); -9999 = spin-up artefact (year < 0); NaN = no z830.
    z830_depth : (n_time,) float32
        Depth of the z830 horizon; NaN where absent.
    """
    FILL        = 9.9e35
    SPINUP_FILL = np.float32(-9999.0)

    dens  = np.where(np.abs(dens)  > FILL, np.nan, dens ).astype(np.float64)
    year  = np.where(np.abs(year)  > FILL, np.nan, year ).astype(np.float64)
    depth = np.where(np.abs(depth) > FILL, np.nan, depth).astype(np.float64)

    n_layers, n_time = dens.shape
    ice_mask = (dens >= 830.0) & np.isfinite(dens) & np.isfinite(year)
    cum_ice  = np.cumprod(ice_mask.astype(np.int8), axis=0)
    z830_layer           = cum_ice.sum(axis=0).astype(int) - 1
    has_continuous_ice   = z830_layer >= 0
    idx                  = np.clip(z830_layer, 0, n_layers - 1)
    yr_at_z830           = year [idx, np.arange(n_time)]
    dep_at_z830          = depth[idx, np.arange(n_time)]
    fm_z830              = (time_sim - yr_at_z830).astype(np.float32)

    valid_hist = (has_continuous_ice & np.isfinite(yr_at_z830)
                  & (fm_z830 > 0) & (yr_at_z830 >= 0.0))
    spinup_era = (has_continuous_ice & np.isfinite(yr_at_z830) & (yr_at_z830 < 0.0))

    fm_out    = np.where(valid_hist, fm_z830,    np.nan).astype(np.float32)
    fm_out    = np.where(spinup_era, SPINUP_FILL, fm_out)
    depth_out = np.where(has_continuous_ice & np.isfinite(dep_at_z830),
                         dep_at_z830, np.nan).astype(np.float32)
    return fm_out, depth_out


# =============================================================================
# Firn memory — public interface
# =============================================================================

def compute_firn_memory_column(nc_file, time_sim, target_depth=10.0):
    """
    Compute firn memory for one per-column 2D output file.

    Parameters
    ----------
    nc_file : Path or str
    time_sim : (n_time,) array — fractional simulation year (0 = model start)
    target_depth : float — depth in metres (default 10 m)

    Returns
    -------
    (n_time,) float32 — firn memory in years; NaN where column < target_depth
    """
    nc_file = Path(nc_file)
    ds      = xr.open_dataset(nc_file, decode_times=False)
    depth   = ds['depth'].values
    year    = ds['year'].values
    dens    = ds['dens'].values
    ds.close()

    nt     = min(depth.shape[1], len(time_sim))
    fm_10m = _compute_firn_memory(depth[:, :nt], year[:, :nt], time_sim[:nt], target_depth)
    _, z830_depth = _compute_firn_memory_z830(
        dens[:, :nt], year[:, :nt], depth[:, :nt], time_sim[:nt])
    return np.where(z830_depth >= target_depth, fm_10m, np.nan).astype(np.float32)


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

    Loops over per-column 2D files, computes firn memory at target_depth and
    z830 for each column, and assembles the results onto the FGRN055 grid.
    Batch operation (~58 265 files) — run as a job.

    Output variables: firn_memory, firn_memory_mean, firn_memory_std,
                      firn_memory_z830, firn_memory_z830_mean, firn_memory_z830_std,
                      lat, lon.
    """
    if output_dir is None:
        output_dir = config.PROCESSED_DIR
    output_path = Path(output_dir) / output_file

    masks    = xr.open_dataset(config.MASKS_FILE)
    n_rlat   = len(masks['rlat'])
    n_rlon   = len(masks['rlon'])
    lat_2d   = masks['lat'].values
    lon_2d   = masks['lon'].values
    masks.close()

    ref      = xr.open_dataset(config.Z830_FILE, decode_times=False)
    time_arr = ref['time'].values
    n_time   = len(time_arr)
    ref.close()
    time_sim = time_arr - time_arr[0]

    grid     = pd.read_csv(config.POINTLIST_FILE, header=None, sep=',')
    n_cols   = len(grid)
    rlat_col = grid.iloc[:, 5].astype(int).values
    rlon_col = grid.iloc[:, 6].astype(int).values

    if file_limit is not None:
        n_cols = min(n_cols, file_limit)
        if verbose:
            print(f'  file_limit={file_limit}: processing first {n_cols} columns only')

    fm_out      = np.full((n_time, n_rlat, n_rlon), np.nan, dtype=np.float32)
    fm_z830_out = np.full((n_time, n_rlat, n_rlon), np.nan, dtype=np.float32)
    input_dir   = Path(config.SCRATCH_DIR) / config.PROJECT_NAME / 'output'
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
            fname = input_dir / f'FGRN055_era055_2D_{nn + 1}.nc'
            if not fname.exists():
                continue
            ri, rj = rlat_col[nn], rlon_col[nn]
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
                fm_out     [:nt, ri, rj] = np.where(z830_depth >= target_depth, fm_10m, np.nan)
                fm_z830_out[:nt, ri, rj] = fm_z830
            except Exception as exc:
                if verbose:
                    print(f'  WARNING: skipping file {nn + 1}: {exc}')
                continue
            n_processed += 1
            if verbose and n_processed % 1000 == 0:
                print(f'  processed {n_processed} / {n_cols}')

    if verbose:
        print(f'  done — {n_processed} columns written', flush=True)

    out = xr.Dataset(
        {
            'firn_memory': xr.DataArray(
                fm_out, dims=('time', 'rlat', 'rlon'),
                attrs={'long_name': f'Firn memory at {target_depth} m depth',
                       'units': 'years',
                       'description': (
                           f'Time for a surface layer to reach {target_depth} m.  '
                           'NaN in ablation zone or where firn < target_depth.')}),
            'firn_memory_mean': xr.DataArray(
                np.nanmean(fm_out, axis=0).astype(np.float32), dims=('rlat', 'rlon'),
                attrs={'long_name': f'Time-mean firn memory at {target_depth} m', 'units': 'years'}),
            'firn_memory_std': xr.DataArray(
                np.nanstd(fm_out, axis=0).astype(np.float32), dims=('rlat', 'rlon'),
                attrs={'long_name': f'Temporal std of firn memory at {target_depth} m', 'units': 'years'}),
            'firn_memory_z830': xr.DataArray(
                fm_z830_out, dims=('time', 'rlat', 'rlon'),
                attrs={'long_name': 'Firn memory at z830 (firn–ice transition)',
                       'units': 'years',
                       'description': (
                           'Positive: z830 formed 1939–2023.  '
                           '-9999: spin-up artefact (z830 pre-dates simulation).  '
                           'NaN: no continuous basal ice.'),
                       'spinup_fill_value': -9999.0}),
            'firn_memory_z830_mean': xr.DataArray(
                np.nanmean(fm_z830_out, axis=0).astype(np.float32), dims=('rlat', 'rlon'),
                attrs={'long_name': 'Time-mean firn memory at z830', 'units': 'years'}),
            'firn_memory_z830_std': xr.DataArray(
                np.nanstd(fm_z830_out, axis=0).astype(np.float32), dims=('rlat', 'rlon'),
                attrs={'long_name': 'Temporal std of firn memory at z830', 'units': 'years'}),
            'lat': xr.DataArray(lat_2d, dims=('rlat', 'rlon'),
                attrs={'units': 'degrees_north', 'standard_name': 'latitude'}),
            'lon': xr.DataArray(lon_2d, dims=('rlat', 'rlon'),
                attrs={'units': 'degrees_east', 'standard_name': 'longitude'}),
        },
        coords={
            'time': xr.DataArray(time_arr,
                attrs={'long_name': 'time', 'units': 'fractional calendar year'}),
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


# =============================================================================
# Firn ice content — parallel workers
# =============================================================================

ICE_DENSITY = 917.0  # kg m⁻³

_FIC_W = {}   # shared state initialised once per worker process


def _fic_worker_init(n_time, input_dir_str):
    _FIC_W['n_time']    = n_time
    _FIC_W['input_dir'] = Path(input_dir_str)


def _fic_worker(nn_ri_rj):
    """Process one per-column 2D file. Returns (ri, rj, nt, fic) or None."""
    nn, ri, rj = nn_ri_rj
    n_time     = _FIC_W['n_time']
    fname      = _FIC_W['input_dir'] / f'FGRN055_era055_2D_{nn + 1}.nc'
    if not fname.exists():
        return None
    try:
        ds    = xr.open_dataset(fname, decode_times=False)
        dens  = ds['dens'].values
        depth = ds['depth'].values
        ds.close()
        nt  = min(dens.shape[1], n_time)
        fic = _compute_firn_ice_content(dens[:, :nt], depth[:, :nt])
        return (ri, rj, nt, fic)
    except Exception as exc:
        print(f'  WARNING col {nn + 1}: {exc}', flush=True)
        return None


# =============================================================================
# Firn ice content — core computation
# =============================================================================

def _compute_firn_ice_content(dens, depth):
    """
    Integrated ice content (kg m⁻²) in the firn column at each timestep.

    Sums ρ·dz over layers above the firn-ice transition (continuous basal ice
    body) that have density >= ICE_DENSITY.

    Parameters
    ----------
    dens  : (n_layers, n_time) float array — layer density (kg m⁻³)
    depth : (n_layers, n_time) float array — layer centre depth (m)

    Returns
    -------
    (n_time,) float32 — kg m⁻²; NaN where no active layers; 0 where no firn ice
    """
    FILL  = 9.9e35
    dens  = np.where(np.abs(dens)  > FILL, np.nan, dens ).astype(np.float64)
    depth = np.where(np.abs(depth) > FILL, np.nan, depth).astype(np.float64)

    n_layers, n_time = dens.shape

    dz = np.empty_like(depth)
    dz[1:-1] = (depth[:-2] - depth[2:]) / 2
    dz[0]    = depth[0]  - depth[1]
    dz[-1]   = depth[-2] - depth[-1]

    is_ice  = (dens >= ICE_DENSITY) & np.isfinite(dens)
    cum_ice = np.cumprod(is_ice.astype(np.int8), axis=0)
    n_ice_base       = cum_ice.sum(axis=0)
    layer_idx        = np.arange(n_layers)[:, np.newaxis]
    above_transition = layer_idx >= n_ice_base[np.newaxis, :]
    is_firn_ice      = above_transition & is_ice & np.isfinite(dz) & (dz > 0)

    ice_content = np.nansum(dens * dz * is_firn_ice, axis=0).astype(np.float32)
    # NaN = no firn layers (fully glaciated column or no active column data)
    # 0   = firn present but no embedded ice lenses
    # >0  = firn with embedded ice content
    n_firn = np.sum((dens < ICE_DENSITY) & np.isfinite(dens), axis=0)
    return np.where(n_firn > 0, ice_content, np.nan).astype(np.float32)


# =============================================================================
# Firn ice content — public interface
# =============================================================================

def compute_firn_thickness_column(nc_file):
    """
    Compute total firn thickness for one per-column 2D output file.

    Returns
    -------
    (n_time,) float32 — firn thickness in metres; NaN where no firn present
    """
    FILL  = 9.9e35
    nc_file = Path(nc_file)
    ds      = xr.open_dataset(nc_file, decode_times=False)
    dens    = ds['dens'].values.astype(np.float64)
    depth   = ds['depth'].values.astype(np.float64)
    ds.close()

    dens  = np.where(np.abs(dens)  > FILL, np.nan, dens)
    depth = np.where(np.abs(depth) > FILL, np.nan, depth)

    dz = np.empty_like(depth)
    dz[1:-1] = (depth[:-2] - depth[2:]) / 2
    dz[0]    = depth[0]  - depth[1]
    dz[-1]   = depth[-2] - depth[-1]

    is_firn = (dens < ICE_DENSITY) & np.isfinite(dens) & np.isfinite(dz) & (dz > 0)
    thickness = np.nansum(np.where(is_firn, dz, 0.0), axis=0).astype(np.float32)
    n_firn    = np.sum(is_firn, axis=0)
    return np.where(n_firn > 0, thickness, np.nan).astype(np.float32)


def compute_firn_ice_content_column(nc_file):
    """
    Compute firn ice content for one per-column 2D output file.

    Returns
    -------
    (n_time,) float32 — integrated firn ice content in kg m⁻²
    """
    nc_file = Path(nc_file)
    ds      = xr.open_dataset(nc_file, decode_times=False)
    dens    = ds['dens'].values
    depth   = ds['depth'].values
    ds.close()
    return _compute_firn_ice_content(dens, depth)


def create_gridded_firn_ice_content(
    output_file,
    output_dir=None,
    file_limit=None,
    n_workers=1,
    verbose=True,
):
    """
    Create a gridded (time, rlat, rlon) NetCDF of integrated firn ice content.

    Loops over per-column 2D files, computes firn ice content for each column,
    and assembles the results onto the FGRN055 grid.
    Batch operation (~58 265 files) — run as a job.

    Output variables: firn_ice_content, firn_ice_content_mean, firn_ice_content_std,
                      lat, lon, rotated_pole, y_FDM, x_FDM.
    """
    if output_dir is None:
        output_dir = config.PROCESSED_DIR
    output_path = Path(output_dir) / output_file

    masks    = xr.open_dataset(config.MASKS_FILE)
    n_rlat   = len(masks['rlat'])
    n_rlon   = len(masks['rlon'])
    lat_2d   = masks['lat'].values
    lon_2d   = masks['lon'].values
    masks.close()

    grid_file = (Path('/home/nld4814/perm/code/IMAU-FDM/reference/FGRN055')
                 / 'FGRN055_grid.nc')
    with xr.open_dataset(grid_file) as gds:
        rlat_deg           = gds['rlat'].values
        rlon_deg           = gds['rlon'].values
        rotated_pole_attrs = dict(gds['rotated_pole'].attrs)

    ref      = xr.open_dataset(config.Z830_FILE, decode_times=False)
    time_arr = ref['time'].values
    n_time   = len(time_arr)
    ref.close()

    grid     = pd.read_csv(config.POINTLIST_FILE, header=None, sep=',')
    n_cols   = len(grid)
    rlat_col = grid.iloc[:, 5].astype(int).values
    rlon_col = grid.iloc[:, 6].astype(int).values

    if file_limit is not None:
        n_cols = min(n_cols, file_limit)
        if verbose:
            print(f'  file_limit={file_limit}: processing first {n_cols} columns only')

    fic_out   = np.full((n_time, n_rlat, n_rlon), np.nan, dtype=np.float32)
    input_dir = Path(config.SCRATCH_DIR) / config.PROJECT_NAME / 'output'
    n_processed = 0

    if n_workers > 1:
        from concurrent.futures import ProcessPoolExecutor
        tasks = [(nn, int(rlat_col[nn]), int(rlon_col[nn])) for nn in range(n_cols)]
        if verbose:
            print(f'  parallel mode: {n_workers} workers, {n_cols} columns', flush=True)
        with ProcessPoolExecutor(
            max_workers=n_workers,
            initializer=_fic_worker_init,
            initargs=(n_time, str(input_dir)),
        ) as pool:
            for result in pool.map(_fic_worker, tasks, chunksize=100):
                if result is not None:
                    ri, rj, nt, fic = result
                    fic_out[:nt, ri, rj] = fic
                    n_processed += 1
                    if verbose and n_processed % 1000 == 0:
                        print(f'  processed {n_processed} / {n_cols}', flush=True)
    else:
        for nn in range(n_cols):
            fname = input_dir / f'FGRN055_era055_2D_{nn + 1}.nc'
            if not fname.exists():
                continue
            ri, rj = rlat_col[nn], rlon_col[nn]
            try:
                ds    = xr.open_dataset(fname, decode_times=False)
                dens  = ds['dens'].values
                depth = ds['depth'].values
                ds.close()
                nt  = min(dens.shape[1], n_time)
                fic = _compute_firn_ice_content(dens[:, :nt], depth[:, :nt])
                fic_out[:nt, ri, rj] = fic
            except Exception as exc:
                if verbose:
                    print(f'  WARNING: skipping file {nn + 1}: {exc}')
                continue
            n_processed += 1
            if verbose and n_processed % 1000 == 0:
                print(f'  processed {n_processed} / {n_cols}')

    if verbose:
        print(f'  done — {n_processed} columns written', flush=True)

    out = xr.Dataset(
        {
            'firn_ice_content': xr.DataArray(
                fic_out, dims=('time', 'rlat', 'rlon'),
                attrs={'long_name'  : 'Integrated ice content in the firn column',
                       'units'      : 'kg m-2',
                       'description': (
                           'Sum of (density × layer_thickness) for layers above '
                           'the firn-ice transition with density >= 917 kg m-3.  '
                           'Divide by 1000 to convert to m w.e.  '
                           'NaN: no active layers (ablation zone).  '
                           'Zero: firn column with no embedded ice.'),
                       'grid_mapping': 'rotated_pole',
                       'coordinates' : 'lon lat'}),
            'firn_ice_content_mean': xr.DataArray(
                np.nanmean(fic_out, axis=0).astype(np.float32), dims=('rlat', 'rlon'),
                attrs={'long_name': 'Time-mean integrated firn ice content', 'units': 'kg m-2',
                       'grid_mapping': 'rotated_pole', 'coordinates': 'lon lat'}),
            'firn_ice_content_std': xr.DataArray(
                np.nanstd(fic_out, axis=0).astype(np.float32), dims=('rlat', 'rlon'),
                attrs={'long_name': 'Temporal std of integrated firn ice content', 'units': 'kg m-2',
                       'grid_mapping': 'rotated_pole', 'coordinates': 'lon lat'}),
            'lat': xr.DataArray(lat_2d, dims=('rlat', 'rlon'),
                attrs={'long_name': 'latitude', 'standard_name': 'latitude',
                       'units': 'degrees_north'}),
            'lon': xr.DataArray(lon_2d, dims=('rlat', 'rlon'),
                attrs={'long_name': 'longitude', 'standard_name': 'longitude',
                       'units': 'degrees_east'}),
            'rotated_pole': xr.DataArray(np.int32(0), attrs=rotated_pole_attrs),
            'y_FDM': xr.DataArray(np.arange(n_rlat, dtype=np.int32), dims=('rlat',),
                attrs={'long_name': 'row index in IMAU-FDM grid (0-based)', 'units': '1'}),
            'x_FDM': xr.DataArray(np.arange(n_rlon, dtype=np.int32), dims=('rlon',),
                attrs={'long_name': 'column index in IMAU-FDM grid (0-based)', 'units': '1'}),
        },
        coords={
            'time': xr.DataArray(time_arr,
                attrs={'long_name': 'Time in fractional calendar years', 'units': 'years'}),
            'rlat': xr.DataArray(rlat_deg, dims=('rlat',),
                attrs={'axis': 'Y', 'long_name': 'latitude in rotated pole grid',
                       'standard_name': 'grid_latitude', 'units': 'degrees'}),
            'rlon': xr.DataArray(rlon_deg, dims=('rlon',),
                attrs={'axis': 'X', 'long_name': 'longitude in rotated pole grid',
                       'standard_name': 'grid_longitude', 'units': 'degrees'}),
        },
        attrs={
            'title'      : 'IMAU-FDM integrated firn ice content',
            'source'     : 'IMAU-FDM version 1.2+',
            'domain'     : 'FGRN055',
            'institution': 'IMAU, Utrecht University',
            'history'    : f'Created on {datetime.utcnow().isoformat()}',
            'Conventions': 'CF-1.6',
        },
    )
    out.to_netcdf(output_path)
    if verbose:
        print(f'Saved: {output_path}')
    return output_path


# =============================================================================
# CLI entry point
# =============================================================================

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Create gridded firn analysis NetCDF')
    parser.add_argument('analysis', choices=['memory', 'ice_content'],
                        help='memory: firn memory at target_depth and z830; '
                             'ice_content: integrated firn ice content')
    parser.add_argument('--output',       default=None,
                        help='Output filename (default depends on analysis)')
    parser.add_argument('--output-dir',   default=None,
                        help='Output directory (default: config.PROCESSED_DIR)')
    parser.add_argument('--target-depth', type=float, default=10.0,
                        help='[memory] Depth in metres (default: 10.0)')
    parser.add_argument('--file-limit',   type=int, default=None,
                        help='Stop after N columns — for test runs')
    parser.add_argument('--n-workers',    type=int, default=1,
                        help='Parallel worker processes (default: 1)')
    args = parser.parse_args()

    if args.analysis == 'memory':
        create_gridded_firn_memory(
            output_file  = args.output or 'FDM_firn_memory_FGRN055_1939-2023_2D.nc',
            output_dir   = args.output_dir,
            target_depth = args.target_depth,
            file_limit   = args.file_limit,
            n_workers    = args.n_workers,
            verbose      = True,
        )
    else:
        create_gridded_firn_ice_content(
            output_file = args.output or 'FDM_firn_ice_content_FGRN055_1939-2023_2D.nc',
            output_dir  = args.output_dir,
            file_limit  = args.file_limit,
            n_workers   = args.n_workers,
            verbose     = True,
        )
