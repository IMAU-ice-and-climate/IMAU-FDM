"""
Integrated ice content within the firn column from IMAU-FDM 2D monthly output.

The firn-ice transition is defined as the top of the continuous ice body that
extends uninterrupted from the column base (layer 0) upward — i.e. the first
layer from the bottom with density < 917 kg/m³.  Ice lenses and perched ice
slabs embedded in the firn above this transition are the target quantity.

Physical quantity
-----------------
For each column at each monthly timestep, the script finds the firn-ice
transition layer and sums

    firn_ice_content = Σ  ρ_k · dz_k    (kg m⁻²)

over all layers k that are (a) above the firn-ice transition and (b) have
density ≥ 917 kg m⁻³ (i.e. ice lenses / refreezing layers within the firn).

Divide by 1000 to convert to m w.e.

Public interface
----------------
compute_firn_ice_content_column(nc_file)
    Process one 2D per-column file.  Returns a (n_time,) float32 array.

create_gridded_firn_ice_content(output_file, ...)
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

_W = {}  # shared state initialised once per worker process


def _fic_worker_init(n_time, input_dir_str):
    _W['n_time']    = n_time
    _W['input_dir'] = Path(input_dir_str)


def _fic_worker(nn_ri_rj):
    """Process one per-column 2D file. Returns (ri, rj, nt, fic) or None."""
    nn, ri, rj = nn_ri_rj
    n_time     = _W['n_time']
    fname      = _W['input_dir'] / f'FGRN055_era055_2D_{nn + 1}.nc'
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
# Core per-column computation
# =============================================================================

ICE_DENSITY = 917.0  # kg m⁻³


def _compute_firn_ice_content(dens, depth):
    """
    Compute the integrated ice content within the firn column at each timestep.

    The firn-ice transition is the top of the *continuous* ice body that
    extends uninterrupted from the column base (layer 0) upward — equivalent
    to the first layer from the bottom with density < 917 kg m⁻³.  Only ice
    embedded above this transition (ice lenses, refreezing layers) contributes.

    Parameters
    ----------
    dens  : (n_layers, n_time) float array
        Layer density (kg m⁻³).  Layer 0 = deepest (oldest); highest active
        index ≈ surface.
    depth : (n_layers, n_time) float array
        Layer centre depth (m below surface).  Decreases from layer 0 to surface.

    Returns
    -------
    (n_time,) float32 array
        Integrated ice content in the firn (kg m⁻²).
        NaN where the column has no active layers.
        Zero where the column has active layers but no ice above the transition.
    """
    FILL = 9.9e35

    dens  = np.where(np.abs(dens)  > FILL, np.nan, dens ).astype(np.float64)
    depth = np.where(np.abs(depth) > FILL, np.nan, depth).astype(np.float64)

    n_layers, n_time = dens.shape

    # Layer thickness via central differences on depth.
    # depth[k] > depth[k+1] (deepest at k=0), so all dz values are positive.
    dz = np.empty_like(depth)
    dz[1:-1] = (depth[:-2] - depth[2:]) / 2  # central diff for interior layers
    dz[0]    = depth[0]  - depth[1]           # one-sided at column base
    dz[-1]   = depth[-2] - depth[-1]          # one-sided at surface

    # Continuous basal ice mask: cumulative AND from base (layer 0) upward.
    # cum_ice[k, t] = 1 iff every layer from 0 to k has density >= ICE_DENSITY.
    is_ice  = (dens >= ICE_DENSITY) & np.isfinite(dens)
    cum_ice = np.cumprod(is_ice.astype(np.int8), axis=0)  # (n_layers, n_time)

    # Number of uninterrupted ice layers from the base = top of basal ice + 1.
    # Layer indices 0 … n_ice_base-1 form the continuous basal ice body.
    # Layers at index >= n_ice_base are in the firn column.
    n_ice_base = cum_ice.sum(axis=0)  # (n_time,)

    # Boolean mask: layer is above the firn-ice transition
    layer_idx        = np.arange(n_layers)[:, np.newaxis]           # (n_layers, 1)
    above_transition = layer_idx >= n_ice_base[np.newaxis, :]       # (n_layers, n_time)

    # Firn ice: above the transition, dens >= ICE_DENSITY, valid dz
    is_firn_ice = above_transition & is_ice & np.isfinite(dz) & (dz > 0)

    # Integrated ice content (kg m⁻²)
    ice_content = np.nansum(dens * dz * is_firn_ice, axis=0).astype(np.float32)

    # Return NaN where the column has no active (non-NaN) layers at all
    n_valid = np.sum(np.isfinite(dens), axis=0)
    return np.where(n_valid > 0, ice_content, np.nan).astype(np.float32)


# =============================================================================
# Single-column convenience wrapper
# =============================================================================

def compute_firn_ice_content_column(nc_file):
    """
    Compute the firn ice content time series for one per-column 2D output file.

    Parameters
    ----------
    nc_file : Path or str
        Path to a ``FGRN055_era055_2D_N.nc`` per-column output file.

    Returns
    -------
    (n_time,) float32 array
        Integrated firn ice content in kg m⁻².
    """
    nc_file = Path(nc_file)
    ds      = xr.open_dataset(nc_file, decode_times=False)
    dens    = ds['dens'].values
    depth   = ds['depth'].values
    ds.close()
    return _compute_firn_ice_content(dens, depth)


# =============================================================================
# Gridded post-processing (batch — run as a job)
# =============================================================================

def create_gridded_firn_ice_content(
    output_file,
    output_dir=None,
    file_limit=None,
    n_workers=1,
    verbose=True,
):
    """
    Create a gridded (time, rlat, rlon) NetCDF of integrated firn ice content.

    Loops over per-column 2D files listed in the FGRN055 pointlist, computes
    firn ice content for each column, and assembles the results onto the 2-D
    model grid.

    This is a batch operation (~58 265 files).  Run it as a job script,
    not interactively.  For a quick test pass file_limit=500 and a test output
    filename.

    Parameters
    ----------
    output_file : str
        Filename of the output NetCDF file (placed in output_dir).
    output_dir : Path or str, optional
        Directory for the output file.  Defaults to config.PROCESSED_DIR.
    file_limit : int, optional
        If given, stop after processing this many columns.  Useful for testing.
    n_workers : int
        Number of parallel worker processes (default 1).
    verbose : bool
        Print progress every 1 000 columns (default True).

    Output variables
    ----------------
    firn_ice_content      (time, rlat, rlon)  kg m⁻²
    firn_ice_content_mean (rlat, rlon)         time-mean
    firn_ice_content_std  (rlat, rlon)         temporal std
    lat / lon             (rlat, rlon)         geographic coordinates
    rotated_pole          scalar               CF grid-mapping variable
    y_FDM / x_FDM         (rlat) / (rlon)     original 0-based FDM grid indices
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

    # --- Rotated pole coordinates and grid-mapping from grid reference file ---
    grid_file = (Path('/home/nld4814/perm/code/IMAU-FDM/reference/FGRN055')
                 / 'FGRN055_grid.nc')
    with xr.open_dataset(grid_file) as gds:
        rlat_deg        = gds['rlat'].values
        rlon_deg        = gds['rlon'].values
        rotated_pole_attrs = dict(gds['rotated_pole'].attrs)

    # --- Time axis from the gridded z830 file (same monthly steps as 2D files) ---
    ref      = xr.open_dataset(config.Z830_FILE, decode_times=False)
    time_arr = ref['time'].values   # fractional calendar year, shape (n_time,)
    n_time   = len(time_arr)
    ref.close()

    # --- Point list ---
    grid     = pd.read_csv(config.POINTLIST_FILE, header=None, sep=',')
    n_cols   = len(grid)
    rlat_col = grid.iloc[:, 5].astype(int).values
    rlon_col = grid.iloc[:, 6].astype(int).values

    if file_limit is not None:
        n_cols = min(n_cols, file_limit)
        if verbose:
            print(f'  file_limit={file_limit}: processing first {n_cols} columns only')

    # --- Allocate output array ---
    fic_out = np.full((n_time, n_rlat, n_rlon), np.nan, dtype=np.float32)

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
            file_num = nn + 1
            fname    = input_dir / f'FGRN055_era055_2D_{file_num}.nc'
            if not fname.exists():
                continue
            ri = rlat_col[nn]
            rj = rlon_col[nn]
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
                    print(f'  WARNING: skipping file {file_num}: {exc}')
                continue
            n_processed += 1
            if verbose and n_processed % 1000 == 0:
                print(f'  processed {n_processed} / {n_cols}')

    if verbose:
        print(f'  done — {n_processed} columns written', flush=True)

    # --- Summary statistics ---
    fic_mean = np.nanmean(fic_out, axis=0).astype(np.float32)
    fic_std  = np.nanstd( fic_out, axis=0).astype(np.float32)

    # --- Build and save dataset ---
    out = xr.Dataset(
        {
            'firn_ice_content': xr.DataArray(
                fic_out, dims=('time', 'rlat', 'rlon'),
                attrs={
                    'long_name'  : 'Integrated ice content in the firn column',
                    'units'      : 'kg m-2',
                    'description': (
                        'Sum of (density × layer_thickness) for all layers above '
                        'the firn-ice transition (top of the continuous basal ice '
                        'body) that have density >= 917 kg m-3.  Represents ice '
                        'lenses and refreezing layers embedded in the firn.  '
                        'Divide by 1000 to convert to m w.e.  '
                        'NaN in columns with no active layers (ablation zone).  '
                        'Zero in firn columns with no ice lenses.'
                    ),
                    'grid_mapping': 'rotated_pole',
                    'coordinates' : 'lon lat',
                },
            ),
            'firn_ice_content_mean': xr.DataArray(
                fic_mean, dims=('rlat', 'rlon'),
                attrs={
                    'long_name': 'Time-mean integrated firn ice content',
                    'units'    : 'kg m-2',
                    'grid_mapping': 'rotated_pole',
                    'coordinates' : 'lon lat',
                },
            ),
            'firn_ice_content_std': xr.DataArray(
                fic_std, dims=('rlat', 'rlon'),
                attrs={
                    'long_name': 'Temporal std of integrated firn ice content',
                    'units'    : 'kg m-2',
                    'grid_mapping': 'rotated_pole',
                    'coordinates' : 'lon lat',
                },
            ),
            'lat': xr.DataArray(
                lat_2d, dims=('rlat', 'rlon'),
                attrs={
                    'long_name'    : 'latitude',
                    'standard_name': 'latitude',
                    'units'        : 'degrees_north',
                },
            ),
            'lon': xr.DataArray(
                lon_2d, dims=('rlat', 'rlon'),
                attrs={
                    'long_name'    : 'longitude',
                    'standard_name': 'longitude',
                    'units'        : 'degrees_east',
                },
            ),
            'rotated_pole': xr.DataArray(
                np.int32(0), attrs=rotated_pole_attrs,
            ),
            'y_FDM': xr.DataArray(
                np.arange(n_rlat, dtype=np.int32), dims=('rlat',),
                attrs={
                    'long_name': 'row index in IMAU-FDM grid (0-based)',
                    'units'    : '1',
                },
            ),
            'x_FDM': xr.DataArray(
                np.arange(n_rlon, dtype=np.int32), dims=('rlon',),
                attrs={
                    'long_name': 'column index in IMAU-FDM grid (0-based)',
                    'units'    : '1',
                },
            ),
        },
        coords={
            'time': xr.DataArray(
                time_arr,
                attrs={
                    'long_name': 'Time in fractional calendar years',
                    'units'    : 'years',
                },
            ),
            'rlat': xr.DataArray(
                rlat_deg,
                attrs={
                    'axis'         : 'Y',
                    'long_name'    : 'latitude in rotated pole grid',
                    'standard_name': 'grid_latitude',
                    'units'        : 'degrees',
                },
            ),
            'rlon': xr.DataArray(
                rlon_deg,
                attrs={
                    'axis'         : 'X',
                    'long_name'    : 'longitude in rotated pole grid',
                    'standard_name': 'grid_longitude',
                    'units'        : 'degrees',
                },
            ),
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


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Create gridded integrated firn ice content NetCDF'
    )
    parser.add_argument(
        '--output',
        default='FDM_firn_ice_content_FGRN055_1939-2023_2D.nc',
        help='Output filename (default: FDM_firn_ice_content_FGRN055_1939-2023_2D.nc)',
    )
    parser.add_argument(
        '--output-dir', default=None,
        help='Output directory (default: config.PROCESSED_DIR)',
    )
    parser.add_argument(
        '--file-limit', type=int, default=None,
        help='Stop after N columns — useful for test runs',
    )
    parser.add_argument(
        '--n-workers', type=int, default=1,
        help='Number of parallel worker processes (default: 1)',
    )
    args = parser.parse_args()

    create_gridded_firn_ice_content(
        output_file = args.output,
        output_dir  = args.output_dir,
        file_limit  = args.file_limit,
        n_workers   = args.n_workers,
        verbose     = True,
    )
