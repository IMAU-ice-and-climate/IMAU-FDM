"""
Ice lens / ice slab detection from IMAU-FDM 2Ddetail subsurface profiles.

An ice lens (or slab) is identified in a firn column wherever a layer reaches
a density >= `min_ice_dens` (default 900 kg/m³) while at least one layer
*below* it remains below that threshold.  This "ice-on-top-of-firn" structure
is the configuration that impedes downward water percolation: meltwater
accumulates above the ice and in reality runs off laterally, but IMAU-FDM's
bucket scheme may route it to the base instead.

Public interface
----------------
detect_ice_lenses_in_column(nc_file)
    Process one 2Ddetail per-column file.  Returns an xr.Dataset with three
    time series: `has_ice_lens`, `depth_top_lens`, `ice_lens_thickness`.

create_gridded_ice_lens_file(output_file)
    Loop over all per-column 2Ddetail files, compute the diagnostics for
    each, assemble the results onto the FGRN055 grid, and save a NetCDF.
    This is a batch operation (58 265 files) and should be run as a job.
"""

import sys
from pathlib import Path

import numpy as np
import xarray as xr
import pandas as pd
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent))
import config


# =============================================================================
# Parallel worker (module-level so it is picklable)
# =============================================================================

_W = {}   # shared state initialised once per worker via ProcessPoolExecutor initializer


def _il_worker_init(n_time, min_ice_dens, require_firn_above, input_dir_str):
    _W['n_time']            = n_time
    _W['min_ice_dens']      = min_ice_dens
    _W['require_firn_above'] = require_firn_above
    _W['input_dir']         = Path(input_dir_str)


def _il_worker(nn_ri_rj):
    """Process one 2Ddetail file. Returns (ri, rj, nt, has_lens, depth_top, thickness) or None."""
    nn, ri, rj         = nn_ri_rj
    n_time             = _W['n_time']
    min_ice_dens       = _W['min_ice_dens']
    require_firn_above = _W['require_firn_above']
    fname              = _W['input_dir'] / f'FGRN055_era055_2Ddetail_{nn + 1}.nc'
    if not fname.exists():
        return None
    try:
        ds   = xr.open_dataset(fname)
        dens = ds['dens'].values
        dep  = ds['depth'].values
        dz   = ds['dz'].values
        ds.close()
        nt   = min(dens.shape[1], n_time)
        diag = _compute_ice_lens_diagnostics(
            dens[:, :nt], dep, dz, min_ice_dens,
            require_firn_above=require_firn_above)
        return (ri, rj, nt,
                diag['has_ice_lens'].astype(np.int8),
                diag['depth_top_lens'],
                diag['ice_lens_thickness'])
    except Exception as exc:
        print(f'  WARNING col {nn + 1}: {exc}', flush=True)
        return None


# =============================================================================
# Core per-column detection (fully vectorised over layers × time)
# =============================================================================

def _compute_ice_lens_diagnostics(dens, depth, dz, min_ice_dens=900.0,
                                   require_firn_above=False):
    """
    Detect ice-on-firn layers in a firn column.

    Parameters
    ----------
    dens  : (n_layers, n_time) float array — layer densities (kg/m³)
    depth : (n_layers,) float array        — layer centre depths (m)
    dz    : (n_layers,) float array        — layer thicknesses (m)
    min_ice_dens : float
        Density threshold classifying a layer as ice (default 900 kg/m³).
    require_firn_above : bool
        If False (default): any ice layer with firn below qualifies, including
        superimposed ice at the very surface.
        If True: also require at least one firn layer *above* the ice layer,
        selecting only ice that is embedded within the firn column (classic
        ice-lens / slab geometry: firn → ice → firn).

    Returns
    -------
    dict with keys:
        has_ice_lens      : (n_time,) bool    — any qualifying ice layer present?
        depth_top_lens    : (n_time,) float   — depth (m) of shallowest ice lens;
                                                NaN when has_ice_lens is False
        ice_lens_thickness: (n_time,) float   — total thickness (m) of all
                                                qualifying ice layers; 0.0 when
                                                has_ice_lens is False
    """
    ice  = dens >= min_ice_dens   # (n_layers, n_time) bool
    firn = ~ice                   # (n_layers, n_time) bool

    # firn_below[i, t] = True if any layer j > i has dens[j, t] < min_ice_dens
    # = total firn count − cumulative firn count through layer i
    cumsum_firn = np.cumsum(firn, axis=0)          # layers are axis-0
    total_firn  = firn.sum(axis=0, keepdims=True)  # (1, n_time)
    firn_below  = (total_firn - cumsum_firn) > 0   # (n_layers, n_time)

    # ice_lens[i, t]: this layer is ice AND has permeable firn below it
    ice_lens = ice & firn_below   # (n_layers, n_time)

    if require_firn_above:
        # firn_above[i, t] = True if any layer j < i has dens[j, t] < min_ice_dens
        # = cumulative firn count through layer i-1 (exclusive)
        firn_above = (cumsum_firn - firn) > 0      # (n_layers, n_time)
        ice_lens   = ice_lens & firn_above

    # --- scalar outputs per timestep ---
    has_ice_lens = ice_lens.any(axis=0)   # (n_time,) bool

    # depth of shallowest qualifying layer (argmax returns 0 if no True → mask)
    first_idx        = ice_lens.argmax(axis=0)             # (n_time,)
    depth_top_lens   = np.where(has_ice_lens,
                                depth[first_idx],
                                np.nan).astype(np.float32)

    # total thickness of all qualifying ice layers
    ice_lens_thickness = (ice_lens * dz[:, np.newaxis]).sum(axis=0)  # (n_time,)
    ice_lens_thickness = np.where(has_ice_lens,
                                  ice_lens_thickness,
                                  0.0).astype(np.float32)

    return {
        'has_ice_lens'      : has_ice_lens,
        'depth_top_lens'    : depth_top_lens,
        'ice_lens_thickness': ice_lens_thickness,
    }


# =============================================================================
# Single-column convenience wrapper
# =============================================================================

def detect_ice_lenses_in_column(nc_file, min_ice_dens=900.0,
                                require_firn_above=False, time_coord=None):
    """
    Detect ice lenses in one 2Ddetail per-column output file.

    Parameters
    ----------
    nc_file : Path or str
        Path to a ``*_2Ddetail_N.nc`` output file.
    min_ice_dens : float
        Density threshold for ice classification (default 900 kg/m³).
    require_firn_above : bool
        If True, only flag ice layers that also have firn above them (embedded
        ice lenses).  If False (default), also flags superimposed ice at the
        very surface.  See ``_compute_ice_lens_diagnostics`` for details.
    time_coord : array-like, optional
        Fractional-year time values to attach as a coordinate.  If None,
        the function uses the reference time axis from config.T10M_FILE
        (same 10-day spacing, same length as the 2Ddetail ind_t dimension).

    Returns
    -------
    xr.Dataset
        Coordinates: time (fractional year, length ind_t)
        Variables:
          has_ice_lens       (time,) bool   — any ice lens in column?
          depth_top_lens     (time,) float  — depth (m) of shallowest lens
          ice_lens_thickness (time,) float  — total ice lens thickness (m)
        Attributes copied from the source file.
    """
    nc_file = Path(nc_file)
    ds = xr.open_dataset(nc_file)

    dens  = ds['dens'].values   # (n_layers, n_time)
    depth = ds['depth'].values  # (n_layers,)
    dz    = ds['dz'].values     # (n_layers,)
    src_attrs = dict(ds.attrs)
    ds.close()

    diag = _compute_ice_lens_diagnostics(dens, depth, dz, min_ice_dens,
                                         require_firn_above=require_firn_above)

    # Build time coordinate
    if time_coord is None:
        ref = xr.open_dataset(config.T10M_FILE, decode_times=False)
        time_coord = ref['time'].values[:dens.shape[1]]
        ref.close()

    out = xr.Dataset(
        {
            'has_ice_lens': xr.DataArray(
                diag['has_ice_lens'].astype(np.int8),
                dims='time',
                attrs={
                    'long_name': 'Ice lens present (ice≥900 on top of firn<900)',
                    'units'    : '1',
                    'flag_values': '0 1',
                    'flag_meanings': 'no_ice_lens ice_lens_present',
                },
            ),
            'depth_top_lens': xr.DataArray(
                diag['depth_top_lens'],
                dims='time',
                attrs={
                    'long_name': 'Depth of shallowest ice lens',
                    'units'    : 'm',
                    'note'     : 'NaN when no ice lens is present',
                },
            ),
            'ice_lens_thickness': xr.DataArray(
                diag['ice_lens_thickness'],
                dims='time',
                attrs={
                    'long_name': 'Total thickness of ice-on-firn layers',
                    'units'    : 'm',
                },
            ),
        },
        coords={'time': time_coord},
    )
    out.attrs.update(src_attrs)
    out.attrs['ice_lens_min_density']   = min_ice_dens
    out.attrs['ice_lens_require_firn_above'] = int(require_firn_above)
    out.attrs['source_file'] = str(nc_file)
    return out


# =============================================================================
# Gridded post-processing (batch — run as a job)
# =============================================================================

def create_gridded_ice_lens_file(
    output_file,
    output_dir=None,
    min_ice_dens=900.0,
    require_firn_above=False,
    file_limit=None,
    n_workers=1,
    verbose=True,
):
    """
    Create a gridded (time, rlat, rlon) NetCDF of ice lens diagnostics.

    Loops over all per-column 2Ddetail files listed in the FGRN055 pointlist,
    computes ice lens diagnostics for each, and places them on the 2-D grid.

    This is a batch operation (~58 265 files).  Run it as a job script, not
    interactively.

    Parameters
    ----------
    output_file : Path or str
        Name of the output NetCDF file.
    output_dir : Path or str, optional
        Directory for the output file.  Defaults to config.PROCESSED_DIR.
    min_ice_dens : float
        Density threshold for ice (default 900 kg/m³).
    file_limit : int, optional
        Stop after processing this many columns — useful for test runs.
    verbose : bool
        Print progress every 1000 columns.

    Output file variables
    ---------------------
    has_ice_lens       (time, rlat, rlon)  int8   — 1 if any lens, 0 otherwise
    depth_top_lens     (time, rlat, rlon)  float  — depth (m) of shallowest lens
    ice_lens_thickness (time, rlat, rlon)  float  — total ice lens thickness (m)
    lat / lon          (rlat, rlon)        float  — geographic coordinates
    """
    if output_dir is None:
        output_dir = config.PROCESSED_DIR
    output_path = Path(output_dir) / output_file

    # --- Load grid metadata ---
    masks    = xr.open_dataset(config.MASKS_FILE)
    n_rlat   = len(masks['rlat'])
    n_rlon   = len(masks['rlon'])
    lat_2d   = masks['lat'].values
    lon_2d   = masks['lon'].values
    masks.close()

    # --- Reference time axis (from existing 2Ddetail gridded file) ---
    ref      = xr.open_dataset(config.T10M_FILE, decode_times=False)
    time_arr = ref['time'].values   # (n_time,)
    n_time   = len(time_arr)
    ref.close()

    # --- Point list: maps file_num → (rlat_idx, rlon_idx) ---
    grid     = pd.read_csv(config.POINTLIST_FILE, header=None, sep=',')
    n_cols   = len(grid)
    rlat_col = grid.iloc[:, 5].astype(int).values
    rlon_col = grid.iloc[:, 6].astype(int).values

    if file_limit is not None:
        n_cols = min(n_cols, file_limit)
        if verbose:
            print(f'  file_limit={file_limit}: processing first {n_cols} columns only')

    # --- Allocate output arrays ---
    has_lens_out  = np.zeros((n_time, n_rlat, n_rlon), dtype=np.int8)
    depth_out     = np.full( (n_time, n_rlat, n_rlon), np.nan, dtype=np.float32)
    thick_out     = np.zeros((n_time, n_rlat, n_rlon), dtype=np.float32)

    output_dir_obj = Path(config.SCRATCH_DIR) / config.PROJECT_NAME / 'output'

    # --- Main loop ---
    n_processed = 0

    if n_workers > 1:
        from concurrent.futures import ProcessPoolExecutor
        tasks = [(nn, int(rlat_col[nn]), int(rlon_col[nn])) for nn in range(n_cols)]
        if verbose:
            print(f'  parallel mode: {n_workers} workers, {n_cols} columns', flush=True)
        with ProcessPoolExecutor(
            max_workers=n_workers,
            initializer=_il_worker_init,
            initargs=(n_time, min_ice_dens, require_firn_above, str(output_dir_obj)),
        ) as pool:
            for result in pool.map(_il_worker, tasks, chunksize=100):
                if result is not None:
                    ri, rj, nt, has_lens, depth_top, thickness = result
                    has_lens_out[:nt, ri, rj] = has_lens
                    depth_out   [:nt, ri, rj] = depth_top
                    thick_out   [:nt, ri, rj] = thickness
                    n_processed += 1
                    if verbose and n_processed % 1000 == 0:
                        print(f'  processed {n_processed} / {n_cols}', flush=True)
    else:
        for nn in range(n_cols):
            file_num = nn + 1
            fname = output_dir_obj / f'FGRN055_era055_2Ddetail_{file_num}.nc'
            if not fname.exists():
                continue
            ri = rlat_col[nn]
            rj = rlon_col[nn]
            try:
                ds   = xr.open_dataset(fname)
                dens = ds['dens'].values
                dep  = ds['depth'].values
                dz   = ds['dz'].values
                ds.close()
                nt = min(dens.shape[1], n_time)
                diag = _compute_ice_lens_diagnostics(
                    dens[:, :nt], dep, dz, min_ice_dens,
                    require_firn_above=require_firn_above)
                has_lens_out[:nt, ri, rj] = diag['has_ice_lens'].astype(np.int8)
                depth_out   [:nt, ri, rj] = diag['depth_top_lens']
                thick_out   [:nt, ri, rj] = diag['ice_lens_thickness']
            except Exception as exc:
                if verbose:
                    print(f'  WARNING: skipping file {file_num}: {exc}')
                continue
            n_processed += 1
            if verbose and n_processed % 1000 == 0:
                print(f'  processed {n_processed} / {n_cols}')

    if verbose:
        print(f'  done — {n_processed} columns written', flush=True)

    # --- Assemble and save ---
    out = xr.Dataset(
        {
            'has_ice_lens': xr.DataArray(
                has_lens_out, dims=('time', 'rlat', 'rlon'),
                attrs={'long_name': 'Ice lens present (1=yes)', 'units': '1'},
            ),
            'depth_top_lens': xr.DataArray(
                depth_out, dims=('time', 'rlat', 'rlon'),
                attrs={'long_name': 'Depth of shallowest ice lens', 'units': 'm'},
            ),
            'ice_lens_thickness': xr.DataArray(
                thick_out, dims=('time', 'rlat', 'rlon'),
                attrs={'long_name': 'Total ice lens thickness', 'units': 'm'},
            ),
            'lat': xr.DataArray(lat_2d, dims=('rlat', 'rlon'),
                                attrs={'units': 'degrees_north'}),
            'lon': xr.DataArray(lon_2d, dims=('rlat', 'rlon'),
                                attrs={'units': 'degrees_east'}),
        },
        coords={
            'time': xr.DataArray(time_arr,
                                 attrs={'long_name': 'time', 'units': 'fractional year'}),
            'rlat': np.arange(n_rlat, dtype=np.float32),
            'rlon': np.arange(n_rlon, dtype=np.float32),
        },
        attrs={
            'title'              : 'IMAU-FDM ice lens diagnostics (top 20 m)',
            'source'             : 'IMAU-FDM version 1.2+',
            'domain'             : 'FGRN055',
            'institution'        : 'IMAU, Utrecht University',
            'history'            : f'Created on {datetime.utcnow().isoformat()}',
            'Conventions'        : 'CF-1.6',
            'ice_lens_threshold' : f'{min_ice_dens} kg/m3',
            'description'        : (
                'Ice lens = layer with density >= threshold that has at least '
                'one layer with density < threshold below it (ice on top of '
                'permeable firn, blocking downward percolation).'
            ),
        },
    )

    out.to_netcdf(output_path)
    if verbose:
        print(f'Saved: {output_path}')
    return output_path


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Create gridded ice lens diagnostics NetCDF')
    parser.add_argument('--output', default='FDM_ice_lens_FGRN055_1939-2023_2Ddetail.nc',
                        help='Output filename')
    parser.add_argument('--output-dir', default=None,
                        help='Output directory (default: config.PROCESSED_DIR)')
    parser.add_argument('--min-ice-dens', type=float, default=900.0,
                        help='Density threshold for ice classification (default: 900.0 kg/m³)')
    parser.add_argument('--require-firn-above', action='store_true',
                        help='Only flag embedded lenses (firn→ice→firn); default includes surface ice')
    parser.add_argument('--file-limit', type=int, default=None,
                        help='Stop after N columns — useful for test runs')
    parser.add_argument('--n-workers', type=int, default=1,
                        help='Number of parallel worker processes (default: 1)')
    args = parser.parse_args()

    create_gridded_ice_lens_file(
        output_file       = args.output,
        output_dir        = args.output_dir,
        min_ice_dens      = args.min_ice_dens,
        require_firn_above= args.require_firn_above,
        file_limit        = args.file_limit,
        n_workers         = args.n_workers,
        verbose           = True,
    )
