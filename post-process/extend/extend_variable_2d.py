#!/usr/bin/env python3
"""
Extend existing 2D-derived gridded variable files to 2025.

Reads merged per-column 2D files (1939-2025, from extend_pointfiles.py output),
applies the same profile reduction as make_2d_files.py (find_depth_at_threshold
or get_value_at_depth), returns only the new timesteps, and appends them to the
existing 1939-2023 gridded file.

Note: Unlike 1D files, no resampling or detrending is needed — each timestep in
the per-column 2D file corresponds directly to one output timestep (30-day).  The
extension portion is simply result[orig_ntime_output:], so only the 2024-2025
per-column data is actually needed (the full merged file is used for convenience).

Usage:
    python extend_variable_2d.py --output-var z830 --var dens --threshold 830
    python extend_variable_2d.py --output-var z550 --var dens --threshold 550
    python extend_variable_2d.py --output-var T10m_2d --var temp --depth 10
"""

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).parent.parent / 'create_1D_2D_2Ddetail_files'))
from utils import (
    load_pointlist,
    load_mask,
    get_grid_dimensions,
    create_output_dataset,
)
from make_2d_files import find_depth_at_threshold, get_value_at_depth
import extend_variable_config as cfg


# ---------------------------------------------------------------------------
# Per-point worker (runs in subprocess)
# ---------------------------------------------------------------------------

def _process_point(args):
    """Read one merged 2D per-column file and return only the extension portion."""
    (point_id, rlat_idx, rlon_idx, var_name, threshold, target_depth,
     secondary_var, orig_ntime_output) = args

    filepath = cfg.EXT_INPUT_DIR / cfg.EXT_2D_FILENAME_PATTERN.format(point_id=point_id)
    if not filepath.exists():
        return (rlat_idx, rlon_idx, None)

    try:
        with xr.open_dataset(filepath) as ds:
            if var_name not in ds:
                return (rlat_idx, rlon_idx, None)

            values = ds[var_name].values  # (layer, ind_t)

            if secondary_var and secondary_var in ds:
                depth = ds[secondary_var].values
            elif 'depth' in ds:
                depth = ds['depth'].values
            else:
                depth = np.arange(values.shape[0]) * 0.15  # 0.15 m/layer: fallback from make_2d_files.py; 2D files should always have 'depth'

            if threshold is not None:
                result = find_depth_at_threshold(depth, values, threshold)
            elif target_depth is not None:
                result = get_value_at_depth(depth, values, target_depth)
            else:
                result = values[0, :]
    except Exception:
        return (rlat_idx, rlon_idx, None)

    if len(result) <= orig_ntime_output:
        return (rlat_idx, rlon_idx, None)

    return (rlat_idx, rlon_idx, result[orig_ntime_output:])


# ---------------------------------------------------------------------------
# Per-variable orchestration
# ---------------------------------------------------------------------------

def extend_variable(output_var, var_name, threshold=None, target_depth=None,
                    secondary_var='depth', workers=None, max_points=None, verbose=True):
    """
    Extend one 2D-derived gridded variable file from 1939-2023 to 1939-2025.

    Parameters
    ----------
    output_var    : str   Output variable name (e.g. 'z830', 'z550')
    var_name      : str   Source variable in per-column file (e.g. 'dens', 'temp')
    threshold     : float Density threshold for find_depth_at_threshold
    target_depth  : float Target depth in metres for get_value_at_depth
    secondary_var : str   Depth variable name in per-column file (default 'depth')
    workers       : int   Parallel workers
    max_points    : int   Cap point count for testing
    verbose       : bool

    Returns
    -------
    Path   path to the written output file
    """
    workers = workers or cfg.NUM_WORKERS

    # Find existing gridded file — infer timestep from name
    candidates = sorted(cfg.ORIG_OUTPUT_DIR.glob(
        f'FDM_{output_var}_{cfg.DOMAIN}_{cfg.ORIG_DATE_TAG}_*.nc'))
    if not candidates:
        raise FileNotFoundError(
            f"No existing gridded file for '{output_var}' in {cfg.ORIG_OUTPUT_DIR}")
    orig_path = candidates[0]
    stem = orig_path.stem   # e.g. FDM_z830_FGRN055_1939-2023_30day
    timestep = stem.split('_')[-1]

    if verbose:
        print(f"\n  Output var: {output_var}  |  source var: {var_name}  "
              f"|  timestep: {timestep}  |  workers: {workers}")
        if threshold is not None:
            print(f"  Mode: depth at threshold={threshold}")
        elif target_depth is not None:
            print(f"  Mode: value at depth={target_depth}m")

    # Read orig_ntime_output and time values from existing file
    with xr.open_dataset(orig_path, decode_times=False) as ds_orig:
        orig_ntime_output = len(ds_orig.time)
        orig_time_values  = ds_orig.time.values
        orig_attrs        = dict(ds_orig.attrs)

    # Determine total timesteps from a sample per-column 2D file
    pointlist = load_pointlist(cfg.POINTLIST_FILE)
    total_ind_t = None
    for _, row in pointlist.iterrows():
        fp = cfg.EXT_INPUT_DIR / cfg.EXT_2D_FILENAME_PATTERN.format(
            point_id=int(row['point_id']))
        if fp.exists():
            with xr.open_dataset(fp) as ds_s:
                total_ind_t = ds_s.sizes.get('ind_t', None)
            break
    if total_ind_t is None:
        raise RuntimeError("Could not determine total timesteps from merged 2D files")

    new_ntime = total_ind_t - orig_ntime_output
    if new_ntime <= 0:
        raise RuntimeError(
            f"No new timesteps: merged file has {total_ind_t} steps, "
            f"existing gridded file has {orig_ntime_output}")

    # New time values — extrapolate from the existing uniform step
    dt = float(orig_time_values[1] - orig_time_values[0])
    new_time_values = orig_time_values[-1] + dt * np.arange(1, new_ntime + 1)

    if verbose:
        print(f"  Original output steps: {orig_ntime_output}  |  "
              f"New output steps (2024-2025): {new_ntime}")
        print(f"  New time range: {new_time_values[0]:.4f} – {new_time_values[-1]:.4f}")

    # Grid dimensions
    n_rlat, n_rlon = get_grid_dimensions(cfg.MASK_FILE)
    output_new = np.full((new_ntime, n_rlat, n_rlon), np.nan, dtype=np.float32)

    if max_points is not None:
        pointlist = pointlist.head(max_points)

    args_list = [
        (int(row['point_id']), int(row['rlat_idx']), int(row['rlon_idx']),
         var_name, threshold, target_depth, secondary_var, orig_ntime_output)
        for _, row in pointlist.iterrows()
    ]

    completed = failed = 0
    total = len(args_list)
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(_process_point, a): a for a in args_list}
        for future in as_completed(futures):
            rlat_idx, rlon_idx, data = future.result()
            if data is not None:
                n = min(len(data), new_ntime)
                output_new[:n, rlat_idx, rlon_idx] = data[:n]
                completed += 1
            else:
                failed += 1
            done = completed + failed
            if verbose and done % 5000 == 0:
                print(f"    {done}/{total} ({100*done/total:.1f}%)")

    if verbose:
        print(f"  Completed: {completed}  Failed/missing: {failed}")

    # Load existing gridded data and concatenate
    with xr.open_dataset(orig_path, decode_times=False) as ds_orig:
        orig_data = ds_orig[output_var].values  # (orig_ntime_output, n_rlat, n_rlon)
        var_metadata = {
            'long_name': ds_orig[output_var].attrs.get('long_name', output_var),
            'units':     ds_orig[output_var].attrs.get('units', ''),
        }

    combined_data = np.concatenate([orig_data, output_new], axis=0)
    combined_time = np.concatenate([orig_time_values, new_time_values])

    mask_ds = load_mask(cfg.MASK_FILE)
    ds_out = create_output_dataset(
        var_name=output_var,
        data=combined_data,
        time_values=combined_time,
        mask_ds=mask_ds,
        var_metadata=var_metadata,
        grid_file=cfg.GRID_FILE,
        detrended=False,
        timestep=timestep,
    )

    ds_out.attrs['history'] = (
        orig_attrs.get('history', '') +
        f'\n{datetime.now().isoformat()}: Extended to {cfg.EXT_MODEL_END.date()} '
        f'(extend_variable_2d.py)'
    )

    cfg.EXT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out_filename = f'FDM_{output_var}_{cfg.DOMAIN}_{cfg.EXT_DATE_TAG}_{timestep}.nc'
    out_path = cfg.EXT_OUTPUT_DIR / out_filename
    ds_out.to_netcdf(out_path, format='NETCDF4', engine='netcdf4',
                     encoding={output_var: {'zlib': True, 'complevel': 4}})

    if verbose:
        print(f"  Written: {out_path}")

    return out_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Extend 2D-derived gridded variable files from 1939-2023 to 1939-2025.',
        epilog="""
Examples:
  python extend_variable_2d.py --output-var z830 --var dens --threshold 830
  python extend_variable_2d.py --output-var z550 --var dens --threshold 550
  python extend_variable_2d.py --output-var T10m_2d --var temp --depth 10
        """,
    )
    parser.add_argument('--output-var', required=True,
                        help='Output variable name (e.g. z830, z550)')
    parser.add_argument('--var', '-v', required=True,
                        help='Source variable in per-column 2D file (e.g. dens, temp)')
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--threshold', '-t', type=float,
                      help='Find depth where source variable reaches this threshold')
    mode.add_argument('--depth', '-d', type=float,
                      help='Extract source variable value at this depth (metres)')
    parser.add_argument('--secondary-var', default='depth',
                        help='Depth variable name in per-column file (default: depth)')
    parser.add_argument('--workers', '-w', type=int, default=None)
    parser.add_argument('--max-points', type=int, default=None)
    parser.add_argument('--quiet', '-q', action='store_true')
    args = parser.parse_args()

    try:
        extend_variable(
            output_var=args.output_var,
            var_name=args.var,
            threshold=args.threshold,
            target_depth=args.depth,
            secondary_var=args.secondary_var,
            workers=args.workers,
            max_points=args.max_points,
            verbose=not args.quiet,
        )
    except Exception as e:
        print(f"ERROR: {e}")
        return 1

    if not args.quiet:
        print(f"\nDone. Output written to {cfg.EXT_OUTPUT_DIR}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
