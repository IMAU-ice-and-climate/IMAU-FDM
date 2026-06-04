#!/usr/bin/env python3
"""
Extend existing 2Ddetail-derived gridded variable files to 2025.

Reads merged per-column 2Ddetail files (1939-2025, from extend_pointfiles.py
output), applies the same profile reduction as make_2Ddetail_files.py (average/
sum over layer range, or get_value_at_depth), returns only the new timesteps, and
appends them to the existing 1939-2023 gridded file.

Note: No resampling or detrending is needed — each timestep in the per-column
2Ddetail file corresponds directly to one output timestep (10-day).  Only the
2024-2025 per-column data is actually used; the full merged file is read for
convenience and sliced at orig_ntime_output.

Usage:
    python extend_variable_2ddetail.py --output-var T10m --var temp --depth 10
    python extend_variable_2ddetail.py --output-var SSN  --var dens --depth-begin 0 --depth-end 0.5
    python extend_variable_2ddetail.py --output-var LWC_surf --var lwc --depth-begin 0 --depth-end 0.2
    python extend_variable_2ddetail.py --output-var rho_2m --var dens --depth-begin 0 --depth-end 2
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
from make_2Ddetail_files import (
    average_over_layers,
    sum_over_layers,
    get_value_at_depth,
)
import extend_variable_config as cfg


# ---------------------------------------------------------------------------
# Per-point worker (runs in subprocess)
# ---------------------------------------------------------------------------

def _process_point(args):
    """Read one merged 2Ddetail per-column file and return only the extension portion."""
    (point_id, rlat_idx, rlon_idx, var_name,
     z_begin, z_end, target_depth, layer_thickness, operation,
     orig_ntime_output) = args

    filepath = cfg.EXT_INPUT_DIR / cfg.EXT_2DDETAIL_FILENAME_PATTERN.format(point_id=point_id)
    if not filepath.exists():
        return (rlat_idx, rlon_idx, None)

    try:
        with xr.open_dataset(filepath) as ds:
            if var_name not in ds:
                return (rlat_idx, rlon_idx, None)

            values = ds[var_name].values  # (layer, ind_t)

            if 'depth' in ds:
                depth = ds['depth'].values
            else:
                n_layers = values.shape[0]
                depth = np.arange(n_layers) * layer_thickness + layer_thickness / 2

            if target_depth is not None:
                result = get_value_at_depth(depth, values, target_depth, layer_thickness)
            elif z_begin is not None and z_end is not None:
                if operation == 'sum':
                    result = sum_over_layers(values, z_begin, z_end)
                else:
                    result = average_over_layers(values, z_begin, z_end)
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

def extend_variable(output_var, var_name,
                    z_begin=None, z_end=None,
                    target_depth=None, layer_thickness=0.04, operation='average',
                    workers=None, max_points=None, verbose=True):
    """
    Extend one 2Ddetail-derived gridded variable file from 1939-2023 to 1939-2025.

    Parameters
    ----------
    output_var      : str   Output variable name (e.g. 'T10m', 'SSN', 'LWC_surf')
    var_name        : str   Source variable in per-column 2Ddetail file
    z_begin         : int   Starting layer index (0-based from surface)
    z_end           : int   Ending layer index (exclusive)
    target_depth    : float Target depth in metres for point extraction
    layer_thickness : float Layer thickness in metres (default 0.04)
    operation       : str   'average' or 'sum' for layer range
    workers         : int   Parallel workers
    max_points      : int   Cap point count for testing
    verbose         : bool

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
    timestep = orig_path.stem.split('_')[-1]

    if verbose:
        mode_str = (f"depth={target_depth}m" if target_depth is not None
                    else f"layers {z_begin}:{z_end} ({operation})")
        print(f"\n  Output var: {output_var}  |  source var: {var_name}  "
              f"|  mode: {mode_str}  |  timestep: {timestep}  |  workers: {workers}")

    # Read orig_ntime_output and time values from existing file
    with xr.open_dataset(orig_path, decode_times=False) as ds_orig:
        orig_ntime_output = len(ds_orig.time)
        orig_time_values  = ds_orig.time.values
        orig_attrs        = dict(ds_orig.attrs)

    # Determine total timesteps from a sample per-column 2Ddetail file
    pointlist = load_pointlist(cfg.POINTLIST_FILE)
    total_ind_t = None
    for _, row in pointlist.iterrows():
        fp = cfg.EXT_INPUT_DIR / cfg.EXT_2DDETAIL_FILENAME_PATTERN.format(
            point_id=int(row['point_id']))
        if fp.exists():
            with xr.open_dataset(fp) as ds_s:
                total_ind_t = ds_s.sizes.get('ind_t', None)
            break
    if total_ind_t is None:
        raise RuntimeError("Could not determine total timesteps from merged 2Ddetail files")

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
         var_name, z_begin, z_end, target_depth, layer_thickness, operation,
         orig_ntime_output)
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
        orig_data = ds_orig[output_var].values
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
        f'(extend_variable_2ddetail.py)'
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
        description='Extend 2Ddetail-derived gridded variable files from 1939-2023 to 1939-2025.',
        epilog="""
Examples:
  python extend_variable_2ddetail.py --output-var T10m --var temp --depth 10
  python extend_variable_2ddetail.py --output-var SSN  --var dens --depth-begin 0 --depth-end 0.5
  python extend_variable_2ddetail.py --output-var LWC_surf --var lwc --depth-begin 0 --depth-end 0.2
  python extend_variable_2ddetail.py --output-var rho_2m --var dens --depth-begin 0 --depth-end 2
        """,
    )
    parser.add_argument('--output-var', required=True,
                        help='Output variable name (e.g. T10m, SSN, LWC_surf)')
    parser.add_argument('--var', '-v', required=True,
                        help='Source variable in per-column 2Ddetail file')

    mode = parser.add_argument_group('Processing mode (choose one)')
    mode.add_argument('--depth', '-d', type=float,
                      help='Extract value at this single depth (metres)')
    mode.add_argument('--z-begin', type=int,
                      help='Starting layer index (0-based from surface)')
    mode.add_argument('--z-end', type=int,
                      help='Ending layer index (exclusive)')
    mode.add_argument('--depth-begin', type=float,
                      help='Starting depth in metres (from surface)')
    mode.add_argument('--depth-end', type=float,
                      help='Ending depth in metres (from surface)')
    mode.add_argument('--operation', choices=['average', 'sum'], default='average',
                      help='Operation for layer/depth range (default: average)')

    parser.add_argument('--layer-thickness', type=float, default=0.04,
                        help='Layer thickness in metres (default: 0.04)')
    parser.add_argument('--workers', '-w', type=int, default=None)
    parser.add_argument('--max-points', type=int, default=None)
    parser.add_argument('--quiet', '-q', action='store_true')
    args = parser.parse_args()

    # Resolve processing mode
    layer_thickness = args.layer_thickness
    z_begin = args.z_begin
    z_end = args.z_end
    target_depth = args.depth

    has_layer = z_begin is not None and z_end is not None
    has_depth_range = args.depth_begin is not None and args.depth_end is not None
    has_single = target_depth is not None

    if not (has_layer or has_depth_range or has_single):
        parser.error("Specify one of: --z-begin/--z-end, --depth-begin/--depth-end, or --depth")
    if sum([has_layer, has_depth_range, has_single]) > 1:
        parser.error("Use only one mode: layer indices, depth range, or single depth")

    if has_depth_range:
        z_begin = int(args.depth_begin / layer_thickness)
        z_end   = int(np.ceil(args.depth_end / layer_thickness))
        if not args.quiet:
            print(f"Depth range {args.depth_begin}m–{args.depth_end}m → layers {z_begin}:{z_end}")

    try:
        extend_variable(
            output_var=args.output_var,
            var_name=args.var,
            z_begin=z_begin,
            z_end=z_end,
            target_depth=target_depth,
            layer_thickness=layer_thickness,
            operation=args.operation,
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
