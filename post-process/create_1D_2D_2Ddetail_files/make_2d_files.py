#!/usr/bin/env python3
"""
Create gridded maps from IMAU-FDM 2D profile output files.

This script processes 2D profile files (time × depth) to create gridded output,
such as the depth of critical density levels (z550, z830, z917).

Usage:
    python make_2d_files.py --output-dir /path/to/output/ --var dens --threshold 830 --output-var z830
    python make_2d_files.py --output-dir /path/to/output/ --var temp --depth 10 --output-var T10m

Examples:
    # Find depth where density reaches 830 kg/m³
    python make_2d_files.py -o /scratch/run/output/ -v dens -t 830 --output-var z830

    # Find depth where density reaches 550 kg/m³ (firn-ice transition)
    python make_2d_files.py -o /scratch/run/output/ -v dens -t 550 --output-var z550

    # Get temperature at 10m depth
    python make_2d_files.py -o /scratch/run/output/ -v temp -d 10 --output-var T10m
"""

import argparse
import sys
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
import xarray as xr
from tqdm import tqdm

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent))
import config as post_config
from run_config import RunConfig, load_pointlist, load_mask, validate_time_aggregation
from utils import create_output_dataset


def find_depth_at_threshold(depth, values, threshold, direction='max'):
    """
    Find the depth where a variable crosses a threshold.

    Parameters
    ----------
    depth : np.ndarray
        Depth values (1D or 2D with shape [layer, time])
    values : np.ndarray
        Variable values with shape [layer, time]
    threshold : float
        Threshold value to find
    direction : str
        'max': find maximum depth where value < threshold (for density)
        'min': find minimum depth where value > threshold

    Returns
    -------
    np.ndarray
        Depth at threshold for each timestep
    """
    n_layers, n_times = values.shape

    # Handle depth array
    if depth.ndim == 1:
        depth_2d = np.tile(depth[:, np.newaxis], (1, n_times))
    else:
        depth_2d = depth

    result = np.full(n_times, np.nan)

    for t in range(n_times):
        vals = values[:, t]
        deps = depth_2d[:, t] if depth_2d.ndim == 2 else depth

        # Find valid (non-NaN) values
        valid = ~np.isnan(vals)
        if not np.any(valid):
            continue

        if direction == 'max':
            # Find maximum depth where value < threshold
            mask = (vals < threshold) & valid
            if np.any(mask):
                result[t] = np.max(deps[mask])
        else:
            # Find minimum depth where value > threshold
            mask = (vals > threshold) & valid
            if np.any(mask):
                result[t] = np.min(deps[mask])

    return result


def get_value_at_depth(depth, values, target_depth):
    """
    Get variable value at a specific depth.

    Parameters
    ----------
    depth : np.ndarray
        Depth values (1D or 2D)
    values : np.ndarray
        Variable values with shape [layer, time]
    target_depth : float
        Target depth in meters

    Returns
    -------
    np.ndarray
        Value at target depth for each timestep
    """
    n_layers, n_times = values.shape

    if depth.ndim == 1:
        # Find nearest layer to target depth
        layer_idx = np.argmin(np.abs(depth - target_depth))
        return values[layer_idx, :]
    else:
        # Depth varies with time
        result = np.full(n_times, np.nan)
        for t in range(n_times):
            deps = depth[:, t]
            valid = ~np.isnan(deps)
            if np.any(valid):
                layer_idx = np.argmin(np.abs(deps[valid] - target_depth))
                result[t] = values[valid, t][layer_idx]
        return result


def process_point(args):
    """
    Process a single grid point.

    Parameters
    ----------
    args : tuple
        (point_num, file_path, var_name, threshold, target_depth, secondary_var)

    Returns
    -------
    tuple
        (point_num, rlat_idx, rlon_idx, result_array) or None if failed
    """
    point_num, file_path, rlat_idx, rlon_idx, var_name, threshold, target_depth, secondary_var = args

    if not file_path.exists():
        return None

    try:
        with xr.open_dataset(file_path) as ds:
            # Get the variable
            if var_name not in ds:
                return None

            values = ds[var_name].values  # [layer, time]

            # Get depth
            if secondary_var and secondary_var in ds:
                depth = ds[secondary_var].values
            elif 'depth' in ds:
                depth = ds['depth'].values
            else:
                # Create depth array from layer index
                n_layers = values.shape[0]
                depth = np.arange(n_layers) * 0.15  # Assume 15cm layers

            # Process based on mode
            if threshold is not None:
                result = find_depth_at_threshold(depth, values, threshold)
            elif target_depth is not None:
                result = get_value_at_depth(depth, values, target_depth)
            else:
                # Default: return surface value
                result = values[0, :]

            return (point_num, rlat_idx, rlon_idx, result)

    except Exception as e:
        print(f"Error processing point {point_num}: {e}")
        return None




def main():
    parser = argparse.ArgumentParser(
        description='Create gridded maps from IMAU-FDM 2D profile output.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Required arguments
    parser.add_argument('-o', '--output-dir', required=True,
                        help='Directory containing model output files')

    # Processing mode (one of these required)
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument('-t', '--threshold', type=float,
                            help='Find depth where variable reaches this threshold')
    mode_group.add_argument('-d', '--depth', type=float,
                            help='Extract value at this depth (meters)')

    # Variable specification
    parser.add_argument('-v', '--var', default='dens',
                        help='Variable to process (default: dens)')
    parser.add_argument('--secondary-var', default='depth',
                        help='Secondary variable for depth (default: depth)')
    parser.add_argument('--output-var', required=True,
                        help='Name for output variable (e.g., z830, T10m)')

    # Optional arguments
    parser.add_argument('--reference-dir',
                        help='Directory containing reference files (mask, pointlist)')
    parser.add_argument('--processed-dir',
                        help='Output directory for processed files')
    parser.add_argument('-n', '--num-workers', type=int, default=None,
                        help='Number of parallel workers (default: all CPUs)')
    parser.add_argument('--start-year', type=float,
                        help='Start year for output (default: auto-detect)')
    parser.add_argument('--end-year', type=float,
                        help='End year for output (default: auto-detect)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print configuration without processing')

    args = parser.parse_args()

    # Load configuration
    print("Loading run configuration...")
    config = RunConfig(
        output_dir=args.output_dir,
        reference_dir=args.reference_dir,
        processed_output_dir=args.processed_dir,
    )

    if args.dry_run:
        print(config.summary())
        return

    print(config.summary())

    # Validate that the configured aggregation is not finer than the model output timestep
    validate_time_aggregation(
        post_config.TIME_AGGREGATION_2D,
        config.timestep_2d or 2592000,
        context='2D files',
    )

    # Check for required files
    mask_path = config.get_mask_path()
    pointlist_path = config.get_pointlist_path()

    if not mask_path or not mask_path.exists():
        print(f"ERROR: Mask file not found. Specify --reference-dir")
        sys.exit(1)

    if not pointlist_path or not pointlist_path.exists():
        print(f"ERROR: Pointlist file not found. Specify --reference-dir")
        sys.exit(1)

    # Load mask and pointlist
    print(f"Loading mask from {mask_path}...")
    mask_data = load_mask(mask_path)

    print(f"Loading pointlist from {pointlist_path}...")
    pointlist = load_pointlist(pointlist_path)

    # Extract grid indices from pointlist (last two columns are rlat, rlon)
    rlat_indices = pointlist[:, -2].astype(int)
    rlon_indices = pointlist[:, -1].astype(int)
    n_points = len(pointlist)

    print(f"Processing {n_points} points...")

    # Determine time array
    n_timesteps = config.n_timesteps_2d
    if n_timesteps is None:
        # Read from first file
        sample_file = config.get_2d_file(1)
        if sample_file.exists():
            with xr.open_dataset(sample_file) as ds:
                time_dim = [d for d in ds.dims if 'ind_t' in d or 'time' in d.lower()]
                if time_dim:
                    n_timesteps = ds.dims[time_dim[0]]

    if n_timesteps is None:
        print("ERROR: Could not determine number of timesteps")
        sys.exit(1)

    # Create time array
    time_array = pd.DatetimeIndex([config.get_datetime(t, '2D') for t in range(n_timesteps)])

    # Filter by year range if specified
    if args.start_year:
        time_mask = time_array.year >= args.start_year
        if args.end_year:
            time_mask &= time_array.year <= args.end_year
        time_indices = np.where(time_mask)[0]
    else:
        time_indices = np.arange(n_timesteps)

    time_array = time_array[time_indices]
    n_times_output = len(time_array)

    print(f"Output time range: {time_array[0]} to {time_array[-1]} ({n_times_output} timesteps)")

    # Initialize output grid
    nlat, nlon = config.grid_shape
    grid_data = np.full((n_times_output, nlat, nlon), np.nan, dtype=np.float32)

    # Prepare tasks
    tasks = []
    for i in range(n_points):
        point_num = i + 1  # 1-based indexing
        file_path = config.get_2d_file(point_num)
        rlat_idx = rlat_indices[i]
        rlon_idx = rlon_indices[i]

        tasks.append((
            point_num, file_path, rlat_idx, rlon_idx,
            args.var, args.threshold, args.depth, args.secondary_var
        ))

    # Process in parallel
    n_workers = args.num_workers or None
    completed = 0
    errors = 0

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(process_point, task): task[0] for task in tasks}

        with tqdm(total=n_points, desc='Processing points',miniters=5000) as pbar:
            for future in as_completed(futures):
                result = future.result()
                if result is not None:
                    point_num, rlat_idx, rlon_idx, data = result
                    # Select time indices
                    if len(time_indices) < len(data):
                        data = data[time_indices]
                    grid_data[:, rlat_idx, rlon_idx] = data
                    completed += 1
                else:
                    errors += 1
                pbar.update(1)

    print(f"Completed: {completed}, Errors/Missing: {errors}")

    # Build variable metadata from args
    if args.threshold is not None:
        var_metadata = {
            'long_name': f'Depth where {args.var} = {args.threshold}',
            'units': 'm',
        }
    elif args.depth is not None:
        var_metadata = {
            'long_name': f'{args.var} at {args.depth}m depth',
            'units': '',
        }
    else:
        var_metadata = {'long_name': args.output_var, 'units': ''}

    # Create output dataset
    print("Creating output dataset...")
    ds = create_output_dataset(
        var_name=args.output_var,
        data=grid_data,
        time_values=time_array,
        mask_ds=mask_data,
        var_metadata=var_metadata,
        timestep=post_config.TIME_AGGREGATION_2D,
        domain=config.domain,
    )

    # Add extra variable-specific attributes (threshold/depth values)
    if args.threshold is not None:
        ds[args.output_var].attrs['threshold'] = args.threshold
    elif args.depth is not None:
        ds[args.output_var].attrs['depth'] = args.depth

    # Save output
    config.processed_output_dir.mkdir(parents=True, exist_ok=True)
    output_file = config.processed_output_dir / config.get_output_filename(args.output_var, '2D')

    print(f"Saving to {output_file}...")
    ds.to_netcdf(output_file, encoding={args.output_var: {'zlib': True, 'complevel': 4}})

    print("Done!")
    print(f"Output: {output_file}")


if __name__ == '__main__':
    main()
