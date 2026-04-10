#!/usr/bin/env python3
"""
Create gridded maps from IMAU-FDM 2Ddetail profile output files.

This script processes 2Ddetail files (high-resolution depth profiles, typically
4cm layers in upper 20m) to create gridded output such as surface snow density
(SSN) or near-surface temperature.

Three modes for specifying depth:
    1. Layer indices: --z-begin and --z-end (0-based layer numbers)
    2. Depth range: --depth-begin and --depth-end (meters from surface)
    3. Single depth: -d/--depth (meters from surface)

Examples:
    # Surface snow density using layer indices (layers 0-12 = upper ~52cm)
    python make_2Ddetail_files.py -o /scratch/run/output/ -v dens --z-begin 0 --z-end 13 --output-var SSN

    # Surface snow density using depth range (0 to 0.5m)
    python make_2Ddetail_files.py -o /scratch/run/output/ -v dens --depth-begin 0 --depth-end 0.5 --output-var SSN

    # Average density in upper 2m (using depth range)
    python make_2Ddetail_files.py -o /scratch/run/output/ -v dens --depth-begin 0 --depth-end 2 --output-var rho_2m

    # Temperature at 10m depth (single point)
    python make_2Ddetail_files.py -o /scratch/run/output/ -v temp -d 10 --output-var T10m

    # Surface liquid water content (upper 20cm)
    python make_2Ddetail_files.py -o /scratch/run/output/ -v lwc --depth-begin 0 --depth-end 0.2 --output-var LWC_surf
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


def average_over_layers(values, z_begin, z_end):
    """
    Average variable over a range of layers.

    Parameters
    ----------
    values : np.ndarray
        Variable values with shape [layer, time]
    z_begin : int
        Starting layer index (0-based, from surface)
    z_end : int
        Ending layer index (exclusive)

    Returns
    -------
    np.ndarray
        Averaged values for each timestep
    """
    n_layers, n_times = values.shape

    # Ensure indices are within bounds
    z_begin = max(0, z_begin)
    z_end = min(n_layers, z_end)

    if z_begin >= z_end:
        return np.full(n_times, np.nan)

    # Average over layers (ignoring NaN)
    layer_data = values[z_begin:z_end, :]
    return np.nanmean(layer_data, axis=0)


def get_value_at_depth(depth, values, target_depth, layer_thickness):
    """
    Get variable value at a specific depth.

    Parameters
    ----------
    depth : np.ndarray
        Depth values (1D array of layer centers)
    values : np.ndarray
        Variable values with shape [layer, time]
    target_depth : float
        Target depth in meters
    layer_thickness : float
        Thickness of each layer in meters

    Returns
    -------
    np.ndarray
        Value at target depth for each timestep
    """
    # For 2Ddetail files, depth is typically fixed
    if depth.ndim == 1:
        # Find nearest layer
        layer_idx = np.argmin(np.abs(depth - target_depth))
        return values[layer_idx, :]
    else:
        # Depth varies (shouldn't happen for 2Ddetail but handle anyway)
        n_layers, n_times = values.shape
        result = np.full(n_times, np.nan)
        for t in range(n_times):
            deps = depth[:, t] if depth.ndim == 2 else depth
            layer_idx = np.argmin(np.abs(deps - target_depth))
            result[t] = values[layer_idx, t]
        return result


def sum_over_layers(values, z_begin, z_end):
    """
    Sum variable over a range of layers.

    Parameters
    ----------
    values : np.ndarray
        Variable values with shape [layer, time]
    z_begin : int
        Starting layer index (0-based, from surface)
    z_end : int
        Ending layer index (exclusive)

    Returns
    -------
    np.ndarray
        Summed values for each timestep
    """
    n_layers, n_times = values.shape

    z_begin = max(0, z_begin)
    z_end = min(n_layers, z_end)

    if z_begin >= z_end:
        return np.full(n_times, np.nan)

    layer_data = values[z_begin:z_end, :]
    return np.nansum(layer_data, axis=0)


def process_point(args):
    """
    Process a single grid point.

    Parameters
    ----------
    args : tuple
        (point_num, file_path, rlat_idx, rlon_idx, var_name, z_begin, z_end,
         target_depth, layer_thickness, operation)

    Returns
    -------
    tuple
        (point_num, rlat_idx, rlon_idx, result_array) or None if failed
    """
    (point_num, file_path, rlat_idx, rlon_idx, var_name,
     z_begin, z_end, target_depth, layer_thickness, operation) = args

    if not file_path.exists():
        return None

    try:
        with xr.open_dataset(file_path) as ds:
            if var_name not in ds:
                return None

            values = ds[var_name].values  # [layer, time]

            # Get depth array if available
            if 'depth' in ds:
                depth = ds['depth'].values
            else:
                # Create depth array from layer thickness
                n_layers = values.shape[0]
                depth = np.arange(n_layers) * layer_thickness + layer_thickness / 2

            # Process based on mode
            if target_depth is not None:
                result = get_value_at_depth(depth, values, target_depth, layer_thickness)
            elif z_begin is not None and z_end is not None:
                if operation == 'sum':
                    result = sum_over_layers(values, z_begin, z_end)
                else:  # average
                    result = average_over_layers(values, z_begin, z_end)
            else:
                # Default: surface layer
                result = values[0, :]

            return (point_num, rlat_idx, rlon_idx, result)

    except Exception as e:
        print(f"Error processing point {point_num}: {e}")
        return None


# Variable metadata for common outputs
VARIABLE_METADATA = {
    'SSN': {
        'long_name': 'Surface snow density',
        'units': 'kg m-3',
        'description': 'Average density in upper firn layers',
    },
    'T10m': {
        'long_name': '10-meter firn temperature',
        'units': 'K',
        'description': 'Temperature at 10m depth',
    },
    'LWC_surf': {
        'long_name': 'Surface liquid water content',
        'units': '',
        'description': 'Liquid water content in upper layers',
    },
    'rho_2m': {
        'long_name': 'Average density in upper 2m',
        'units': 'kg m-3',
        'description': 'Mean density in upper 2 meters of firn',
    },
}


def main():
    parser = argparse.ArgumentParser(
        description='Create gridded maps from IMAU-FDM 2Ddetail output.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Required arguments
    parser.add_argument('-o', '--output-dir', required=True,
                        help='Directory containing model output files')
    parser.add_argument('-v', '--var', required=True,
                        help='Variable to process (dens, temp, lwc, refreeze)')
    parser.add_argument('--output-var', required=True,
                        help='Name for output variable (e.g., SSN, T10m)')

    # Processing mode - three options:
    # 1. Layer indices (--z-begin, --z-end)
    # 2. Depth range in meters (--depth-begin, --depth-end)
    # 3. Single depth point (-d/--depth)
    mode_group = parser.add_argument_group('Processing mode (choose one set)')
    mode_group.add_argument('--z-begin', type=int,
                            help='Starting layer index (0-based from surface)')
    mode_group.add_argument('--z-end', type=int,
                            help='Ending layer index (exclusive)')
    mode_group.add_argument('--depth-begin', type=float,
                            help='Starting depth in meters (from surface)')
    mode_group.add_argument('--depth-end', type=float,
                            help='Ending depth in meters (from surface)')
    mode_group.add_argument('-d', '--depth', type=float,
                            help='Extract value at this single depth (meters)')
    mode_group.add_argument('--operation', choices=['average', 'sum'], default='average',
                            help='Operation for layer/depth range (default: average)')

    # Optional arguments
    parser.add_argument('--reference-dir',
                        help='Directory containing reference files')
    parser.add_argument('--processed-dir',
                        help='Output directory for processed files')
    parser.add_argument('--layer-thickness', type=float,
                        help='Layer thickness in meters (default: auto-detect)')
    parser.add_argument('-n', '--num-workers', type=int, default=post_config.NUM_WORKERS,
                        help='Number of parallel workers (default: from config/SLURM_CPUS_PER_TASK)')
    parser.add_argument('--start-year', type=float,
                        help='Start year for output')
    parser.add_argument('--end-year', type=float,
                        help='End year for output')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print configuration without processing')

    args = parser.parse_args()

    # Validate arguments - need one of:
    # 1. Both --z-begin and --z-end (layer indices)
    # 2. Both --depth-begin and --depth-end (depth range in meters)
    # 3. --depth (single depth point)
    has_layer_range = args.z_begin is not None and args.z_end is not None
    has_depth_range = args.depth_begin is not None and args.depth_end is not None
    has_single_depth = args.depth is not None

    if not (has_layer_range or has_depth_range or has_single_depth):
        parser.error("Specify one of: --z-begin/--z-end, --depth-begin/--depth-end, or --depth")

    if sum([has_layer_range, has_depth_range, has_single_depth]) > 1:
        parser.error("Use only one mode: layer indices, depth range, or single depth")

    # Load configuration
    print("Loading run configuration...")
    config = RunConfig(
        output_dir=args.output_dir,
        reference_dir=args.reference_dir,
        processed_output_dir=args.processed_dir,
    )

    print(config.summary())

    # Validate that the configured aggregation is not finer than the model output timestep
    validate_time_aggregation(
        post_config.TIME_AGGREGATION_2Ddetail,
        config.timestep_2ddetail or 864000,
        context='2Ddetail files',
    )

    # Get layer thickness
    layer_thickness = args.layer_thickness or config.detail_thickness or 0.04
    print(f"Layer thickness: {layer_thickness}m")

    # Convert depth range to layer indices if specified
    z_begin = args.z_begin
    z_end = args.z_end

    if args.depth_begin is not None and args.depth_end is not None:
        # Convert depth (m) to layer index
        # Layer centers are at: layer_thickness/2, 3*layer_thickness/2, ...
        # Layer i covers depth from i*layer_thickness to (i+1)*layer_thickness
        z_begin = int(args.depth_begin / layer_thickness)
        z_end = int(np.ceil(args.depth_end / layer_thickness))
        print(f"Depth range {args.depth_begin}m - {args.depth_end}m -> layers {z_begin} - {z_end}")

    if args.dry_run:
        return

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

    rlat_indices = pointlist[:, -2].astype(int)
    rlon_indices = pointlist[:, -1].astype(int)
    n_points = len(pointlist)

    print(f"Processing {n_points} points...")

    # Determine time array
    n_timesteps = config.n_timesteps_2ddetail
    if n_timesteps is None:
        sample_file = config.get_2ddetail_file(1)
        if sample_file.exists():
            with xr.open_dataset(sample_file) as ds:
                time_dim = [d for d in ds.dims if 'ind_t' in d or 'time' in d.lower()]
                if time_dim:
                    n_timesteps = ds.dims[time_dim[0]]

    if n_timesteps is None:
        print("ERROR: Could not determine number of timesteps")
        sys.exit(1)

    # Create time array
    time_array = pd.DatetimeIndex([config.get_datetime(t, '2Ddetail') for t in range(n_timesteps)])

    # Filter by year range
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
        point_num = i + 1
        file_path = config.get_2ddetail_file(point_num)
        rlat_idx = rlat_indices[i]
        rlon_idx = rlon_indices[i]

        tasks.append((
            point_num, file_path, rlat_idx, rlon_idx,
            args.var, z_begin, z_end, args.depth,
            layer_thickness, args.operation
        ))

    # Process in parallel
    n_workers = args.num_workers or None
    completed = 0
    errors = 0

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(process_point, task): task[0] for task in tasks}

        with tqdm(total=n_points, desc='Processing points', miniters=5000) as pbar:
            for future in as_completed(futures):
                result = future.result()
                if result is not None:
                    point_num, rlat_idx, rlon_idx, data = result
                    if len(time_indices) < len(data):
                        data = data[time_indices]
                    grid_data[:, rlat_idx, rlon_idx] = data
                    completed += 1
                else:
                    errors += 1
                pbar.update(1)

    print(f"Completed: {completed}, Errors/Missing: {errors}")

    # Build variable metadata
    if args.output_var in VARIABLE_METADATA:
        var_metadata = VARIABLE_METADATA[args.output_var]
        extra_attrs = {}
    elif args.depth is not None:
        var_metadata = {'long_name': f'{args.var} at {args.depth}m depth', 'units': ''}
        extra_attrs = {'depth': args.depth}
    elif z_begin is not None and z_end is not None:
        depth_begin_m = z_begin * layer_thickness
        depth_end_m = z_end * layer_thickness
        var_metadata = {
            'long_name': f'{args.operation} of {args.var} from {depth_begin_m:.2f}m to {depth_end_m:.2f}m',
            'units': '',
        }
        extra_attrs = {
            'depth_begin': depth_begin_m,
            'depth_end': depth_end_m,
            'z_begin': z_begin,
            'z_end': z_end,
        }
    else:
        var_metadata = {'long_name': args.output_var, 'units': ''}
        extra_attrs = {}

    # Create output dataset
    print("Creating output dataset...")
    ds = create_output_dataset(
        var_name=args.output_var,
        data=grid_data,
        time_values=time_array,
        mask_ds=mask_data,
        var_metadata=var_metadata,
        timestep=post_config.TIME_AGGREGATION_2Ddetail,
        domain=config.domain,
    )

    # Add extra variable-specific attributes
    ds[args.output_var].attrs.update(extra_attrs)

    # Slice to output period if configured
    if post_config.OUTPUT_START is not None or post_config.OUTPUT_END is not None:
        t_start = str(post_config.OUTPUT_START.date()) if post_config.OUTPUT_START else None
        t_end   = str(post_config.OUTPUT_END.date())   if post_config.OUTPUT_END   else None
        ds = ds.sel(time=slice(t_start, t_end))
        print(f"Output period: {ds.time.values[0]} to {ds.time.values[-1]}")

    # Save output
    config.processed_output_dir.mkdir(parents=True, exist_ok=True)
    output_file = config.processed_output_dir / config.get_output_filename(args.output_var, '2Ddetail')

    print(f"Saving to {output_file}...")
    ds.to_netcdf(output_file, encoding={args.output_var: {'zlib': True, 'complevel': 4}})

    print("Done!")
    print(f"Output: {output_file}")


if __name__ == '__main__':
    main()
