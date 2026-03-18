#!/usr/bin/env python3
"""
IMAU-FDM 1D to Gridded Maps Post-Processing

This script converts individual 1D netCDF files (one per grid point) into
gridded maps (3D: time, rlat, rlon) for visualization and analysis.

Features:
- Configurable time aggregation (daily, 10-day, monthly)
- Spinup-aware detrending for h_surf and FirnAir
- Parallel processing for efficient handling of ~58,000 files
- Compatible with interactive use or batch job submission

Usage:
    # Process single variable
    python make_1d_maps.py --var h_surf

    # Process all variables with parallel processing
    python make_1d_maps.py --var all --workers 16

    # Customize options
    python make_1d_maps.py --var h_surf --timestep monthly --spinup-end 1975

    # List available variables
    python make_1d_maps.py --list-vars

Author: Generated for IMAU-FDM post-processing
"""

import argparse
import os
import sys
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

import numpy as np
import xarray as xr
from tqdm import tqdm

# Import local modules
import config
from run_config import validate_time_aggregation
from utils import (
    load_pointlist,
    load_mask,
    get_grid_dimensions,
    create_time_array,
    detrend_spinup_aware,
    find_spinup_end_index,
    resample_timeseries,
    get_output_time_axis,
    create_output_dataset,
    read_1d_file,
)


def process_single_point(args):
    """
    Process a single point (for parallel execution).

    Parameters
    ----------
    args : tuple
        (point_id, rlat_idx, rlon_idx, var_name, input_dir, filename_pattern,
         ntime_daily, spinup_end_idx, needs_detrend, timestep, model_start, aggregation)

    Returns
    -------
    tuple
        (rlat_idx, rlon_idx, resampled_data) or (rlat_idx, rlon_idx, None) if failed
    """
    (point_id, rlat_idx, rlon_idx, var_name, input_dir, filename_pattern,
     ntime_daily, spinup_end_idx, needs_detrend, timestep, model_start, aggregation) = args

    # Construct filepath
    filename = filename_pattern.format(point_id=point_id)
    filepath = input_dir / filename

    # Read data
    data = read_1d_file(filepath, var_name)
    if data is None:
        return (rlat_idx, rlon_idx, None)

    # Trim to expected length if needed
    if len(data) > ntime_daily:
        data = data[:ntime_daily]
    elif len(data) < ntime_daily:
        # Pad with NaN if shorter (shouldn't happen normally)
        padded = np.full(ntime_daily, np.nan)
        padded[:len(data)] = data
        data = padded

    # Apply detrending if needed
    if needs_detrend and spinup_end_idx > 0:
        data = detrend_spinup_aware(data, np.arange(len(data)), spinup_end_idx)

    # Resample to target timestep
    if timestep != 'daily':
        time_indices = np.arange(len(data))
        resampled_data, _ = resample_timeseries(
            data, time_indices, method=timestep, start_date=model_start,
            aggregation=aggregation
        )
    else:
        resampled_data = data

    return (rlat_idx, rlon_idx, resampled_data)


def process_variable(var_name, timestep=None, spinup_start=None, spinup_end=None,
                     workers=None, output_dir=None, max_points=None, verbose=True):
    """
    Process a single variable from all 1D files to gridded output.

    Parameters
    ----------
    var_name : str
        Variable name to process
    timestep : str, optional
        Time aggregation ('daily', '10day', 'monthly'). Default: from config
    spinup_start : datetime, optional
        Start of spinup period. Default: from config
    spinup_end : datetime, optional
        End of spinup period. Default: from config
    workers : int, optional
        Number of parallel workers. Default: from config or all CPUs
    output_dir : Path, optional
        Output directory. Default: from config
    max_points : int, optional
        Maximum number of points to process (for testing). Default: all points
    verbose : bool
        Print progress information

    Returns
    -------
    Path
        Path to the output file
    """
    # Use defaults from config if not specified
    timestep = timestep or config.TIME_AGGREGATION_1D
    spinup_start = spinup_start or config.SPINUP_START
    spinup_end = spinup_end or config.SPINUP_END
    output_dir = Path(output_dir) if output_dir else config.OUTPUT_DIR
    workers = workers or config.NUM_WORKERS or mp.cpu_count()

    # Validate that the requested aggregation is not finer than the input timestep
    validate_time_aggregation(timestep, config.INPUT_TIMESTEP_SECONDS, context='1D files')

    # Check if variable exists
    if var_name not in config.VARIABLES:
        raise ValueError(f"Unknown variable: {var_name}. Use --list-vars to see available variables.")

    var_metadata = config.VARIABLES[var_name]
    needs_detrend = var_name in config.DETREND_VARIABLES

    if verbose:
        print(f"\nProcessing variable: {var_name}")
        print(f"  Time aggregation: {timestep}")
        print(f"  Detrending: {'Yes (spinup-aware)' if needs_detrend else 'No'}")
        print(f"  Workers: {workers}")

    # Load pointlist and mask
    if verbose:
        print("  Loading pointlist and mask...")
    pointlist = load_pointlist(config.POINTLIST_FILE)
    mask_ds = load_mask(config.MASK_FILE)

    # Calculate time dimensions
    ntime_daily = len(create_time_array(config.MODEL_START, config.MODEL_END))
    time_values = get_output_time_axis(
        config.MODEL_START, config.MODEL_END, method=timestep
    )
    ntime_output = len(time_values)

    if verbose:
        print(f"  Input timesteps (daily): {ntime_daily}")
        print(f"  Output timesteps ({timestep}): {ntime_output}")

    # Calculate spinup end index
    daily_times = create_time_array(config.MODEL_START, config.MODEL_END)
    spinup_end_idx = find_spinup_end_index(daily_times, spinup_end) if needs_detrend else 0

    if verbose and needs_detrend:
        print(f"  Spinup period: {spinup_start.date()} to {spinup_end.date()}")
        print(f"  Spinup end index: {spinup_end_idx}")

    # Get grid dimensions from mask file
    n_rlat, n_rlon = get_grid_dimensions(config.MASK_FILE)
    if verbose:
        print(f"  Grid dimensions: {n_rlat} x {n_rlon}")

    # Initialize output array
    output_data = np.full((ntime_output, n_rlat, n_rlon), np.nan, dtype=np.float32)

    # Limit points if max_points is specified (for testing)
    if max_points is not None:
        pointlist = pointlist.head(max_points)

    # Get aggregation method for this variable
    aggregation = config.get_aggregation_method(var_name)
    if verbose:
        print(f"  Aggregation method: {aggregation}")

    # Prepare arguments for parallel processing
    # Note: explicitly convert to int to avoid pandas float conversion issue
    args_list = [
        (
            int(row['point_id']),
            int(row['rlat_idx']),
            int(row['rlon_idx']),
            var_name,
            config.INPUT_DIR,
            config.INPUT_FILENAME_PATTERN,
            ntime_daily,
            spinup_end_idx,
            needs_detrend,
            timestep,
            config.MODEL_START,
            aggregation,
        )
        for _, row in pointlist.iterrows()
    ]

    # Process files in parallel
    if verbose:
        print(f"  Processing {len(args_list)} points...")

    completed = 0
    failed = 0

    with ProcessPoolExecutor(max_workers=workers) as executor:
        # Submit all tasks
        futures = {executor.submit(process_single_point, args): args for args in args_list}

        # Collect results with periodic progress updates
        total_points = len(args_list)
        processed = 0
        for future in as_completed(futures):
            rlat_idx, rlon_idx, data = future.result()
            if data is not None:
                output_data[:, rlat_idx, rlon_idx] = data
                completed += 1
            else:
                failed += 1
            processed += 1
            if verbose and processed % 5000 == 0:
                print(f"  {var_name}: {processed}/{total_points} points ({100*processed/total_points:.1f}%)")

    if verbose:
        print(f"  Completed: {completed}, Failed/Missing: {failed}")

    # Create output dataset
    ds = create_output_dataset(
        var_name=var_name,
        data=output_data,
        time_values=time_values,
        mask_ds=mask_ds,
        var_metadata=var_metadata,
        detrended=needs_detrend,
        timestep=timestep,
    )

    # Create output directory if needed
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate output filename
    output_filename = config.get_output_filename(var_name, timestep, needs_detrend)
    output_path = output_dir / output_filename

    # Save to NetCDF
    if verbose:
        print(f"  Saving to: {output_path}")

    ds.to_netcdf(output_path, format='NETCDF4', engine='netcdf4')

    if verbose:
        print(f"  Done!")

    return output_path


def process_all_variables(variables=None, timestep=None, spinup_start=None,
                          spinup_end=None, workers=None, output_dir=None,
                          max_points=None, verbose=True):
    """
    Process multiple variables.

    Parameters
    ----------
    variables : list, optional
        List of variable names. Default: all variables
    timestep : str, optional
        Time aggregation
    spinup_start : datetime, optional
        Start of spinup period
    spinup_end : datetime, optional
        End of spinup period
    workers : int, optional
        Number of parallel workers
    output_dir : Path, optional
        Output directory
    max_points : int, optional
        Maximum number of points to process (for testing)
    verbose : bool
        Print progress information

    Returns
    -------
    list
        List of output file paths
    """
    if variables is None:
        variables = config.get_variable_names()

    output_paths = []
    for i, var_name in enumerate(variables, 1):
        if verbose:
            print(f"\n{'='*60}")
            print(f"Variable {i}/{len(variables)}: {var_name}")
            print('='*60)

        try:
            path = process_variable(
                var_name=var_name,
                timestep=timestep,
                spinup_start=spinup_start,
                spinup_end=spinup_end,
                workers=workers,
                output_dir=output_dir,
                max_points=max_points,
                verbose=verbose,
            )
            output_paths.append(path)
        except Exception as e:
            print(f"  ERROR processing {var_name}: {e}")
            continue

    return output_paths


def list_variables():
    """Print available variables and their metadata."""
    print("\nAvailable variables:")
    print("-" * 60)
    for var_name, meta in config.VARIABLES.items():
        detrend_note = " [DETREND]" if meta.get('needs_detrend') else ""
        print(f"  {var_name:12s} - {meta['long_name']} ({meta['units']}){detrend_note}")
    print("-" * 60)
    print("\nVariables marked [DETREND] will have spinup-aware detrending applied.")


def main():
    """Main entry point for command-line usage."""
    parser = argparse.ArgumentParser(
        description='Convert IMAU-FDM 1D output files to gridded maps.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python make_1d_maps.py --var h_surf
  python make_1d_maps.py --var all --workers 16
  python make_1d_maps.py --var h_surf FirnAir Runoff --timestep monthly
  python make_1d_maps.py --list-vars
        """
    )

    parser.add_argument(
        '--var', '-v',
        nargs='+',
        help='Variable(s) to process. Use "all" for all variables.'
    )
    parser.add_argument(
        '--timestep', '-t',
        choices=['daily', '10day', 'monthly'],
        default=None,
        help=f'Time aggregation for 1D resampling (default: {config.TIME_AGGREGATION_1D})'
    )
    parser.add_argument(
        '--spinup-start',
        type=int,
        default=None,
        help=f'Spinup start year (default: {config.SPINUP_START.year})'
    )
    parser.add_argument(
        '--spinup-end',
        type=int,
        default=None,
        help=f'Spinup end year (default: {config.SPINUP_END.year})'
    )
    parser.add_argument(
        '--workers', '-w',
        type=int,
        default=None,
        help='Number of parallel workers (default: all CPUs)'
    )
    parser.add_argument(
        '--output-dir', '-o',
        type=str,
        default=None,
        help=f'Output directory (default: {config.OUTPUT_DIR})'
    )
    parser.add_argument(
        '--max-points',
        type=int,
        default=None,
        help='Maximum number of points to process (for testing)'
    )
    parser.add_argument(
        '--list-vars',
        action='store_true',
        help='List available variables and exit'
    )
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Suppress progress output'
    )

    args = parser.parse_args()

    # Handle --list-vars
    if args.list_vars:
        list_variables()
        return 0

    # Require --var unless listing variables
    if not args.var:
        parser.print_help()
        print("\nError: --var is required. Use --list-vars to see available variables.")
        return 1

    # Parse spinup dates
    spinup_start = datetime(args.spinup_start, 1, 1) if args.spinup_start else None
    spinup_end = datetime(args.spinup_end, 1, 1) if args.spinup_end else None

    # Parse output directory
    output_dir = Path(args.output_dir) if args.output_dir else None

    # Process variables
    verbose = not args.quiet

    if 'all' in args.var:
        variables = None  # All variables
    else:
        variables = args.var

    if variables is None:
        # Process all variables
        output_paths = process_all_variables(
            variables=None,
            timestep=args.timestep,
            spinup_start=spinup_start,
            spinup_end=spinup_end,
            workers=args.workers,
            output_dir=output_dir,
            max_points=args.max_points,
            verbose=verbose,
        )
    elif len(variables) == 1:
        # Process single variable
        output_path = process_variable(
            var_name=variables[0],
            timestep=args.timestep,
            spinup_start=spinup_start,
            spinup_end=spinup_end,
            workers=args.workers,
            output_dir=output_dir,
            max_points=args.max_points,
            verbose=verbose,
        )
        output_paths = [output_path]
    else:
        # Process multiple specific variables
        output_paths = process_all_variables(
            variables=variables,
            timestep=args.timestep,
            spinup_start=spinup_start,
            spinup_end=spinup_end,
            workers=args.workers,
            output_dir=output_dir,
            max_points=args.max_points,
            verbose=verbose,
        )

    if verbose:
        print(f"\n{'='*60}")
        print("Processing complete!")
        print(f"Output files: {len(output_paths)}")
        for path in output_paths:
            print(f"  {path}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
