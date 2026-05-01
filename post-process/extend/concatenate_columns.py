"""
Merge per-column FDM output files from an original run and an extended run.

For each point the extended file covers the full model period (e.g. 1939-2025)
but has NaN for all timesteps already present in the original run (1939-2023).
The merge simply uses the extended file as-is and overwrites its NaN region with
the original run data, producing a combined file with no NaN gaps.

Usage
-----
    python3 concatenate_columns.py [--start N] [--end N] [--workers W]

    --start   First point ID to process (default: 1)
    --end     Last point ID to process inclusive (default: N_POINTS from config)
    --workers Number of parallel workers (default: from config)

If a point has no extended file it is skipped and the original is not copied.
"""

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import netCDF4 as nc
import numpy as np

import config


def find_ext_file(point_id: int, ftype: str) -> Path | None:
    """Return path to first extended run file found for this point, or None."""
    for d in config.EXT_OUTPUT_DIRS:
        p = d / f'{config.DOMAIN}_{config.FORCING}_{ftype}_{point_id}.nc'
        if p.exists():
            return p
    return None


def _merge_file(orig_path: Path, ext_path: Path, out_path: Path) -> None:
    """Write a merged file: extended file with its NaN region replaced by orig data."""
    with nc.Dataset(orig_path, 'r') as orig, nc.Dataset(ext_path, 'r') as ext:
        orig_nt = len(orig.dimensions['ind_t'])

        with nc.Dataset(out_path, 'w') as out:
            # Global attributes — keep orig metadata but update model_end_datetime
            for attr in orig.ncattrs():
                val = ext.getncattr(attr) if attr == 'model_end_datetime' else orig.getncattr(attr)
                out.setncattr(attr, val)

            # Dimensions from extended file (has the full combined length)
            for dname, dim in ext.dimensions.items():
                out.createDimension(dname, None if dim.isunlimited() else len(dim))

            # Variables
            for vname, var in ext.variables.items():
                out_var = out.createVariable(vname, var.dtype, var.dimensions)
                out_var.setncatts({a: var.getncattr(a) for a in var.ncattrs()})

                if 'ind_t' not in var.dimensions:
                    # Time-invariant (e.g. depth, dz in 2Ddetail)
                    out_var[:] = var[:]
                    continue

                # Copy extended data, then overwrite NaN region with original
                data = var[:]
                t_axis = var.dimensions.index('ind_t')
                if t_axis == 0:
                    data[:orig_nt] = orig.variables[vname][:]
                else:
                    data[:, :orig_nt] = orig.variables[vname][:]
                out_var[:] = data


def process_point(point_id: int) -> str:
    """Concatenate all file types for one point. Returns a status string."""
    # Check all file types are available before writing anything
    orig_files = {}
    ext_files = {}
    for ftype in config.FILE_TYPES:
        orig = config.ORIG_OUTPUT_DIR / f'{config.DOMAIN}_{config.FORCING}_{ftype}_{point_id}.nc'
        if not orig.exists():
            return f'{point_id}: SKIP (missing original {ftype})'
        ext = find_ext_file(point_id, ftype)
        if ext is None:
            return f'{point_id}: SKIP (missing extended {ftype})'
        orig_files[ftype] = orig
        ext_files[ftype] = ext

    config.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for ftype in config.FILE_TYPES:
        out_path = config.OUTPUT_DIR / f'{config.DOMAIN}_{config.FORCING}_{ftype}_{point_id}.nc'
        _merge_file(orig_files[ftype], ext_files[ftype], out_path)

    return f'{point_id}: OK'


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--start', type=int, default=1,
                        help='First point ID (default: 1)')
    parser.add_argument('--end', type=int, default=config.N_POINTS,
                        help=f'Last point ID inclusive (default: {config.N_POINTS})')
    parser.add_argument('--workers', type=int, default=config.NUM_WORKERS,
                        help='Parallel workers (default: from config)')
    args = parser.parse_args()

    point_ids = range(args.start, args.end + 1)
    n_total = len(point_ids)
    n_ok = n_skip = 0

    print(f'Processing points {args.start}–{args.end} with {args.workers or "all"} workers')
    print(f'Original : {config.ORIG_OUTPUT_DIR}')
    print(f'Extended : {[str(d) for d in config.EXT_OUTPUT_DIRS]}')
    print(f'Output   : {config.OUTPUT_DIR}')

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(process_point, pid): pid for pid in point_ids}
        for i, future in enumerate(as_completed(futures), 1):
            result = future.result()
            if result.endswith('OK'):
                n_ok += 1
            else:
                n_skip += 1
                print(result)
            if i % 1000 == 0 or i == n_total:
                print(f'  {i}/{n_total} done — {n_ok} OK, {n_skip} skipped')

    print(f'\nDone: {n_ok} combined, {n_skip} skipped')


if __name__ == '__main__':
    main()
