#!/usr/bin/env python3
"""
Slice decadal RACMO files into yearly files for a single variable.
Replaces makeFDMyears_FGRN055.sc

Run from within the domain subdirectory (which contains config.py):

    cd FGRN055_era055
    python3 ../make_fdm_years.py <varname>
"""

import argparse
import subprocess
import sys
import numpy as np
import xarray as xr
from pathlib import Path

import config


def get_year_indices(infile):
    """Return {year: (start_idx, end_idx)} by reading the time variable."""
    ds = xr.open_dataset(infile, decode_times=True)
    years = ds["time"].dt.year.values
    ds.close()
    result = {}
    for year in np.unique(years):
        idx = np.where(years == year)[0]
        result[int(year)] = (int(idx[0]), int(idx[-1]))
    return result


def make_years(varname):
    years_dir = config.YEARS_DIR
    years_dir.mkdir(parents=True, exist_ok=True)

    decade_files = sorted(config.RAW_DIR.glob(f"{varname}.KNMI-*{config.FNAME_SUFFIX}"))
    if not decade_files:
        print(f"ERROR: no raw files found for '{varname}' in {config.RAW_DIR}")
        sys.exit(1)

    for infile in decade_files:
        print(f"\nProcessing {infile.name} ...")
        year_indices = get_year_indices(infile)

        for year, (start_t, end_t) in sorted(year_indices.items()):
            if year < config.TS_START_YEAR or year > config.TS_END_YEAR:
                continue

            outfile = years_dir / f"{varname}_{config.PROJECT_NAME}_forFDM_Year{year}.nc"

            if outfile.exists():
                print(f"  {year}: already present, skipping")
                continue

            print(f"  {varname} {year}: t={start_t}..{end_t}")
            subprocess.run(
                ["ncks", "-d", f"time,{start_t},{end_t}", str(infile), str(outfile)],
                check=True,
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Slice RACMO decade files into yearly files."
    )
    parser.add_argument("varname", help="Variable name (e.g. precip)")
    args = parser.parse_args()
    make_years(args.varname)
