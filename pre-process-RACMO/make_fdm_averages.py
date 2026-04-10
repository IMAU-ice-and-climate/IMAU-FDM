#!/usr/bin/env python3
"""
Compute spinup-period averages for a single variable.
Replaces makeFDMaverages_FGRN055.sc

Requires that make_fdm_years.py has been run first so that yearly files exist.

Run from within the domain subdirectory (which contains config.py):

    cd FGRN055_era055
    python3 ../make_fdm_averages.py <varname>
"""

import argparse
import subprocess
from pathlib import Path

import config


def make_averages(varname):
    years_dir = config.YEARS_DIR
    ave_dir   = config.AVE_DIR
    ave_dir.mkdir(parents=True, exist_ok=True)

    fname_out = (
        ave_dir
        / f"{varname}_{config.PROJECT_NAME}-{config.TS_START_YEAR}"
          f"_{config.AVE_START_YEAR}-{config.AVE_END_YEAR}_ave.nc"
    )
    if fname_out.exists():
        print(f"Average file already present: {fname_out.name}, skipping")
        return

    # Per-year temporal averages
    yave_files = []
    for year in range(config.AVE_START_YEAR, config.AVE_END_YEAR + 1):
        fname = years_dir / f"{varname}_{config.PROJECT_NAME}_forFDM_Year{year}.nc"
        if not fname.exists():
            print(f"  Warning: {fname.name} not found, skipping year {year}")
            continue
        fname_yave = years_dir / f"{varname}_Yave_{year}.nc"
        print(f"  Averaging year {year}")
        subprocess.run(["ncra", str(fname), str(fname_yave)], check=True)
        yave_files.append(fname_yave)

    if not yave_files:
        print("ERROR: no yearly average files created; aborting")
        return

    # Concatenate yearly averages, then average across years
    temp1 = years_dir / f"{varname}_ave_temp1.nc"
    temp2 = years_dir / f"{varname}_ave_temp2.nc"

    print(f"  Concatenating {len(yave_files)} yearly averages")
    try:
        subprocess.run(
            ["ncrcat"] + [str(f) for f in yave_files] + [str(temp1)],
            check=True,
        )
        subprocess.run(["ncra", str(temp1), str(temp2)], check=True)
        subprocess.run(["nccopy", "-k", "classic", str(temp2), str(fname_out)], check=True)
    finally:
        for p in [temp1, temp2]:
            if p.exists():
                p.unlink()
        for f in yave_files:
            if f.exists():
                f.unlink()

    print(f"  Done: {fname_out.name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute spinup-period averages.")
    parser.add_argument("varname", help="Variable name (e.g. precip)")
    args = parser.parse_args()
    make_averages(args.varname)
