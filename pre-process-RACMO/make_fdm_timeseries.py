#!/usr/bin/env python3
"""
Create longitude-band timeseries files for a single variable.
Replaces makeFDMtimeseries_FGRN055.sc

Steps:
  1. (optional) Slice each yearly file into longitude-band parts.
  2. Concatenate all parts across years into a single timeseries per band.

Run from within the domain subdirectory (which contains config.py):

    cd FGRN055_era055
    python3 ../make_fdm_timeseries.py <varname> [--no-parts]

Extending an existing timeseries with new RACMO years:
    python3 ../make_fdm_timeseries.py <varname> --extend-start 2024 --extend-end 2025 --extend-mode parts-only
    python3 ../make_fdm_timeseries.py <varname> --extend-start 2024 --extend-end 2025 --extend-mode append

--extend-mode parts-only
    Slice the new yearly files into longitude-band parts and save them in a
    separate parts_extend/ subfolder. Existing timeseries files are untouched.
    Use this when you want to inspect or archive the new parts before committing.

--extend-mode append
    Same as parts-only, but also appends the new parts to the existing
    timeseries files. The existing timeseries file is replaced with a new file
    that spans the original years plus the new years. The old file is renamed
    to <name>.bak before being overwritten.
"""

import argparse
import subprocess
from pathlib import Path

import config


def _slice_years_into_parts(varname, year_range, parts_dir):
    """Slice yearly files for years in year_range into parts_dir."""
    years_dir = config.YEARS_DIR
    parts_dir.mkdir(parents=True, exist_ok=True)

    for year in year_range:
        fname_in = years_dir / f"{varname}_{config.PROJECT_NAME}_forFDM_Year{year}.nc"
        if not fname_in.exists():
            print(f"  Warning: {fname_in.name} not found, skipping year {year}")
            continue

        start_t = 0
        for part in range(1, config.NUM_LONG_BANDS + 1):
            end_t = start_t + config.CELL_WIDTH - 1
            fname_out = parts_dir / f"{varname}_{year}_part{part}.nc"

            if fname_out.exists():
                print(f"  {year} part {part}: already present, skipping")
            else:
                print(f"  Slicing {year} part {part} (rlon {start_t}..{end_t})")
                subprocess.run(
                    ["ncks", "-d", f"rlon,{start_t},{end_t}", str(fname_in), str(fname_out)],
                    check=True,
                )
            start_t = end_t + 1


def make_timeseries(varname, do_parts=True):
    """Build full timeseries files from scratch (normal mode)."""
    years_dir = config.YEARS_DIR
    parts_dir = years_dir / "parts"
    files_dir = config.FILES_DIR
    files_dir.mkdir(parents=True, exist_ok=True)

    if do_parts:
        _slice_years_into_parts(
            varname,
            range(config.TS_START_YEAR, config.TS_END_YEAR + 1),
            parts_dir,
        )

    tmp_path = years_dir / f"{varname}_temp_timeseries.nc"

    for part in range(1, config.NUM_LONG_BANDS + 1):
        fname_final = (
            files_dir
            / f"{varname}_{config.PROJECT_NAME}_{config.TS_START_YEAR}-{config.TS_END_YEAR}_p{part}.nc"
        )

        if fname_final.exists():
            print(f"  Part {part} timeseries already present, skipping")
            continue

        part_files = sorted(parts_dir.glob(f"{varname}_*_part{part}.nc"))
        if not part_files:
            print(f"  Warning: no part {part} files found for {varname}, skipping")
            continue

        print(f"  Concatenating {len(part_files)} year files → part {part} timeseries")
        try:
            subprocess.run(
                ["ncrcat"] + [str(f) for f in part_files] + [str(tmp_path)],
                check=True,
            )
            subprocess.run(
                ["nccopy", "-k", "classic", str(tmp_path), str(fname_final)],
                check=True,
            )
        finally:
            if tmp_path.exists():
                tmp_path.unlink()

    print(f"\nDone: timeseries files written to {files_dir}")


def extend_timeseries(varname, extend_start, extend_end, mode):
    """
    Extend an existing timeseries with new yearly files.

    mode='parts-only': create parts in parts_extend/, leave timeseries untouched.
    mode='append':     create parts, then append them to existing timeseries files.
    """
    years_dir  = config.YEARS_DIR
    parts_dir  = years_dir / "parts_extend"
    files_dir  = config.FILES_DIR

    print(f"\nExtending {varname} timeseries with years {extend_start}-{extend_end} "
          f"(mode: {mode})")

    _slice_years_into_parts(varname, range(extend_start, extend_end + 1), parts_dir)

    if mode == "parts-only":
        print(f"\nDone: new parts written to {parts_dir}")
        return

    # mode == "append": find each existing timeseries file and append new parts
    tmp_path = years_dir / f"{varname}_temp_extend.nc"

    # Existing timeseries files may have any end-year in their name; find them.
    existing_files = sorted(files_dir.glob(
        f"{varname}_{config.PROJECT_NAME}_{config.TS_START_YEAR}-*_p*.nc"
    ))
    if not existing_files:
        print(f"  Warning: no existing timeseries files found in {files_dir} — "
              "run normal make_timeseries first")
        return

    for existing in existing_files:
        # Extract part number from filename (e.g. _p7.nc → 7)
        stem = existing.stem
        part = int(stem.rsplit("_p", 1)[1])

        new_parts = sorted(parts_dir.glob(f"{varname}_*_part{part}.nc"))
        if not new_parts:
            print(f"  Part {part}: no new parts found, skipping")
            continue

        # Derive new filename with updated end year
        new_fname = (
            files_dir
            / f"{varname}_{config.PROJECT_NAME}_{config.TS_START_YEAR}-{extend_end}_p{part}.nc"
        )

        print(f"  Appending {len(new_parts)} new year(s) to part {part} → {new_fname.name}")
        try:
            subprocess.run(
                ["ncrcat", str(existing)] + [str(f) for f in new_parts] + [str(tmp_path)],
                check=True,
            )
            subprocess.run(
                ["nccopy", "-k", "classic", str(tmp_path), str(new_fname)],
                check=True,
            )
        finally:
            if tmp_path.exists():
                tmp_path.unlink()

        # Rename old file to .bak (don't silently delete)
        bak = existing.with_suffix(".nc.bak")
        existing.rename(bak)
        print(f"  Old file backed up as {bak.name}")

    print(f"\nDone: extended timeseries files written to {files_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create or extend longitude-band timeseries from yearly files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("varname", help="Variable name (e.g. precip)")
    parser.add_argument(
        "--no-parts",
        action="store_true",
        help="(Normal mode) Skip slicing into parts; assume parts already exist",
    )
    parser.add_argument(
        "--extend-start",
        type=int,
        metavar="YEAR",
        help="First year of new RACMO data to add (triggers extend mode)",
    )
    parser.add_argument(
        "--extend-end",
        type=int,
        metavar="YEAR",
        help="Last year of new RACMO data to add",
    )
    parser.add_argument(
        "--extend-mode",
        choices=["parts-only", "append"],
        default="parts-only",
        help=(
            "parts-only: create new parts in parts_extend/ only; "
            "append: also append new parts to existing timeseries files "
            "(default: parts-only)"
        ),
    )
    args = parser.parse_args()

    if args.extend_start or args.extend_end:
        if not (args.extend_start and args.extend_end):
            parser.error("--extend-start and --extend-end must both be provided")
        extend_timeseries(args.varname, args.extend_start, args.extend_end, args.extend_mode)
    else:
        make_timeseries(args.varname, do_parts=not args.no_parts)
