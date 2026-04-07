"""
make_racmo_monthly.py
---------------------
For one RACMO variable, merge all longitude-band files, resample to monthly
means, and write a single (time, rlat, rlon) NetCDF file.

Usage:
    python3 make_racmo_monthly.py --var ff10m
    python3 make_racmo_monthly.py --var precip --output-dir /path/to/output
"""

import argparse
import glob
import os
import sys

import numpy as np
import xarray as xr

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
TS_DIR     = "/home/nld4814/scratch/FGRN055_era055/input/timeseries"
OUTPUT_DIR = "/home/nld4814/scratch/FGRN055_era055/input/monthly/"
VARIABLES  = ["ff10m", "precip", "evap", "sndiv", "snowfall", "snowmelt", "tskin"]


def make_monthly(var: str, ts_dir: str, output_dir: str) -> None:
    files = sorted(glob.glob(f"{ts_dir}/{var}_*.nc"))
    if not files:
        print(f"ERROR: no files found for variable '{var}' in {ts_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Variable : {var}")
    print(f"Bands    : {len(files)}")
    print(f"Output   : {output_dir}")

    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{var}_FGRN055_era055_1939-2023_monthly.nc")

    if os.path.exists(output_file):
        print(f"Output already exists, skipping: {output_file}")
        return

    monthly_bands = []
    for i, f in enumerate(files):
        print(f"  [{i+1}/{len(files)}] {os.path.basename(f)}", flush=True)
        with xr.open_dataset(f) as ds:
            da = ds[var].squeeze("height", drop=True)
            monthly = da.resample(time="ME").mean()
            monthly_bands.append(monthly.load())

    print("Concatenating bands along rlon...", flush=True)
    result = xr.concat(monthly_bands, dim="rlon")

    # Restore lat/lon coordinates (take from last opened file — same rlat grid)
    # They were dropped because each band has different lat/lon values;
    # re-attach by concatenating them too.
    lat_bands = []
    lon_bands = []
    for f in files:
        with xr.open_dataset(f) as ds:
            lat_bands.append(ds["lat"].values)  # (rlat, rlon_band)
            lon_bands.append(ds["lon"].values)

    lat_full = np.concatenate(lat_bands, axis=1)  # (rlat, rlon_full)
    lon_full = np.concatenate(lon_bands, axis=1)

    result = result.assign_coords(
        lat=(("rlat", "rlon"), lat_full),
        lon=(("rlat", "rlon"), lon_full),
    )

    # Encoding: use float32 + zlib compression to keep file size small
    encoding = {var: {"dtype": "float32", "zlib": True, "complevel": 4}}

    print(f"Writing {output_file} ...", flush=True)
    result.to_netcdf(output_file, encoding=encoding)
    print("Done.")


def main():
    parser = argparse.ArgumentParser(description="Create monthly RACMO means")
    parser.add_argument("--var", required=True, choices=VARIABLES,
                        help="Variable to process")
    parser.add_argument("--ts-dir", default=TS_DIR,
                        help="Directory containing the lon-band timeseries files")
    parser.add_argument("--output-dir", default=OUTPUT_DIR,
                        help="Directory to write monthly output files")
    args = parser.parse_args()

    make_monthly(args.var, args.ts_dir, args.output_dir)


if __name__ == "__main__":
    main()
