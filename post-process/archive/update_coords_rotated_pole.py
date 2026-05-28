#!/usr/bin/env python3
"""
Update existing post-processed FDM output files to use proper rotated pole CRS.

Changes applied to each file:
  - rlat/rlon coordinate values replaced with actual degrees (from FGRN055_grid.nc)
  - rlat/rlon attributes updated to CF-1.6 convention (standard_name, units=degrees)
  - rotated_pole grid mapping variable added
  - y_FDM / x_FDM variables added (original 0-based FDM indices)
  - Data variable grid_mapping updated to 'rotated_pole'

Uses xarray assign_coords + r/w (temp file then rename) required on NFS.
"""

import os
import glob
import numpy as np
import xarray as xr

GRID_FILE = '/home/nld4814/perm/code/IMAU-FDM/reference/FGRN055/FGRN055_grid.nc'
POST_DIR  = '/home/nld4814/scratch/run_FGRN055-era055_1939-2023/post-process'


def update_file(src_path, grid_ds):
    tmp_path = src_path + '.tmp'
    rp_attrs = dict(grid_ds['rotated_pole'].attrs)
    try:
        ds = xr.open_dataset(src_path)

        # Replace rlat/rlon with actual degree coordinates (attrs come from grid_ds)
        ds.coords.update({"rlat": grid_ds["rlat"], "rlon": grid_ds["rlon"]})

        # Drop EPSG:3413 variables added by the previous CRS processing step
        ds = ds.drop_vars(['x', 'y', 'crs'], errors='ignore')

        # Add rotated pole grid mapping variable
        ds['rotated_pole'] = xr.Variable([], np.int32(0), attrs=rp_attrs)

        # Add FDM index variables
        ds['y_FDM'] = xr.Variable(
            ['rlat'], np.arange(len(grid_ds['rlat']), dtype=np.int32),
            attrs={'long_name': 'row index in IMAU-FDM grid (0-based)', 'units': '1'}
        )
        ds['x_FDM'] = xr.Variable(
            ['rlon'], np.arange(len(grid_ds['rlon']), dtype=np.int32),
            attrs={'long_name': 'column index in IMAU-FDM grid (0-based)', 'units': '1'}
        )

        # Update data variables: set grid_mapping, clear encoding conflicts
        for vname in ds.data_vars:
            if 'rlat' in ds[vname].dims and 'rlon' in ds[vname].dims:
                ds[vname].attrs['grid_mapping'] = 'rotated_pole'
                for key in ('coordinates', 'missing_value', '_FillValue'):
                    ds[vname].encoding.pop(key, None)
                ds[vname].attrs.pop('missing_value', None)

        ds.to_netcdf(tmp_path)
        ds.close()

        os.rename(tmp_path, src_path)

    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise


def main():
    grid_ds = xr.open_dataset(GRID_FILE)

    nc_files = sorted(f for f in glob.glob(os.path.join(POST_DIR, '*.nc'))
                      if not os.path.basename(f).startswith('_'))
    print(f"Found {len(nc_files)} files to update.")

    for i, fpath in enumerate(nc_files):
        fname = os.path.basename(fpath)
        print(f"[{i+1}/{len(nc_files)}] {fname} ...", flush=True)
        update_file(fpath, grid_ds)
        print(f"  done.", flush=True)

    print("\nAll files updated successfully.")


if __name__ == '__main__':
    main()
