"""Add EPSG:3413 CRS metadata to IMAU-FDM post-processed output files.

Adds to every .nc file in the post-process directory:
  - a 'crs' scalar variable (CF grid mapping + WKT for GDAL/QGIS)
  - 2D 'x' and 'y' coordinate arrays in metres (EPSG:3413)
  - 'grid_mapping' and 'coordinates' attributes on data variables
  - Conventions updated to CF-1.8

Each file is rewritten using 'r' (read) + 'w' (create) modes — never 'r+'.
HDF5 requires file locking for 'r+', which fails on NFS.  'r' and 'w' need
no locks and work reliably on NFS.
"""

import os
os.environ.setdefault('HDF5_USE_FILE_LOCKING', 'FALSE')

import shutil
import argparse
import netCDF4 as nc
import numpy as np
from pathlib import Path
from datetime import datetime

# ---------------------------------------------------------------------------
# EPSG:3413 — WGS 84 / NSIDC Sea Ice Polar Stereographic North
# Verified: mask lat/lon projected through EPSG:3413 matches stored X/Y to
# floating-point precision.  The mask global attribute claiming lon_0=-39°
# and lat_ts=71° is incorrect; the true CRS is lon_0=-45°, lat_ts=70°.
# ---------------------------------------------------------------------------

_CRS_WKT = (
    'PROJCS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",'
    'GEOGCS["WGS 84",'
    'DATUM["WGS_1984",'
    'SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],'
    'AUTHORITY["EPSG","6326"]],'
    'PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],'
    'UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],'
    'AUTHORITY["EPSG","4326"]],'
    'PROJECTION["Polar_Stereographic"],'
    'PARAMETER["latitude_of_origin",70],'
    'PARAMETER["central_meridian",-45],'
    'PARAMETER["scale_factor",1],'
    'PARAMETER["false_easting",0],'
    'PARAMETER["false_northing",0],'
    'UNIT["metre",1,AUTHORITY["EPSG","9001"]],'
    'AXIS["Easting",SOUTH],'
    'AXIS["Northing",SOUTH],'
    'AUTHORITY["EPSG","3413"]]'
)

_CRS_ATTRS = {
    'grid_mapping_name':                     'polar_stereographic',
    'epsg_code':                             'EPSG:3413',
    'latitude_of_projection_origin':         90.0,
    'straight_vertical_longitude_from_pole': -45.0,
    'standard_parallel':                     70.0,
    'false_easting':                         0.0,
    'false_northing':                        0.0,
    'semi_major_axis':                       6378137.0,
    'inverse_flattening':                    298.257223563,
    'crs_wkt':                               _CRS_WKT,
    'proj4_params':                          (
        '+proj=stere +lat_0=90 +lon_0=-45 +lat_ts=70 '
        '+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
    ),
}

# ---------------------------------------------------------------------------


def add_crs_to_all(output_dir, mask_file, dry_run=False):
    """Add CRS metadata to every .nc file in output_dir."""
    output_dir = Path(output_dir)

    # Clean up any temp files left by a previously interrupted run.
    for stale in output_dir.glob('*.tmp.nc'):
        stale.unlink()
        print(f'Removed stale temp file: {stale.name}')

    nc_files = sorted(
        f for f in output_dir.glob('*.nc')
        if not f.name.endswith('.tmp.nc')
    )
    if not nc_files:
        print(f'No .nc files found in {output_dir}')
        return

    print('Loading X/Y from mask file ...')
    with nc.Dataset(mask_file) as ds:
        x_m = np.array(ds.variables['X'][:] * 1000.0, dtype=np.float32)  # km → m
        y_m = np.array(ds.variables['Y'][:] * 1000.0, dtype=np.float32)

    print(f'Processing {len(nc_files)} file(s) ...\n')
    for f in nc_files:
        if dry_run:
            print(f'[dry-run] {f.name}')
        else:
            print(f'  {f.name}', flush=True)
            _update_file(f, x_m, y_m)
            print('    done.')

    print('\nFinished.')


def _update_file(filepath, x_m, y_m):
    """Rewrite filepath with CRS metadata added.

    Opens the source in 'r' mode and writes a new file in 'w' mode —
    never 'r+'.  This avoids HDF5 file-locking errors on NFS.
    The original is replaced only after a clean, complete write.
    """
    tmp = filepath.with_suffix('.tmp.nc')
    try:
        with nc.Dataset(filepath, 'r') as src:
            with nc.Dataset(tmp, 'w', format=src.file_format) as dst:
                _copy_nc(src, dst)
                _add_crs_metadata(dst, x_m, y_m)
        shutil.move(str(tmp), str(filepath))
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise


def _copy_nc(src, dst):
    """Copy all global attributes, dimensions, and variables from src to dst."""
    dst.setncatts({k: src.getncattr(k) for k in src.ncattrs()})

    for dname, dim in src.dimensions.items():
        dst.createDimension(dname, None if dim.isunlimited() else len(dim))

    for vname, sv in src.variables.items():
        filt   = sv.filters() or {}
        chunks = sv.chunking()
        fill   = sv.getncattr('_FillValue') if '_FillValue' in sv.ncattrs() else False

        dv = dst.createVariable(
            vname, sv.datatype, sv.dimensions,
            zlib=filt.get('zlib', False),
            complevel=filt.get('complevel', 4),
            chunksizes=None if chunks in (None, 'contiguous') else list(chunks),
            fill_value=fill,
        )
        dv.setncatts({k: sv.getncattr(k) for k in sv.ncattrs() if k != '_FillValue'})

        # Stream large time-indexed arrays in chunks to limit memory use.
        if 'time' in sv.dimensions and sv.ndim > 1:
            t_ax = sv.dimensions.index('time')
            for i in range(0, sv.shape[t_ax], 100):
                slc = tuple(
                    slice(i, min(i + 100, sv.shape[t_ax])) if j == t_ax else slice(None)
                    for j in range(sv.ndim)
                )
                dv[slc] = sv[slc]
        elif sv.ndim == 0:
            dv.assignValue(sv.getValue())
        else:
            dv[:] = sv[:]


def _add_crs_metadata(ds, x_m, y_m):
    """Add crs variable, x/y coordinates, and grid_mapping attributes."""
    # CRS scalar variable
    if 'crs' not in ds.variables:
        crs_var = ds.createVariable('crs', 'i4')
    else:
        crs_var = ds.variables['crs']
    for k, v in _CRS_ATTRS.items():
        crs_var.setncattr(k, v)

    # x / y 2-D coordinate arrays (rlat × rlon, metres)
    for name, data, sname in [
        ('x', x_m, 'projection_x_coordinate'),
        ('y', y_m, 'projection_y_coordinate'),
    ]:
        if name not in ds.variables:
            coord = ds.createVariable(name, 'f4', ('rlat', 'rlon'), fill_value=False)
        else:
            coord = ds.variables[name]
        coord[:] = data
        coord.standard_name = sname
        coord.units = 'm'
        coord.grid_mapping = 'crs'

    # Link data variables to the CRS and lat/lon auxiliary coordinates
    for vname, var in ds.variables.items():
        if 'time' in var.dimensions and var.ndim > 1:
            var.grid_mapping = 'crs'
            var.coordinates  = 'lat lon'

    ds.setncattr('Conventions', 'CF-1.8')
    old_hist = ds.getncattr('history') if 'history' in ds.ncattrs() else ''
    ds.setncattr(
        'history',
        old_hist + f'\n{datetime.now().isoformat()}: Added EPSG:3413 CRS (add_crs_metadata.py)',
    )


# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import sys
    parser = argparse.ArgumentParser(
        description='Add EPSG:3413 CRS to IMAU-FDM post-process files.'
    )
    parser.add_argument('--dry-run', action='store_true',
                        help='List files without modifying them.')
    args = parser.parse_args()

    sys.path.insert(0, str(Path(__file__).parent))
    from create_1D_2D_2Ddetail_files.config import MASK_FILE, OUTPUT_DIR

    add_crs_to_all(OUTPUT_DIR, MASK_FILE, dry_run=args.dry_run)
