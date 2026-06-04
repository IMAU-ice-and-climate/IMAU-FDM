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


def _load_rotated_pole_attrs(grid_file):
    """Read rotated_pole attributes from the reference grid NetCDF file."""
    with nc.Dataset(grid_file, 'r') as gds:
        rp = gds.variables['rotated_pole']
        return {k: rp.getncattr(k) for k in rp.ncattrs()}


def _build_crs_attrs(epsg_code):
    """Generate CF-compliant CRS attributes from an EPSG code using pyproj."""
    from pyproj import CRS
    crs = CRS.from_user_input(epsg_code)
    attrs = crs.to_cf()   # includes crs_wkt, grid_mapping_name, standard_parallel, etc.
    attrs['epsg_code'] = epsg_code
    return attrs


def _load_crs_attrs(grid_file):
    """Read crs attributes from the reference grid NetCDF file."""
    with nc.Dataset(grid_file, 'r') as gds:
        crs_var = gds.variables['crs']
        return {k: crs_var.getncattr(k) for k in crs_var.ncattrs()}


def _load_update_metadata(mask_file, grid_file):
    """Load X/Y arrays and CRS attrs needed to update a file."""
    with nc.Dataset(mask_file) as ds:
        x_m = np.array(ds.variables['X'][:] * 1000.0, dtype=np.float32)  # km → m
        y_m = np.array(ds.variables['Y'][:] * 1000.0, dtype=np.float32)
    rotated_pole_attrs = _load_rotated_pole_attrs(grid_file)
    crs_attrs          = _load_crs_attrs(grid_file)
    return x_m, y_m, rotated_pole_attrs, crs_attrs


def add_crs_to_file(filepath, mask_file, grid_file, dry_run=False):
    """Add CRS metadata to a single .nc file."""
    filepath = Path(filepath)
    print('Loading X/Y and CRS attrs ...')
    x_m, y_m, rotated_pole_attrs, crs_attrs = _load_update_metadata(mask_file, grid_file)
    if dry_run:
        print(f'[dry-run] {filepath.name}')
    else:
        print(f'  {filepath.name}', flush=True)
        _update_file(filepath, x_m, y_m, rotated_pole_attrs, crs_attrs)
        print('    done.')


def add_crs_to_all(output_dir, mask_file, grid_file, dry_run=False):
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

    print('Loading X/Y and CRS attrs ...')
    x_m, y_m, rotated_pole_attrs, crs_attrs = _load_update_metadata(mask_file, grid_file)

    print(f'Processing {len(nc_files)} file(s) ...\n')
    for f in nc_files:
        if dry_run:
            print(f'[dry-run] {f.name}')
        else:
            print(f'  {f.name}', flush=True)
            _update_file(f, x_m, y_m, rotated_pole_attrs, crs_attrs)
            print('    done.')

    print('\nFinished.')


def _update_file(filepath, x_m, y_m, rotated_pole_attrs, crs_attrs):
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
                _add_crs_metadata(dst, x_m, y_m, rotated_pole_attrs, crs_attrs)
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


def _add_crs_metadata(ds, x_m, y_m, rotated_pole_attrs, crs_attrs):
    """Add/update rotated_pole, crs variable, x/y coordinates, and grid_mapping attributes."""
    # Rotated pole grid-mapping variable (attrs loaded from reference grid file)
    if 'rotated_pole' not in ds.variables:
        rp_var = ds.createVariable('rotated_pole', 'i4')
    else:
        rp_var = ds.variables['rotated_pole']
    for k, v in rotated_pole_attrs.items():
        rp_var.setncattr(k, v)

    # CRS scalar variable (attrs loaded from reference grid file)
    if 'crs' not in ds.variables:
        crs_var = ds.createVariable('crs', 'i4')
    else:
        crs_var = ds.variables['crs']
    for k, v in crs_attrs.items():
        crs_var.setncattr(k, v)

    epsg = crs_attrs.get('epsg_code', 'crs')
    # x / y 2-D coordinate arrays (rlat × rlon, metres)
    for name, data, sname, lname in [
        ('x', x_m, 'projection_x_coordinate', f'x coordinate ({epsg})'),
        ('y', y_m, 'projection_y_coordinate', f'y coordinate ({epsg})'),
    ]:
        if name not in ds.variables:
            coord = ds.createVariable(name, 'f4', ('rlat', 'rlon'), fill_value=False)
        else:
            coord = ds.variables[name]
        coord[:] = data
        coord.standard_name = sname
        coord.long_name     = lname
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


def create_grid_file(source_file, grid_file, epsg_code='EPSG:3413'):
    """Create (or overwrite) a reference grid file from a RACMO or FDM mask file.

    Extracts the rotated_pole variable and the rlat/rlon coordinate arrays
    from source_file and writes them to grid_file, together with a crs variable
    whose attributes are generated from epsg_code via pyproj.  The result is the
    single source of truth for projection metadata used by add_crs_to_all() and
    utils.add_crs_to_dataset().

    source_file must contain:
      - a 'rotated_pole' scalar variable with CF grid-mapping attributes
      - 1-D 'rlat' and 'rlon' coordinate variables (degrees)

    Parameters
    ----------
    source_file : str or Path
        RACMO NetCDF file (or any CF file with rotated_pole + rlat/rlon).
    grid_file : str or Path
        Output path for the new grid file.  Written via a temp file so
        NFS partial writes never corrupt an existing grid file.
    epsg_code : str
        EPSG code for the projected CRS stored in the output X/Y coordinates
        (default: 'EPSG:3413' — WGS 84 / NSIDC Sea Ice Polar Stereographic North).
    """
    source_file = Path(source_file)
    grid_file   = Path(grid_file)
    tmp         = grid_file.with_suffix('.tmp.nc')

    crs_attrs = _build_crs_attrs(epsg_code)

    try:
        with nc.Dataset(source_file, 'r') as src:
            # Validate required variables
            for required in ('rotated_pole', 'rlat', 'rlon'):
                if required not in src.variables:
                    raise ValueError(
                        f"'{required}' not found in {source_file.name}. "
                        "Source must be a CF file with rotated_pole, rlat, and rlon."
                    )

            rlat = src.variables['rlat']
            rlon = src.variables['rlon']
            rp   = src.variables['rotated_pole']

            with nc.Dataset(tmp, 'w', format='NETCDF4') as dst:
                # Global attributes
                dst.title   = f'{grid_file.stem} rotated pole grid coordinates'
                dst.source  = f'Extracted from {source_file.name}'
                dst.history = f'{datetime.now().isoformat()}: Created by create_grid_file() in add_crs_metadata.py'

                # Dimensions
                dst.createDimension('rlat', len(rlat))
                dst.createDimension('rlon', len(rlon))

                # rlat coordinate
                rlat_var = dst.createVariable('rlat', rlat.datatype, ('rlat',))
                rlat_var.setncatts({k: rlat.getncattr(k) for k in rlat.ncattrs()})
                rlat_var[:] = rlat[:]

                # rlon coordinate
                rlon_var = dst.createVariable('rlon', rlon.datatype, ('rlon',))
                rlon_var.setncatts({k: rlon.getncattr(k) for k in rlon.ncattrs()})
                rlon_var[:] = rlon[:]

                # rotated_pole scalar — copy all attrs, override long_name
                rp_var = dst.createVariable('rotated_pole', 'i4')
                rp_var.setncatts({k: rp.getncattr(k) for k in rp.ncattrs()})
                rp_var.long_name = 'RACMO Rotated Pole'

                # crs scalar — attrs generated from epsg_code via pyproj
                crs_var = dst.createVariable('crs', 'i4')
                for k, v in crs_attrs.items():
                    crs_var.setncattr(k, v)

        shutil.move(str(tmp), str(grid_file))
        print(f'Written: {grid_file}')

    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise


# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import sys
    parser = argparse.ArgumentParser(
        description='Add EPSG:3413 CRS to IMAU-FDM post-process files.'
    )
    parser.add_argument('--dry-run', action='store_true',
                        help='List files without modifying them (add-crs mode only).')
    parser.add_argument('--file', metavar='FILE',
                        help='Process a single .nc file instead of all files in OUTPUT_DIR.')
    parser.add_argument('--create-grid', metavar='SOURCE_FILE',
                        help='Create/overwrite the reference grid file from SOURCE_FILE '
                             '(a RACMO NetCDF with rotated_pole + rlat/rlon) instead of '
                             'running the CRS-update pass.')
    args = parser.parse_args()

    sys.path.insert(0, str(Path(__file__).parent))
    from create_1D_2D_2Ddetail_files.config import MASK_FILE, GRID_FILE, OUTPUT_DIR

    if args.create_grid:
        create_grid_file(args.create_grid, GRID_FILE)
    elif args.file:
        add_crs_to_file(args.file, MASK_FILE, GRID_FILE, dry_run=args.dry_run)
    else:
        add_crs_to_all(OUTPUT_DIR, MASK_FILE, GRID_FILE, dry_run=args.dry_run)
