"""
Utility functions for IMAU-FDM 1D to gridded maps post-processing.

Contains functions for:
- Loading pointlist and mask files
- Detrending timeseries (spinup-aware)
- Resampling timeseries to different aggregations
- Creating output datasets
"""

import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta
from pathlib import Path

# ---------------------------------------------------------------------------
# EPSG:3413 — WGS 84 / NSIDC Sea Ice Polar Stereographic North
# Verified against FGRN055_Masks.nc: projecting lat/lon through EPSG:3413
# gives zero residuals vs the stored X/Y arrays.
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
}


def load_pointlist(pointlist_file):
    """
    Load and parse the pointlist file.

    The pointlist maps point IDs (1-based, corresponding to 1D file numbers)
    to grid coordinates.

    Parameters
    ----------
    pointlist_file : str or Path
        Path to the pointlist file (IN_ll_FGRN055.txt)

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: point_id, lon, lat, rlat_idx, rlon_idx
    """
    # Read the comma-separated file
    df = pd.read_csv(
        pointlist_file,
        header=None,
        names=['lon', 'lat', 'field2', 'field3', 'field4', 'rlat_idx', 'rlon_idx'],
        skipinitialspace=True
    )

    # Add point_id (1-based index matching file naming)
    df['point_id'] = (df.index + 1).astype(int)

    # Convert grid indices to integers
    df['rlat_idx'] = df['rlat_idx'].astype(int)
    df['rlon_idx'] = df['rlon_idx'].astype(int)

    return df[['point_id', 'lon', 'lat', 'rlat_idx', 'rlon_idx']]


def load_mask(mask_file):
    """
    Load the mask file containing grid coordinates and ice mask.

    Parameters
    ----------
    mask_file : str or Path
        Path to the mask NetCDF file

    Returns
    -------
    xr.Dataset
        Dataset containing IceMask, lat, lon, rlat, rlon, and rotated_pole
    """
    ds = xr.open_dataset(mask_file)
    return ds


def get_grid_dimensions(mask_file):
    """
    Get grid dimensions from the mask file.

    Parameters
    ----------
    mask_file : str or Path
        Path to the mask NetCDF file

    Returns
    -------
    tuple
        (n_rlat, n_rlon) grid dimensions
    """
    with xr.open_dataset(mask_file) as ds:
        n_rlat = len(ds['rlat'])
        n_rlon = len(ds['rlon'])
    return n_rlat, n_rlon


def create_time_array(start_date, end_date, timestep_seconds=86400):
    """
    Create a time array for the model run period.

    Parameters
    ----------
    start_date : datetime
        Model start date
    end_date : datetime
        Model end date
    timestep_seconds : int
        Timestep in seconds (default: 86400 for daily)

    Returns
    -------
    np.ndarray
        Array of datetime objects
    """
    timestep = timedelta(seconds=timestep_seconds)
    times = []
    current = start_date
    while current <= end_date:
        times.append(current)
        current += timestep
    return np.array(times)


def datetime_to_fractional_year(dt):
    """
    Convert datetime to fractional year.

    Parameters
    ----------
    dt : datetime or np.ndarray of datetime
        Datetime(s) to convert

    Returns
    -------
    float or np.ndarray
        Fractional year(s)
    """
    if isinstance(dt, np.ndarray):
        return np.array([datetime_to_fractional_year(d) for d in dt])

    year = dt.year
    year_start = datetime(year, 1, 1)
    year_end = datetime(year + 1, 1, 1)
    year_fraction = (dt - year_start).total_seconds() / (year_end - year_start).total_seconds()
    return year + year_fraction


def detrend_spinup_aware(data, time_indices, spinup_end_idx):
    """
    Apply spinup-aware detrending to a timeseries.

    This removes the linear trend during the spinup period while preserving
    the trend after spinup. A continuity jump is applied to prevent
    discontinuity at the spinup boundary.

    Based on the Antarctica implementation in ANT27_make_1D_ERA5_FDMv12A.py.

    Parameters
    ----------
    data : np.ndarray
        1D array of data values
    time_indices : np.ndarray
        Array of time indices (0, 1, 2, ...)
    spinup_end_idx : int
        Index where spinup period ends

    Returns
    -------
    np.ndarray
        Detrended data
    """
    detrended = data.copy()
    ntime = len(data)

    if spinup_end_idx <= 0 or spinup_end_idx >= ntime:
        # No valid spinup period, return original data
        return detrended

    # Calculate trend during spinup period
    # rc0 = (data[spinup_end] - data[0]) / spinup_end
    rc0 = (data[spinup_end_idx] - data[0]) / spinup_end_idx

    # Remove trend from spinup period
    for i in range(spinup_end_idx):
        detrended[i] = data[i] - rc0 * i

    # Calculate jump at spinup boundary to maintain continuity
    jump = detrended[spinup_end_idx - 1] - data[spinup_end_idx - 1]

    # Apply jump to post-spinup data (preserves original trend after spinup)
    for i in range(spinup_end_idx, ntime):
        detrended[i] = data[i] + jump

    return detrended


def find_spinup_end_index(times, spinup_end_date, timestep_seconds=86400):
    """
    Find the index corresponding to the spinup end date.

    Parameters
    ----------
    times : np.ndarray
        Array of datetime objects or fractional years
    spinup_end_date : datetime
        Date when spinup period ends
    timestep_seconds : int
        Timestep in seconds

    Returns
    -------
    int
        Index of the spinup end
    """
    if isinstance(times[0], datetime):
        # Find index of first time >= spinup_end_date
        for i, t in enumerate(times):
            if t >= spinup_end_date:
                return i
        return len(times)
    else:
        # Assume fractional years
        spinup_end_year = datetime_to_fractional_year(spinup_end_date)
        for i, t in enumerate(times):
            if t >= spinup_end_year:
                return i
        return len(times)


def resample_timeseries(data, times, method='10day', start_date=None, aggregation='mean'):
    """
    Resample a timeseries to a different time aggregation.

    Parameters
    ----------
    data : np.ndarray
        1D array of data values (daily)
    times : np.ndarray
        Array of datetime objects or indices
    method : str
        Resampling method: 'daily', '10day', or 'monthly'
    start_date : datetime, optional
        Start date for time axis (used if times are indices)
    aggregation : str
        Aggregation method: 'mean' for rates/state variables, 'sum' for flux totals

    Returns
    -------
    tuple
        (resampled_data, resampled_times) as numpy arrays
    """
    if method == 'daily':
        return data, times

    # Convert times to datetime if needed
    if isinstance(times[0], (int, np.integer)):
        if start_date is None:
            raise ValueError("start_date required when times are indices")
        times = np.array([start_date + timedelta(days=int(t)) for t in times])

    # Create pandas Series for easy resampling
    series = pd.Series(data, index=pd.DatetimeIndex(times))

    # Choose aggregation function
    if aggregation == 'sum':
        agg_func = 'sum'
    else:
        agg_func = 'mean'

    if method == '10day':
        # Resample to 10-day periods
        resampled = series.resample('10D').agg(agg_func)
    elif method == 'monthly':
        resampled = series.resample('ME').agg(agg_func)
    else:
        raise ValueError(f"Unknown resampling method: {method}")

    # Remove any NaN values from incomplete periods at the end
    resampled = resampled.dropna()

    return resampled.values, resampled.index.to_pydatetime()


def get_output_time_axis(start_date, end_date, method='10day'):
    """
    Create the output time axis for the gridded data.

    Parameters
    ----------
    start_date : datetime
        Model start date
    end_date : datetime
        Model end date
    method : str
        Time aggregation method: 'daily', '10day', or 'monthly'

    Returns
    -------
    tuple
        (time_values, time_units, calendar)
        time_values as fractional years for compatibility with existing tools
    """
    # Create full daily time array
    daily_times = create_time_array(start_date, end_date)

    if method == 'daily':
        output_times = daily_times
    elif method == '10day':
        # Create 10-day means starting from start_date
        output_times = []
        current = start_date
        while current <= end_date:
            output_times.append(current + timedelta(days=5))  # Center of 10-day window
            current += timedelta(days=10)
        output_times = np.array(output_times)
    elif method == 'monthly':
        # Create monthly means (use mid-month)
        output_times = []
        current = start_date.replace(day=1)
        while current <= end_date:
            # Middle of month (roughly day 15)
            mid_month = current.replace(day=15)
            if mid_month >= start_date and mid_month <= end_date:
                output_times.append(mid_month)
            # Move to next month
            if current.month == 12:
                current = current.replace(year=current.year + 1, month=1)
            else:
                current = current.replace(month=current.month + 1)
        output_times = np.array(output_times)
    else:
        raise ValueError(f"Unknown method: {method}")

    # Convert to fractional years for compatibility
    time_values = datetime_to_fractional_year(output_times)

    return time_values, output_times


def add_crs_to_dataset(ds, mask_source=None):
    """
    Add EPSG:3413 CRS metadata to an xarray Dataset in place.

    Adds a 'crs' grid-mapping variable, 2D 'x'/'y' coordinate arrays in
    metres, and 'grid_mapping'/'coordinates' attributes on every data
    variable that has a 'time' dimension.

    Parameters
    ----------
    ds : xr.Dataset
        Output dataset to modify.
    mask_source : str, Path, xr.Dataset, or None
        Source for the X/Y coordinate arrays (must contain 'X' and 'Y' in km).
        - Path / str  : path to the mask .nc file
        - xr.Dataset  : already-loaded mask dataset
        - None (default): look up automatically from config.py using DOMAIN:
              <BASE_DIR>/reference/<DOMAIN>/<DOMAIN>_Masks.nc
    """
    if mask_source is None:
        from config import BASE_DIR, DOMAIN
        mask_source = BASE_DIR / 'reference' / DOMAIN / f'{DOMAIN}_Masks.nc'

    if isinstance(mask_source, (str, Path)):
        mask_path = Path(mask_source)
        if not mask_path.exists():
            try:
                from config import BASE_DIR, DOMAIN
                hint = (
                    f"  DOMAIN   = '{DOMAIN}'  (config.py)\n"
                    f"  BASE_DIR = '{BASE_DIR}'  (config.py)\n"
                    f"Change DOMAIN or BASE_DIR in config.py, or pass mask_source explicitly."
                )
            except ImportError:
                hint = "Pass mask_source explicitly (config.py could not be imported)."
            raise FileNotFoundError(
                f"Mask file not found: {mask_path}\n"
                f"Expected pattern: <BASE_DIR>/reference/<DOMAIN>/<DOMAIN>_Masks.nc\n"
                f"{hint}"
            )
        mask_ds = xr.open_dataset(mask_path)
        close_after = True
    else:
        mask_ds = mask_source
        close_after = False

    try:
        x_m = (mask_ds['X'].values * 1000.0).astype(np.float32)  # km → m
        y_m = (mask_ds['Y'].values * 1000.0).astype(np.float32)
    finally:
        if close_after:
            mask_ds.close()

    # CRS scalar variable
    ds['crs'] = xr.DataArray(np.int32(0), attrs=_CRS_ATTRS)

    # x / y 2-D coordinate arrays
    ds['x'] = xr.DataArray(
        x_m, dims=['rlat', 'rlon'],
        attrs={'standard_name': 'projection_x_coordinate',
               'units': 'm', 'grid_mapping': 'crs'},
    )
    ds['y'] = xr.DataArray(
        y_m, dims=['rlat', 'rlon'],
        attrs={'standard_name': 'projection_y_coordinate',
               'units': 'm', 'grid_mapping': 'crs'},
    )

    # Link data variables to the CRS and lat/lon auxiliary coordinates
    for vname, var in ds.data_vars.items():
        if 'time' in var.dims and var.ndim > 1:
            ds[vname].attrs['grid_mapping'] = 'crs'
            ds[vname].attrs['coordinates']  = 'lat lon'

    ds.attrs['Conventions'] = 'CF-1.8'


def create_output_dataset(var_name, data, time_values, mask_ds, var_metadata,
                          detrended=False, timestep='10day'):
    """
    Create an xarray Dataset for output.

    Parameters
    ----------
    var_name : str
        Variable name
    data : np.ndarray
        3D data array (time, rlat, rlon)
    time_values : np.ndarray
        Time coordinate values (fractional years)
    mask_ds : xr.Dataset
        Mask dataset containing lat, lon, rlat, rlon coordinates
    var_metadata : dict
        Variable metadata (long_name, units)
    detrended : bool
        Whether detrending was applied
    timestep : str
        Time aggregation used

    Returns
    -------
    xr.Dataset
        Output dataset ready for saving
    """
    # Create the dataset
    ds = xr.Dataset(
        data_vars={
            var_name: (['time', 'rlat', 'rlon'], data.astype(np.float32)),
            'lat': (['rlat', 'rlon'], mask_ds['lat'].values),
            'lon': (['rlat', 'rlon'], mask_ds['lon'].values),
        },
        coords={
            'time': time_values,
            'rlat': mask_ds['rlat'].values,
            'rlon': mask_ds['rlon'].values,
        }
    )

    # Add rotated pole if available
    if 'rotated_pole' in mask_ds:
        ds['rotated_pole'] = mask_ds['rotated_pole']

    # Set variable attributes
    ds[var_name].attrs = {
        'long_name': var_metadata.get('long_name', var_name),
        'units': var_metadata.get('units', ''),
        'detrended': str(detrended),
        '_FillValue': np.float32(np.nan),
    }

    # Set coordinate attributes
    ds['time'].attrs = {
        'long_name': 'Time in fractional years',
        'units': 'years',
    }
    ds['lat'].attrs = {
        'long_name': 'latitude',
        'units': 'degrees_north',
    }
    ds['lon'].attrs = {
        'long_name': 'longitude',
        'units': 'degrees_east',
    }

    # Set global attributes
    ds.attrs = {
        'title': f'IMAU-FDM gridded output: {var_metadata.get("long_name", var_name)}',
        'source': 'IMAU-FDM version 1.2+',
        'domain': 'FGRN055',
        'institution': 'IMAU, Utrecht University',
        'history': f'Created on {datetime.now().isoformat()}',
        'time_aggregation': timestep,
        'Conventions': 'CF-1.8',
    }

    # Add CRS metadata (requires X/Y arrays in the mask dataset)
    if 'X' in mask_ds and 'Y' in mask_ds:
        add_crs_to_dataset(ds, mask_ds)

    return ds


def read_1d_file(filepath, var_name):
    """
    Read a single variable from a 1D output file.

    Parameters
    ----------
    filepath : str or Path
        Path to the 1D NetCDF file
    var_name : str
        Variable name to read

    Returns
    -------
    np.ndarray or None
        1D array of data values, or None if file doesn't exist or error
    """
    try:
        with xr.open_dataset(filepath) as ds:
            data = ds[var_name].values
        return data
    except (FileNotFoundError, KeyError, OSError):
        return None
