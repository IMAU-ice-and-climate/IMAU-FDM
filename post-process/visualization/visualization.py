"""
Visualization functions for IMAU-FDM gridded output.

Table of Contents
-----------------
Data Loading
  load_gridded_data       Load a gridded variable from post-processed output

Plotting
  plot_map                Spatial map of a variable at one timestep (or mean)
  plot_timeseries         Time series at a specific lon/lat location
  plot_annual_mean        Spatial map of annual mean for a given year
  plot_temporal_mean      Spatial map averaged over a year range
  plot_temporal_mean_difference  Difference map between two period means

Multi-Variable Comparison
  plot_multi_var_map      Side-by-side maps of multiple variables at one time

Ice Sheet Facies Classification
  classify_zones          Classify ablation / percolation / dry-snow zones
  plot_zones              Map of zone classification for a given year
  plot_zones_timeseries   Panel of zone maps for multiple years

Statistics and Export
  print_stats             Print min/max/mean/std/NaN counts for a variable
  export_annual_means     Save annual means to a NetCDF file

Animation
  create_animation_gif    Create animated GIF stepping through timesteps

Masked Difference Plotting
  plot_masked_difference  Difference/ratio map filtered by a mask condition
  plot_process_comparison 2x2 overview of melt vs runoff process partitioning

Mask Overlays (add to existing axes — pass the same ds used for the map)
  add_elevation_contours  Elevation contour lines at a fixed interval (m)
  add_land_outline        Land/ice-sheet coastline from LSM_GR
  add_basin_outlines      Mouginot basin outlines + labels + coastline (all-in-one)

Spatial Mean Time Series
  plot_timeseries_gis           Area-weighted mean over the whole GIS
  plot_timeseries_by_basin      Area-weighted mean per Mouginot basin
  plot_timeseries_elevation_band  Area-weighted mean above/below an elevation
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap
from pathlib import Path

# Import local config
import config

# Okabe-Ito colour-blind safe palette (black last so it acts as a fallback)
OKABE_ITO_COLORS = [
    '#E69F00',  # orange
    '#56B4E9',  # sky blue
    '#009E73',  # bluish green
    '#F0E442',  # yellow
    '#0072B2',  # blue
    '#D55E00',  # vermillion
    '#CC79A7',  # reddish purple
    '#000000',  # black
]

# Diverging colormap: red (negative) → white → blue (positive).
# Use for any map centred on zero: dFAC anomalies, runoff change, difference maps, etc.
DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    'fdm_diverging', ['#D92B2B', '#FFFFFF', '#3275DA']
)

# Canonical colours for ice-sheet facies classification (Okabe-Ito safe)
FACIES_COLORS = [
    '#D55E00',  # ablation      — vermillion
    '#0072B2',  # percolation   — blue
    '#999999',  # dry snow      — grey
]

def _infer_var_name(ds):
    """Return the single FDM data variable in ds by matching against config.VARIABLES."""
    candidates = [v for v in ds.data_vars if v in config.VARIABLES]
    if len(candidates) == 1:
        return candidates[0]
    raise ValueError(
        f"Cannot infer variable: dataset has {len(candidates)} FDM candidates "
        f"({candidates}). Pass var_name explicitly."
    )


# =============================================================================
# Data Loading
# =============================================================================

def load_gridded_data(var_name, timestep='10day', detrended=None, output_dir=None,
                      date_tag=None):
    """
    Load a gridded variable from the processed output.

    Parameters
    ----------
    var_name : str
        Variable name (e.g., 'h_surf', 'FirnAir')
    timestep : str
        Time aggregation used ('daily', '10day', 'monthly')
    detrended : bool, optional
        Whether to load detrended version. If None, auto-detect.
    output_dir : Path, optional
        Output directory. If None, uses config.OUTPUT_DIR.
    date_tag : str, optional
        Date range string in the filename (e.g. '1939-2025').
        Defaults to config.DATE_TAG ('1939-2023').

    Returns
    -------
    xr.Dataset
        Loaded dataset
    """
    if output_dir is None:
        output_dir = config.OUTPUT_DIR

    if detrended is None:
        # Check if variable should be detrended
        detrended = var_name in config.DETREND_VARIABLES

    filename = config.get_output_filename(var_name, timestep, detrended, date_tag)
    filepath = output_dir / filename

    if not filepath.exists():
        if not detrended:
            # Never silently substitute detrended data when raw was requested —
            # that would produce misleading results.
            raise FileNotFoundError(
                f"Non-detrended file not found: {filepath}\n"
                f"Only the detrended version exists. Re-run post-processing "
                f"to generate the raw file, or load with detrended=True."
            )
        # detrended requested but missing → fall back to raw with a note
        alt_filename = config.get_output_filename(var_name, timestep, False, date_tag)
        alt_filepath = output_dir / alt_filename
        if alt_filepath.exists():
            print(f"Note: Loading {alt_filename} (detrended version not found)")
            filepath = alt_filepath
        else:
            raise FileNotFoundError(f"File not found: {filepath}")

    #print(f"Loading: {filepath.name}")
    return xr.open_dataset(filepath, chunks={})


# =============================================================================
# Plotting Functions
# =============================================================================

def plot_map(ds, var_name=None, time_idx=None, time_value=None, ax=None,
             cmap='viridis', vmin=None, vmax=None,
             add_colorbar=True, use_latlon=False):
    """
    Plot a spatial map of a variable using xarray's native plotting.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing the variable
    var_name : str
        Variable name to plot
    time_idx : int, optional
        Time index to plot
    time_value : float, optional
        Time value (fractional year) to plot. Takes precedence over time_idx.
    ax : matplotlib axes, optional
        Axes to plot on
    cmap : str
        Colormap name
    vmin, vmax : float, optional
        Color scale limits
    add_colorbar : bool
        Whether to add a colorbar
    use_latlon : bool
        If True, use lat/lon coordinates for axes; otherwise use rlat/rlon indices

    Returns
    -------
    matplotlib axes
    """
    if var_name is None:
        var_name = _infer_var_name(ds)
    # Select time slice and build title
    if time_value is not None:
        data = ds[var_name].sel(time=time_value, method='nearest')
        time_str = f"{float(data.time):.2f}"
    elif time_idx is not None:
        data = ds[var_name].isel(time=time_idx)
        time_str = f"{float(data.time):.2f}"
    else:
        data = ds[var_name].mean(dim='time')
        time_str = "mean"

    # Build title from variable metadata
    long_name = ds[var_name].attrs.get('long_name', var_name)
    units = ds[var_name].attrs.get('units', '')
    title = f"{long_name} ({time_str})"

    # Create figure if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 10))

    # Colorbar settings - shrink to match plot height
    cbar_kwargs = {'shrink': 0.6, 'aspect': 20, 'pad': 0.02}
    if units:
        cbar_kwargs['label'] = units

    # Plot using xarray
    if use_latlon:
        # Use lat/lon coordinates
        data_plot = data.assign_coords(
            x=(['rlat', 'rlon'], ds['lon'].values),
            y=(['rlat', 'rlon'], ds['lat'].values)
        )
        p = data_plot.plot(
            ax=ax, x='x', y='y',
            cmap=cmap, vmin=vmin, vmax=vmax,
            add_colorbar=add_colorbar,
            cbar_kwargs=cbar_kwargs if add_colorbar else None
        )
        ax.set_xlabel('Longitude', fontsize=16)
        ax.set_ylabel('Latitude', fontsize=16)
    else:
        # Use native rlat/rlon grid indices
        p = data.plot(
            ax=ax,
            cmap=cmap, vmin=vmin, vmax=vmax,
            add_colorbar=add_colorbar,
            cbar_kwargs=cbar_kwargs if add_colorbar else None
        )
        ax.set_xlabel('rlon (grid index)', fontsize=16)
        ax.set_ylabel('rlat (grid index)', fontsize=16)

    if add_colorbar and hasattr(p, 'colorbar') and p.colorbar is not None:
        if units:
            p.colorbar.set_label(units, fontsize=16)
        p.colorbar.ax.tick_params(labelsize=12)

    ax.set_aspect('equal')
    ax.set_title(title, fontsize=20)

    return ax


def plot_timeseries(ds, lon_point, lat_point, var_name=None, ax=None, label=None):
    """
    Plot time series at a specific location.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing the variable
    var_name : str
        Variable name to plot
    lon_point, lat_point : float
        Coordinates of the point to extract
    ax : matplotlib axes, optional
        Axes to plot on
    label : str, optional
        Line label

    Returns
    -------
    matplotlib axes
    """
    if var_name is None:
        var_name = _infer_var_name(ds)
    # Find nearest grid point
    lon_2d = ds['lon'].values
    lat_2d = ds['lat'].values

    # Calculate distance to all points
    dist = np.sqrt((lon_2d - lon_point)**2 + (lat_2d - lat_point)**2)
    min_idx = np.unravel_index(np.nanargmin(dist), dist.shape)

    # Extract time series
    data = ds[var_name].isel(rlat=min_idx[0], rlon=min_idx[1])
    actual_lon = lon_2d[min_idx]
    actual_lat = lat_2d[min_idx]

    # Create figure if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=(14, 5))
    ax.set_prop_cycle(color=OKABE_ITO_COLORS)

    # Plot
    if label is None:
        label = f"({actual_lon:.2f}°, {actual_lat:.2f}°)"

    ax.plot(ds['time'], data, label=label, linewidth=0.8)

    # Labels
    units = ds[var_name].attrs.get('units', '')
    long_name = ds[var_name].attrs.get('long_name', var_name)
    ax.set_xlabel('Year', fontsize=16)
    ax.set_ylabel(f"{var_name} ({units})", fontsize=16)
    ax.set_title(f"{long_name} at ({actual_lon:.2f}°, {actual_lat:.2f}°)", fontsize=20)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    return ax


def plot_annual_mean(ds, year, var_name=None, ax=None, **kwargs):
    """
    Plot annual mean for a specific year.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset
    var_name : str
        Variable name
    year : int
        Year to average
    ax : matplotlib axes, optional
        Axes to plot on
    **kwargs :
        Additional arguments passed to plot_map
    """
    if var_name is None:
        var_name = _infer_var_name(ds)
    # Select data for the year
    year_data = ds.sel(time=slice(year, year + 1))
    year_mean = year_data.mean(dim='time')

    # Create temporary dataset for plotting
    ds_plot = xr.Dataset()
    ds_plot[var_name] = year_mean[var_name].expand_dims('time')
    ds_plot['time'] = [year + 0.5]  # Mid-year for title display
    ds_plot['lat'] = ds['lat']
    ds_plot['lon'] = ds['lon']
    ds_plot = ds_plot.assign_coords(rlat=ds['rlat'], rlon=ds['rlon'])

    # Copy attributes
    ds_plot[var_name].attrs = ds[var_name].attrs

    return plot_map(ds_plot, var_name, time_idx=0, ax=ax, **kwargs)


def plot_temporal_mean(ds, year_start, year_end, var_name=None, ax=None, **kwargs):
    """
    Plot temporal mean of a variable over a given year range.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing the variable
    var_name : str
        Variable name to plot
    year_start : int or float
        Start of averaging period (inclusive)
    year_end : int or float
        End of averaging period (inclusive)
    ax : matplotlib axes, optional
        Axes to plot on
    **kwargs :
        Additional arguments passed to plot_map (cmap, vmin, vmax, etc.)

    Returns
    -------
    (matplotlib axes, xr.DataArray)
        Axes and the computed period-mean DataArray.

    Examples
    --------
    >>> ds = load_gridded_data('FirnAir')
    >>> ax, data = plot_temporal_mean(ds, 1940, 1970)
    >>> ax, data = plot_temporal_mean(ds, 1940, 1970, 'FirnAir')
    """
    if var_name is None:
        var_name = _infer_var_name(ds)
    period_data = ds.sel(time=slice(year_start, year_end))
    if len(period_data.time) == 0:
        raise ValueError(f"No data found between {year_start} and {year_end}")

    period_mean = period_data.mean(dim='time')

    ds_plot = xr.Dataset()
    ds_plot[var_name] = period_mean[var_name].expand_dims('time')
    ds_plot['time'] = [(year_start + year_end) / 2.0]
    ds_plot['lat'] = ds['lat']
    ds_plot['lon'] = ds['lon']
    ds_plot = ds_plot.assign_coords(rlat=ds['rlat'], rlon=ds['rlon'])
    ds_plot[var_name].attrs = ds[var_name].attrs

    ax = plot_map(ds_plot, var_name, time_idx=0, ax=ax, **kwargs)

    long_name = ds[var_name].attrs.get('long_name', var_name)
    ax.set_title(f"{long_name} ({year_start}–{year_end} mean)", fontsize=20)

    return ax, period_mean[var_name]


def plot_temporal_mean_difference(ds, period1, period2, var_name=None, ax=None,
                                  cmap=DIVERGING_CMAP, vcenter=0.0, cax=None, **kwargs):
    """
    Plot the difference between temporal means of two periods (period2 - period1).

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing the variable
    var_name : str
        Variable name to plot
    period1 : tuple of (int or float, int or float)
        (year_start, year_end) for the first (reference) period
    period2 : tuple of (int or float, int or float)
        (year_start, year_end) for the second period
    ax : matplotlib axes, optional
        Axes to plot on
    cmap : str
        Colormap. Default 'RdBu_r' (diverging).
    vcenter : float
        Value to center the colormap on. Default 0.0.
    **kwargs :
        Additional arguments passed to plot_map (vmin, vmax, add_colorbar, etc.)

    Returns
    -------
    (matplotlib axes, xr.DataArray)
        Axes and the computed difference DataArray (period2_mean − period1_mean).

    Examples
    --------
    >>> ds = load_gridded_data('FirnAir')
    >>> ax, diff = plot_temporal_mean_difference(ds, (1940, 1970), (1995, 2023))
    >>> ax, diff = plot_temporal_mean_difference(ds, (1940, 1970), (1995, 2023),
    ...                                          vmin=-5, vmax=5)
    """
    if var_name is None:
        var_name = _infer_var_name(ds)
    from matplotlib.colors import TwoSlopeNorm

    def _period_mean(year_start, year_end):
        data = ds.sel(time=slice(year_start, year_end))
        if len(data.time) == 0:
            raise ValueError(f"No data found between {year_start} and {year_end}")
        return data[var_name].mean(dim='time')

    mean1 = _period_mean(*period1)
    mean2 = _period_mean(*period2)
    diff = mean2 - mean1

    # Build a minimal single-timestep dataset so plot_map handles axes/colorbar
    mid_time = (period1[0] + period1[1] + period2[0] + period2[1]) / 4.0
    ds_plot = xr.Dataset()
    ds_plot[var_name] = diff.expand_dims('time')
    ds_plot['time'] = [mid_time]
    ds_plot['lat'] = ds['lat']
    ds_plot['lon'] = ds['lon']
    ds_plot = ds_plot.assign_coords(rlat=ds['rlat'], rlon=ds['rlon'])
    ds_plot[var_name].attrs = ds[var_name].attrs

    # Symmetric color limits around vcenter if not overridden
    vmin = kwargs.pop('vmin', None)
    vmax = kwargs.pop('vmax', None)
    if vmin is None and vmax is None:
        abs_max = float(np.nanmax(np.abs(diff.values)))
        vmin, vmax = vcenter - abs_max, vcenter + abs_max

    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)

    # plot_map doesn't accept norm directly, so plot manually
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 10))

    units = ds[var_name].attrs.get('units', '')
    add_colorbar = kwargs.pop('add_colorbar', True)

    if add_colorbar:
        # If no cax provided, create one attached to ax so the colorbar steals
        # space from the map axes only (not from other axes in the figure).
        if cax is None:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.08)
        cbar_kwargs = {'cax': cax}
        if units:
            cbar_kwargs['label'] = units
    else:
        cbar_kwargs = None

    diff.plot(ax=ax, cmap=cmap, norm=norm,
              add_colorbar=add_colorbar,
              cbar_kwargs=cbar_kwargs)

    long_name = ds[var_name].attrs.get('long_name', var_name)
    ax.set_title(
        f"{long_name}\n({period2[0]}–{period2[1]} mean) − ({period1[0]}–{period1[1]} mean)",
        fontsize=20
    )
    ax.set_aspect('equal')
    ax.set_xlabel('rlon (grid index)', fontsize=16)
    ax.set_ylabel('rlat (grid index)', fontsize=16)
    if add_colorbar and cax is not None:
        cax.yaxis.label.set_size(16)
        cax.tick_params(labelsize=12)

    return ax, diff


# =============================================================================
# Multi-Variable Comparison
# =============================================================================

def plot_multi_var_map(var_names, time_value, cmaps=None, vlims=None, figsize=None):
    """
    Plot multiple variables side by side at the same time.

    Parameters
    ----------
    var_names : list of str
        List of variable names to plot (will load each automatically)
    time_value : float
        Time value (fractional year) to plot
    cmaps : list of str, optional
        List of colormaps for each variable. If None, uses 'viridis' for all.
    vlims : list of tuple, optional
        List of (vmin, vmax) tuples for each variable. If None, uses auto-scaling.
    figsize : tuple, optional
        Figure size. If None, auto-calculated based on number of variables.

    Returns
    -------
    fig, axes
        Figure and axes objects
    """
    n_vars = len(var_names)

    # Set defaults
    if cmaps is None:
        cmaps = ['viridis'] * n_vars
    if vlims is None:
        vlims = [(None, None)] * n_vars
    if figsize is None:
        figsize = (5 * n_vars, 8)

    # Create figure
    fig, axes = plt.subplots(1, n_vars, figsize=figsize)
    if n_vars == 1:
        axes = [axes]

    # Plot each variable
    for ax, var_name, cmap, (vmin, vmax) in zip(axes, var_names, cmaps, vlims):
        ds = load_gridded_data(var_name)
        plot_map(ds, var_name, time_value=time_value, ax=ax,
                 cmap=cmap, vmin=vmin, vmax=vmax)
        ds.close()

    plt.tight_layout()
    return fig, axes


# =============================================================================
# Ice Sheet Facies Classification
# =============================================================================

def classify_zones(year, melt_threshold=5.0):
    """
    Classify ice sheet zones for a given year based on annual totals.

    Parameters
    ----------
    year : int
        Year to classify.
    melt_threshold : float
        Cumulative annual surface melt (mm w.e.) above which a pixel is
        considered to have significant melt. Default 5.0.

    Returns
    -------
    xr.DataArray
        Zone classification:
        0 = No data (NaN)
        1 = Ablation zone  (annual melt > annual accumulation)
        2 = Percolation zone (melt > threshold AND NOT ablation)
        3 = Dry snow zone  (melt <= threshold AND NOT ablation)
    """
    ds_melt = load_gridded_data('surfmelt')
    ds_vacc = load_gridded_data('vacc')

    yr = slice(year, year + 1)
    # Annual cumulative melt [mm w.e.]
    melt_year = ds_melt['surfmelt'].sel(time=yr).sum(dim='time')
    # Annual mean accumulation rate [m/yr] → convert to mm w.e./yr
    accum_year = ds_vacc['vacc'].sel(time=yr).mean(dim='time') * 1000.0

    zones = xr.zeros_like(melt_year)

    # Ablation: annual melt exceeds annual accumulation
    ablation_mask = melt_year > accum_year
    zones = zones.where(~ablation_mask, 1)

    # Percolation: significant melt but not losing mass overall
    percolation_mask = (melt_year > melt_threshold) & ~ablation_mask
    zones = zones.where(~percolation_mask, 2)

    # Dry snow: little or no melt and not ablating
    dry_snow_mask = (melt_year <= melt_threshold) & ~ablation_mask
    zones = zones.where(~dry_snow_mask, 3)

    # Preserve NaN mask: off-ice pixels are NaN in vacc even if melt=0
    valid = ~np.isnan(melt_year) & ~np.isnan(accum_year)
    zones = zones.where(valid, np.nan)

    zones.attrs['long_name'] = f'Ice sheet facies classification ({year})'
    zones.attrs['flag_values'] = [1, 2, 3]
    zones.attrs['flag_meanings'] = 'ablation_zone percolation_zone dry_snow_zone'

    zones = zones.assign_coords(
        lat=ds_melt['lat'],
        lon=ds_melt['lon'],
    )

    ds_melt.close()
    ds_vacc.close()

    return zones


def plot_zones(year, ax=None, melt_threshold=5.0, add_colorbar=True):
    """
    Plot ablation and percolation zones for a given year.

    Parameters
    ----------
    year : int
        Year to plot
    ax : matplotlib axes, optional
        Axes to plot on
    melt_threshold : float
        Cumulative annual surface melt (mm w.e.) above which a pixel is
        classified as percolation zone. Default 5.0.
    add_colorbar : bool
        Whether to add a colorbar. Set to False when using
        plot_zones_timeseries, which adds one shared colorbar. Default True.

    Returns
    -------
    matplotlib axes
    """
    # Get zone classification
    zones = classify_zones(year, melt_threshold)

    # Create figure if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 10))

    cmap = ListedColormap(FACIES_COLORS)
    norm = BoundaryNorm([0.5, 1.5, 2.5, 3.5], cmap.N)

    # Plot
    im = zones.plot(
        ax=ax,
        cmap=cmap,
        norm=norm,
        add_colorbar=False
    )

    if add_colorbar:
        cbar = plt.colorbar(im, ax=ax, shrink=0.6, aspect=20, pad=0.02)
        cbar.set_ticks([1, 2, 3])
        cbar.set_ticklabels(['Ablation', 'Percolation', 'Dry Snow'])
        cbar.ax.tick_params(labelsize=12)

    ax.set_aspect('equal')
    ax.set_xlabel('rlon (grid index)', fontsize=16)
    ax.set_ylabel('rlat (grid index)', fontsize=16)
    ax.set_title(f'Ice Sheet Facies ({year})', fontsize=20)

    return ax


def plot_zones_timeseries(years, melt_threshold=5.0, ncols=4, figsize=None):
    """
    Plot zone maps for multiple years.

    Parameters
    ----------
    years : list of int
        Years to plot
    melt_threshold : float
        Cumulative annual surface melt (mm w.e.) above which a pixel is
        classified as percolation zone. Default 5.0.
    ncols : int
        Number of columns in subplot grid. Default 4.
    figsize : tuple, optional
        Figure size. If None, auto-calculated.

    Returns
    -------
    fig, axes
        Figure and axes objects
    """
    n_years = len(years)
    nrows = int(np.ceil(n_years / ncols))

    if figsize is None:
        figsize = (4 * ncols, 5 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = np.atleast_2d(axes).flatten()

    for i, year in enumerate(years):
        plot_zones(year, ax=axes[i],
                   melt_threshold=melt_threshold,
                   add_colorbar=False)
        axes[i].set_title(str(year), fontsize=20)

    # Hide unused axes
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(f'Ice Sheet Facies Timeseries (\u2265{melt_threshold} mm w.e.)',
                 fontsize=24)

    # Single shared colorbar — placed in reserved space below all subplots
    zones_cmap = ListedColormap(FACIES_COLORS)
    zones_norm = BoundaryNorm([0.5, 1.5, 2.5, 3.5], zones_cmap.N)
    sm = plt.cm.ScalarMappable(cmap=zones_cmap, norm=zones_norm)
    sm.set_array([])

    # tight_layout rect leaves space at top (suptitle) and bottom (colorbar)
    plt.tight_layout(rect=[0, 0.10, 1, 0.95])
    # Colorbar axes: centred horizontally, [left, bottom, width, height]
    cbar_width = 0.30
    cax = fig.add_axes([(1 - cbar_width) / 2, 0.03, cbar_width, 0.025])
    cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')
    cbar.set_ticks([1, 2, 3])
    cbar.set_ticklabels(['Ablation', 'Percolation', 'Dry Snow'])
    cbar.ax.tick_params(labelsize=12)

    return fig, axes


# =============================================================================
# Statistics and Export
# =============================================================================

def print_stats(ds, var_name=None):
    """Print basic statistics for a variable."""
    if var_name is None:
        var_name = _infer_var_name(ds)
    data = ds[var_name]

    print(f"\nStatistics for {var_name}:")
    print(f"  Shape: {data.shape}")
    print(f"  Min: {float(data.min()):.4f}")
    print(f"  Max: {float(data.max()):.4f}")
    print(f"  Mean: {float(data.mean()):.4f}")
    print(f"  Std: {float(data.std()):.4f}")
    print(f"  NaN count: {int(np.isnan(data.values).sum())}")
    print(f"  Valid count: {int((~np.isnan(data.values)).sum())}")


def export_annual_means(ds, var_name=None, output_file=None):
    """
    Export annual means to a separate file.

    Parameters
    ----------
    ds : xr.Dataset
        Input dataset
    var_name : str
        Variable name
    output_file : str
        Output file path
    """
    if var_name is None:
        var_name = _infer_var_name(ds)
    # Calculate annual means
    ds_annual = ds.groupby(np.floor(ds.time)).mean(dim='time')
    ds_annual = ds_annual.rename({'floor': 'year'})

    # Save
    ds_annual.to_netcdf(output_file)
    print(f"Saved annual means to: {output_file}")


# =============================================================================
# Animation/GIF Creation
# =============================================================================

def create_animation_gif(var_name, output_path, time_slice=None, fps=10,
                         cmap='viridis', vmin=None, vmax=None, dpi=100,
                         figsize=(8, 10), use_latlon=False, title_fmt=None):
    """
    Create an animated GIF showing a variable at each timestep.

    Parameters
    ----------
    var_name : str
        Variable name to animate (will load automatically)
    output_path : str or Path
        Output path for the GIF file
    time_slice : slice, optional
        Time slice to animate (e.g., slice(2010, 2020)). If None, uses all times.
    fps : int
        Frames per second for the animation. Default 10.
    cmap : str
        Colormap name. Default 'viridis'.
    vmin, vmax : float, optional
        Color scale limits. If None, computed from the data range.
    dpi : int
        Resolution of each frame. Default 100.
    figsize : tuple
        Figure size. Default (8, 10).
    use_latlon : bool
        If True, use lat/lon coordinates for axes. Default False.
    title_fmt : str, optional
        Custom title format string. Use {var_name}, {long_name}, {time}, {units}.
        Default: "{long_name} ({time:.2f})"

    Returns
    -------
    Path
        Path to the created GIF file

    Examples
    --------
    >>> create_animation_gif('Runoff', 'runoff_animation.gif',
    ...                      time_slice=slice(2010, 2020), cmap='Blues')
    >>> create_animation_gif('surfmelt', 'melt_2019.gif',
    ...                      time_slice=slice(2019, 2020), fps=5)
    """
    import matplotlib.animation as animation
    from PIL import Image
    import io

    # Load data
    ds = load_gridded_data(var_name)

    # Select time range
    if time_slice is not None:
        ds = ds.sel(time=time_slice)

    data = ds[var_name]
    times = ds['time'].values

    # Get metadata
    long_name = data.attrs.get('long_name', var_name)
    units = data.attrs.get('units', '')

    # Compute color limits if not provided
    if vmin is None:
        vmin = float(data.min())
    if vmax is None:
        vmax = float(data.max())

    # Default title format
    if title_fmt is None:
        title_fmt = "{long_name} ({time:.2f})"

    # Create frames
    frames = []
    n_times = len(times)

    print(f"Creating animation with {n_times} frames...")

    for i, t in enumerate(times):
        if (i + 1) % 50 == 0 or i == 0:
            print(f"  Processing frame {i + 1}/{n_times}")

        fig, ax = plt.subplots(figsize=figsize)

        # Get data for this timestep
        frame_data = data.isel(time=i)

        # Colorbar settings
        cbar_kwargs = {'shrink': 0.6, 'aspect': 20, 'pad': 0.02}
        if units:
            cbar_kwargs['label'] = units

        # Plot
        if use_latlon:
            frame_plot = frame_data.assign_coords(
                x=(['rlat', 'rlon'], ds['lon'].values),
                y=(['rlat', 'rlon'], ds['lat'].values)
            )
            p = frame_plot.plot(
                ax=ax, x='x', y='y',
                cmap=cmap, vmin=vmin, vmax=vmax,
                add_colorbar=True,
                cbar_kwargs=cbar_kwargs
            )
            ax.set_xlabel('Longitude', fontsize=16)
            ax.set_ylabel('Latitude', fontsize=16)
        else:
            p = frame_data.plot(
                ax=ax,
                cmap=cmap, vmin=vmin, vmax=vmax,
                add_colorbar=True,
                cbar_kwargs=cbar_kwargs
            )
            ax.set_xlabel('rlon (grid index)', fontsize=16)
            ax.set_ylabel('rlat (grid index)', fontsize=16)

        if hasattr(p, 'colorbar') and p.colorbar is not None:
            if units:
                p.colorbar.set_label(units, fontsize=16)
            p.colorbar.ax.tick_params(labelsize=12)

        ax.set_aspect('equal')

        # Set title
        title = title_fmt.format(
            var_name=var_name,
            long_name=long_name,
            time=float(t),
            units=units
        )
        ax.set_title(title, fontsize=20)

        plt.tight_layout()

        # Save frame to buffer
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=dpi)
        buf.seek(0)
        frames.append(Image.open(buf).copy())
        buf.close()
        plt.close(fig)

    # Save as GIF
    output_path = Path(output_path)
    duration = int(1000 / fps)  # milliseconds per frame

    print(f"Saving GIF to {output_path}...")
    frames[0].save(
        output_path,
        save_all=True,
        append_images=frames[1:],
        duration=duration,
        loop=0
    )

    ds.close()
    print(f"Animation saved: {output_path}")

    return output_path


# =============================================================================
# Masked Difference Plotting
# =============================================================================

def plot_masked_difference(var1_name, var2_name, time_value=None, time_slice=None,
                           threshold1=0.0, threshold2=0.0,
                           comparison='var1_only', ax=None,
                           cmap=DIVERGING_CMAP, figsize=(10, 12)):
    """
    Plot the difference between two variables, filtered by a mask condition.

    This is useful for identifying where one process occurs but another doesn't,
    e.g., where surface melt occurs but runoff doesn't (indicating refreezing).

    Parameters
    ----------
    var1_name : str
        First variable name (will load automatically)
    var2_name : str
        Second variable name (will load automatically)
    time_value : float, optional
        Specific time (fractional year) to plot. If None with no time_slice,
        uses annual mean over the entire dataset.
    time_slice : slice, optional
        Time range to aggregate (e.g., slice(2010, 2020)). Takes precedence
        over time_value if both are None for computing means.
    threshold1 : float
        Threshold for var1. Default 0.0.
    threshold2 : float
        Threshold for var2. Default 0.0.
    comparison : str
        Type of comparison:
        - 'var1_only': Show where var1 > threshold1 AND var2 <= threshold2
        - 'var2_only': Show where var2 > threshold2 AND var1 <= threshold1
        - 'both': Show where both var1 > threshold1 AND var2 > threshold2
        - 'neither': Show where both var1 <= threshold1 AND var2 <= threshold2
        - 'difference': Show var1 - var2 everywhere (no masking)
        - 'ratio': Show var1 / var2 where var2 > threshold2
    ax : matplotlib axes, optional
        Axes to plot on
    cmap : str
        Colormap. Default 'RdBu_r'.
    figsize : tuple
        Figure size. Default (10, 12).

    Returns
    -------
    fig, ax
        Figure, axes

    Examples
    --------
    >>> # Where does melt occur but not runoff? (indicates refreezing)
    >>> fig, ax = plot_masked_difference('surfmelt', 'Runoff',
    ...                                        time_slice=slice(2019, 2020),
    ...                                        comparison='var1_only')

    >>> # Where does runoff occur without melt? (rain-driven runoff)
    >>> fig, ax = plot_masked_difference('Runoff', 'surfmelt',
    ...                                        time_value=2012.5,
    ...                                        comparison='var1_only')

    >>> # Show difference between melt and runoff (positive = retained)
    >>> fig, ax = plot_masked_difference('surfmelt', 'Runoff',
    ...                                        time_slice=slice(2010, 2020),
    ...                                        comparison='difference')
    """
    # Load data
    ds1 = load_gridded_data(var1_name)
    ds2 = load_gridded_data(var2_name)

    # Select/aggregate time
    if time_value is not None:
        data1 = ds1[var1_name].sel(time=time_value, method='nearest')
        data2 = ds2[var2_name].sel(time=time_value, method='nearest')
        time_str = f"{float(data1.time):.2f}"
    elif time_slice is not None:
        data1 = ds1[var1_name].sel(time=time_slice).sum(dim='time')
        data2 = ds2[var2_name].sel(time=time_slice).sum(dim='time')
        time_str = f"{time_slice.start}-{time_slice.stop}"
    else:
        data1 = ds1[var1_name].sum(dim='time')
        data2 = ds2[var2_name].sum(dim='time')
        time_str = "total"

    # Get metadata
    long_name1 = ds1[var1_name].attrs.get('long_name', var1_name)
    long_name2 = ds2[var2_name].attrs.get('long_name', var2_name)
    units1 = ds1[var1_name].attrs.get('units', '')
    units2 = ds2[var2_name].attrs.get('units', '')

    # Compute result based on comparison type
    if comparison == 'var1_only':
        mask = (data1 > threshold1) & (data2 <= threshold2)
        result = data1.where(mask)
        title = f"{long_name1} where {var2_name} absent ({time_str})"
        units = units1
        cmap = 'Reds'
    elif comparison == 'var2_only':
        mask = (data2 > threshold2) & (data1 <= threshold1)
        result = data2.where(mask)
        title = f"{long_name2} where {var1_name} absent ({time_str})"
        units = units2
        cmap = 'Blues'
    elif comparison == 'both':
        mask = (data1 > threshold1) & (data2 > threshold2)
        result = (data1 + data2).where(mask) / 2  # Average of both
        title = f"Both {var1_name} and {var2_name} present ({time_str})"
        units = units1
        cmap = 'Purples'
    elif comparison == 'neither':
        mask = (data1 <= threshold1) & (data2 <= threshold2)
        result = xr.ones_like(data1).where(mask)
        title = f"Neither {var1_name} nor {var2_name} ({time_str})"
        units = ''
        cmap = 'Greys'
    elif comparison == 'difference':
        result = data1 - data2
        title = f"{long_name1} - {long_name2} ({time_str})"
        units = units1 if units1 == units2 else f"{units1} - {units2}"
    elif comparison == 'ratio':
        result = data1 / data2.where(data2 > threshold2)
        title = f"{var1_name} / {var2_name} ({time_str})"
        units = ''
        cmap = DIVERGING_CMAP
    else:
        raise ValueError(f"Unknown comparison type: {comparison}")

    # Create figure if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    # Colorbar settings
    cbar_kwargs = {'shrink': 0.6, 'aspect': 20, 'pad': 0.02}
    if units:
        cbar_kwargs['label'] = units

    # Plot
    p = result.plot(
        ax=ax,
        cmap=cmap,
        add_colorbar=True,
        cbar_kwargs=cbar_kwargs
    )

    if hasattr(p, 'colorbar') and p.colorbar is not None:
        if units:
            p.colorbar.set_label(units, fontsize=16)
        p.colorbar.ax.tick_params(labelsize=12)

    ax.set_aspect('equal')
    ax.set_xlabel('rlon (grid index)', fontsize=16)
    ax.set_ylabel('rlat (grid index)', fontsize=16)
    ax.set_title(title, fontsize=20)

    plt.tight_layout()

    ds1.close()
    ds2.close()

    return fig, ax#, result


def plot_process_comparison(time_slice=None, time_value=None, figsize=(16, 10)):
    """
    Create a 2x2 comparison plot showing where different melt/runoff processes occur.

    This provides a quick overview of:
    - Where melt occurs but runoff doesn't (refreezing)
    - Where runoff occurs but melt doesn't (rain-driven)
    - Where both occur (melt-driven runoff)
    - The difference between melt and runoff

    Parameters
    ----------
    time_slice : slice, optional
        Time range to aggregate. Default uses full dataset.
    time_value : float, optional
        Specific time to plot. Overrides time_slice.
    figsize : tuple
        Figure size. Default (16, 10).

    Returns
    -------
    fig, axes
        Figure and axes objects
    """
    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # Melt but no runoff (refreezing)
    plot_masked_difference('surfmelt', 'Runoff',
                           time_value=time_value, time_slice=time_slice,
                           comparison='var1_only', ax=axes[0, 0])
    axes[0, 0].set_title('Melt but no Runoff\n(Refreezing)', fontsize=20)

    # Runoff but no melt (rain-driven)
    plot_masked_difference('Runoff', 'surfmelt',
                           time_value=time_value, time_slice=time_slice,
                           comparison='var1_only', ax=axes[0, 1])
    axes[0, 1].set_title('Runoff but no Melt\n(Rain-driven)', fontsize=20)

    # Both present
    plot_masked_difference('surfmelt', 'Runoff',
                           time_value=time_value, time_slice=time_slice,
                           comparison='both', ax=axes[1, 0])
    axes[1, 0].set_title('Both Melt and Runoff\n(Melt-driven runoff)', fontsize=20)

    # Difference
    plot_masked_difference('surfmelt', 'Runoff',
                           time_value=time_value, time_slice=time_slice,
                           comparison='difference', ax=axes[1, 1])
    axes[1, 1].set_title('Melt - Runoff\n(Positive = retained)', fontsize=20)

    plt.tight_layout()
    return fig, axes


# =============================================================================
# Mask Overlays
# =============================================================================

# Default basin names for Mouginot_basins (IDs 1-7).
# Check the mask file if you are unsure which ID maps to which basin.
MOUGINOT_BASIN_NAMES = {
    1: 'NO', 2: 'NE', 3: 'CE', 4: 'SE', 5: 'SW', 6: 'CW', 7: 'NW',
}


def add_elevation_contours(ds, ax, interval=500, ds_mask=None,
                           color='black', linewidth=0.5, alpha=0.6,
                           label_contours=True, fontsize=12):
    """
    Add elevation contour lines to an existing map axes.

    Parameters
    ----------
    ax : matplotlib axes
    ds : xr.Dataset
        The same dataset used to make the map (from load_gridded_data).
        Its rlat/rlon coordinate values (rotated lat/lon degrees) are used
        to place the contours in the correct coordinate system.
    interval : int
        Contour interval in metres. Default 500.
    ds_mask : xr.Dataset, optional
        Pre-loaded mask dataset. If None, loads from config.MASK_FILE.
    color : str
        Contour colour. Default 'black'.
    linewidth : float
        Contour line width. Default 0.5.
    alpha : float
        Opacity. Default 0.6.
    label_contours : bool
        Whether to label contour lines with elevation values. Default True.
    fontsize : int
        Font size for contour labels. Default 7.

    Returns
    -------
    QuadContourSet

    Examples
    --------
    >>> ds = load_gridded_data('FirnAir')
    >>> ax = plot_map(ds, 'FirnAir', time_idx=0)
    >>> add_elevation_contours(ds, ax)
    """
    rlat = ds['rlat'].values
    rlon = ds['rlon'].values

    opened = ds_mask is None
    if opened:
        ds_mask = xr.open_dataset(config.MASK_FILE).squeeze('time', drop=True)

    try:
        topo = ds_mask['Topography'].values.copy()
        topo[ds_mask['LSM_GR'].values == 0] = np.nan
        levels = np.arange(interval, float(np.nanmax(topo)) + interval, interval)
        cs = ax.contour(rlon, rlat, topo,
                        levels=levels, colors=color,
                        linewidths=linewidth, alpha=alpha)
        if label_contours:
            ax.clabel(cs, fmt='%d m', fontsize=fontsize, inline=True)
    finally:
        if opened:
            ds_mask.close()

    return cs


def add_land_outline(ds, ax, ds_mask=None,
                     color='black', linewidth=0.8, alpha=0.7):
    """
    Add the land/ice-sheet outline (LSM_GR coastline) to an existing map axes.

    Parameters
    ----------
    ds : xr.Dataset
        The same dataset used to make the map (from load_gridded_data).
    ax : matplotlib axes
    ds_mask : xr.Dataset, optional
        Pre-loaded mask dataset. If None, loads from config.MASK_FILE.
    color : str
        Outline colour. Default 'black'.
    linewidth : float
        Line width. Default 0.8.
    alpha : float
        Opacity. Default 0.7.

    Returns
    -------
    QuadContourSet
    """
    rlat = ds['rlat'].values
    rlon = ds['rlon'].values

    opened = ds_mask is None
    if opened:
        ds_mask = xr.open_dataset(config.MASK_FILE).squeeze('time', drop=True)

    try:
        lsm = ds_mask['LSM_GR'].values
        cs = ax.contour(rlon, rlat, lsm, levels=[0.5],
                        colors=color, linewidths=linewidth, alpha=alpha)
    finally:
        if opened:
            ds_mask.close()

    return cs


def add_basin_outlines(ds, ax, basin_var='Mouginot_basins', ds_mask=None,
                       color='black', linewidth=1.2, alpha=0.9,
                       add_labels=True, label_names=None, label_fontsize=12,
                       add_coast=True, coast_color='black',
                       coast_linewidth=0.8, coast_alpha=0.7):
    """
    Add basin boundary outlines to an existing map axes.

    Also draws the land/sea outline and labels each basin by default.

    Parameters
    ----------
    ds : xr.Dataset
        The same dataset used to make the map (from load_gridded_data).
    ax : matplotlib axes
    basin_var : str
        Variable in the mask file. Options: 'Mouginot_basins', 'GRACE_basins',
        'Basins'. Default 'Mouginot_basins'.
    ds_mask : xr.Dataset, optional
        Pre-loaded mask dataset. If None, loads from config.MASK_FILE.
    color : str
        Basin outline colour. Default 'black'.
    linewidth : float
        Basin outline line width. Default 1.2.
    alpha : float
        Basin outline opacity. Default 0.9.
    add_labels : bool
        If True, place basin name at the centroid of each basin. Default True.
    label_names : dict, optional
        Mapping {basin_id (int): name (str)}. If None and basin_var is
        'Mouginot_basins', uses MOUGINOT_BASIN_NAMES.
    label_fontsize : int
        Font size for basin labels. Default 9.
    add_coast : bool
        If True, also draw the land/sea outline from LSM_GR. Default True.
    coast_color : str
        Coastline colour. Default 'black'.
    coast_linewidth : float
        Coastline line width. Default 0.8.
    coast_alpha : float
        Coastline opacity. Default 0.7.

    Returns
    -------
    QuadContourSet
        The basin boundary contour set.

    Examples
    --------
    >>> ds = load_gridded_data('FirnAir')
    >>> ax = plot_map(ds, 'FirnAir', time_idx=0)
    >>> add_basin_outlines(ds, ax)
    """
    rlat = ds['rlat'].values
    rlon = ds['rlon'].values

    opened = ds_mask is None
    if opened:
        ds_mask = xr.open_dataset(config.MASK_FILE).squeeze('time', drop=True)

    try:
        basins = ds_mask[basin_var].values.astype(float)
        basin_ids = np.unique(basins[~np.isnan(basins)])
        levels = basin_ids[:-1] + 0.5
        cs = ax.contour(rlon, rlat, basins,
                        levels=levels, colors=color,
                        linewidths=linewidth, alpha=alpha)

        if add_coast:
            lsm = ds_mask['LSM'].values
            ax.contour(rlon, rlat, lsm, levels=[0.5],
                       colors=coast_color, linewidths=coast_linewidth,
                       alpha=coast_alpha)

        if add_labels:
            names = label_names
            if names is None:
                names = MOUGINOT_BASIN_NAMES if basin_var == 'Mouginot_basins' \
                    else {int(i): str(int(i)) for i in basin_ids if i > 0}

            # Build 2D coordinate grids to compute per-basin centroids
            rlat_2d, rlon_2d = np.meshgrid(rlat, rlon, indexing='ij')
            for basin_id in basin_ids:
                if basin_id == 0:
                    continue
                mask = basins == basin_id
                if not np.any(mask):
                    continue
                centroid_rlon = rlon_2d[mask].mean()
                centroid_rlat = rlat_2d[mask].mean()
                name = names.get(int(basin_id), str(int(basin_id)))
                if name=="CW":
                    centroid_rlon = centroid_rlon+1
                elif name=="SE":
                    centroid_rlat = centroid_rlat+1
                ax.text(centroid_rlon, centroid_rlat, name,
                        ha='center', va='center',
                        fontsize=label_fontsize, fontweight='bold',
                        color=color,
                        bbox=dict(boxstyle='round,pad=0.15', fc='white',
                                  ec='none', alpha=0.6))
                

    finally:
        if opened:
            ds_mask.close()

    return cs


# =============================================================================
# Spatial Mean Time Series
# =============================================================================

def _load_mask_aligned(ds):
    """Load the mask file with rlat/rlon coordinates aligned to ds (rotated degrees)."""
    ds_mask = xr.open_dataset(config.MASK_FILE).squeeze('time', drop=True)
    return ds_mask.assign_coords(rlat=ds['rlat'].values, rlon=ds['rlon'].values)


def _spatial_mean_timeseries(data_var, mask_2d, agg='mean', annual=True,
                              year_start=None, drop_last_year=True):
    """
    Compute spatial mean time series.

    All grid cells have equal area (5.5 km²), so this is a simple mean
    over the valid pixels — no area weighting needed.

    Parameters
    ----------
    data_var : xr.DataArray
        Variable with dimensions (time, rlat, rlon).
    mask_2d : xr.DataArray or None
        Boolean mask (rlat, rlon). True = include. If None, all non-NaN
        pixels are included (off-ice pixels are already NaN in the data).
    agg : str
        'sum' or 'mean' for annual aggregation.
    annual : bool
        If True, aggregate to annual values.
    year_start : float or None
        If given, crop data to times >= year_start before aggregating.
    drop_last_year : bool
        If True (default), drop the last annual bin, which is typically
        partial (the run start/end rarely aligns with Jan 1).

    Returns
    -------
    xr.DataArray
        Time series with dimension 'year' (if annual) or 'time'.
    """
    data_sliced = data_var.sel(time=slice(year_start, None)) if year_start is not None else data_var

    if mask_2d is not None:
        mean_ts = data_sliced.where(mask_2d).mean(dim=['rlat', 'rlon'])
    else:
        mean_ts = data_sliced.mean(dim=['rlat', 'rlon'])

    if annual:
        years = np.floor(mean_ts['time'].values)
        mean_ts = mean_ts.assign_coords(year=('time', years))
        if agg == 'sum':
            mean_ts = mean_ts.groupby('year').sum()
        else:
            mean_ts = mean_ts.groupby('year').mean()
        if drop_last_year:
            mean_ts = mean_ts.isel(year=slice(None, -1))

    return mean_ts


def _plot_variability_and_rolling(ax, x_sub, y_sub, x_ann, y_ann,
                                   rolling_window, label, kwargs,
                                   rolling_center=True):
    """
    Plot sub-annual values as a thin transparent line and a rolling mean as a
    bold line.  Both lines share the same colour (auto-cycled unless 'color'
    is supplied in kwargs).

    rolling_center=True  → centered window (peak timing preserved)
    rolling_center=False → trailing window (first value = first annual value)
    """
    import pandas as pd
    y_roll = (pd.Series(y_ann)
              .rolling(rolling_window, center=rolling_center, min_periods=1)
              .mean().values)

    # Strip formatting keys the caller controls so we can set them explicitly.
    base_kw = {k: v for k, v in kwargs.items()
               if k not in ('linewidth', 'lw', 'alpha', 'label')}
    line, = ax.plot(x_sub, y_sub, linewidth=0.5, alpha=0.2, **base_kw)

    roll_kw = {k: v for k, v in base_kw.items() if k != 'color'}
    ax.plot(x_ann, y_roll, linewidth=2, label=label,
            color=line.get_color(), **roll_kw)


def plot_timeseries_gis(ds, var_name=None, ax=None, label='GIS',
                        year_start=1940, rolling_window=10, **kwargs):
    """
    Plot spatial mean time series over the whole Greenland Ice Sheet.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing the variable (from load_gridded_data).
    var_name : str
        Variable name to plot.
    ax : matplotlib axes, optional
    label : str
        Line label. Default 'GIS'.
    year_start : float
        Start year for the time series. Default 1940 (first full year).
    rolling_window : int or None
        If set, plot sub-annual values as a thin background line and overlay
        a bold N-year rolling mean. Set to None to plot annual values only.
        Default 10.
    **kwargs
        Passed to ax.plot().

    Returns
    -------
    (matplotlib axes, dict)
        Axes and ``{'sub': ts_sub, 'ann': ts_ann}`` DataArrays (sub-annual
        and annual series). Use these to re-plot without recomputing.
    """
    if var_name is None:
        var_name = _infer_var_name(ds)
    agg = config.get_aggregation_method(var_name)

    if ax is None:
        _, ax = plt.subplots(figsize=(14, 5))
    ax.set_prop_cycle(color=OKABE_ITO_COLORS)

    ts_sub = _spatial_mean_timeseries(ds[var_name], None, agg=agg,
                                      annual=False, year_start=year_start)
    ts_ann = _spatial_mean_timeseries(ds[var_name], None, agg=agg,
                                      annual=True, year_start=year_start)

    if rolling_window is not None:
        _plot_variability_and_rolling(ax,
                                      ts_sub['time'].values, ts_sub.values,
                                      ts_ann['year'].values, ts_ann.values,
                                      rolling_window, label, kwargs)
    else:
        ax.plot(ts_ann['year'].values, ts_ann.values, label=label, **kwargs)

    units = ds[var_name].attrs.get('units', '')
    long_name = ds[var_name].attrs.get('long_name', var_name)
    ax.set_xlabel('Year', fontsize=16)
    ax.set_ylabel(f"{var_name} ({units})", fontsize=16)
    ax.set_title(f"{long_name} — GIS mean", fontsize=20)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    return ax#, {'sub': ts_sub, 'ann': ts_ann}


def plot_timeseries_by_basin(ds, var_name=None, basin_var='Mouginot_basins',
                              basin_names=None, ax=None,
                              year_start=1940, rolling_window=1,
                              difference=False, **kwargs):
    """
    Plot spatial mean time series for each basin.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing the variable.
    var_name : str
        Variable name to plot.
    basin_var : str
        Basin variable in the mask file. Default 'Mouginot_basins'.
    basin_names : dict, optional
        Mapping {basin_id (int): label (str)}. If None and basin_var is
        'Mouginot_basins', uses MOUGINOT_BASIN_NAMES.
    ax : matplotlib axes, optional
    year_start : float
        Start year for the time series. Default 1940 (first full year).
    rolling_window : int or None
        Bold N-year rolling mean overlaid on the thin sub-annual background.
        Set to None to plot annual values with no rolling smoothing.
        Default 10.
    difference : bool
        If True, subtract each basin's first annual value (at year_start) so
        the series shows the anomaly (e.g. dFAC). Default False.
    **kwargs
        Passed to ax.plot() for all lines.

    Returns
    -------
    (matplotlib axes, dict)
        Axes and ``{basin_name: {'sub': ts_sub, 'ann': ts_ann}}`` DataArrays
        for each basin. Use these to re-plot without recomputing.
    """
    if var_name is None:
        var_name = _infer_var_name(ds)
    ds_mask = _load_mask_aligned(ds)
    basins = ds_mask[basin_var]

    basin_ids = sorted(
        int(v) for v in np.unique(basins.values[~np.isnan(basins.values)]) if v > 0
    )

    if basin_names is None:
        basin_names = MOUGINOT_BASIN_NAMES if basin_var == 'Mouginot_basins' \
            else {i: str(i) for i in basin_ids}

    if ax is None:
        _, ax = plt.subplots(figsize=(14, 5))
    ax.set_prop_cycle(color=OKABE_ITO_COLORS)

    agg = config.get_aggregation_method(var_name)
    data_out = {}
    for basin_id in basin_ids:
        mask = basins == basin_id
        label = basin_names.get(basin_id, str(basin_id))

        # Always compute sub-annual and annual series
        ts_sub = _spatial_mean_timeseries(ds[var_name], mask, agg=agg,
                                          annual=False, year_start=year_start)
        ts_ann = _spatial_mean_timeseries(ds[var_name], mask, agg=agg,
                                          annual=True, year_start=year_start)

        if difference:
            ts_sub = ts_sub - float(ts_sub.values[0])
            ts_ann = ts_ann - float(ts_ann.values[0])

        data_out[label] = {'sub': ts_sub, 'ann': ts_ann}

        if rolling_window is not None:
            # Use a trailing rolling mean for difference plots so the bold line
            # also starts at 0 (first value = first annual value = 0).
            rolling_center = not difference
            _plot_variability_and_rolling(ax,
                                          ts_sub['time'].values, ts_sub.values,
                                          ts_ann['year'].values, ts_ann.values,
                                          rolling_window, label, kwargs,
                                          rolling_center=rolling_center)
        else:
            ax.plot(ts_ann['year'].values, ts_ann.values, label=label, **kwargs)

    ds_mask.close()

    units = ds[var_name].attrs.get('units', '')
    long_name = ds[var_name].attrs.get('long_name', var_name)
    prefix = 'Δ' if difference else ''
    ax.set_xlabel('Year', fontsize=16)
    ax.set_ylabel(f"{prefix}{var_name} ({units})", fontsize=16)
    ax.set_title(f"{prefix}{long_name} by basin", fontsize=20)
    if difference:
        ax.axhline(0, color='k', linewidth=0.8, linestyle='--', alpha=0.4)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    return ax#, data_out


def plot_timeseries_elevation_band(ds, elevation_threshold, var_name=None,
                                   ax=None, year_start=1940, rolling_window=1,
                                   label_above=None, label_below=None, **kwargs):
    """
    Plot area-weighted mean time series above and below an elevation threshold.

    Both lines are plotted on the same axes. Only ice-sheet pixels are included
    (Icemask_GR == 1).

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing the variable.
    var_name : str
        Variable name to plot.
    elevation_threshold : float
        Elevation in metres at which to split the domain.
    ax : matplotlib axes, optional
    year_start : float
        Start year for the time series. Default 1940 (first full year).
    rolling_window : int or None
        If set, plot sub-annual values as a thin background line and overlay
        a bold N-year rolling mean. Set to None to plot annual values only.
        Default 10.
    label_above : str, optional
        Label for the above-threshold line. Default '>XXXX m'.
    label_below : str, optional
        Label for the below-threshold line. Default '≤XXXX m'.
    **kwargs
        Passed to ax.plot() for both lines.

    Returns
    -------
    matplotlib axes
    """
    if var_name is None:
        var_name = _infer_var_name(ds)
    ds_mask = _load_mask_aligned(ds)
    ice_mask = ds_mask['Icemask_GR'] == 1
    topo = ds_mask['Topography']

    mask_above = ice_mask & (topo > elevation_threshold)
    mask_below = ice_mask & (topo <= elevation_threshold)

    agg = config.get_aggregation_method(var_name)
    ds_mask.close()

    if label_above is None:
        label_above = f'>{elevation_threshold:.0f} m'
    if label_below is None:
        label_below = f'≤{elevation_threshold:.0f} m'

    if ax is None:
        _, ax = plt.subplots(figsize=(14, 5))
    ax.set_prop_cycle(color=OKABE_ITO_COLORS)

    for mask, label in ((mask_above, label_above), (mask_below, label_below)):
        if rolling_window is not None:
            ts_sub = _spatial_mean_timeseries(ds[var_name], mask, agg=agg,
                                              annual=False, year_start=year_start)
            ts_ann = _spatial_mean_timeseries(ds[var_name], mask, agg=agg,
                                              annual=True, year_start=year_start)
            _plot_variability_and_rolling(ax,
                                          ts_sub['time'].values, ts_sub.values,
                                          ts_ann['year'].values, ts_ann.values,
                                          rolling_window, label, kwargs)
        else:
            ts = _spatial_mean_timeseries(ds[var_name], mask, agg=agg,
                                          annual=True, year_start=year_start)
            ax.plot(ts['year'].values, ts.values, label=label, **kwargs)

    units = ds[var_name].attrs.get('units', '')
    long_name = ds[var_name].attrs.get('long_name', var_name)
    ax.set_xlabel('Year', fontsize=16)
    ax.set_ylabel(f"{var_name} ({units})", fontsize=16)
    annual = rolling_window is None
    ax.set_title(
        f"{long_name} — above/below {elevation_threshold:.0f} m"
        f"{'  (annual)' if annual else ''}",
        fontsize=20
    )
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    return ax
