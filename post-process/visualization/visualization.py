"""
Visualization functions for IMAU-FDM gridded output.

This module provides tools for:
- Loading and plotting spatial maps
- Time series extraction and visualization
- Multi-variable comparison
- Ice sheet facies classification
- Animation/GIF creation
- Masked difference plotting
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from pathlib import Path

# Import local config
import config


# =============================================================================
# Data Loading
# =============================================================================

def load_gridded_data(var_name, timestep='10day', detrended=None, output_dir=None):
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

    filename = config.get_output_filename(var_name, timestep, detrended)
    filepath = output_dir / filename

    if not filepath.exists():
        # Try alternative (non-detrended if detrended not found, or vice versa)
        alt_filename = config.get_output_filename(var_name, timestep, not detrended)
        alt_filepath = output_dir / alt_filename
        if alt_filepath.exists():
            print(f"Note: Loading {alt_filename} instead")
            filepath = alt_filepath
        else:
            raise FileNotFoundError(f"File not found: {filepath}")

    #print(f"Loading: {filepath.name}")
    return xr.open_dataset(filepath)


# =============================================================================
# Plotting Functions
# =============================================================================

def plot_map(ds, var_name, time_idx=None, time_value=None, ax=None,
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
        data_plot.plot(
            ax=ax, x='x', y='y',
            cmap=cmap, vmin=vmin, vmax=vmax,
            add_colorbar=add_colorbar,
            cbar_kwargs=cbar_kwargs if add_colorbar else None
        )
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
    else:
        # Use native rlat/rlon grid indices
        data.plot(
            ax=ax,
            cmap=cmap, vmin=vmin, vmax=vmax,
            add_colorbar=add_colorbar,
            cbar_kwargs=cbar_kwargs if add_colorbar else None
        )
        ax.set_xlabel('rlon (grid index)')
        ax.set_ylabel('rlat (grid index)')

    ax.set_aspect('equal')
    ax.set_title(title)

    return ax


def plot_timeseries(ds, var_name, lon_point, lat_point, ax=None, label=None):
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

    # Plot
    if label is None:
        label = f"({actual_lon:.2f}°, {actual_lat:.2f}°)"

    ax.plot(ds['time'], data, label=label, linewidth=0.8)

    # Labels
    units = ds[var_name].attrs.get('units', '')
    long_name = ds[var_name].attrs.get('long_name', var_name)
    ax.set_xlabel('Year')
    ax.set_ylabel(f"{var_name} ({units})")
    ax.set_title(f"{long_name} at ({actual_lon:.2f}°, {actual_lat:.2f}°)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    return ax


def plot_annual_mean(ds, var_name, year, ax=None, **kwargs):
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

def classify_zones(year, firn_air_threshold=0.5, melt_threshold=5.0):
    """
    Classify ice sheet zones for a given year based on annual mean values.

    Parameters
    ----------
    year : int
        Year to classify
    firn_air_threshold : float
        Threshold for firn air content (m). Below this = ablation zone.
        Default 0.5 m.
    melt_threshold : float
        Threshold for cumulative annual surface melt (mm w.e.).
        Above this = melt present. Default 0.0.

    Returns
    -------
    xr.DataArray
        Zone classification:
        0 = No data (NaN)
        1 = Ablation zone (FirnAir < threshold)
        2 = Percolation zone (melt > 0 AND FirnAir >= threshold)
        3 = Dry snow zone (no melt AND FirnAir >= threshold)
    """
    # Load FirnAir and surfmelt data
    ds_firn = load_gridded_data('FirnAir')
    ds_melt = load_gridded_data('surfmelt')

    # Get annual means
    firn_year = ds_firn['FirnAir'].sel(time=slice(year, year + 1)).mean(dim='time')
    melt_year = ds_melt['surfmelt'].sel(time=slice(year, year + 1)).sum(dim='time')  # Cumulative melt

    # Create zone array (start with NaN = 0)
    zones = xr.zeros_like(firn_year)

    # Classify zones
    # Ablation: low firn air content
    ablation_mask = firn_year < firn_air_threshold
    zones = zones.where(~ablation_mask, 1)

    # Percolation: has melt AND has firn air
    percolation_mask = (melt_year > melt_threshold) & (firn_year >= firn_air_threshold)
    zones = zones.where(~percolation_mask, 2)

    # Dry snow: no significant melt AND has firn air
    dry_snow_mask = (melt_year <= melt_threshold) & (firn_year >= firn_air_threshold)
    zones = zones.where(~dry_snow_mask, 3)

    # Keep NaN where original data was NaN
    zones = zones.where(~np.isnan(firn_year), np.nan)

    # Add metadata
    zones.attrs['long_name'] = f'Ice sheet facies classification ({year})'
    zones.attrs['flag_values'] = [1, 2, 3]
    zones.attrs['flag_meanings'] = 'ablation_zone percolation_zone dry_snow_zone'

    # Store coordinates for plotting
    zones = zones.assign_coords(
        lat=ds_firn['lat'],
        lon=ds_firn['lon']
    )

    ds_firn.close()
    ds_melt.close()

    return zones


def plot_zones(year, ax=None, firn_air_threshold=0.5, melt_threshold=5.0):
    """
    Plot ablation and percolation zones for a given year.

    Parameters
    ----------
    year : int
        Year to plot
    ax : matplotlib axes, optional
        Axes to plot on
    firn_air_threshold : float
        Threshold for firn air content (m). Default 0.5 m.
    melt_threshold : float
        Threshold for surface melt. Default 0.0.

    Returns
    -------
    matplotlib axes
    """
    # Get zone classification
    zones = classify_zones(year, firn_air_threshold, melt_threshold)

    # Create figure if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 10))

    # Custom colormap: ablation=red, percolation=orange, dry snow=blue
    colors = ['#d62728', '#ff7f0e', '#1f77b4']  # red, orange, blue
    cmap = ListedColormap(colors)
    bounds = [0.5, 1.5, 2.5, 3.5]
    norm = BoundaryNorm(bounds, cmap.N)

    # Plot
    im = zones.plot(
        ax=ax,
        cmap=cmap,
        norm=norm,
        add_colorbar=False
    )

    # Add colorbar with labels
    cbar = plt.colorbar(im, ax=ax, shrink=0.6, aspect=20, pad=0.02)
    cbar.set_ticks([1, 2, 3])
    cbar.set_ticklabels(['Ablation', 'Percolation', 'Dry Snow'])

    ax.set_aspect('equal')
    ax.set_xlabel('rlon (grid index)')
    ax.set_ylabel('rlat (grid index)')
    ax.set_title(f'Ice Sheet Facies ({year})')

    return ax


def plot_zones_timeseries(years, firn_air_threshold=0.5, melt_threshold=5.0,
                          ncols=4, figsize=None):
    """
    Plot zone maps for multiple years.

    Parameters
    ----------
    years : list of int
        Years to plot
    firn_air_threshold : float
        Threshold for firn air content (m). Default 0.5 m.
    melt_threshold : float
        Threshold for surface melt. Default 0.0.
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
                   firn_air_threshold=firn_air_threshold,
                   melt_threshold=melt_threshold)

    # Hide unused axes
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()
    return fig, axes


# =============================================================================
# Statistics and Export
# =============================================================================

def print_stats(ds, var_name):
    """Print basic statistics for a variable."""
    data = ds[var_name]

    print(f"\nStatistics for {var_name}:")
    print(f"  Shape: {data.shape}")
    print(f"  Min: {float(data.min()):.4f}")
    print(f"  Max: {float(data.max()):.4f}")
    print(f"  Mean: {float(data.mean()):.4f}")
    print(f"  Std: {float(data.std()):.4f}")
    print(f"  NaN count: {int(np.isnan(data.values).sum())}")
    print(f"  Valid count: {int((~np.isnan(data.values)).sum())}")


def export_annual_means(ds, var_name, output_file):
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
            frame_plot.plot(
                ax=ax, x='x', y='y',
                cmap=cmap, vmin=vmin, vmax=vmax,
                add_colorbar=True,
                cbar_kwargs=cbar_kwargs
            )
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')
        else:
            frame_data.plot(
                ax=ax,
                cmap=cmap, vmin=vmin, vmax=vmax,
                add_colorbar=True,
                cbar_kwargs=cbar_kwargs
            )
            ax.set_xlabel('rlon (grid index)')
            ax.set_ylabel('rlat (grid index)')

        ax.set_aspect('equal')

        # Set title
        title = title_fmt.format(
            var_name=var_name,
            long_name=long_name,
            time=float(t),
            units=units
        )
        ax.set_title(title)

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
                           cmap='RdBu_r', figsize=(10, 12)):
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
    fig, ax, result_data
        Figure, axes, and the computed data array

    Examples
    --------
    >>> # Where does melt occur but not runoff? (indicates refreezing)
    >>> fig, ax, data = plot_masked_difference('surfmelt', 'Runoff',
    ...                                        time_slice=slice(2019, 2020),
    ...                                        comparison='var1_only')

    >>> # Where does runoff occur without melt? (rain-driven runoff)
    >>> fig, ax, data = plot_masked_difference('Runoff', 'surfmelt',
    ...                                        time_value=2012.5,
    ...                                        comparison='var1_only')

    >>> # Show difference between melt and runoff (positive = retained)
    >>> fig, ax, data = plot_masked_difference('surfmelt', 'Runoff',
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
        cmap = 'RdYlBu_r'
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
    result.plot(
        ax=ax,
        cmap=cmap,
        add_colorbar=True,
        cbar_kwargs=cbar_kwargs
    )

    ax.set_aspect('equal')
    ax.set_xlabel('rlon (grid index)')
    ax.set_ylabel('rlat (grid index)')
    ax.set_title(title)

    plt.tight_layout()

    ds1.close()
    ds2.close()

    return fig, ax, result


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
    axes[0, 0].set_title('Melt but no Runoff\n(Refreezing)')

    # Runoff but no melt (rain-driven)
    plot_masked_difference('Runoff', 'surfmelt',
                           time_value=time_value, time_slice=time_slice,
                           comparison='var1_only', ax=axes[0, 1])
    axes[0, 1].set_title('Runoff but no Melt\n(Rain-driven)')

    # Both present
    plot_masked_difference('surfmelt', 'Runoff',
                           time_value=time_value, time_slice=time_slice,
                           comparison='both', ax=axes[1, 0])
    axes[1, 0].set_title('Both Melt and Runoff\n(Melt-driven runoff)')

    # Difference
    plot_masked_difference('surfmelt', 'Runoff',
                           time_value=time_value, time_slice=time_slice,
                           comparison='difference', ax=axes[1, 1])
    axes[1, 1].set_title('Melt - Runoff\n(Positive = retained)')

    plt.tight_layout()
    return fig, axes
