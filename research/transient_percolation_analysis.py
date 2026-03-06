"""
Transient percolation zone analysis for IMAU-FDM output.

Three analyses, each supporting any gridded variable registered in
config.VARIABLE_FILES (currently 'surfmelt' and 'Runoff'):

  1. map_melt_frequency
     Map of how many years each pixel exceeded an annual event threshold.
     Can cover the full simulation or be restricted to the last N years.

  2. plot_melt_histogram
     Side-by-side bar chart comparing the distribution of per-pixel event
     counts between two user-specified time periods.

  3. plot_percolation_migration
     Map classifying whether pixels stayed in the low-event category, stayed
     in the high-event category, or transitioned between the two when
     comparing an early and a late window.

All functions default to variable='surfmelt'.  Pass variable='Runoff' (or any
other key in config.VARIABLE_FILES) to run the same analysis on a different
field.  All other defaults come from config.py and can be overridden inline.
"""

import sys
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap, BoundaryNorm

sys.path.insert(0, str(Path(__file__).parent))
import config


# =============================================================================
# Internal helper
# =============================================================================

def _load_annual_totals(variable='surfmelt', nc_file=None):
    """
    Load a gridded variable file and return annual totals per pixel.

    The source file has a fractional-year 'time' coordinate (e.g., 1990.27).
    This function truncates each value to its calendar year, groups all
    timesteps within a year together, and sums them to give the annual total
    at each grid point.

    Off-ice pixels (NaN across all timesteps) are preserved as NaN.
    The 2-D lat/lon coordinates and variable metadata (long_name, units) are
    carried through so that xarray's .plot() and plot titles work correctly.

    Parameters
    ----------
    variable : str
        Variable name as it appears in the NetCDF file and in
        config.VARIABLE_FILES.  Default: 'surfmelt'.
    nc_file : Path or str, optional
        Override the file path.  If None, uses config.VARIABLE_FILES[variable].

    Returns
    -------
    xr.DataArray
        Shape (year, rlat, rlon).  Dimension 'year' is integer calendar years.
        Attributes 'long_name' and 'units' are copied from the source file.
    """
    if nc_file is None:
        if variable not in config.VARIABLE_FILES:
            raise ValueError(
                f"Unknown variable '{variable}'. "
                f"Add it to config.VARIABLE_FILES. "
                f"Currently available: {list(config.VARIABLE_FILES)}"
            )
        nc_file = config.VARIABLE_FILES[variable]

    # chunks splits only the time axis so peak memory is ~chunk_size x spatial x 4 bytes
    # instead of the full ~1.8 GB time series — safe to call twice in the same cell.
    ds = xr.open_dataset(nc_file, chunks={'time': 30})

    # Identify off-ice pixels before any arithmetic (sum over all-NaN gives 0).
    # Compute immediately to hold only a tiny 2-D bool array rather than a dask
    # graph that would re-read the file again during the later groupby.
    all_nan_mask = ds[variable].isnull().all('time').compute()

    # Assign integer year coordinate and group-sum within each year
    year_vals = ds['time'].values.astype(int)
    annual = (
        ds[variable]
        .assign_coords(year=('time', year_vals))
        .groupby('year')
        .sum('time')
        .compute()   # materialise small annual array; releases full time series
    )

    # Restore the off-ice NaN mask
    annual = annual.where(~all_nan_mask)

    # Carry lat/lon and variable metadata for plotting
    annual = annual.assign_coords(lat=ds['lat'], lon=ds['lon'])
    annual.attrs['long_name'] = ds[variable].attrs.get('long_name', variable)
    annual.attrs['units']     = ds[variable].attrs.get('units', 'mm w.e.')

    ds.close()
    return annual


# =============================================================================
# Internal colourmap helper
# =============================================================================

def _make_zero_grey_cmap(base_cmap_name, max_count, zero_color='#aaaaaa'):
    """
    Return (cmap, norm) where count=0 maps to `zero_color` and counts
    1 through max_count are sampled uniformly from base_cmap_name.
    """
    base = plt.get_cmap(base_cmap_name)
    n = max(int(max_count), 1)
    pos_colors = base(np.linspace(0.05, 1.0, n))
    all_colors = np.vstack([mcolors.to_rgba(zero_color), pos_colors])
    cmap = ListedColormap(all_colors)
    boundaries = np.arange(-0.5, n + 1.5)
    norm = BoundaryNorm(boundaries, cmap.N)
    return cmap, norm


# =============================================================================
# Analysis 1 – Event frequency map
# =============================================================================

def map_melt_frequency(
    variable='surfmelt',
    threshold=None,
    years_back=None,
    ax=None,
    cmap='YlOrRd',
    figsize=(8, 10),
):
    """
    Map the number of years each pixel exceeded an annual event threshold.

    Parameters
    ----------
    variable : str
        Variable to analyse.  Must be a key in config.VARIABLE_FILES.
        Default: 'surfmelt'.
    threshold : float, optional
        Annual total (mm w.e.) required to count as an event year.
        Default: config.EVENT_THRESHOLDS[variable].
    years_back : int, optional
        If given, restrict the count to the last `years_back` calendar years
        (e.g., years_back=10 for the 10 most recent years in the dataset).
        If None, counts over the full simulation.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on.  A new figure is created if None.
    cmap : str
        Colormap.  Default 'YlOrRd'.
    figsize : tuple
        Figure size when creating a new figure.

    Returns
    -------
    count : xr.DataArray
        Number of event years at each grid point.
    ax : matplotlib.axes.Axes
    """
    if threshold is None:
        threshold = config.EVENT_THRESHOLDS[variable]

    annual = _load_annual_totals(variable)
    long_name = annual.attrs['long_name']
    units     = annual.attrs['units']

    # Optionally restrict to last N years
    if years_back is not None:
        last_year  = int(annual['year'].values[-1])
        start_year = last_year - years_back + 1
        annual     = annual.sel(year=slice(start_year, last_year))
        period_str = f'last {years_back} years ({start_year}–{last_year})'
    else:
        first_year = int(annual['year'].values[0])
        last_year  = int(annual['year'].values[-1])
        period_str = f'{first_year}–{last_year}'

    # Count event years; preserve the off-ice NaN mask
    valid = ~annual.isnull().all('year')
    count = (annual > threshold).sum('year').where(valid)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    max_count = int(count.max().item())
    custom_cmap, norm = _make_zero_grey_cmap(cmap, max_count)

    # Colorbar ticks: always include 0, then every 10/5/1 depending on range
    step = 10 if max_count > 30 else (5 if max_count > 10 else 1)
    ticks = [0] + list(range(step, max_count, step)) + [max_count]
    ticks = sorted(set(ticks))  # deduplicate if step lands on max_count

    count.plot(
        ax=ax,
        cmap=custom_cmap,
        norm=norm,
        cbar_kwargs={
            'label': 'Number of event years',
            'shrink': 0.6,
            'aspect': 20,
            'pad': 0.02,
            'ticks': ticks,
        },
    )

    # Coastline: contour at the land/ocean boundary (topo > 0 = land)
    _masks = xr.open_dataset(config.MASKS_FILE)
    _topo  = _masks['Topography'].squeeze().values.astype(np.float32)
    _masks.close()
    _coast = np.where(np.isfinite(_topo) & (_topo > 0), 1.0, 0.0)
    ax.contour(_coast, levels=[0.5], colors='k', linewidths=0.7, alpha=0.85)

    ax.set_aspect('equal')
    ax.set_title(
        f'{long_name} event frequency  ({period_str})\n'
        f'threshold: {threshold} {units} yr⁻¹'
    )
    ax.set_xlabel('rlon (grid index)')
    ax.set_ylabel('rlat (grid index)')

    return count, ax


# =============================================================================
# Analysis 2 – Event frequency histogram
# =============================================================================

def plot_melt_histogram(
    variable='surfmelt',
    period1=None,
    period2=None,
    threshold=None,
    ax=None,
    period1_label=None,
    period2_label=None,
    figsize=(10, 6),
):
    """
    Side-by-side bar chart comparing per-pixel event counts between two periods.

    The x-axis shows how many years within the period a pixel exceeded the
    threshold.  The y-axis shows how many ice-covered pixels fall into each
    count bin.  Two grouped bars per bin allow direct visual comparison.

    Parameters
    ----------
    variable : str
        Variable to analyse.  Default: 'surfmelt'.
    period1 : tuple of (int, int), optional
        (start_year, end_year) for the first period, both inclusive.
        Default: config.PERIOD_1.
    period2 : tuple of (int, int), optional
        (start_year, end_year) for the second period, both inclusive.
        Default: config.PERIOD_2.
    threshold : float, optional
        Annual event threshold (mm w.e.).
        Default: config.EVENT_THRESHOLDS[variable].
    ax : matplotlib.axes.Axes, optional
    period1_label, period2_label : str, optional
        Legend labels.  Defaults to the year ranges.
    figsize : tuple

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    if period1 is None:
        period1 = config.PERIOD_1
    if period2 is None:
        period2 = config.PERIOD_2
    if threshold is None:
        threshold = config.EVENT_THRESHOLDS[variable]
    if period1_label is None:
        period1_label = f'{period1[0]}–{period1[1]}'
    if period2_label is None:
        period2_label = f'{period2[0]}–{period2[1]}'

    annual = _load_annual_totals(variable)
    long_name = annual.attrs['long_name']
    units     = annual.attrs['units']

    # Keep only valid (ice-covered) pixels
    valid = ~annual.isnull().all('year')

    def _event_counts(period):
        subset = annual.sel(year=slice(period[0], period[1]))
        count  = (subset > threshold).sum('year').where(valid)
        return count.values[valid.values].astype(int)

    counts1 = _event_counts(period1)
    counts2 = _event_counts(period2)

    # Grouped bar chart with one bar per integer count value
    max_count = max(counts1.max(), counts2.max())
    bins  = np.arange(0, max_count + 2)
    h1, _ = np.histogram(counts1, bins=bins)
    h2, _ = np.histogram(counts2, bins=bins)

    x     = bins[:-1]
    width = 0.4

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    ax.bar(x - width / 2, h1, width, label=period1_label, color='steelblue', alpha=0.85)
    ax.bar(x + width / 2, h2, width, label=period2_label, color='tomato',    alpha=0.85)

    ax.set_xlabel(f'{long_name} event years in period')
    ax.set_ylabel('Number of pixels')
    ax.set_title(
        f'{long_name} event count distribution per pixel\n'
        f'threshold: {threshold} {units} yr⁻¹'
    )
    ax.set_xticks(x)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    return ax


# =============================================================================
# Analysis 3 – Zone migration map
# =============================================================================

def plot_percolation_migration(
    variable='surfmelt',
    mode='events',
    early_window=None,
    late_window=None,
    # --- events mode ---
    dry_events_max=None,
    percolation_events_min=None,
    threshold=None,
    # --- amount mode ---
    melt_threshold=10,
    ax=None,
    figsize=(8, 10),
    add_legend=True,
    topo_contour_interval=200,
):
    """
    Map pixels that shifted from a low to a high regime between two windows.

    Two classification modes are available via the `mode` parameter:

    mode='events'  (default)
        A pixel is classified as low  if annual melt exceeded `threshold`
        in ≤ dry_events_max years within the window.
        A pixel is classified as high if it exceeded the threshold in
        ≥ percolation_events_min years.
        Pixels with counts in between are left ambiguous (NaN).

    mode='amount'
        A pixel is classified as high if its *mean* annual melt over the
        window is ≥ melt_threshold, otherwise as low.  There is no
        ambiguous zone.

    In both modes the resulting map has four categories:
      1  Remained low
      2  Remained high
      3  Expanded: low → high   ← main signal
      4  Contracted: high → low
      NaN  Off-ice (always), or ambiguous count (events mode only)

    Parameters
    ----------
    variable : str
        Variable to analyse.  Default: 'surfmelt'.
    mode : {'events', 'amount'}
        Classification mode.  Default: 'events'.
    early_window : tuple of (int, int), optional
        (start_year, end_year) for the early comparison window.
        Default: config.EARLY_WINDOW.
    late_window : tuple of (int, int), optional
        (start_year, end_year) for the late comparison window.
        Default: config.LATE_WINDOW.
    dry_events_max : int, optional
        [events mode] Max event years to classify as low.
        Default: config.DRY_EVENTS_MAX.
    percolation_events_min : int, optional
        [events mode] Min event years to classify as high.
        Default: config.PERCOLATION_EVENTS_MIN.
    threshold : float, optional
        [events mode] Annual total (mm w.e.) required to count as an event.
        Default: config.EVENT_THRESHOLDS[variable].
    melt_threshold : float
        [amount mode] Mean annual melt (mm w.e.) separating low from high.
        Required when mode='amount'.
    ax : matplotlib.axes.Axes, optional
    figsize : tuple
    add_legend : bool
        Whether to draw the category legend on the axes (default True).
        Set to False when embedding in a panel with a shared legend.
    topo_contour_interval : int or None
        Draw labeled topographic contour lines every this many metres so you
        can read the elevation at the orange/red→blue transition boundary.
        Default: 200 m.  Set to None to suppress contours.

    Returns
    -------
    classification : xr.DataArray
        Integer category map (float dtype; NaN for off-ice / ambiguous).
    ax : matplotlib.axes.Axes
    """
    if mode not in ('events', 'amount'):
        raise ValueError(f"mode must be 'events' or 'amount', got '{mode}'")
    if mode == 'amount' and melt_threshold is None:
        raise ValueError("melt_threshold must be provided when mode='amount'")

    if early_window is None:
        early_window = config.EARLY_WINDOW
    if late_window is None:
        late_window = config.LATE_WINDOW
    if dry_events_max is None:
        dry_events_max = config.DRY_EVENTS_MAX
    if percolation_events_min is None:
        percolation_events_min = config.PERCOLATION_EVENTS_MIN
    if threshold is None:
        threshold = config.EVENT_THRESHOLDS[variable]

    annual = _load_annual_totals(variable)
    long_name = annual.attrs['long_name']
    units     = annual.attrs['units']
    valid = ~annual.isnull().all('year')

    # --- classify_window: returns (is_low, is_high) boolean DataArrays ---
    if mode == 'events':
        def _classify_window(window):
            subset = annual.sel(year=slice(window[0], window[1]))
            count  = (subset > threshold).sum('year').where(valid)
            return count <= dry_events_max, count >= percolation_events_min
    else:  # amount
        def _classify_window(window):
            subset    = annual.sel(year=slice(window[0], window[1]))
            mean_melt = subset.mean('year').where(valid)
            is_high   = mean_melt >= melt_threshold
            return ~is_high, is_high

    was_low, was_high = _classify_window(early_window)
    now_low, now_high = _classify_window(late_window)

    # Build classification: start all-NaN, fill each category in turn.
    # DataArray.where(cond, other) keeps original where cond is True,
    # replaces with `other` where cond is False.
    base = xr.full_like(was_low.astype(float), fill_value=np.nan)
    classification = base.where(~(was_low  & now_low),  1.0)  # stayed low
    classification = classification.where(~(was_high & now_high), 2.0)  # stayed high
    classification = classification.where(~(was_low  & now_high), 3.0)  # expanded
    classification = classification.where(~(was_high & now_low),  4.0)  # contracted

    # Final safety mask: off-ice pixels stay NaN
    classification = classification.where(valid)

    # --- Plot ---
    colors = ['#AED6F1', '#F39C12', '#E74C3C', '#8E44AD']  # blue, orange, red, purple
    cat_cmap = ListedColormap(colors)
    norm     = BoundaryNorm([0.5, 1.5, 2.5, 3.5, 4.5], cat_cmap.N)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    classification.plot(ax=ax, cmap=cat_cmap, norm=norm, add_colorbar=False)

    # Load topography once for both coastline and elevation contours
    masks = xr.open_dataset(config.MASKS_FILE)
    topo  = masks['Topography'].squeeze().values.astype(np.float32)
    masks.close()

    # Coastline: contour at the land/ocean boundary (topo > 0 = land)
    coast = np.where(np.isfinite(topo) & (topo > 0), 1.0, 0.0)
    ax.contour(coast, levels=[0.5], colors='k', linewidths=0.7, alpha=0.85)

    # Topographic contours so the elevation at the transition boundary is readable
    if topo_contour_interval is not None:
        topo_min = int(np.nanmin(topo[topo > 0]) // topo_contour_interval) * topo_contour_interval
        topo_max = int(np.nanmax(topo) // topo_contour_interval + 1) * topo_contour_interval
        levels = np.arange(topo_min, topo_max + 1, topo_contour_interval)
        cs = ax.contour(topo, levels=levels, colors='k',
                        linewidths=0.4, alpha=0.5)
        ax.clabel(cs, fmt='%d m', fontsize=6, inline=True, inline_spacing=2)

    # Physically meaningful category names
    if variable == 'surfmelt':
        low_name  = 'dry facies'
        high_name = 'percolation zone'
    else:
        low_name  = f'low {long_name.lower()}'
        high_name = f'high {long_name.lower()}'

    # Legend labels depend on mode
    if mode == 'events':
        low_label  = f'Remained {low_name}  (≤{dry_events_max} event yrs in both)'
        high_label = f'Remained {high_name}  (≥{percolation_events_min} event yrs in both)'
    else:
        low_label  = f'Remained {low_name}  (mean < {melt_threshold} {units})'
        high_label = f'Remained {high_name}  (mean ≥ {melt_threshold} {units})'

    patches = [
        mpatches.Patch(color='#AED6F1', label=low_label),
        mpatches.Patch(color='#F39C12', label=high_label),
        mpatches.Patch(color='#E74C3C', label=f'Expanded: {low_name} → {high_name}'),
        mpatches.Patch(color='#8E44AD', label=f'Contracted: {high_name} → {low_name}'),
    ]
    if add_legend:
        ax.legend(handles=patches, loc='lower right', fontsize=8, framealpha=0.9)

    n_early   = early_window[1] - early_window[0] + 1
    n_late    = late_window[1]  - late_window[0]  + 1
    early_str = f'{early_window[0]}–{early_window[1]}  ({n_early} yr)'
    late_str  = f'{late_window[0]}–{late_window[1]}  ({n_late} yr)'

    if mode == 'events':
        criterion_str = (
            f'Low: ≤{dry_events_max} event yrs,  '
            f'High: ≥{percolation_events_min} event yrs  '
            f'(event threshold: {threshold} {units} yr⁻¹)'
        )
    else:
        criterion_str = (
            f'Low: mean annual {long_name.lower()} < {melt_threshold} {units}  '
            f'High: ≥ {melt_threshold} {units}'
        )

    ax.set_aspect('equal')
    ax.set_title(
        f'{long_name} zone migration\n'
        f'Early: {early_str}   Late: {late_str}\n'
        f'{criterion_str}'
    )
    ax.set_xlabel('rlon (grid index)')
    ax.set_ylabel('rlat (grid index)')

    return classification, ax


# =============================================================================
# Analysis 3b – Multi-year migration panel
# =============================================================================

def plot_migration_panel(
    year_pairs,
    variable='surfmelt',
    mode='events',
    window_width=1,
    dry_events_max=None,
    percolation_events_min=None,
    threshold=None,
    melt_threshold=10,
    ncols=3,
    figsize=None,
    title=None,
    topo_contour_interval=200,
):
    """
    Plot a grid of zone-migration maps, one per year-pair, with a shared legend.

    Each element of `year_pairs` is a ``(y1, y2)`` tuple.  The early window
    runs from ``y1`` to ``y1 + window_width - 1`` and the late window from
    ``y2`` to ``y2 + window_width - 1``.  For single-year comparisons
    (the default) each window is just that one year.

    Parameters
    ----------
    year_pairs : list of (int, int)
        Each tuple ``(y1, y2)`` gives the start year of the early and late
        windows, e.g. ``[(2019, 2020), (2020, 2021), (2021, 2022)]``.
    variable : str
        Variable to analyse.  Default: 'surfmelt'.
    mode : {'events', 'amount'}
        Classification mode (see plot_percolation_migration).
    window_width : int
        Number of years in each window.  Default: 1 (single-year snapshots).
    dry_events_max : int, optional
        [events mode] Max event years to classify as low.
        Default: config.DRY_EVENTS_MAX.
    percolation_events_min : int, optional
        [events mode] Min event years to classify as high.
        Default: config.PERCOLATION_EVENTS_MIN.
    threshold : float, optional
        [events mode] Annual event threshold (mm w.e.).
        Default: config.EVENT_THRESHOLDS[variable].
    melt_threshold : float
        [amount mode] Mean annual melt (mm w.e.) separating low from high.
    ncols : int
        Number of columns in the subplot grid.  Default: 3.
    figsize : tuple, optional
        Figure size.  Auto-sized if None.
    title : str, optional
        Overall figure title (suptitle).

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : (n_panels,) ndarray of matplotlib.axes.Axes
    """
    # Resolve config defaults needed for the shared legend
    if dry_events_max is None:
        dry_events_max = config.DRY_EVENTS_MAX
    if percolation_events_min is None:
        percolation_events_min = config.PERCOLATION_EVENTS_MIN
    if threshold is None:
        threshold = config.EVENT_THRESHOLDS.get(variable, 10.0)

    n            = len(year_pairs)
    nrows        = int(np.ceil(n / ncols))
    ncols_actual = min(n, ncols)

    if figsize is None:
        figsize = (ncols_actual * 4, nrows * 5 + 1.5)

    from matplotlib.gridspec import GridSpec as _GS
    top_val = 0.90 if title else 0.97
    fig     = plt.figure(figsize=figsize)
    _gs     = _GS(nrows, ncols_actual, figure=fig,
                  left=0.02, right=0.98, bottom=0.13, top=top_val,
                  wspace=0.05, hspace=0.15)
    axes_flat = [fig.add_subplot(_gs[i // ncols_actual, i % ncols_actual])
                 for i in range(n)]

    for i, (y1, y2) in enumerate(year_pairs):
        ax = axes_flat[i]
        early_window = (y1, y1 + window_width - 1)
        late_window  = (y2, y2 + window_width - 1)

        plot_percolation_migration(
            variable=variable,
            mode=mode,
            early_window=early_window,
            late_window=late_window,
            dry_events_max=dry_events_max,
            percolation_events_min=percolation_events_min,
            threshold=threshold,
            melt_threshold=melt_threshold,
            ax=ax,
            add_legend=False,
            topo_contour_interval=topo_contour_interval,
        )
        ax.set_title(f'{y1}\u2192{y2}', fontsize=10)
        ax.set_xlabel('')
        ax.set_ylabel('')

    # Build a shared legend from the known category colours
    if variable == 'surfmelt':
        low_name  = 'dry facies'
        high_name = 'percolation zone'
    else:
        low_name  = f'low {variable}'
        high_name = f'high {variable}'

    units = 'mm w.e.'
    if mode == 'events':
        low_label  = f'Remained {low_name}  (\u2264{dry_events_max} event yrs in both)'
        high_label = f'Remained {high_name}  (\u2265{percolation_events_min} event yrs in both)'
    else:
        low_label  = f'Remained {low_name}  (mean < {melt_threshold} {units})'
        high_label = f'Remained {high_name}  (mean \u2265 {melt_threshold} {units})'

    legend_patches = [
        mpatches.Patch(color='#AED6F1', label=low_label),
        mpatches.Patch(color='#F39C12', label=high_label),
        mpatches.Patch(color='#E74C3C', label=f'Expanded: {low_name} \u2192 {high_name}'),
        mpatches.Patch(color='#8E44AD', label=f'Contracted: {high_name} \u2192 {low_name}'),
    ]
    # Legend placed in the bottom strip reserved by the GridSpec bottom margin
    fig.legend(
        handles=legend_patches,
        loc='lower center',
        ncol=2,
        fontsize=9,
        framealpha=0.9,
        bbox_to_anchor=(0.5, 0.01),
    )

    if title:
        fig.suptitle(title, fontsize=12)

    return fig, axes_flat
