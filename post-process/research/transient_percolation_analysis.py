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

    ds = xr.open_dataset(nc_file)

    # Identify off-ice pixels before any arithmetic (sum over all-NaN gives 0)
    all_nan_mask = ds[variable].isnull().all('time')

    # Assign integer year coordinate and group-sum within each year
    year_vals = ds['time'].values.astype(int)
    annual = (
        ds[variable]
        .assign_coords(year=('time', year_vals))
        .groupby('year')
        .sum('time')
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

    count.plot(
        ax=ax,
        cmap=cmap,
        cbar_kwargs={
            'label': 'Number of event years',
            'shrink': 0.6,
            'aspect': 20,
            'pad': 0.02,
        },
    )
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
    early_window=None,
    late_window=None,
    dry_events_max=None,
    percolation_events_min=None,
    threshold=None,
    ax=None,
    figsize=(8, 10),
):
    """
    Map pixels that shifted from a low-event to a high-event regime.

    Each pixel is classified in an early and a late window based on how many
    event years it accumulated.  The resulting map has four categories:

      1  Remained low   (≤ dry_events_max in both windows)
      2  Remained high  (≥ percolation_events_min in both windows)
      3  Expanded: low → high                    ← main signal
      4  Contracted: high → low
      NaN  Off-ice, or ambiguous (count between the two thresholds)

    When variable='surfmelt' the low/high categories correspond to dry facies
    and percolation zone respectively; for other variables the same thresholds
    define analogous low-event and high-event regimes.

    Parameters
    ----------
    variable : str
        Variable to analyse.  Default: 'surfmelt'.
    early_window : tuple of (int, int), optional
        (start_year, end_year) for the early comparison window.
        Default: config.EARLY_WINDOW.
    late_window : tuple of (int, int), optional
        (start_year, end_year) for the late comparison window.
        Default: config.LATE_WINDOW.
    dry_events_max : int, optional
        Maximum event years in a window to classify as low-event.
        Default: config.DRY_EVENTS_MAX (2).
    percolation_events_min : int, optional
        Minimum event years in a window to classify as high-event.
        Default: config.PERCOLATION_EVENTS_MIN (3).
    threshold : float, optional
        Annual event threshold (mm w.e.).
        Default: config.EVENT_THRESHOLDS[variable].
    ax : matplotlib.axes.Axes, optional
    figsize : tuple

    Returns
    -------
    classification : xr.DataArray
        Integer category map (float dtype; NaN for off-ice / ambiguous).
    ax : matplotlib.axes.Axes
    """
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
    valid = ~annual.isnull().all('year')

    def _count_events(window):
        subset = annual.sel(year=slice(window[0], window[1]))
        return (subset > threshold).sum('year').where(valid)

    events_early = _count_events(early_window)
    events_late  = _count_events(late_window)

    # Classification masks (NaN pixels compare as False and stay unclassified)
    was_low  = events_early <= dry_events_max
    now_low  = events_late  <= dry_events_max
    was_high = events_early >= percolation_events_min
    now_high = events_late  >= percolation_events_min

    # Build classification: start all-NaN, fill each category in turn.
    # DataArray.where(cond, other) keeps original where cond is True,
    # replaces with `other` where cond is False.
    classification = xr.full_like(events_early, fill_value=np.nan, dtype=float)
    classification = classification.where(~(was_low  & now_low),  1.0)  # stayed low
    classification = classification.where(~(was_high & now_high), 2.0)  # stayed high
    classification = classification.where(~(was_low  & now_high), 3.0)  # expanded
    classification = classification.where(~(was_high & now_low),  4.0)  # contracted

    # Final safety mask: off-ice pixels stay NaN
    classification = classification.where(valid)

    # --- Plot ---
    colors = ['#AED6F1', '#F39C12', '#E74C3C', '#8E44AD']  # blue, orange, red, purple
    cmap   = ListedColormap(colors)
    norm   = BoundaryNorm([0.5, 1.5, 2.5, 3.5, 4.5], cmap.N)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    classification.plot(ax=ax, cmap=cmap, norm=norm, add_colorbar=False)

    # Use physically meaningful labels for surfmelt; generic labels otherwise
    if variable == 'surfmelt':
        low_name  = 'dry facies'
        high_name = 'percolation zone'
    else:
        low_name  = f'low {long_name.lower()}'
        high_name = f'high {long_name.lower()}'

    n_early = early_window[1] - early_window[0] + 1
    n_late  = late_window[1]  - late_window[0]  + 1
    patches = [
        mpatches.Patch(color='#AED6F1',
                       label=f'Remained {low_name}  (≤{dry_events_max} event yrs in both)'),
        mpatches.Patch(color='#F39C12',
                       label=f'Remained {high_name}  (≥{percolation_events_min} event yrs in both)'),
        mpatches.Patch(color='#E74C3C',
                       label=f'Expanded: {low_name} → {high_name}'),
        mpatches.Patch(color='#8E44AD',
                       label=f'Contracted: {high_name} → {low_name}'),
    ]
    ax.legend(handles=patches, loc='lower right', fontsize=8, framealpha=0.9)

    early_str = f'{early_window[0]}–{early_window[1]}  ({n_early} yr)'
    late_str  = f'{late_window[0]}–{late_window[1]}  ({n_late} yr)'
    ax.set_aspect('equal')
    ax.set_title(
        f'{long_name} zone migration\n'
        f'Early: {early_str}   Late: {late_str}\n'
        f'Low: ≤{dry_events_max} event yrs,  High: ≥{percolation_events_min} event yrs'
    )
    ax.set_xlabel('rlon (grid index)')
    ax.set_ylabel('rlat (grid index)')

    return classification, ax
