"""
Zone transition analysis for IMAU-FDM gridded output.

Identifies pixels that switched between the dry-snow and percolation zones
between consecutive years, based on whether annual surface melt exceeded a
threshold.  Also flags pixels that stayed in each zone for context.

Public interface
----------------
classify_zone_transitions(year, melt_threshold=10.0)
    Compare year-1 vs year.  Returns an xr.DataArray with four categories:
      1 = dry both years       (unchanged)
      2 = percolation both     (unchanged)
      3 = dry  → percolation   (newly melt-active)
      4 = percolation → dry    (newly melt-inactive)
    NaN where data are missing or the pixel is off-ice.

plot_zone_transitions(year, melt_threshold=10.0, ax=None, show_unchanged=True)
    Map of the four (or two) transition categories for a given year pair.
"""

import sys
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

sys.path.insert(0, str(Path(__file__).parent))
import config


# =============================================================================
# Internal helpers
# =============================================================================

def _annual_melt(year):
    """Annual cumulative surface melt [mm w.e.] for *year*."""
    ds = xr.open_dataset(config.VARIABLE_FILES['surfmelt'], decode_times=False)
    out = ds['surfmelt'].sel(time=slice(year, year + 1)).sum(dim='time')
    ds.close()
    return out


def _ice_mask(year):
    """
    Boolean DataArray — True where the pixel is on the ice sheet in *year*.

    vacc is NaN off-ice, so any pixel where vacc is finite is on-ice.
    """
    ds = xr.open_dataset(config.VARIABLE_FILES['vacc'], decode_times=False)
    vacc = ds['vacc'].sel(time=slice(year, year + 1)).mean(dim='time')
    ds.close()
    return vacc.notnull()


# =============================================================================
# Classification
# =============================================================================

def classify_zone_transitions(year, melt_threshold=10.0):
    """
    Classify each grid pixel by its melt-zone transition between year-1 and year.

    Parameters
    ----------
    year : int
        The *current* year.  The *prior* year is year-1.
    melt_threshold : float
        Annual cumulative surface melt (mm w.e.) that separates the dry-snow
        zone (below threshold) from the percolation zone (at or above).
        Default 10.0 mm w.e.

    Returns
    -------
    xr.DataArray  (rlat, rlon)
        Zone classification:
          1 = dry  both years          (unchanged dry-snow)
          2 = perc both years          (unchanged percolation)
          3 = dry  → percolation       (newly melt-active)
          4 = percolation → dry        (newly melt-inactive)
          NaN = off-ice in either year, or no data
    """
    melt_prev = _annual_melt(year - 1)
    melt_curr = _annual_melt(year)

    # Off-ice mask: must be on-ice in both years
    on_ice = _ice_mask(year - 1) & _ice_mask(year)

    perc_prev = melt_prev >= melt_threshold
    perc_curr = melt_curr >= melt_threshold

    transitions = xr.full_like(melt_curr, np.nan)
    transitions = transitions.where( perc_prev |  perc_curr, 1)   # dry  both
    transitions = transitions.where(~perc_prev | ~perc_curr, 2)   # perc both
    transitions = transitions.where( perc_prev | ~perc_curr, 3)   # dry → perc
    transitions = transitions.where(~perc_prev |  perc_curr, 4)   # perc → dry

    transitions = transitions.where(on_ice, np.nan)

    transitions.attrs.update({
        'long_name'     : f'Melt-zone transition {year - 1}\u2192{year}',
        'melt_threshold': melt_threshold,
        'flag_values'   : [1, 2, 3, 4],
        'flag_meanings' : 'dry_unchanged perc_unchanged dry_to_perc perc_to_dry',
    })
    return transitions


# =============================================================================
# Plotting
# =============================================================================

def plot_zone_transitions(year, melt_threshold=10.0, ax=None,
                          show_unchanged=True, show_colorbar=False,
                          add_title=True):
    """
    Map of melt-zone transitions between year-1 and year.

    Parameters
    ----------
    year : int
        Current year; compared against year-1.
    melt_threshold : float
        Annual melt threshold separating dry from percolation (mm w.e.).
    ax : matplotlib Axes, optional
        Axes to draw on.  Created if None.
    show_unchanged : bool
        If True (default), show all four categories.
        If False, show only the two transition categories (3 and 4); unchanged
        pixels are masked to white so transitions stand out.
    add_title : bool
        If True (default), set a descriptive axes title.  Set to False when
        embedding in a panel so the caller can set a short title instead.

    Returns
    -------
    ax : matplotlib Axes
    """
    transitions = classify_zone_transitions(year, melt_threshold)

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 10))

    if show_unchanged:
        colors = ['#c8c8c8', '#42628F', '#2ca02c', '#d62728']
        bounds = [0.5, 1.5, 2.5, 3.5, 4.5]
        labels = [
            'Dry (unchanged)',
            'Percolation (unchanged)',
            'Dry \u2192 Percolation',
            'Percolation \u2192 Dry',
        ]
        plot_data = transitions
    else:
        colors = ['#2ca02c', '#d62728']
        bounds = [2.5, 3.5, 4.5]
        labels = ['Dry \u2192 Percolation', 'Percolation \u2192 Dry']
        plot_data = transitions.where(transitions >= 3)

    cmap = ListedColormap(colors)
    norm = BoundaryNorm(bounds, cmap.N)

    im = plot_data.plot(ax=ax, cmap=cmap, norm=norm, add_colorbar=False)

    if show_colorbar:
    
        cbar = plt.colorbar(im, ax=ax, shrink=0.6, aspect=20, pad=0.02)
        ticks = [1, 2, 3, 4] if show_unchanged else [3, 4]
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(labels)

    ax.set_aspect('equal')
    ax.set_xlabel('rlon (grid index)')
    ax.set_ylabel('rlat (grid index)')
    if add_title:
        ax.set_title(
            f'Zone transitions {year - 1}\u2192{year}\n'
            f'(threshold: {melt_threshold} mm w.e.)'
        )
    return ax


# =============================================================================
# Multi-year panel with shared horizontal colorbar
# =============================================================================

def plot_zone_transition_panel(
    years,
    melt_threshold=10.0,
    show_unchanged=True,
    ncols=3,
    figsize=None,
    title=None,
):
    """
    Multi-year grid of zone-transition maps with a single shared horizontal
    colorbar at the bottom.

    Each element of `years` is compared against `year - 1`.  Subplot titles
    show only the year pair (e.g. "2019→2020").

    Parameters
    ----------
    years : list of int
        Each year is compared against year-1.
    melt_threshold : float
        Annual melt threshold separating dry from percolation (mm w.e.).
    show_unchanged : bool
        If True (default), show all four categories.
        If False, show only the two transition categories (3 and 4).
    ncols : int
        Number of subplot columns.  Default 3.
    figsize : tuple, optional
        Figure size.  Auto-sized if None.
    title : str, optional
        Overall figure suptitle.

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : (n_panels,) ndarray of matplotlib.axes.Axes
    """
    n            = len(years)
    nrows        = int(np.ceil(n / ncols))
    ncols_actual = min(n, ncols)

    if figsize is None:
        figsize = (ncols_actual * 4, nrows * 5 + 1.2)

    from matplotlib.gridspec import GridSpec as _GS
    top_val = 0.90 if title else 0.97
    fig     = plt.figure(figsize=figsize)
    _gs     = _GS(nrows, ncols_actual, figure=fig,
                  left=0.02, right=0.98, bottom=0.10, top=top_val,
                  wspace=0.05, hspace=0.15)
    axes_flat = [fig.add_subplot(_gs[i // ncols_actual, i % ncols_actual])
                 for i in range(n)]

    for i, year in enumerate(years):
        ax = axes_flat[i]
        plot_zone_transitions(
            year,
            melt_threshold=melt_threshold,
            ax=ax,
            show_unchanged=show_unchanged,
            show_colorbar=False,
            add_title=False,
        )
        ax.set_title(f'{year - 1}\u2192{year}', fontsize=10)
        ax.set_xlabel('')
        ax.set_ylabel('')

    # Shared horizontal colorbar using a ScalarMappable so we don't need
    # a reference to any individual subplot's image artist.
    if show_unchanged:
        colors = ['#c8c8c8', '#42628F', '#2ca02c', '#d62728']
        bounds = [0.5, 1.5, 2.5, 3.5, 4.5]
        labels = [
            'Dry (unchanged)',
            'Percolation (unchanged)',
            'Dry \u2192 Percolation',
            'Percolation \u2192 Dry',
        ]
        ticks = [1, 2, 3, 4]
    else:
        colors = ['#2ca02c', '#d62728']
        bounds = [2.5, 3.5, 4.5]
        labels = ['Dry \u2192 Percolation', 'Percolation \u2192 Dry']
        ticks = [3, 4]

    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import BoundaryNorm as _BN
    sm = ScalarMappable(cmap=ListedColormap(colors), norm=_BN(bounds, len(colors)))
    sm.set_array([])

    # Colorbar in a fixed strip below the GridSpec rows
    cbar_ax = fig.add_axes([0.15, 0.03, 0.70, 0.022])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(labels, fontsize=9)

    if title:
        fig.suptitle(title, fontsize=12)

    return fig, axes_flat
