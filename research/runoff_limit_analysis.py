"""
Runoff limit analysis for IMAU-FDM output.

The "runoff limit" (RL) in a given year is the maximum topographic elevation
of any Greenland ice-sheet cell that produces annual runoff >= a threshold.
Above the RL, essentially all melt refreezes in the firn.

Public interface
----------------
compute_runoff_limit_timeseries(runoff_threshold=1.0)
    Return (years, rl_elevation) for each simulation year.

map_runoff_limit(year_range=None, ...)
    Map of mean annual runoff with the mean RL elevation contoured.

map_melt_above_limit_ice_lens(year_range=None, ice_lens_file=None, ...)
    Melt at cells above the RL where ice lenses are present.

map_ice_lens_summary(ice_lens_file=None, year_range=None, ...)
    Two-panel map: ice lens frequency and mean thickness.

plot_lens_vs_melt(ice_lens_file=None, year_range=None, ...)
    Per-year scatter of mean melt vs mean ice lens thickness above the RL.
"""

import sys
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

sys.path.insert(0, str(Path(__file__).parent))
import config


# =============================================================================
# Internal helpers
# =============================================================================

def _load_grid():
    """Return (topo, ice) — both (rlat, rlon) arrays from the masks file."""
    masks   = xr.open_dataset(config.MASKS_FILE)
    topo    = masks['Topography'].squeeze().values.astype(np.float32)  # m a.s.l.
    icemask = masks['Icemask_GR'].squeeze().values
    masks.close()
    return topo, (icemask > 0)


def _annual_totals(data, time_frac):
    """
    Aggregate (n_time, rlat, rlon) 10-day data to per-year totals.

    Returns
    -------
    years  : (n_years,) int
    annual : (n_years, rlat, rlon) float32
    """
    yr = np.floor(time_frac).astype(int)
    years = np.unique(yr)
    annual = np.stack([data[yr == y].sum(axis=0) for y in years]).astype(np.float32)
    return years, annual


def _annual_means(data, time_frac):
    """Aggregate (n_time, rlat, rlon) to per-year means."""
    yr = np.floor(time_frac).astype(int)
    years = np.unique(yr)
    annual = np.stack([data[yr == y].mean(axis=0) for y in years]).astype(np.float32)
    return years, annual


def _load_annual_variable(variable, nc_file=None, stat='sum'):
    """Load a 10-day gridded variable and return (years, annual_values).

    stat : 'sum' (default) or 'mean' — whether to aggregate timesteps within
           each year by summing or averaging.
    """
    if nc_file is None:
        nc_file = config.VARIABLE_FILES[variable]
    ds   = xr.open_dataset(nc_file, decode_times=False)
    data = ds[variable].values
    time = ds['time'].values
    ds.close()
    agg = _annual_totals if stat == 'sum' else _annual_means
    return agg(data, time)


def _rl_per_year(annual_runoff, topo, ice, runoff_threshold):
    """Return per-year runoff limit elevation array (NaN if no runoff that year)."""
    rl = np.full(len(annual_runoff), np.nan, dtype=np.float32)
    for i, runoff_yr in enumerate(annual_runoff):
        mask = (runoff_yr >= runoff_threshold) & ice
        if mask.any():
            rl[i] = np.nanmax(topo[mask])
    return rl


# =============================================================================
# a) Runoff limit time series
# =============================================================================

def compute_runoff_limit_timeseries(runoff_threshold=10.0, nc_file=None):
    """
    Compute the annual runoff limit elevation time series.

    The runoff limit for year Y is the maximum topographic elevation of any
    Greenland ice-sheet cell with annual total runoff >= ``runoff_threshold``.

    Parameters
    ----------
    runoff_threshold : float
        Minimum annual runoff (mm w.e.) to count a cell as a runoff producer.
        Default: 10.0 mm w.e.
    nc_file : Path or str, optional
        Runoff NetCDF.  Default: config.VARIABLE_FILES['Runoff'].

    Returns
    -------
    years        : (n_years,) int array
    rl_elevation : (n_years,) float array — runoff limit (m a.s.l.), NaN if
                   no runoff produced that year
    """
    topo, ice = _load_grid()
    years, annual_runoff = _load_annual_variable('Runoff', nc_file)
    rl = _rl_per_year(annual_runoff, topo, ice, runoff_threshold)
    return years, rl


# =============================================================================
# a) Map: runoff extent + RL contour
# =============================================================================

def map_runoff_limit(year_range=None, runoff_threshold=10.0,
                     ax=None, figsize=(7, 9), cmap='YlOrRd'):
    """
    Map of mean annual runoff with the mean runoff limit elevation contoured.

    Parameters
    ----------
    year_range : (int, int), optional
        Inclusive year range to average over.  Default: all years.
    runoff_threshold : float
        Runoff threshold for RL definition (mm w.e./yr).  Default: 10.0.
    ax : matplotlib.axes.Axes, optional
    figsize : tuple
    cmap : str

    Returns
    -------
    (fig, ax), mean_rl_elevation : float
    """
    from matplotlib.colors import LogNorm

    topo, ice = _load_grid()
    years, annual_runoff = _load_annual_variable('Runoff')

    if year_range is not None:
        sel = (years >= year_range[0]) & (years <= year_range[1])
        years, annual_runoff = years[sel], annual_runoff[sel]

    mean_runoff = annual_runoff.mean(axis=0)
    rl_arr = _rl_per_year(annual_runoff, topo, ice, runoff_threshold)
    mean_rl = float(np.nanmean(rl_arr))

    # Only display ice-sheet cells with runoff > 0
    plot_data = np.where(ice & (mean_runoff > 0), mean_runoff, np.nan)

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    im = ax.imshow(plot_data, origin='lower',
                   norm=LogNorm(vmin=1, vmax=5000), cmap=cmap,
                   interpolation='nearest')

    # Coastline
    coast = np.where(np.isfinite(topo) & (topo > 0), 1.0, 0.0)
    ax.contour(coast, levels=[0.5], colors='k', linewidths=0.7, alpha=0.85)

    # RL contour: boundary around cells that produce mean annual runoff > 0
    # (contouring topography at mean_rl elevation is misleading because the
    # topographic contour extends to interior regions where no runoff occurs)
    runoff_mask = np.where(ice, (mean_runoff > 0).astype(np.float32), np.nan)
    ax.contour(runoff_mask, levels=[0.5], colors=['navy'],
               linewidths=1.8, linestyles='--')

    if standalone:
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        cb.set_label('Mean annual runoff (mm w.e.)')

    yr_str = f'{years[0]}–{years[-1]}' if len(years) > 1 else str(years[0])
    ax.set_title(f'Runoff & runoff limit  ({yr_str})\n'
                 f'Mean RL = {mean_rl:.0f} m a.s.l.', fontsize=10)
    ax.set_aspect('equal')
    ax.axis('off')

    return (fig, ax), mean_rl


# =============================================================================
# b) Map: melt above RL in ice-lens cells
# =============================================================================

def map_melt_above_limit_ice_lens(year_range=None, year=None, runoff_threshold=10.0,
                                   ice_lens_file=None, topo_contour_interval=500,
                                   ax=None, figsize=(7, 9)):
    """
    Map of mean annual surfmelt above the runoff limit where ice lenses exist.

    Only melt that co-occurs with an ice lens **in the same 10-day timestep** is
    included.  For each year the annual RL elevation is computed from runoff, then
    only cells above that RL and on the ice sheet are considered.

    Parameters
    ----------
    year_range : (int, int), optional
        Inclusive year range.  Default: all years.
    year : int, optional
        Single year — shorthand for year_range=(year, year).
        If both year and year_range are given, year takes precedence.
    runoff_threshold : float
        Minimum annual runoff (mm w.e.) defining the RL.  Default: 10.0.
    ice_lens_file : Path or str, optional
        Default: config.ICE_LENS_FILE.
    topo_contour_interval : int or None
        Draw topographic contour lines every this many metres.  Default: 500 m.
        Set to None to suppress contours.
    ax : matplotlib.axes.Axes, optional
    figsize : tuple

    Returns
    -------
    fig, ax
    """
    if year is not None:
        year_range = (year, year)

    if ice_lens_file is None:
        ice_lens_file = config.ICE_LENS_FILE
    ice_lens_file = Path(ice_lens_file)
    if not ice_lens_file.exists():
        raise FileNotFoundError(
            f'Ice lens file not found: {ice_lens_file}\n'
            'Run create_gridded_ice_lens_file() first.')

    topo, ice = _load_grid()
    years_r, annual_runoff = _load_annual_variable('Runoff')

    # Read only the 1-D time axes — no full data arrays in memory yet
    melt_ds   = xr.open_dataset(config.VARIABLE_FILES['surfmelt'], decode_times=False)
    melt_time = melt_ds['time'].values   # (n_time_m,) float — tiny
    il_ds     = xr.open_dataset(ice_lens_file, decode_times=False)
    il_time   = il_ds['time'].values     # (n_time_il,) float — tiny

    # Match each melt timestep to the nearest ice-lens timestep (tolerance 0.02 yr ≈ 7 days).
    # Both files use the same 10-day cadence so matches should be exact.
    TOL = 0.02
    melt_to_il = np.full(len(melt_time), -1, dtype=int)
    for ti, t in enumerate(melt_time):
        diffs = np.abs(il_time - t)
        idx   = int(np.argmin(diffs))
        if diffs[idx] < TOL:
            melt_to_il[ti] = idx

    melt_yr_arr = np.floor(melt_time).astype(int)
    common = sorted(set(years_r.tolist()) & set(melt_yr_arr.tolist()))
    if year_range is not None:
        common = [y for y in common if year_range[0] <= y <= year_range[1]]
    if not common:
        melt_ds.close(); il_ds.close()
        raise ValueError('No years in common between runoff, melt, and ice lens files.')

    accum        = np.zeros_like(topo, dtype=np.float64)
    n_years_used = 0

    for y in common:
        ri      = np.searchsorted(years_r, y)
        rl_mask = (annual_runoff[ri] >= runoff_threshold) & ice
        if not rl_mask.any():
            continue
        rl_y       = float(np.nanmax(topo[rl_mask]))
        above_mask = ice & (topo > rl_y)

        # Indices for this year's melt timesteps and their matched ice-lens indices
        t_sel   = np.where(melt_yr_arr == y)[0]
        il_idxs = melt_to_il[t_sel]
        valid   = t_sel[il_idxs >= 0]
        valid_il = il_idxs[il_idxs >= 0]
        if len(valid) == 0:
            continue

        # Load only this year's slices from disk — avoids holding full arrays in RAM
        melt_yr = melt_ds['surfmelt'].isel(time=valid.tolist()).values     # (n_steps, rlat, rlon)
        il_yr   = il_ds['has_ice_lens'].isel(time=valid_il.tolist()).values # (n_steps, rlat, rlon)

        # Accumulate melt only where ice lens is present in the same step
        for i in range(len(valid)):
            combined = above_mask & (il_yr[i] > 0)
            accum[combined] += melt_yr[i][combined]

        n_years_used += 1

    melt_ds.close()
    il_ds.close()

    if n_years_used == 0:
        raise ValueError('No valid years found.')

    mean_melt = (accum / n_years_used).astype(np.float32)
    mean_melt = np.where(ice & (mean_melt > 0), mean_melt, np.nan)
    vmax = float(np.nanpercentile(mean_melt, 98)) if np.any(np.isfinite(mean_melt)) else 500.0

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    im = ax.imshow(mean_melt, origin='lower', cmap='plasma',
                   vmin=0, vmax=vmax, interpolation='nearest')

    # Coastline
    coast = np.where(np.isfinite(topo) & (topo > 0), 1.0, 0.0)
    ax.contour(coast, levels=[0.5], colors='k', linewidths=0.7, alpha=0.85)

    # Topographic contours
    if topo_contour_interval is not None:
        topo_min = int(np.nanmin(topo[topo > 0]) // topo_contour_interval) * topo_contour_interval
        topo_max = int(np.nanmax(topo) // topo_contour_interval + 1) * topo_contour_interval
        topo_levels = np.arange(topo_min, topo_max + 1, topo_contour_interval)
        cs = ax.contour(topo, levels=topo_levels, colors='k', linewidths=0.4, alpha=0.5)
        ax.clabel(cs, fmt='%d m', fontsize=6, inline=True, inline_spacing=2)

    if standalone:
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        single_yr = len(common) == 1
        cb.set_label('Annual surfmelt (mm w.e.)' if single_yr
                     else 'Mean annual surfmelt (mm w.e.)')

    yr_str = str(common[0]) if len(common) == 1 else f'{common[0]}–{common[-1]}'
    ax.set_title(f'Melt above runoff limit — ice-lens cells (same timestep)\n{yr_str}',
                 fontsize=10)
    ax.set_aspect('equal')
    ax.axis('off')

    return fig, ax


# =============================================================================
# c) Map: ice lens presence and thickness
# =============================================================================

def map_ice_lens_summary(ice_lens_file=None, year_range=None,
                         topo_contour_interval=500, figsize=(14, 8)):
    """
    Two-panel map: ice lens frequency (fraction of 10-day steps with a lens)
    and mean ice lens thickness when a lens is present.

    Parameters
    ----------
    ice_lens_file : Path or str, optional
        Default: config.ICE_LENS_FILE.
    year_range : (int, int), optional
    topo_contour_interval : int or None
        Draw topographic contour lines every this many metres.  Default: 500 m.
        Set to None to suppress contours.
    figsize : tuple

    Returns
    -------
    fig, (ax_freq, ax_thick)
    """
    if ice_lens_file is None:
        ice_lens_file = config.ICE_LENS_FILE
    ice_lens_file = Path(ice_lens_file)
    if not ice_lens_file.exists():
        raise FileNotFoundError(
            f'Ice lens file not found: {ice_lens_file}\n'
            'Run create_gridded_ice_lens_file() first.')

    topo, ice = _load_grid()

    il_ds    = xr.open_dataset(ice_lens_file, decode_times=False)
    il_has   = il_ds['has_ice_lens'].values        # keep as int8: ~4× smaller than float32
    il_thick = il_ds['ice_lens_thickness'].values  # float32
    il_time  = il_ds['time'].values
    il_ds.close()

    if year_range is not None:
        yr = np.floor(il_time).astype(int)
        sel = (yr >= year_range[0]) & (yr <= year_range[1])
        il_has, il_thick = il_has[sel], il_thick[sel]

    # Frequency: sum bool array (1 byte each) then divide — avoids a float32 3D copy
    freq = (il_has > 0).sum(axis=0, dtype=np.float32) / len(il_has)
    freq = np.where(ice, freq, np.nan)

    # Mean thickness when lens present — avoid a full 3D copy by zeroing in-place
    # then dividing by the count of lens-present steps.
    n_lens = (il_has > 0).sum(axis=0).astype(np.float32)   # (rlat, rlon)
    il_thick[il_has == 0] = 0.0                              # in-place: no extra copy
    sum_thick = il_thick.sum(axis=0).astype(np.float32)
    mean_thick = np.where(n_lens > 0, sum_thick / n_lens, np.nan)
    mean_thick = np.where(ice, mean_thick, np.nan).astype(np.float32)
    del n_lens, sum_thick

    vmax_thick = float(np.nanpercentile(mean_thick, 98)) if np.any(np.isfinite(mean_thick)) else 1.0

    fig = plt.figure(figsize=figsize)
    gs  = GridSpec(1, 2, figure=fig,
                   left=0.03, right=0.95, bottom=0.05, top=0.88, wspace=0.12)
    ax_freq  = fig.add_subplot(gs[0, 0])
    ax_thick = fig.add_subplot(gs[0, 1])

    coast = np.where(np.isfinite(topo) & (topo > 0), 1.0, 0.0)

    # Topographic contours (computed once, applied to both panels)
    if topo_contour_interval is not None:
        topo_min = int(np.nanmin(topo[topo > 0]) // topo_contour_interval) * topo_contour_interval
        topo_max = int(np.nanmax(topo) // topo_contour_interval + 1) * topo_contour_interval
        topo_levels = np.arange(topo_min, topo_max + 1, topo_contour_interval)

    im1 = ax_freq.imshow(freq, origin='lower', cmap='Blues',
                         vmin=0, vmax=1, interpolation='nearest')
    fig.colorbar(im1, ax=ax_freq, fraction=0.04, pad=0.02,
                 label='Fraction of 10-day steps with ice lens')
    ax_freq.contour(coast, levels=[0.5], colors='k', linewidths=0.7, alpha=0.85)
    if topo_contour_interval is not None:
        cs1 = ax_freq.contour(topo, levels=topo_levels, colors='k', linewidths=0.4, alpha=0.5)
        ax_freq.clabel(cs1, fmt='%d m', fontsize=6, inline=True, inline_spacing=2)
    ax_freq.set_title('Ice lens frequency', fontsize=11)
    ax_freq.set_aspect('equal')
    ax_freq.axis('off')

    im2 = ax_thick.imshow(mean_thick, origin='lower', cmap='hot_r',
                           vmin=0, vmax=vmax_thick, interpolation='nearest')
    fig.colorbar(im2, ax=ax_thick, fraction=0.04, pad=0.02,
                 label='Mean ice lens thickness when present (m)')
    ax_thick.contour(coast, levels=[0.5], colors='k', linewidths=0.7, alpha=0.85)
    if topo_contour_interval is not None:
        cs2 = ax_thick.contour(topo, levels=topo_levels, colors='k', linewidths=0.4, alpha=0.5)
        ax_thick.clabel(cs2, fmt='%d m', fontsize=6, inline=True, inline_spacing=2)
    ax_thick.set_title('Mean ice lens thickness', fontsize=11)
    ax_thick.set_aspect('equal')
    ax_thick.axis('off')

    yr_str = (f' ({year_range[0]}–{year_range[1]})' if year_range else '')
    fig.suptitle(f'Ice lens diagnostics{yr_str}', fontsize=13)

    return fig, (ax_freq, ax_thick)


# =============================================================================
# d) Scatter: annual mean ice lens thickness vs annual mean melt above RL
# =============================================================================

def plot_lens_vs_melt(ice_lens_file=None, year_range=None,
                      runoff_threshold=10.0, melt_stat='sum',
                      ax=None, figsize=(7, 5)):
    """
    Per-year scatter: mean surfmelt vs mean ice lens thickness above the RL.

    For each year, cells are restricted to those above that year's runoff limit
    on the Greenland ice sheet.  Each scatter point is one calendar year,
    coloured by year.  A linear trend line is added when n >= 3 years.

    Parameters
    ----------
    ice_lens_file : Path or str, optional
        Default: config.ICE_LENS_FILE.
    year_range : (int, int), optional
    runoff_threshold : float
    melt_stat : {'sum', 'mean'}
        Whether to show the annual total ('sum', mm w.e./yr, default) or the
        mean per 10-day timestep ('mean', mm w.e./step) on the x-axis.
    ax : matplotlib.axes.Axes, optional
    figsize : tuple

    Returns
    -------
    fig, ax
    """
    if ice_lens_file is None:
        ice_lens_file = config.ICE_LENS_FILE
    ice_lens_file = Path(ice_lens_file)
    if not ice_lens_file.exists():
        raise FileNotFoundError(
            f'Ice lens file not found: {ice_lens_file}\n'
            'Run create_gridded_ice_lens_file() first.')

    if melt_stat not in ('sum', 'mean'):
        raise ValueError(f"melt_stat must be 'sum' or 'mean', got '{melt_stat}'")

    topo, ice = _load_grid()
    years_r, annual_runoff = _load_annual_variable('Runoff')
    years_m, annual_melt   = _load_annual_variable('surfmelt', stat=melt_stat)

    # Ice lens thickness: annual mean of (thickness when lens present)
    il_ds    = xr.open_dataset(ice_lens_file, decode_times=False)
    il_has   = il_ds['has_ice_lens'].values    # (n_time, rlat, rlon) int8
    il_thick = il_ds['ice_lens_thickness'].values  # (n_time, rlat, rlon) float32
    il_time  = il_ds['time'].values
    il_ds.close()

    il_yr   = np.floor(il_time).astype(int)
    il_yrs  = np.unique(il_yr)
    # Annual mean thickness (over steps where a lens exists)
    with np.errstate(all='ignore'):
        il_thick_ann = np.stack([
            np.nanmean(np.where(il_has[il_yr == y] > 0,
                                il_thick[il_yr == y], np.nan), axis=0)
            for y in il_yrs
        ]).astype(np.float32)   # (n_years, rlat, rlon)

    common = sorted(set(years_r.tolist()) & set(years_m.tolist()) & set(il_yrs.tolist()))
    if year_range is not None:
        common = [y for y in common if year_range[0] <= y <= year_range[1]]
    if not common:
        raise ValueError('No years in common between runoff, melt, and ice lens files.')

    pt_years, pt_melt, pt_thick = [], [], []
    for y in common:
        ri = np.searchsorted(years_r, y)
        mi = np.searchsorted(years_m, y)
        ii = np.searchsorted(il_yrs,  y)

        rl_mask = (annual_runoff[ri] >= runoff_threshold) & ice
        if not rl_mask.any():
            continue
        rl_y = float(np.nanmax(topo[rl_mask]))

        above = ice & (topo > rl_y)
        if not above.any():
            continue

        pt_years.append(y)
        pt_melt.append(float(np.nanmean(annual_melt[mi][above])))
        pt_thick.append(float(np.nanmean(il_thick_ann[ii][above])))

    pt_years = np.array(pt_years)
    pt_melt  = np.array(pt_melt,  dtype=float)
    pt_thick = np.array(pt_thick, dtype=float)

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    sc = ax.scatter(pt_melt, pt_thick, c=pt_years, cmap='viridis',
                    edgecolors='k', linewidths=0.4, s=55, zorder=3)
    if standalone:
        cb = fig.colorbar(sc, ax=ax)
        cb.set_label('Year')

    # if len(pt_melt) >= 3:
    #     p = np.polyfit(pt_melt, pt_thick, 1)
    #     xfit = np.linspace(pt_melt.min(), pt_melt.max(), 100)
    #     ax.plot(xfit, np.polyval(p, xfit), 'r--', linewidth=1.3,
    #             label=f'slope={p[0]:.3f} m/(mm w.e.)')
    #     ax.legend(fontsize=9)

    _melt_xlabel = ('Mean annual surfmelt above RL (mm w.e./yr)'
                    if melt_stat == 'sum' else
                    'Mean surfmelt per 10-day step above RL (mm w.e.)')
    ax.set_xlabel(_melt_xlabel, fontsize=10)
    ax.set_ylabel('Mean ice lens thickness above RL (m)', fontsize=10)
    ax.set_title('Ice lens thickness vs melt above runoff limit', fontsize=11)
    ax.grid(True, alpha=0.3)

    return fig, ax
