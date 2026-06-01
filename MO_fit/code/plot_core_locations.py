import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import xarray as xr
from pyproj import Transformer
from datetime import datetime
import os

MASK_FILE = os.path.join(os.path.dirname(__file__), "../../../reference/FGRN055/FGRN055_Masks.nc")
DATA_DIR  = os.path.join(os.path.dirname(__file__), "../data/")


def _load_merged(data_dir):
    year = datetime.now().year
    path = os.path.join(data_dir, f"merged/MERGED_CORE_LIST_greenland_{year}.csv")
    if not os.path.exists(path):
        # fall back to most recent available year
        candidates = sorted([
            f for f in os.listdir(os.path.join(data_dir, "merged"))
            if f.startswith("MERGED_CORE_LIST_greenland_") and f.endswith(".csv")
        ])
        if not candidates:
            raise FileNotFoundError(f"No merged core list found in {data_dir}/merged/")
        path = os.path.join(data_dir, "merged", candidates[-1])
    return pd.read_csv(path)


def plot_core_locations(merged_df=None, mask_file=None, data_dir=None, save_path=None):
    """
    Two-panel map of firn core locations on Greenland:
      left  — colored by depth to 550 kg/m³
      right — colored by depth to 830 kg/m³
    PKM standardized cores: squares; SUMUP cores: circles.
    500 m topography contours shown on ice sheet.
    """
    if data_dir is None:
        data_dir = DATA_DIR
    if mask_file is None:
        mask_file = MASK_FILE
    if merged_df is None:
        merged_df = _load_merged(data_dir)

    # --- load mask ---
    ds = xr.open_dataset(mask_file)
    X    = ds["X"].values        # (rlat, rlon), km
    Y    = ds["Y"].values        # (rlat, rlon), km
    topo   = ds["Topography"].values
    ice    = ds["Icemask_GR"].values   # Greenland ice only
    lsm    = ds["LSM_GR"].values       # Greenland land only

    # --- project core lat/lon → EPSG:3413 km ---
    tr = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
    x_m, y_m = tr.transform(merged_df["longitude"].values, merged_df["latitude"].values)
    xc = x_m / 1000.0
    yc = y_m / 1000.0

    pkm_mask   = (merged_df["source"] == "PKM standardized cores").values
    sumup_mask = (merged_df["source"] == "SUMUP 2024").values

    topo_levels = np.arange(500, 3500, 500)

    depth_vars = ["depth_to_550", "depth_to_830"]
    panel_titles = ["Depth to 550 kg m$^{-3}$ (m)", "Depth to 830 kg m$^{-3}$ (m)"]
    cmap = "plasma"

    fig, axes = plt.subplots(1, 2, figsize=(13,11))

    for ax, depth_var, title in zip(axes, depth_vars, panel_titles):

        ax.set_facecolor("white")

        # land (non-ice Greenland)
        land = np.where((lsm > 0) & (ice < 0.5), 1.0, np.nan)
        ax.pcolormesh(X, Y, land, cmap="Greys", vmin=0, vmax=3,
                      shading="auto", rasterized=True)

        # ice sheet
        ice_show = np.where(ice >= 0.5, 1.0, np.nan)
        ax.pcolormesh(X, Y, ice_show, cmap="Blues", vmin=0.5, vmax=2,
                      shading="auto", alpha=0.4, rasterized=True)

        # topography contours (ice only)
        topo_ice = np.where(ice >= 0.5, topo, np.nan)
        cs = ax.contour(X, Y, topo_ice, levels=topo_levels,
                        colors="grey", linewidths=0.6, alpha=0.8)
        ax.clabel(cs, fmt="%d m", fontsize=6, inline=True)

        # core depths
        depths = merged_df[depth_var].values
        valid  = ~np.isnan(depths)
        vmin   = np.nanmin(depths)
        vmax   = np.nanmax(depths)
        norm   = mcolors.Normalize(vmin=vmin, vmax=vmax)

        # coloured symbols for valid values (SUMUP first so PKM sits on top)
        sc = None
        for mask, marker, size in [(sumup_mask, "o", 35), (pkm_mask, "s", 50)]:
            has_data = mask & valid
            if has_data.any():
                sc = ax.scatter(xc[has_data], yc[has_data],
                                c=depths[has_data], marker=marker,
                                s=size, cmap=cmap, norm=norm, zorder=6,
                                edgecolors="k", linewidths=0.5)

        if sc is not None:
            cbar = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.03, shrink=1)
            cbar.set_label("Depth (m)", fontsize=10)

        ax.set_xlim(-750, 1000)
        ax.set_ylim(-3500, -500)
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=12, pad=8)
        ax.set_xlabel("X (km)", fontsize=10)
        ax.set_ylabel("Y (km)", fontsize=10)
        ax.tick_params(labelsize=8)
        ax.grid(False)

    legend_elements = [
        Line2D([0], [0], marker="s", color="w", markerfacecolor="grey",
               markersize=9, markeredgecolor="k", label="PKM standardized"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="grey",
               markersize=9, markeredgecolor="k", label="SUMUP 2024"),
    ]
    for ax in axes:
        ax.legend(handles=legend_elements, loc="lower right", fontsize=9, framealpha=0.8)

    fig.suptitle("Greenland Firn Core Locations", fontsize=20, y=0.93)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved to {save_path}")

    return fig, axes
