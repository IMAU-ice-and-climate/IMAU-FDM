"""
Create FGRN055 pointlists for the K-transect study area.

Outputs (in run_dir/pointlists/):
  IN_ll_FGRN055_k_transect.txt      — all grid points within the study area polygon
  IN_ll_FGRN055_k_transect_line.txt — continuous west→east transect near 67°N

Both files contain only the 1-based file index (row number in the master pointlist),
one integer per line — the format expected by readpointlist.x.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd

# ---------------------------------------------------------------------------
RESEARCH_DIR   = Path('/home/nld4814/perm/code/IMAU-FDM/research')
POINTLIST_FILE = Path('/home/nld4814/perm/code/IMAU-FDM/reference/FGRN055/IN_ll_FGRN055.txt')
K_TRANSECT_SHP = RESEARCH_DIR / 'shapefiles' / 'k_transect_study_area.shp'
OUTPUT_DIR     = Path('/home/nld4814/perm/code/IMAU-FDM/rundir/pointlists')
TARGET_LAT     = 67.0   # classic K-transect latitude (°N)
# ---------------------------------------------------------------------------

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# --- Load master pointlist --------------------------------------------------
grid = pd.read_csv(POINTLIST_FILE, header=None, sep=',',
                   names=['lon', 'lat', 'c2', 'c3', 'c4', 'rlat_idx', 'rlon_idx'])
grid['file_num'] = grid.index + 1    # 1-based: file numbers start at 1

# --- Spatial filter: points within the study area polygon ------------------
k_area = gpd.read_file(K_TRANSECT_SHP)   # EPSG:4326
gdf    = gpd.GeoDataFrame(
    grid,
    geometry=gpd.points_from_xy(grid['lon'], grid['lat']),
    crs='EPSG:4326',
)
in_area = gpd.sjoin(gdf, k_area[['geometry']], how='inner', predicate='within')
in_area = in_area.drop(columns=['geometry', 'index_right'])

print(f'Points within K-transect study area: {len(in_area)}')

# Save full area pointlist (file indices only)
in_area['file_num'].to_csv(OUTPUT_DIR / 'pointlist_k_transect_area.txt',
                           header=False, index=False)
print(f'Saved: {OUTPUT_DIR}/pointlist_k_transect_area.txt')

# --- Continuous transect: single rlat row nearest TARGET_LAT ---------------
# For each rlat_idx, compute the mean latitude of its points in the area.
by_rlat = (in_area.groupby('rlat_idx')
           .agg(n_pts=('rlon_idx', 'count'), mean_lat=('lat', 'mean'))
           .reset_index())

# Filter to rows within 0.5° of target, pick the one with the most coverage
near = by_rlat[abs(by_rlat['mean_lat'] - TARGET_LAT) < 0.5].sort_values('n_pts', ascending=False)
if near.empty:
    # Relax tolerance
    near = by_rlat.sort_values(by=lambda df: abs(df['mean_lat'] - TARGET_LAT))
best_rlat = int(near.iloc[0]['rlat_idx'])
print(f'Transect rlat_idx = {best_rlat}  (mean lat ≈ {near.iloc[0]["mean_lat"]:.2f}°N, '
      f'{int(near.iloc[0]["n_pts"])} pts)')

# Select and sort west→east (increasing rlon_idx)
row = in_area[in_area['rlat_idx'] == best_rlat].sort_values('rlon_idx').reset_index(drop=True)

# Keep only the longest continuous segment (no gap > 1 in rlon_idx)
rlon_vals = row['rlon_idx'].values
gaps      = np.where(np.diff(rlon_vals) > 1)[0]
if gaps.size:
    # Find the longest contiguous run
    boundaries = np.concatenate([[- 1], gaps, [len(rlon_vals) - 1]])
    lengths    = np.diff(boundaries)
    best_seg   = int(np.argmax(lengths))
    start      = boundaries[best_seg] + 1
    end        = boundaries[best_seg + 1] + 1
    row        = row.iloc[start:end].reset_index(drop=True)

print(f'Continuous transect: {len(row)} points  '
      f'lon {row["lon"].min():.2f}° → {row["lon"].max():.2f}°  '
      f'lat {row["lat"].min():.3f}°–{row["lat"].max():.3f}°')

row['file_num'].to_csv(OUTPUT_DIR / 'pointlist_k_transect_line.txt',
                       header=False, index=False)
print(f'Saved: {OUTPUT_DIR}/pointlist_k_transect_line.txt')
