# Physics

IMAU-FDM is a one-dimensional firn densification model. Each column evolves
independently through:

1. **Accumulation** — new snow is added at the surface with density ρ₀ (initial
   surface snow density, parameterised from temperature and wind).
2. **Firn compaction** — subsurface layers compact via two rate equations
   (dry and wet) calibrated against firn-core observations (MO fitting).
3. **Melt and refreeze** — surface melt is computed from the surface energy
   balance; meltwater percolates downward and can refreeze into ice lenses
   or recharge the firn pore space.
4. **Liquid water transport** — water moves vertically through available pore
   space; capillary retention is accounted for.

## Firn compaction

Compaction follows a two-stage Herron-Langway-type scheme. The densification
rate for a layer is:

```
dρ/dt = c · A · exp(−E / (R T)) · (ρ_ice − ρ)
```

where:
- `c` — empirical coefficient (different for Stage 1, ρ < 550 kg m⁻³, and Stage 2)
- `A` — annual accumulation rate (m w.e. yr⁻¹)
- `E` — activation energy
- `R` — gas constant
- `T` — temperature (K)
- `ρ_ice = 917 kg m⁻³`

The coefficients `c` are calibrated by the MO-fitting procedure (see
[MO Fitting](../mo_fitting/overview)).

## Output layers

| Output type | Layers | Thickness | Interval | Depth |
|-------------|--------|-----------|----------|-------|
| 1D | — | — | Daily | Surface scalars |
| 2D | 3000 | ~4 cm → m | 30-day | ~122 m full column |
| 2Ddetail | 500 | 4 cm | 10-day | Top 20 m |

## Key variables

| Variable | Description | Units |
|----------|-------------|-------|
| `h_surf` | Surface height anomaly | m |
| `FirnAir` | Firn air content (pore space above ρ = 917) | m |
| `Runoff` | Meltwater that leaves the column | mm w.e. |
| `surfmelt` | Surface melt rate | mm w.e. |
| `refreeze` | Refreezing in firn | mm w.e. |
| `TotLwc` | Total liquid water content in column | mm |
| `dens` | Density profile | kg m⁻³ |
| `temp` | Temperature profile | K |
| `z830` | Depth of firn–ice transition (ρ = 830) | m |
| `T10m` | Temperature at 10 m depth | K |

## Coordinate system

Output grids use a **rotated-pole projection** (native RACMO grid). Absolute
coordinates are provided as auxiliary `lat`/`lon` variables. The physical
x/y coordinates are in **EPSG:3413** (WGS 84 / NSIDC Sea Ice Polar
Stereographic North, lon₀ = −45°, lat_ts = 70°).
