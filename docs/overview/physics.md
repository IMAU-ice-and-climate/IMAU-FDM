# Physics

IMAU-FDM is a one-dimensional firn densification model. Each column evolves
independently through:

1. **Accumulation** — new snow is added at the surface with density ρ₀ (initial
   surface snow density).
2. **Firn compaction** — subsurface layers compact via density-dependent rate equations
 calibrated against firn-core observations (see [MO fitting](../mo_fitting/overview)).
3. **Melt and refreeze** — surface melt is computed from the surface energy
   balance; meltwater percolates downward and can refreeze, creating denser firn or ice lenses, or becomes runoff if it reaaches the bottom of the firn column
4. **Liquid water transport** — water moves vertically and instantaneoulsy through
   available pore space

## Firn compaction

Compaction follows the two-stage Arthern (2010) / Ligtenberg (2011)
semi-empirical scheme (with the domain-specific MO tuning). The densification
rate for a layer is:

```
dρ/dt = c · MO_fit · A · g · exp(−Ec / (R T) + Eg/(R T_base)) · (ρ_ice − ρ)
```

where:
- `c` — empirical coefficient (0.07 for Stage 1 (ρ < 550 kg m⁻³), and 0.03 for Stage 2 (ρ > 550 kg m⁻³))
- `MO_fit` - scaling factor based on average accumulation ([MO Fitting](../mo_fitting/overview))
- `A` — average accumulation rate (mm w.e. yr⁻¹)
- `g` - gravitational acceleration (9.81 m s⁻²)
- `Ec` — creep activation energy
- `Eg` - grain growth activation energy
- `R` — gas constant
- `T` — temperature (K) at layer
- `T_base` - temperature at the base of the firn column
- `ρ_ice = 917 kg m⁻³`

## Output layers

| Output type | Layers | Thickness | Interval | Depth |
|-------------|--------|-----------|----------|-------|
| 1D | — | — | Daily* | integrated or surface value |
| 2D | 3000* | variable | 30-day* | full column |
| 2Ddetail | 500* | 4 cm* | 10-day* | top 20 m |

\* defaults, can be adjusted in model_settings

## Coordinate system

Output grids use a **rotated-pole projection** (native RACMO grid). Absolute
coordinates are provided as auxiliary `lat`/`lon` variables. The physical
x/y coordinates are in **EPSG:3413** (WGS 84 / NSIDC Sea Ice Polar
Stereographic North, lon₀ = −45°, lat_ts = 70°).

## Some example variables

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
