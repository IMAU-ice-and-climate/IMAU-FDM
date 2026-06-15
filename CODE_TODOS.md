# Code TODOs (from `TKTKTK` markers)

Tracks the `TKTKTK` placeholders in `source/`. Each entry lists the location, the
in-code note, and a recommended action. Remove the `TKTKTK` marker from the source
(and update this file) when its item is done.

> **Re-scanned 2026-06-15: 10 markers remain.** Of the original 15 (verified
> 2026-06-08), 5 have been cleared — the dead `Find_Grid` subroutine was deleted,
> and 4 stale "is this needed / already set?" markers (`output.f90:474`,
> `output.f90:574`, `model_settings.f90:347`, `:580`) were confirmed non-issues
> and stripped. The remaining 10 below are real-work items; none are bugs.

---

## 1. Make physics choices explicit in `model.toml`

The surface-snow density and thermal-conductivity parameterisations are hardcoded
and branch on `domain` (`FGRN*` vs everything-else). They should be selectable
physics options in `model.toml` so each run records which scheme it used, and so
new domains don't silently fall into the "else" branch.

- **[source/initialise_model.f90:91](source/initialise_model.f90#L91)** — *"make choice of surface snow explicit in model.toml"* (the surface-density branch).
- **[source/initialise_model.f90:108](source/initialise_model.f90#L108)** — Greenland surface density hardcoded as *Fausto et al. 2018* (`Rho0 = 362.1 + 2.78*(T − Tmelt)`).
- **[source/initialise_model.f90:109](source/initialise_model.f90#L109)** — also wants an option for constant surface density (315 kg m⁻³).
- **[source/initialise_model.f90:120](source/initialise_model.f90#L120)** — non-Greenland surface density hardcoded as *Veldhuijsen 2023 (refitting Lenaerts 2012)* (`Rho0 = 82.97 + 0.769*T + 11.67*ff10`).
- **[source/initialise_model.f90:242](source/initialise_model.f90#L242)** — thermal-conductivity scheme in `Init_Temp_Prof` is hardcoded.
- **[source/firn_physics.f90:248](source/firn_physics.f90#L248)** — the same conductivity scheme, duplicated in the `Thermal_Cond` function: *"all these model choices should be set explicitly in model.toml"*.

**Recommended:** add a `[model_physics]`-style block in `model.toml` (e.g.
`surface_density_scheme`, `thermal_conductivity_scheme`) and dispatch on it instead
of on `domain`. **Note the duplication:** the thermal-conductivity formulas appear
in *both* `firn_physics.f90:Thermal_Cond` and `initialise_model.f90:Init_Temp_Prof`
— consolidate into the single `Thermal_Cond` function and call it from both.

---

## 2. Validate `run.toml` fields

- **[source/model_settings.f90:454](source/model_settings.f90#L454)** — `Read_Job`: *"check on domain, forcing, restart_type"*. **Recommended:** validate these against allowed values (`domain ∈ {FGRN055, ANT27}`, `restart_type ∈ {none, spinup, run}`, …) and fail early with a clear message rather than letting a typo propagate.

---

## 3. Constants that should be centralised

- **[source/initialise_model.f90:238](source/initialise_model.f90#L238)** — `om = 2.*const%pi/const%seconds_per_year` *"set in constants; rads per second"*. **Recommended:** if `om` (annual angular frequency) is reused, store it as a derived constant rather than recomputing inline.
- **[source/model_settings.f90:25](source/model_settings.f90#L25)** — `days_per_year` *"remove once timestamps incorporated"*. **Recommended:** tied to the real-CF-timestamps refactor — once per-step timestamps come from the forcing netCDF, `days_per_year` (and `seconds_per_year` derived from it) can be dropped.

---

## 4. Config completeness

- **[source/initialise_variables.f90:50](source/initialise_variables.f90#L50)** — `Calc_Output_Freq`: *"add these to the output_dimensions config?"*. **Recommended:** surface the resulting `numOutput*` / `output*` counts (or their inputs) under `[output_dimensions]` so output cadence is configured in one place.

---

## Summary

| Theme | Markers | Effort |
|-------|---------|--------|
| 1. Physics choices → `model.toml` (+ conductivity dedupe) | 6 — `initialise_model.f90:91,108,109,120,242`, `firn_physics.f90:248` | 🔴 real work |
| 2. Validate `run.toml` fields | 1 — `model_settings.f90:454` | 🟡 small |
| 3. Centralise constants | 2 — `initialise_model.f90:238`, `model_settings.f90:25` | 🟡 trivial / tied to timestamps refactor |
| 4. Output-dimensions config | 1 — `initialise_variables.f90:50` | 🟡 small |

**Cleared since 2026-06-08:** `Find_Grid` (deleted); markers stripped from `output.f90:474`, `output.f90:574`, `model_settings.f90:347`, `model_settings.f90:580` (confirmed non-issues).
