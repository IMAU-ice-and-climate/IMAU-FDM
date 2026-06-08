# Code TODOs (from `TKTKTK` markers)

Generated from the 15 `TKTKTK` placeholders currently in `source/`. Each entry lists
the location, the in-code note, and a recommended action. Remove the `TKTKTK`
marker from the source when its item is done.

Grouped by theme, roughly highest-value first.

> **Verification pass (2026-06-08):** checked each marker against the code. 4 were
> stale non-issues (question answerable, no code change) and 1 was confirmed dead
> code.
>
> **Actioned 2026-06-08 (build passes — `Project compiled successfully`):**
> deleted the dead `Find_Grid` subroutine (+ its `public ::` entry) and stripped
> the 4 stale markers (`output.f90:474`/`:574` got a one-line clarifying comment;
> `model_settings.f90:347`/`:580` comment lines removed). **5 of 15 markers cleared.**
> The remaining 10 are the "real work" items below.

---

## 1. Make physics choices explicit in `model.toml`

The surface-snow density and thermal-conductivity parameterisations are hardcoded
and branch on `domain` (`FGRN*` vs everything-else). They should be selectable
physics options in `model.toml` so a run records exactly which scheme it used, and
so new domains don't silently fall into the "else" branch.

- **[source/initialise_model.f90:135](source/initialise_model.f90#L135)** — *"make choice of surface snow explicit in model.toml"* (the `numSnow == 1` instantaneous-density branch).
- **[source/initialise_model.f90:152-153](source/initialise_model.f90#L152)** — Greenland surface density `Rho0 = 362.1 + 2.78*(T - Tmelt)` is hardcoded as *Fausto et al. 2018*; also wants an option for constant density (315 kg m⁻³).
- **[source/initialise_model.f90:164](source/initialise_model.f90#L164)** — non-Greenland surface density `Rho0 = 82.97 + 0.769*T + 11.67*ff10` hardcoded as *Veldhuijsen 2023 / Lenaerts 2012*.
- **[source/initialise_model.f90:286](source/initialise_model.f90#L286)** — thermal-conductivity scheme in `Init_Temp_Prof` (Paterson 1994 / Calonne 2019 / Reid 1966 blend) is hardcoded.
- **[source/firn_physics.f90:248](source/firn_physics.f90#L248)** — same conductivity scheme duplicated in the `Thermal_Cond` function: *"all these model choices should be set explicitly in model.toml"*.

**Recommended:** add a `[model_physics]`-style block in `model.toml` (e.g.
`surface_density_scheme`, `thermal_conductivity_scheme`) and dispatch on it instead
of on `domain`. **Note the duplication:** the thermal-conductivity formulas appear
in *both* `firn_physics.f90:Thermal_Cond` and `initialise_model.f90:Init_Temp_Prof`
— consolidate into the single `Thermal_Cond` function and call it from both so the
scheme is defined once.

---

## 2. Remove redundant / verify dead code

- **[source/initialise_model.f90:16](source/initialise_model.f90#L16)** — `Find_Grid`: *"should be removed since lat/lon already found in distributor.f90"*. ✅ **CONFIRMED DEAD** — grep shows only the `public ::` export + the definition, **zero call sites** in `source/`. **Simple win:** delete the subroutine and remove it from the `public ::` list. Pure deletion, no logic to rework.
- **[source/output.f90:474](source/output.f90#L474)** — `outputProf = size(out_2D_dens, 1)` *"should already be set, check"*. ✅ **RESOLVED (non-issue).** `outputProf` is the **time** dimension (`ind_t`), and `size(out_2D_dens, 1)` correctly returns it. `Save_out_2D` doesn't receive the count as an argument, so recomputing from the array extent is valid and correct. Threading `outputProf` in as an arg is possible but has no functional benefit. → just drop the marker.
- **[source/output.f90:574](source/output.f90#L574)** — `outputDetail = size(out_2D_det_dens, 1)` *"should already be defined - check"*. ✅ **RESOLVED (non-issue)** — same as above (time dimension; `Save_out_2Ddetail` doesn't receive the count; recompute is correct). → drop the marker.
- **[source/model_settings.f90:580](source/model_settings.f90#L580)** — `Load_Constants` reads `constants.toml`: *"is this needed?"*. ✅ **RESOLVED — yes, it's needed.** Called from `model_main.f90:74`; `settings/FGRN055/constants.toml` exists; it opens that file and calls `Read_Constants` to parse it (same `Load_*`(open)/`Read_*`(parse) pattern as the rest). → drop the marker.

---

## 3. Settings-loading structure

- **[source/model_settings.f90:347](source/model_settings.f90#L347)** — *"why have this be separate from subroutine Load_Model_Settings()"*. ✅ **RESOLVED — intentional pattern.** `Load_Model_Settings` opens `model.toml` → table and calls `Read_Settings` to parse it; this mirrors `Load_Constants`/`Read_Constants` and `Read_Job`. Consistent by design. Merging is optional and low-value. → drop the marker.
- **[source/model_settings.f90:455](source/model_settings.f90#L455)** — `Read_Job`: *"check on domain, forcing, restart_type"*. **Recommended:** validate these run.toml fields against allowed values (`domain ∈ {FGRN055, ANT27}`, `restart_type ∈ {none, spinup, run}`, etc.) and fail early with a clear message rather than letting a typo propagate.

---

## 4. Constants that should be centralised

- **[source/initialise_model.f90:282](source/initialise_model.f90#L282)** — `om = 2.*const%pi/const%seconds_per_year` *"set in constants; rads per second"*. **Recommended:** if `om` (annual angular frequency) is reused, store it as a derived constant in the constants type rather than recomputing inline.
- **[source/model_settings.f90:25](source/model_settings.f90#L25)** — `days_per_year` *"remove once timestamps incorporated"*. **Recommended:** tied to the real-CF-timestamps refactor — once per-step timestamps come from the forcing netCDF, `days_per_year` (and `seconds_per_year` derived from it) can be dropped. See the timestamps work in the output/netCDF-metadata roadmap.

---

## 5. Config completeness

- **[source/initialise_variables.f90:50](source/initialise_variables.f90#L50)** — `Calc_Output_Freq`: *"add these to the output_dimensions config?"*. **Recommended:** the output-frequency derivation (prof/speed/detail) is computed from `writein*` + `dtobs`; consider surfacing the resulting `numOutput*` / `output*` counts (or their inputs) explicitly under `[output_dimensions]` so output cadence is configured in one place.

---

## Summary

After the verification pass, the 15 markers sort into:

| Status | Items | Action |
|--------|-------|--------|
| ✅ Resolved (stale, no code change) | 4 — `output.f90:474`, `output.f90:574`, `model_settings.f90:580`, `model_settings.f90:347` | just delete the markers |
| 🟢 Simple win (confirmed dead code) | 1 — `Find_Grid` (`initialise_model.f90:16`) | delete subroutine + `public ::` entry |
| 🟡 Trivial, low value | 1 — `om` (`initialise_model.f90:282`) | optional: move to a derived constant |
| 🔴 Real work | 5 physics → `model.toml` (+ conductivity dedupe), `455` validation, `50` output-dims config, `25` `days_per_year` (timestamps refactor) | scope separately |

So **none are bugs / already-broken**, 4 are answerable non-issues, and the only
clear quick implementation is removing the dead `Find_Grid`.
