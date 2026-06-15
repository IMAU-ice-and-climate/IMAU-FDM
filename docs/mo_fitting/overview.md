# MO Fitting

Source: `MO_fit/`

The **MO fitting** procedure calibrates the empirical densification-rate
coefficients against observed firn-core density profiles. The `MO-values` are the
accumulation-dependent correction factors applied to the firn densification rate
(the `MO_low` / `MO_high` multipliers in `source/firn_physics.f90`): they scale
the Stage-1 (ρ < 550 kg m⁻³) and Stage-2 (ρ > 550 kg m⁻³) compaction rates as a
function of mean accumulation `acav`, following the Arthern/Ligtenberg-type
parameterisation (the Greenland MO-values were most recently set for the 1939-2025 era5 run; the Antarctic fit follows Veldhuijsen et
al. 2023). 

## In detail

The model should be run at every core location/cell with an MO value of 1 (this can be set in the settings files). A curve is fit to each observed density profile and used to find the depths to $\rho=550$ and $\rho=830$ for each firn core in the observation dataset. 

The optimal MO-factor for each core is the ratio needed to bring the modeled depths into agreement with the observed (core) depths.

$$MO_{550} = d_{model, 550}/d_{obs, 550}$$
$$MO_{830} = (d_{model, 830} - d_{model, 550})/(d_{obs, 830} - d_{obs, 550})$$

These MO values are plotted against the average accumulation rate from the spin-up period $\bar{b}$, and the best-fit parameters are determined for \ref{eq:MO-fits} in python using the the \textit{curve\_fit} function from \textit{scipy.optimize}, which minimises the sum of squared residuals between the fitted curve and the MO values. The fits are of the form `MO = a + b·log(acav)` or
`MO = a·acav^b + c`.

The calibration is global: one set of coefficients is fitted to all cores
simultaneously, weighted by region.

## Notebooks and scripts

All code lives in `MO_fit/code/`:

| File | Description |
|------|-------------|
| `MO_fit.ipynb` | Main fitting notebook — runs optimisation, plots residuals |
| `MO_fit-visualizations.ipynb` | Plots and figures from the fit results |
| `mo_fit_functions.py` | Core fitting functions (depth extraction, log/power fits) |
| `process_firn_cores.ipynb` | Pre-process raw core data to standard format |
| `process_standardized_cores.ipynb` / `.py` | Ingest standardised core dataset |
| `process_sumup_cores.py` / `processing_sumup_cores.ipynb` | Process SUMup dataset; defines core regions (`def_models`) |
| `merge_datasets.ipynb` / `.py` | Merge multiple core datasets |
| `plot_core_locations.py` | Map of firn-core locations |

## Data

`MO_fit/data/` contains processed firn-core density profiles used for
calibration. Raw core data is not redistributed here; see the notebook headers
for data sources.

## Figures

`MO_fit/figures/` contains output plots from the fitting run.

````{admonition} Walkthrough notebook
:class: tip
To embed the MO fitting notebook here, copy `MO_fit/MO_fit.ipynb` to
`docs/notebooks/mo_fitting.ipynb` and add it to `_toc.yml`.
````
