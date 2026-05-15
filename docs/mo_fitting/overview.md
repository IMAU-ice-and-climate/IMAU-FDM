# MO Fitting

Source: `MO_fit/`

The **MO fitting** procedure calibrates the empirical densification-rate
coefficients against observed firn-core density profiles. The acronym MO
refers to Medley-Olathe (or the specific parameterisation family used here).

## What it does

For each firn core in the training dataset, the model is run with a range of
coefficient values. The best-fit coefficients minimise the RMSE between
simulated and observed density–depth profiles.

The calibration is global: one set of coefficients is fitted to all cores
simultaneously, weighted by region.

## Notebooks and scripts

| File | Description |
|------|-------------|
| `MO_fit.ipynb` | Main fitting notebook — runs optimisation, plots residuals |
| `mo_fit_functions.py` | Core fitting functions |
| `process_firn_cores.ipynb` | Pre-process raw core data to standard format |
| `process_standardized_cores.ipynb` | Ingest standardised core dataset |
| `merge_datasets.ipynb` | Merge multiple core datasets |
| `processing_sumup_cores.ipynb` | Process SUMup community dataset |

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
