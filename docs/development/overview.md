# Run Structure & Design Changes

The IMAU-FDM run structure was substantially reworked around the MPI
[distributor](distributor). The key changes:
- distributor.f90 is now the main program. MPI is used to allocate the points efficiently. The previous main.f90 is now main_model.f90.
- there is now a single launch_job script (launch_job/launch_job.sh). this creates the output directory structure. only one thing it only needs to be updated if the settings file structures are updated, otherwise it reads everything from the settings files.
- compilation works a bit differently; this has already been integrated with the old architecture but good to keep in mind that the Makefile is no longer sufficient for fully compiling the model. setting recompile = true in the settings file ensures the model is recompiled each time a new job is launched
- all settings now live in the `settings/` folder. these are written in toml. for each substantially new model run, a new folder with specific settings should be created to keep things clean. if new variables are added to the toml files, care needs to be taken that they are read correctly into the model (in the model_settings.f90 file).
- all settings, model code, and executable are now copied to each model run to ensure anyone in the future can really see what choices were made, what physics was used, etc. increases transparency at the expense of a tiny bit of file bloat.
- the output directories are simplified so the structure is updated:
```
PROJECT_NAME/
├── output/                       # point output files
├── post-process/                 # post-processed output (e.g. gridded timeseries of single variables, created during post-processing)
├── restart/
  ├── spinup                      # for running again, starting after spinup
  ├── run                         # for running again, from the end of this run
├── logfiles
  ├── model_logfiles/             # per-point logfiles
  ├── distributor_logfiles/       # distributor logfiles
├── localcode
  ├── settings/                   # settings files
  ├── source/                     # model code
  ├── imau-fdm                    # executable
  ├── job                         # .env, .sc, submission_iteration files for relaunching job
├── pointlists 
```
