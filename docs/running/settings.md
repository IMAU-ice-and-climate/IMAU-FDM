# Settings Files

There are currently three locations that need to be updated or checked before a run. This will be updated eventually with a set of toml files in `settings`. However currently these are:

1. `rundir/launch_new_job.sc` - This is were you set the model domain and forcing, the pointlist, the SLURM settings, and the restart type and location of the restart files (if relevant)
2. `rundir/start_model_ccab.sc` - Each domain has its own settings. Check to ensure the number of timesteps, number of years, number of spinup years, and the domain dimensions match the input. Also check the output timesteps and dimensions are what you want. 
3. `source/model_settings.f90` - This is where output metadata and paths are set, along with model constants constants. 