# The example file contains the input for one point for Greenland to test the FDM offline and/or on github

To run the FDM using the example input (in `example/input`), run `rundir/launch_example_job.sc`. The output goes to `example/output/`. 

The file structure looks like:

- example
    - input
        - timeseries
        _timeseries output for just one point (point 1 for FGRN055), created `example/pre-process/extract\_single\_point.py`_
        - averages
        _averages for all vars at all points between timesteps in filename_
    -LocalCode
        _new executable copied here from main IMAU-FDM every time example run is launched_
        - imau-fdm.x
    - logfiles
        - `log_IMAU-FDM_cca_1.out` contains the logfile from the model run using the example data
    - ms_files
        _settings file usually created during `./launch_new_job.sc` from `start_model_ccab.sc`_
        - model_settings_FGRN055_1.txt
    - output
        _standard output from FDM_
        - 1D.nc
        - 2D.nc
        - 2Ddet.nc
    - pre-process
        - `extract_single_point.py` contains the code to extract a single timeseries datapoint. to run it, you need to update the fnumb and have the rlat/rlon for the point. you can get this from the logfile created by running the point on the ECMWF
    - restart
        _files allow run to be started from after spinup or after last timestep of previous run_
        - restart_from_2023
        - restart_from_spinup



    

