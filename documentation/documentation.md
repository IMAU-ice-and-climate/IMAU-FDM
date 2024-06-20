# Overview of running IMAU-FDM:
_Sanne Veldhuijsen – Version FDM v1.2AD (Antarctica)_

ATOS:
File system: Home, PERM  code, SCRATCH  input/output files, will be deleted after 4 weeks

Documentation atos: 
https://confluence.ecmwf.int/display/UDOC/HPC2020+User+Guide

Here you can see you SBU usage: 
https://hpc-usage.ecmwf.int/sbu-accounting/my-accounts

1.	Model code  
Directory: /perm/rusv/MSc_thesis/Program_FDMv12AD_SSP585/

Fortran files which call subroutines: 
-	mainprogram_firnmodel.f90: calls initializations subroutines  
-	subprogram.f90: calls spin-up, time-integration and writing subroutines 

Fortran files containing subroutines: 
-	ini_model.f90: reading settings and constants, finding grid point, interpolation of forcing data, initializion of firn pack. 
-	openNetCDF.f90: opens atmospheric forcing files.  
-	output.f90: writes output files 1D , 2D and 2Ddetail (see below for details) 
-	time_loop.f90: spin-up, time-integration firn processes, layer splitting and merging 

Makefile: run this file by typing “Make” to compile fortran 
 IMAU-FDM.x (this is the compiled file)

2.	Files: 
Atmospheric forcing files: 
Location: /ec/res4/scratch/rusv/MSc_thesis/data/SSP585/input/files/
3-hourly: skin temperature, snowmelt, precipitation, snowfall, windspeed, evaporation, snowdrift 

Setting files: 
Location: /perm/rusv/MSc_thesis/data/
-	grid file: IN_ll_ANT27_CB1979-2020_nor.txt 	             (e.g. x,y coordinates)
-	input settings file: input_settings_ANT27_plus.txt    (numLons, numLats, numTimes)
Location: /perm/rusv/MSc_thesis/NPNF_Run/
-	start_model_ccab.sc: model settings for FDM  e.g. amount of layers, timestep, writing output speed, amount of years 
-	launc_new_job.sc  FDM settings: domain, runname, location to executable etc. (further explained in section 3).

Output files:
Location: /ec/res4/scratch/rusv/MSc_thesis/data/SSP585/output/
Ini, 1D (10-days), 2D (30-days), 2Ddetail (10-days): 

Restart files:
Location: /ec/res4/scratch/rusv/MSc_thesis/data/restart_files/
restart_1950, restart_2014,

3.	Bash scripts that launch the model
Location: /perm/rusv/MSc_thesis/NPNF_run/
From here you can run NP jobs (using a full node – all cores) or NF jobs (using a fraction of a node)
For now, you are going to use NF jobs only 

How to launch a new job? 
./launch_new_job.sc 
This script reads the points which are needed to be run (pointlist.txt)
And then it calls: submit_job.sc > npnf_outer_script.sc > npnf_inner_scripts.sc > start_model_ccab.sc 
Thereby the simulations get started. 
All the details of the scripts are for now not relevant, if something doesn’t work you can ask me.  

4.	Processing raw data to gridded data 
Example script: /perm/rusv/MSc_thesis/python_scripts/ANT27_make_gridded_output_SSP585.py
From raw data (1 grid point), you can make 2D files (Antarctic wide or ice shelf wide maps)

5.	Working on ATOS: 
Several: 
squeue  -u “username”  job information 
quota  see how much storage is used 
scancel “jobid”  cancel job 

Touch all files recursively:
In a scratch directory, open a “screen”:
then execute: “find . -type f -exec touch {} +”

How to download files to your own pc, example:
scp -r -J rusv@jump.ecmwf.int rusv@hpc-login:'/ec/res4/scratch/rusv/backup_HPC/output/ant27_exp_MO_2022/*_12922.nc' . 



