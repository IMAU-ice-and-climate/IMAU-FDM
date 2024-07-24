#!/bin/ksh -x

./run_sc_job.sc ts_evap makeFDMfiles_FGRN055-era.sc evap
./run_sc_job.sc ts_snowfall makeFDMfiles_FGRN055-era.sc snowfall
./run_sc_job.sc ts_snowmelt makeFDMfiles_FGRN055-era.sc snowmelt
./run_sc_job.sc ts_precip makeFDMfiles_FGRN055-era.sc precip
./run_sc_job.sc ts_sndiv makeFDMfiles_FGRN055-era.sc sndiv
./run_sc_job.sc ts_tskin makeFDMfiles_FGRN055-era.sc tskin
./run_sc_job.sc ts_ff10m makeFDMfiles_FGRN055-era.sc ff10m
