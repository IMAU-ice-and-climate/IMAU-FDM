#!/bin/ksh -x

./run_sc_job.sc years_evap makeFDMyears_FGRN11.sc evap
./run_sc_job.sc years_snowfall makeFDMyears_FGRN11.sc snowfall
./run_sc_job.sc years_snowmelt makeFDMyears_FGRN11.sc snowmelt
./run_sc_job.sc years_precip makeFDMyears_FGRN11.sc precip
./run_sc_job.sc years_sndiv makeFDMyears_FGRN11.sc sndiv
./run_sc_job.sc years_tskin makeFDMyears_FGRN11.sc tskin
./run_sc_job.sc years_ff10m makeFDMyears_FGRN11.sc ff10m
