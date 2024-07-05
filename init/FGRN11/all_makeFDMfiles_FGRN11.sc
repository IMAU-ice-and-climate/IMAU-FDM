#!/bin/ksh -x

./run_sc_job.sc fil_evap        makeFDMfiles_FGRN11.sc evap
./run_sc_job.sc fil_snowfall    makeFDMfiles_FGRN11.sc snowfall
./run_sc_job.sc fil_snowmelt    makeFDMfiles_FGRN11.sc snowmelt
./run_sc_job.sc fil_precip      makeFDMfiles_FGRN11.sc precip
./run_sc_job.sc fil_sndiv       makeFDMfiles_FGRN11.sc sndiv
./run_sc_job.sc fil_tskin       makeFDMfiles_FGRN11.sc tskin
./run_sc_job.sc fil_ff10m       makeFDMfiles_FGRN11.sc ff10m
