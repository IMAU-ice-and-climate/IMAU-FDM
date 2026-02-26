#!/bin/ksh -x

fname="NCL_1d_"

# submit script    jobname    variab dt atzero
./submit_1d_job.sc ${fname}zs zs 0 0
./submit_1d_job.sc ${fname}vice vice 0 0
./submit_1d_job.sc ${fname}vacc vacc 0 0
./submit_1d_job.sc ${fname}vfc vfc 0 0
./submit_1d_job.sc ${fname}vmelt vmelt 0 0
./submit_1d_job.sc ${fname}vbuoy vbuoy 0 0
./submit_1d_job.sc ${fname}vsub vsub 0 0
./submit_1d_job.sc ${fname}vsnd vsnd 0 0
./submit_1d_job.sc ${fname}vtot vtotal 0 0
./submit_1d_job.sc ${fname}run Runoff 0 0
./submit_1d_job.sc ${fname}fac FirnAir 0 0
./submit_1d_job.sc ${fname}lwc TotLwc 0 0
./submit_1d_job.sc ${fname}refr refreeze 0 0
./submit_1d_job.sc ${fname}rain rain 0 0
./submit_1d_job.sc ${fname}melt surfmelt 0 0
./submit_1d_job.sc ${fname}solin solin 0 0
./submit_1d_job.sc ${fname}ice icemass 0 0
./submit_1d_job.sc ${fname}rho0 Rho0 0 0
