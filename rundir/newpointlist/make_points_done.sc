#!/bin/bash

project_name="FGRN055-era055_1939-2023"
outdir="$SCRATCH/${project_name}/output/"
restartdir="$SCRATCH/restart/${project_name}/"

name_ini="$restartdir*spinup*"
name_1d="$outdir*1D*"
name_2d="$outdir*2D_*"
name_2ddet="$outdir*2Dd*"

ls -1 $name_ini > points_ini_done.txt &
ls -1 $name_1d > points_1D_done.txt &
ls -1 $name_2d > points_2D_done.txt &
ls -1 $name_2ddet > points_2Ddet_done.txt &
