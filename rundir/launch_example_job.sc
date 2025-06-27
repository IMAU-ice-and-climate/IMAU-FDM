#!/bin/bash

username="nld4814" # can be anything, but must be set and may need to update code_dir below
point_numb="1" # must be 1 unless new example files made
domain="FGRN055" # must be FGRN055 until ANT27 example files are made
project_name="example" # must be "example" for example mode to run in IMAU-FDM
forcing="era055" #used to set prefix, should match data
prefix_output="${domain}_${forcing}"
restart_type="none"

code_dir="/perm/${username}/code/IMAU-FDM/" #directory of IMAU-FDM
FDM_executable="${code_dir}imau-fdm.x" #makefile
logfile="${code_dir}example/logfiles/log_IMAU-FDM_cca_1.out" #logfile in example directory

localexe="${code_dir}example/LocalCode/imau-fdm.x" #to copy executable to example directory each time; allows for modifications to main executable during run
cp $FDM_executable $localexe #copies exeutable to example directory

echo "Running example FDM from $localexe; logfile at ${logfile}"
$localexe $username $point_numb $domain $prefix_output $project_name $restart_type $&> ${logfile} # launches example FDM run