#!/bin/bash

paths_file="settings/paths.toml" # can be anything, but must be set and may need to update code_dir below
point_numb="1" # must be 1 unless new example files made
domain="FGRN055" # must be FGRN055 until ANT27 example files are made
project_name="example" # must be "example" for example mode to run in IMAU-FDM
forcing="era055" #used to set prefix, should match data
prefix_output="${domain}_${forcing}"
restart_type="none"

logfile="example/logfiles/log_IMAU-FDM_cca_1.out" #logfile in example directory

localexe="/home/jelle/.local/bin/imau-fdm"

echo "Running example FDM from $localexe; logfile at ${logfile}"
$localexe "$paths_file" $point_numb $prefix_output $project_name $restart_type $&> ${logfile} # launches example FDM run
