#!/bin/bash
shopt -s expand_aliases  # Enables alias expansion.
echo "----------------------------------------------------------------------------------"
echo "$(date +%c): ns_script started"

# source environment file
myenvfile=$1
if [[ "$myenvfile" == "" ]]; then
  echo "ns_script: No environmentfile provided through the command line"
  exit 1
fi  
source $myenvfile

# read out gridpoint to do
gridpoint_c=$2
if [[ "$gridpoint_c" == "" ]]; then
  echo "ns_script: no gridpoint provided through the command line."
  exit 2
fi  
# this let ... is (possibly) needed to remove the leading space in gridpoint_c
let "gridpoint=$gridpoint_c*1"

# create readydir
export readydir="${readydir_base}/iter_${submission_iteration}"
mkdir -p $readydir
export EC_FARMID="ns"

# FDM executable called
echo "$(date +%c): ns_script: start running gridpoint $gridpoint"
$homedir/start_model_ccab.sc $hostname $p2exe/$FDM_executable $gridpoint 

echo "$(date +%c): ns_script: grid point ready, terminate."
echo "$(date +%c): ns_script ends"
echo "--------------------------------------------------------------------------------"
exit 0
# the script is ready
