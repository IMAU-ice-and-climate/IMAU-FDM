#!/bin/bash
shopt -s expand_aliases  # Enables alias expansion.

# FDM settings # copied from run_make_loadscript
export domain="Greenland"
project_name="grainsize"

export outputname="FGRN055_${myname}" #"FGRN055_noCB_nor"
export runname="FGRN055_${myname}" #"FGRN055_noCB_nor_run8"
export p2input="$PERM/code/DATA/IN_ll_FGRN055_paper1_imp.txt" #IN_ll_FGRN055_GrIS_GIC_implicit.txt"
export p2exe="$PERM/imau-fdm-develop/"
#export FDM_executable="IMAU-FDM_np_${myname}.x"
export FDM_executable="imau-fdm.x"
export homedir=`pwd`
gridpointlist="$homedir/pointlists/pointlist_paper1.txt" #pointlist_highres_todo.txt"    #pointlist_paper1_imp_newspin_todo.txt"
export ini_filename=""
export filename_part1="ECMWF_${outputname}"

# hardcoded FDM input
expname="era055/${myname}"
outputdir="$SCRATCH/data/output/$expname/" 
restartdir="$SCRATCH/IMAU-FDM_RACMO23p2/RESTART/"
export p2ms="$SCRATCH/data/ms_files/" # hardcoded in IMAU-FDM
# not hardcoded, FDM output
export p2logs="$SCRATCH/data/logfile/$expname/$runname"
# 
export walltime="48:00:00"   # (hms) walltime of the job 
export cooldown="00:00:30"   # (hms) how long prior end should focus shift to completing running jobs?
export workdir="$SCRATCH/FDMtests/${domain}/${runname}"
export hostname="cca"
export relaunch="no"        # with "no", no new iteration will be launched

# other FDM input
export usern=$USER

# likely not to change
export account_no="spnlberg"
export jobname_base="FDM_${myname}_i"
export nnodes_max=1 #8
export tasks_per_node=64 #9 #128 # this is not to be changed
export FDMs_per_node=128 #9 #128 # play around for the optimal performance 
export EC_hyperthreads=1
export memory_per_task="999Mb"
export taskfactor="1."                  # prior launch at least 5. task per core must be available    
export EC_ecfs=0 			# number of parallel ECFS calls 

# script misc
export workpointlist="$workdir/pointlist"
export readydir_base="$workdir/readypoints"
export readpointexe="$homedir/readpointlist/readpointlist.x"
export distributor="$homedir/readpointlist/distribute_points.x"
export requestdir="$workdir/requests"
export nplogdir="$workdir/nplogs"
export workexe="$workdir/LocalCode"
export submission_iteration=1


echo "We remove all current files in the workdir!"
echo "I'll sleep 5 seconds to allow you to abort this script if starting it was not what you want."
#sleep 5 

# -------- start of the launching procedure -------
echo "Start cleaning working directory"
# prepare launch
mkdir -p $workdir
rm -rf $workdir/*
mkdir -p $nplogdir
mkdir -p $workexe
mkdir -p $requestdir
if [[ ! -r $gridpointlist ]]; then
  echo "The grid point list is missing!"
  exit 1
fi  
cp $gridpointlist ${workpointlist}.txt

# create output & restart directory
mkdir -p $p2ms
mkdir -p $outputdir
mkdir -p $restartdir
mkdir -p $p2logs
  
let "ncpu_tot=$nnodes_max*$tasks_per_node"
echo "Run at max on $tasks_per_node cores on $nnodes_max node(s)."  
echo "So max parallel jobs is ${ncpu_tot}."

echo "Create environment file and submittable script."
./submit_job.sc


exit 0  
# ready

# worklist
# test
# improve performance distributor
