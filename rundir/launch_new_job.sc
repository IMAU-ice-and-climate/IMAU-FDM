#!/bin/bash
shopt -s expand_aliases  # Enables alias expansion.

##### SETS PATHS AND SBATCH OPTIONS FOR FDM ####
##### 
##### project_name = unique name for run; must match pointlist
##### run_location = paths set to defautl for ECMWF, otherwise need to specify
##### domain = FGRN055 for greenland, ANT27 for antarcitca, or custom
##### forcing = era055 for greenalnd, era027 for antarctica, or custom
##### restart_type = none, spinup, run
##### 
##### nnodes = 1 if only running a few points; 8 if doing a continent-sized run
##### taskfactor set automatically depending on nnodes
#####
##### all other vars/paths assume run is on ECMWF & IMAU-FDM is structured as on github


export project_name="FGRN055-era055_1939-2023" # set unique project_name; pointlist must have matching name e.g. pointlist_PROJECT_NAME.txt

if [[ -z "$project_name" ]]; then
  echo "project_name is empty; set before continuing"
  exit 1
fi

# set domain and forcing
export domain="FGRN055" #or ANT27

if [[ "${domain}" == "FGRN055" ]]; then
  forcing="era055"
elif [[ "${domain}" == "ANT27" ]]; then
  forcing="era027"
else
  echo "must specify forcing for unrecognized domain"
  exit 1
fi

export restart_type="none" # none - do spinup; spinup - restart from spinup; (testing -> run - restart from run)

export outputname="${domain}_${forcing}"
export runname="${domain}_${project_name}"

export p2exe="${PERM}/code/IMAU-FDM"
export workdir="${SCRATCH}/${project_name}"
export restartdir="${SCRATCH}/restart/${project_name}" #"$SCRATCH/IMAU-FDM_RACMO23p2/RESTART/"
export outputdir="${workdir}/output" #"$SCRATCH/data/output/$expname/" 

export p2input="$p2exe/reference/${domain}/IN_ll_${domain}.txt"
export FDM_executable="imau-fdm.x"
export homedir=`pwd`
#gridpointlist="$p2exe/rundir/pointlists/pointlist_${project_name}.txt"
#gridpointlist="$p2exe/rundir/pointlists/pointlist_run_1939-2023_FGRN055-era055.txt"
gridpointlist="$p2exe/rundir/pointlists/pointlist_FGRN055-era055_1939-2023.txt"

export p2ms="${workdir}/ms_files" #"$SCRATCH/data/ms_files/" # hardcoded in IMAU-FDM
export p2logs="${workdir}/logfiles" #"$SCRATCH/data/logfile/$expname/$runname"
export filename_part1="${outputname}"
# 
export walltime="48:00:00"   # (hms) walltime of the job 
export cooldown="00:00:30"   # (hms) how long prior end should focus shift to completing running jobs?
export hostname="cca"
export relaunch="yes"        # with "no", no new iteration will be launched

# other FDM input
export usern=$USER

# SBATCH options
export nnodes_max=16 #update to 8 if doing full run, otherwise use 1 for smaller runs
export account_no="spnlberg"
export jobname_base="FDM_${project_name}_i"
export FDMs_per_node=128 #128 # play around for the optimal performance 
export EC_hyperthreads=1
export memory_per_task="999Mb"
export EC_ecfs=0      # number of parallel ECFS calls 
if [[ nnodes_max -gt 2 ]]; then
  export taskfactor="5."   # prior launch at least #taskfactor (3-5) tasks per core must be available, change to 1 if just running test point
  export tasks_per_node=128
else
  export tasks_per_node=64
  export taskfactor="1."
fi

# script misc
export workpointlist="$workdir/pointlist"
export readydir_base="$workdir/readypoints"
export readpointexe="$p2exe/rundir/readpointlist/readpointlist.x"
export distributor="$p2exe/rundir/readpointlist/distribute_points.x"
export requestdir="$workdir/requests"
export nplogdir="$workdir/nplogs"
export workexe="$workdir/LocalCode"
export submission_iteration=1

# echo all filepaths to double check
echo "outputdir: ${outputdir}"
echo "restartdir: ${restartdir}"
echo "p2ms: ${p2ms}"
echo "p2logs: ${p2logs}"
echo "outputname: ${outputname}"
echo "p2input: ${p2input}"
echo "p2exe: ${p2exe}"
echo "homedir: ${homedir}"
echo "gridpointlist: ${gridpointlist}"
echo "workdir: ${workdir}"
echo "workpointlist: ${workpointlist}"
echo "readydir_base: ${readydir_base}"
echo "readpointexe: ${readpointexe}"
echo "distributor: ${distributor}"
echo "requestdir: ${requestdir}"
echo "nplogdir: ${nplogdir}"
echo "workexe: ${workexe}"
echo "restart_type=${restart_type}"

echo "Continuing will submit the job."
echo
read -p "Check paths. Want to continue? (y/n)" -n 1 -r
echo

if [[ $REPLY =~ ^[Yy]$ ]]; then

  # -------- start of the launching procedure -------
  echo "Start cleaning working directory"
  # prepare launch
  mkdir -p $workdir
#  rm -rf $workdir/*
  mkdir -p $nplogdir
  mkdir -p $workexe
  mkdir -p $requestdir
  mkdir -p $readydir_base
  mkdir -p $p2ms
  mkdir -p $outputdir
  mkdir -p $restartdir
  mkdir -p $p2logs

  cp "$p2exe/$FDM_executable" "$workexe/$FDM_executable"

  if [[ ! -r $gridpointlist ]]; then
    echo "The grid point list is missing!"
    exit 1
  fi  
  cp $gridpointlist ${workpointlist}.txt
    
  let "ncpu_tot=$nnodes_max*$tasks_per_node"
  echo "Run at max on $tasks_per_node cores on $nnodes_max node(s)."  
  echo "So max parallel jobs is ${ncpu_tot}."

  echo "Create environment file and submittable script."
  ./submit_job.sc

fi


exit 0  
