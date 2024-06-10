#!/bin/bash
shopt -s expand_aliases  # Enables alias expansion.
echo "----------------------------------------------------------------------------------"
echo "$(date +%c): npnf_outer: started"

# source environment file
myenvfile=$1
if [ "$myenvfile" == "" ]; then
  echo "No environmentfile provided through the command line"
  exit 1
fi  
source $myenvfile

# load this module to allow for numberin of the ranks in a nf job.
#if [[ "$jobtype" == "nf" ]]; then
  #module load cray-snplauncher 
#fi  

echo "npnf_outer: Run at max $maxFDMs jobs at the same time (cooldown = ${coolFDMs} jobs)"

# set ready-dir for this iteration
export readydir="${readydir_base}/iter_${submission_iteration}"
mkdir -p $readydir
rm -f $readydir/*

# clean requestdir, backup workpointlist and start the distributor
rm -f $requestdir/*
cp ${workpointlist}.txt ${workpointlist}_${submission_iteration}.txt
$distributor ${workpointlist}.txt $maxFDMs ${requestdir} -t $p2input & #>& ${workdir}/DP_log.txt &

# make that every code works individually
export PMI_NO_FORK=1

# make logdirectory for inner for this iteration
export innerlogdir="$nplogdir/Threads_iter_${submission_iteration}"
mkdir -p $innerlogdir

#------------------- launch parallel task script -----------
if [[ "$jobtype" == "np" ]]; then
  srun -n ${maxFDMs} --ntasks-per-node 72 --threads-per-core 1 $homedir/npnf_inner_script.sc $myenvfile
elif [[ "$jobtype" == "nf" ]]; then
  srun -n ${maxFDMs} $homedir/npnf_inner_script.sc $myenvfile
else
  echo "The npnf_outer: unknown jobtype $jobtype"
  exit 1
fi  


#------------------- post run work ---------------
echo "$(date +%c): All npnf_inner_scripts are done, wait for distributor"
nwait=0
waits=2
while [[ ! -r "${requestdir}/DP" && $nwait -lt 120 ]]; do
  sleep $waits
  let "waits+=$waits"
done

#echo "$(date +%c): -----------distributor timer info----------"
#cat ${workdir}/timer
#echo "$(date +%c): -------------------------------------------"
#rm ${workdir}/timer

if [[ ! -r "${requestdir}/DP" ]]; then
  echo "$(date +%c): Distributor is not responding in time, terminate job"
  nextaction="stop"
else    
  nextaction=`head -1 "${requestdir}/DP"`
  echo "$(date +%c): Distributor asked for ${nextaction}."
fi

if [[ "$nextaction" == "continue" && "$relaunch" != "no" ]]; then
# run another time 
# update iteration & environmentfile
  let "submission_iteration+=1"
  echo "$(date +%c): Launch iteration ${submission_iteration}."
  
  ${homedir}/submit_job.sc 
fi

# combine list of completed tasks into one file
ls -1 $readydir >> "${readydir_base}/AllCompletedPoints.txt"

# clean work executables
rm $workexe/*

echo "$(date +%c): np_script ended"
echo "--------------------------------------------------------------------------------"
exit 0
# the script is ready
