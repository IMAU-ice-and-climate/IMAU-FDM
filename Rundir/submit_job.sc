#!/bin/bash
# a) This script writes the environmentfile of the new launch
# As we use environment variables, this file needs to be sourced
#   only once a np-job
# b) writes the submittable job file

echo "--------------------------------------------------------------------------------"

# figure out how many points needs to be done
ntodo_c=`$readpointexe ${workpointlist}.txt -count`
wronginput=0
let "ntodo=${ntodo_c}*1" || wronginput=1
if [[ $wronginput -eq 1 || "$ntodo" == "0" ]]; then
  echo "incorrect output of $readpointexe ${workpointlist}.txt -count"
  echo "  namely ${ntodo_c}, abort submitting new jobs."
  exit 1
fi  
echo "submit_job: $ntodo points still to do."

# calculate minimum number of task required to allow usage of all nodes
# first as a real (using) bc, then convert it to integer using printf
mintasks_real=`echo "${FDMs_per_node}*${nnodes_max}*${taskfactor}" | bc`
printf -v mintasks "%.0f" "$mintasks_real"
  
if [[ $mintasks -le $ntodo ]]; then
  echo "submit_job: Enough tasks to use all ${nnodes_max} node(s)."
  nnodes=$nnodes_max
  jobtype="np"   #Change Veldhuijsen
else
# derive needed number of nodes (as a real) and take int value while rounding down
  nnodes_real=`echo "${ntodo}.0/(${FDMs_per_node}.0*${taskfactor})" | bc`
# somehow bc provides already a rounded down integer, nevertheless...  
  nnodes=${nnodes_real%.*} 
  echo "nnodes" $nnodes
  if [[ $nnodes -gt 0 ]]; then
    echo "submit_job: Only ${nnodes} node(s) needed now."
    jobtype="np"  #Change Veldhuijsen
  else
    ntasks_real=`echo "${ntodo}/${taskfactor}" | bc`
    ntasks=${ntasks_real%.*}
      let "halfnode=${tasks_per_node}/2"
      if [[ $ntasks -gt $halfnode ]]; then
        ntasks=$halfnode
      fi
      jobtype="nf"
      echo "submit_job: Only ${ntasks} tasks in a nf job needed now."
      maxFDMs=${ntasks}  
  fi    	 
fi

echo "ntasks = ${ntasks}" 

# get maximum number of parallel FDM runs in a job
if [[ "$jobtype" == "np" ]]; then
# tasks are defined by the number of tasks per node, cool down at the end
  let "maxFDMs=${nnodes}*${FDMs_per_node}"
  let "coolFDMs=${maxFDMs}/2"
elif [[ "$jobtype" == "nf" ]]; then
# a FDM run for every task, no cool down
  maxFDMs=$ntasks
  coolFDMs=$ntasks
else
# ns
  maxFDMs=1
  coolFDMs=1
fi    

myname="${jobname_base}${submission_iteration}_${jobtype}"
runlog="${jobname_base}${submission_iteration}_runlog"    
envfile="${myname}.env"

echo "submit_job: write environment file"
cat << EOE > $workdir/$envfile
# Date  creating Env file: $(date +%c)
# At this point, ${ntodo} point were still to do

# script essentials
export workdir=$workdir
export homedir=$homedir
export readydir_base=$readydir_base
export nplogdir=$nplogdir
export workexe=$workexe
export requestdir=$requestdir
export hostname="$hostname"

export readpointexe=$readpointexe
export distributor=$distributor

# file names
export jobname_base=$jobname_base

# FDM essentials
export domain="$domain"
export p2input="$p2input"
export p2exe="$p2exe"
export FDM_executable="$FDM_executable"

# other FDM input
export p2ms="$p2ms"
export p2logs="$p2logs"
export usern="$usern"
export filename_part1="$filename_part1"
export ini_filename="$ini_filename"
export workpointlist="$workpointlist"

# job management
export submission_iteration=$submission_iteration
export relaunch=$relaunch
export walltime="$walltime"
export cooldown="$cooldown"
export memory_per_task="$memory_per_task"
maxFDMs=$maxFDMs
coolFDMs=$coolFDMs
export tasks_per_node=$tasks_per_node
export FDMs_per_node=$FDMs_per_node
export nnodes_max=$nnodes_max
export taskfactor=$taskfactor
export account_no=$account_no
export EC_ecfs=$EC_ecfs
export jobtype=$jobtype
EOE

if [[ "$jobtype" == "np" ]]; then
  echo "submit_job: Create submittable script for a np job."
  jobname="${jobname_base}${submission_iteration}_np"
  echo "submit_job: Jobname: ${jobname}"
   # set queue script
cat << EOFp > $workdir/${jobname}.sc
#!/bin/bash

#SBATCH -q np 
#SBATCH --cpus-per-task=1
#SBATCH --nodes=$nnodes
#SBATCH -J $jobname 
#SBATCH --time=$walltime
#SBATCH -o $nplogdir/${myname}.log
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=$tasks_per_node
#SBATCH --mem-per-cpu=$memory_per_task
#SBATCH --threads-per-core=$EC_hyperthreads


# launch script
$homedir/npnf_outer_script.sc $workdir/$envfile >& $nplogdir/${runlog}.log

exit 0
EOFp

  # make it executable
  chmod u+x $workdir/${jobname}.sc
  
  sbatch $workdir/${jobname}.sc

elif [[ "$jobtype" == "nf" ]]; then
  echo "submit_job: Create submittable script for a nf job"
  jobname="${jobname_base}${submission_iteration}_nf"
  echo "submit_job: Jobname: ${jobname}"
   # set queue script
cat << EOFf > $workdir/${jobname}.sc
#!/bin/bash

#SBATCH -q nf 
#SBATCH --cpus-per-task=1
#SBATCH -J $jobname 
#SBATCH --time=$walltime
#SBATCH -o $nplogdir/${myname}.log
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=$tasks_per_node      
#SBATCH --mem-per-cpu=$memory_per_task
#SBATCH --threads-per-core=$EC_hyperthreads

# launch script
$homedir/npnf_outer_script.sc $workdir/$envfile >& $nplogdir/${runlog}.log

exit 0
EOFf
  echo "try1"
  
  # make it executable
  chmod u+x $workdir/${jobname}.sc
  
  sbatch $workdir/${jobname}.sc
fi

# create a nicer abort script (so that distribute point creates an updated to-do list
if [[ $submission_iteration == 1 ]]; then
  if [[ "$jobtype" == "np" || "$jobtype" == "nf" ]]; then
    echo "Use CancelMyJob.sc to cancel this job nicely, so that you have an updated pointlist"

    let "threadmx=$maxFDMs-1"

cat << EOC > $homedir/CancelMyJob.sc
#!/bin/bash
echo "Cancel job ${runname}. If that is not what you want, you have five seconds left."
sleep 5

echo "Start canceling"
# post to all thread request files "fatal". To be sure that the distributor picks it up
for it in {0..${threadmx}}; do
  requestfile=${requestdir}/\`printf "%.5d\\n" \$it\`
  echo "fatal" > \$requestfile
done  
echo "Done"
exit 0
EOC
  chmod u+x $homedir/CancelMyJob.sc
  fi
fi

echo "Leave submit_job."
echo "--------------------------------------------------------------------------------"
exit 0
# ready
