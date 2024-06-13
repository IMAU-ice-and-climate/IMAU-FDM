#!/bin/bash
shopt -s expand_aliases  # Enables alias expansion.

# as the script terminates once one terminates, we wait upon stop call of the distributor

# messages to distributor are
# provide  -> i 2 d: a new grid point is requested 
# nonew    -> i 2 d: time is almost done, no new point
# abort    -> i 2 d: time is done, abort calculations / or job is not running properly
# fatal    -> i 2 d: error occured, abort without restart
# wait     -> d 2 i: points are done, (or too little time left), wait for stop
# stop     -> d 2 i: stop calculations
# <number> -> d 2 i: a grid point to do

# script settions
sleeptime=10 # seconds, defines time between different checks
exittime=600 # how many seconds before the end do we terminate
nlistmax=600 # how long (s) do we wait for a new job? Maybe the distributor "died"?

# identify "my rank"
# ranks go from 0 to ntasks-1
if [[ "$jobtype" == "np" ]]; then
  rank=$SLURM_PROCID
else # we assume nf
  rank=$SLURM_PROCID
fi    
rankmss="inner #${rank}:"

echo "rank = ${rank}"
echo "rank = ${rankmss}"

# source environment file
myenvfile=$1
if [[ "$myenvfile" == "" ]]; then
  echo "$(date +%c): ${rankmss} No environmentfile provided through the command line"
  exit 1
fi  
source $myenvfile

# set logfile
logfile="${innerlogdir}/Thread_$rank"
echo "----------------------------------------------------------------------------------" > $logfile
echo "$(date +%c): ${rankmss} started" >> $logfile



# define end second 
hours=`echo $walltime | awk '{print substr($1,1,2)}'`
mins=`echo $walltime | awk '{print substr($1,4,2)}'`
secs=`echo $walltime | awk '{print substr($1,7,2)}'`
endsec=`echo "($hours*3600)+($mins*60)+$secs-$exittime" | bc`

if [[ $rank -ge $coolFDMs ]]; then
  # define cool down point
  hours=`echo $cooldown | awk '{print substr($1,1,2)}'`
  mins=`echo $cooldown | awk '{print substr($1,4,2)}'`
  secs=`echo $cooldown | awk '{print substr($1,7,2)}'`
  coolsec=`echo "$endsec-($hours*3600)-($mins*60)-$secs" | bc`
else
  coolsec=$endsec  
fi
echo "$(date +%c): ${rankmss} Terminate script after $endsec seconds, no new jobs after $coolsec seconds." >> $logfile

# define ready-dir for this iteration
export readydir="${readydir_base}/iter_${submission_iteration}"

# define request file
requestfile=${requestdir}"/"`printf "%.5d\n" $rank`
oldrequest=${requestfile}"_old"
# echo "$(date +%c): ${rankmss} Use requestfile $requestfile"

# change compilor
#prgenvswitchto gnu >> $logfile 2>> $logfile

## hopefully this helps to prevent pmi errors
#export PMI_MMAP_SYNC_WAIT_TIME=360
# Well, it didn't help! 

# start loop
action="request"
newaction="NULL"
gridpoint="NULL"
duration=$SECONDS
nlisten=0
npimerr=0


while [[ "$action" != "stop" ]]; do

 if [[ "$action" == "request" ]]; then
# put request into communicating file
   let "timeleft=${endsec}-${duration}"
   echo "$(date +%c): ${rankmss} Request for a grid point (${timeleft} s left), wait 2 seconds." >> $logfile
   echo "provide" > $requestfile
   echo ${timeleft} >> $requestfile
   action="listen"
   nlisten=0
# wait 2 seconds for response         
   sleep 2
 fi
 
 if [[ "$action" == "listen" ]]; then   
# read out reply
   newaction=`head -1 $requestfile`
# adjust action
   if [[ "$newaction" == "provide" ]]; then
     echo "$(date +%c): ${rankmss} The distributor is slow in responding" >> $logfile
     let "nlisten+=$sleeptime"
     if [[ $nlisten -gt $nlistmax ]]; then
       echo "$(date +%c): ${rankmss} The distributor did not respond within $nlistmax s. abort" >> $logfile
       echo "$(date +%c): ${rankmss} The distributor did not respond within $nlistmax s. abort" 
       action="stop"     
     fi
# it's up to the distributor to remove this message     
   elif [[ "$newaction" == "stop" ]]; then
# we are running out of time, leave loop (but won't be used I presume) 
     echo "$(date +%c): ${rankmss} Leave script on request of distributor" >> $logfile
     action="stop"
   elif [[ "$newaction" == "wait" ]]; then
# the distributor has no new task for us, but let others run until completed   
     echo "$(date +%c): ${rankmss} Wait until called to stop" >> $logfile
     action="wait"        
     mv $requestfile $oldrequest
   else
# we got a gridpoint number  
     expruntime=`head -2 $requestfile | tail -1`
     if [[ "$expruntime" != "" ]]; then
       let "expruntime=${expruntime}*1" # remove leading spaces
     else
       expruntime="<not given>"
     fi
             
     action="run"
     mv $requestfile $oldrequest
     # this let ... is needed to remove the leading space in gridpoint_c  
     let "gridpoint=${newaction}*1"
     
     echo "$(date +%c): ${rankmss} Start running gridpoint ${gridpoint}, expected runtime ${expruntime} minutes" >> $logfile
   fi
 fi
 
 if [[ "$action" == "run" ]]; then  
# make local copy
   localexe="$workexe/IMAU-FDM_np_gp_${gridpoint}.x"
   localscp="$workexe/${gridpoint}_start_model.sc"
   cp $p2exe/$FDM_executable $localexe
   cp $homedir/start_model_ccab.sc $localscp
# launch new task
   $localscp $hostname $localexe $gridpoint >> $logfile 2>> $logfile &
   scriptpid=$!
   rstartsec=$SECONDS
   rprintsec=$rstartsec
# this is the duplicate of the name in start_model_ccab    
   log_fname=${p2logs}/log_IMAU-FDM_${hostname}_${gridpoint}.out
   action="running"
   
# check if everyting goes fine, but first wait 
   sleep 300
   
   if [[ -r $log_fname ]]; then
     mpierr=`grep pmi_mmap_tmp $log_fname`
     mpiokm=`grep _pmi_init $log_fname`
     if [[ "$mpierr" != "" && "$mpiokm" == "" ]]; then
       echo "$(date +%c): ${rankmss} FATAL pmi_mmap_tmp ERROR" >> $logfile
       echo "$(date +%c): ${rankmss} FATAL pmi_mmap_tmp ERROR"
       echo "$(date +%c): ${rankmss} job id = ${scriptpid}" >> $logfile
       kill -9 ${scriptpid} >> $logfile 2>> $logfile
       exepid=`ps T | grep $localexe | grep -v grep`
       echo "$(date +%c): ${rankmss} exe id = ${exepid}" >> $logfile
       kill -9 ${exepid:0:6} >> $logfile 2>> $logfile
       exepid=`ps T | grep $localexe | grep -v grep`
       echo "$(date +%c): ${rankmss} exe id (check) = ${exepid}" >> $logfile
       echo "$(date +%c): ${rankmss} Stop this task or retry." >> $logfile
       let "npimerr+=1"
       
       if [[ $npimerr -ge 5 ]]; then
         echo "$(date +%c): ${rankmss} We give up now" >> $logfile
	 echo "$(date +%c): ${rankmss} We give up now"
         echo "abort" > $requestfile
         action="wait"
       else
         echo "$(date +%c): ${rankmss} Another try" >> $logfile
	 echo "$(date +%c): ${rankmss} `ps T | grep ${gridpoint}_start_model | grep -v grep`" >> $logfile
	 action="run"
       fi
       	 	
     elif [[ "$mpierr" != "" ]]; then
       echo "$(date +%c): ${rankmss} non-fatal pmi_mmap_tmp error" >> $logfile
#       echo "$(date +%c): ${rankmss} non-fatal pmi_mmap_tmp error" # overflows main log file
       npimerr=0
     else
       echo "$(date +%c): ${rankmss} proper launch of job" >> $logfile
       npimerr=0
     fi  
   else
     echo "$(date +%c): ${rankmss} Cannot read logfile" >> $logfile	 
     echo "$(date +%c): ${rankmss} $log_fname" >> $logfile
   fi	 
   
      
 fi
 
# now sleep for a while
 sleep $sleeptime 
   
 if [[ "$action" == "running" ]]; then
# check if task is ready 
   if [ -r $readydir/$gridpoint ]; then
     rendsec=$SECONDS
     runtime=`echo "($rendsec-$rstartsec)/60" | bc` 
# echo to both inner and outer log     
     echo "$(date +%c): ${rankmss} Grid point $gridpoint is ready! (runtime $runtime minutes, expected ${expruntime})" >> $logfile
     echo "$(date +%c): ${rankmss} # $gridpoint ready! (rt $runtime minutes, exprt ${expruntime})" 
# ask for next
     if [[ $duration -lt $coolsec ]]; then
       action="request"
     else
       action="nonew"
     fi    
     echo "$(date +%c): ${rankmss} New action: $action (${duration} wrt ${coolsec})" >> $logfile

   else
     let "runsec=$duration-$rprintsec"
     rundff=`echo "(($duration-$rstartsec)/60)-$expruntime" | bc`
     if [[ $runsec -gt 3600 ]]; then
       echo "$(date +%c): ${rankmss} ... still running ...  $rundff m" >> $logfile
       if [[ $rundff -gt 0 ]]; then
         echo "$(date +%c): ${rankmss} ... $gridpoint still running ... +${rundff} m on ${expruntime} m" 
       fi
       rprintsec=$duration
     fi  
   fi
   
# check if actually the distributor would like to stop
   if [[ -r $requestfile ]]; then
     newaction=`head -1 $requestfile`          
     if [[ "$newaction" == "stop" ]]; then       
       rendsec=$SECONDS
       runtime=`echo "($rendsec-$rstartsec)/60" | bc` 
       echo "$(date +%c): ${rankmss} The distributor asks to stop, which we do" >> $logfile
       echo "$(date +%c): ${rankmss} We stop, but grid point $gridpoint is not yet ready (runtime $runtime minutes, expected ${expruntime})." >> $logfile
       echo "$(date +%c): ${rankmss} Distributor request to stop, but grid point $gridpoint is not yet ready (runtime $runtime minutes, expected ${expruntime})." 
       action="stop"
     else
       echo "$(date +%c): ${rankmss} Strange message $newaction" >> $logfile
       mv $requestfile $oldrequest
     fi
   fi
 fi
 
 if [[ "$action" == "nonew" ]]; then
# notify to distributor that we don't want a new point.
   echo "$(date +%c): ${rankmss} Ask for termination while last point is completed" >> $logfile
   echo "nonew" > $requestfile
   action="wait"
 fi  

 if [[ "$action" == "wait" ]]; then
# wait until the distributor says it's time to stop
   if [[ -r $requestfile ]]; then
     newaction=`head -1 $requestfile`
     if [[ "$newaction" == "stop" ]]; then
       echo "$(date +%c): ${rankmss} The distributor asks to stop, which we do" >> $logfile
       action="stop"
     else
       echo "$(date +%c): ${rankmss} Strange message $newaction" >> $logfile
       mv $requestfile $oldrequest
     fi
   fi      
 fi
# update duration measure
 duration=$SECONDS

# force termination of the script 
 if [[ $duration -ge $endsec && "$action" == "running" ]]; then
   rendsec=$SECONDS
   runtime=`echo "($rendsec-$rstartsec)/60" | bc` 
   echo "$(date +%c): ${rankmss} We stop, but grid point $gridpoint is not yet ready (runtime $runtime minutes, expected ${expruntime})." >> $logfile
   echo "$(date +%c): ${rankmss} We stop, but grid point $gridpoint is not yet ready (runtime $runtime minutes, expected ${expruntime})." 
   echo "abort" > $requestfile
   action="wait"
 fi 
done

#------------------- post run work ---------------

echo "$(date +%c): ${rankmss} Last action is ${action}." >> $logfile
echo "$(date +%c): ${rankmss} ends." >> $logfile
echo "$(date +%c): ${rankmss} ends."
exit 0
# the script is ready
