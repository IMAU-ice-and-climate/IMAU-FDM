# draft file for testing launch of new mpi distributor
# launches distributor.f90

# call fpm run to make up-to-date executable

run_type="hpc" #or offline or other option
path_to_FDM="/home/nld4814/perm/code/IMAU-FDM"
cd $path_to_FDM
if [run_type="hpc"]; then 
    ./compile_hpc.sh
elif [run_type="offline"]; then
    fpm run
else
    echo "Run type -- $run_type -- not recognized"
fi




# 
# call srun w/ new distributor (?)
# something like srun -n ${maxFDMs} --ntasks-per-node $tasks_per_node --threads-per-core 1 $homedir/npnf_inner_script.sc $myenvfile
# or srun -n ${maxFDMs} $homedir/npnf_inner_script.sc $myenvfile

# need to get lat lons (or read lat lons into fortran directly?)

# starting the model
# $exe_id $usern $cpoint $domain $filename_part1 $project_name $restart_type &> ${log_fname} 