# draft file for testing launch of new mpi distributor
# import file paths from toml
# create file structure where data output should go
# if running online, import sbatch info 
# launch distributor.f90

# load the yq module to read toml files

PATHS_TOML="${path_to_FDM}/settings/paths.toml"
RUN_TYPE="ECMWF" #or offline, or other options as added

if [ "$RUN_TYPE" = "ECMWF" ]; then
    module load yq/4.45.2
else
    if ! command -v tomlq &> /dev/null; then
        echo "Error: tomlq not found. Install yq: https://github.com/mikefarah/yq"
        exit 1
    fi
fi


#read in file structure from toml file

# [directories]
code_dir=$(tomlq -r '.directories.code_dir' "$PATHS_TOML")
input_dir=$(tomlq -r '.directories.input_dir' "$PATHS_TOML")
output_dir=$(tomlq -r '.directories.output_dir' "$PATHS_TOML")

# [code-subdirectories]
ref_dir=$(tomlq -r '."code-subdirectories".ref_dir' "$PATHS_TOML")
source_dir=$(tomlq -r '."code-subdirectories".source' "$PATHS_TOML")
executable=$(tomlq -r '."code-subdirectories".executable' "$PATHS_TOML")

# Build full paths
path_to_executable="${code_dir}${source_dir}${executable}"



# compiles

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