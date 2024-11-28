#!/bin/ksh -x

## submits job to turn decadal variable outputs into sliced timeseries	##
##																		##
## 1:job name 2:script name 3: variable				 					##
##																		##
## --------------------------------------------------------------------	##

domain="FGRN055"
forcing="era055"

jobname=$1
scriptname="${2}_${domain}.sc"
vars=$3

project_name="${domain}_${forcing}"
base_dir="$SCRATCH/${project_name}"
script_dir="$PERM/IMAU-FDM/init/${domain}_${forcing}"

ts_start_year=1957
ts_end_year=2023
ave_start_year=1960
ave_end_year=1980
num_long_bands=74 #TODO: update so that files creates num_lon_bands = instead of <
cell_width=5 #TODO: update so that files slices for cell_width, not cell_width + 1

years_dir="${base_dir}/process-RACMO/years-${ts_start_year}/"
parts_dir="${years_dir}parts/"
files_dir="${base_dir}/input/timeseries-${ts_start_year}/"
ave_dir="${base_dir}/input/averages-${ts_start_year}_${ave_start_year}-${ave_end_year}/"
jobfile_dir="${base_dir}/process-RACMO/jobs/"
logfile_dir="${base_dir}/logfiles/process-RACMO/"

mkdir -p ${years_dir}
mkdir -p ${parts_dir}
mkdir -p ${files_dir}
mkdir -p ${ave_dir}
mkdir -p ${logfile_dir}
mkdir -p ${jobfile_dir}

## creates and writes initRacmoFile, then submits with sbatch			##
##																		##
## 1: variable 2: project name 3: base directory 						##
## 4: years_dir 5: files_dir 6: ave_dir									##
## 7: number of longitudinal bands 8: cell width 						##
## 9: ts start year 10: ts end year 									##
## 11: avg start year 12: avg end year 									##
##																		##
## -------------------------------------------------------------------	##

initRacmoFile="${jobfile_dir}preprocess-RACMO_job_${jobname}-${ts_start_year}"

cat <<EOS1 > ${initRacmoFile} 
#!/bin/ksh -f

#SBATCH -q nf
#SBATCH -J ${jobname}_preprocess-${ts_start_year}
#SBATCH --time=48:00:00
#SBATCH -o ${logfile_dir}${jobname}_preprocess-${ts_start_year}.log
#SBATCH --mem-per-cpu=4000mb

module load nco

echo "Start year: " ${ts_start_year}
echo "End year: " ${ts_end_year}
echo "Spinup start year: " ${ave_start_year}
echo "Spinup end year: " ${ave_end_year}
echo "Number of lon bands: " ${num_long_bands}
echo "Cell width: " ${cell_width}

cd ${script_dir}

./${scriptname} ${vars} ${project_name} ${base_dir} ${years_dir} ${files_dir} ${ave_dir} ${num_long_bands} ${cell_width} ${ts_start_year} ${ts_end_year} ${ave_start_year} ${ave_end_year}


EOS1

sbatch ${initRacmoFile}
