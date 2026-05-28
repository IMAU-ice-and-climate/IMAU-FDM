### LAUNCH FDM MODEL RUN ###
############################
#
# To run: update SETTINGS_FOLDER with your settings folder name
# Then ./launch_job.sh in the terminal.
#
# Created by: Elizabeth Case, May 2026
############################




SETTINGS_FOLDER="FGRN055" # For TOML files, create a new folder for each project that updates any settings and speciy via $PATHS_TO_TOML




#python helper function to read in toml files
toml_get() {
    local file="$1" section="$2" key="$3"
    python3 -c "
import tomllib
with open('$file', 'rb') as f:
    d = tomllib.load(f)
print(d['$section']['$key'])
"
}


LAUNCH_DIR="$(pwd)/"
FDM_BASE="$(dirname "$(pwd)")/"
SETTINGS_DIR="${FDM_BASE}/settings/${SETTINGS_FOLDER}/"

PROJECT_NAME=$(toml_get "${SETTINGS_DIR}run.toml" job project_name)
POINTLIST_NAME=$(toml_get "${SETTINGS_DIR}run.toml" job pointlist_name)
DOMAIN=$(toml_get "${SETTINGS_DIR}run.toml" job domain)
RUN_TYPE=$(toml_get "${SETTINGS_DIR}run.toml" job run_type)

#TKTKTK figure out how to get job to resubmit
SUBMISSION_ITERATION=$(toml_get "${SETTINGS_DIR}run.toml" job submission_iteration) #needed for restarting if job doesn't finish

FDM_EXECUTABLE=$(toml_get "${SETTINGS_DIR}run.toml" code-subdirectories fdm_executable)

POINTLIST_PATH="${LAUNCH_DIR}pointlists/${POINTLIST_NAME}.txt"
FDM_EXECECUTABLE_PATH="${FDM_BASE}bin/${FDM_EXECUTABLE}"

# WORKING DIRECTORY PATHS

DATA_DIR=$(toml_get "${SETTINGS_DIR}run.toml" directories data_dir)
WORK_DIR="${DATA_DIR}${PROJECT_NAME}/"
POINT_LOG_DIR="${WORK_DIR}logfiles/model_logfiles/"
DISTRIBUTOR_LOG_DIR="${WORK_DIR}logfiles/distributor_logfiles/"
WORKEXE_DIR="${WORK_DIR}localcode/"
RESTART_DIR="${WORK_DIR}restart/" 
OUTPUT_DIR="${WORK_DIR}output/"
POSTPROCESS_DIR="${WORK_DIR}post-process/"

#also gets updated when submission iteration >1
WORK_POINTLIST_PATH="${WORK_DIR}pointlist_${SUBMISSION_ITERATION}.txt"

DISTRLOG_FILE="${WORK_DIR}launch_${PROJECT_NAME}_${SUBMISSION_ITERATION}.log"

RECOMPILE=$(toml_get "${SETTINGS_DIR}run.toml" job recompile)

# compiles excecutable again to make sure everything is up to date if recompile = true
if [ "$RECOMPILE" == "True" ]; then
  if [ ${RUN_TYPE} == "ECMWF" ]; then
      bash "${FDM_BASE}/compile_hpc.sh" || { echo "ERROR: Compilation failed. Exiting."; exit 1; }
  elif [ ${RUN_TYPE} == "offline" ]; then
      fpm run || { echo "ERROR: Compilation failed. Exiting."; exit 1; }
  else
      echo "ERROR: Run type -- ${RUN_TYPE} -- not recognized"
      exit 1
  fi
fi

echo "Continuing will submit the job."
echo "Working directory: $WORK_DIR"
echo "Pointlist path: $POINTLIST_PATH"

read -p "Check paths. Want to continue? (y/n)" -n 1 -r
echo " "

if [[ $REPLY =~ ^[Yy]$ ]]; then

  # -------- start of the launching procedure ------- #

  # -------- create working directory ------- #
  echo "Make working directory"
  # prepare launch
  mkdir -p $WORK_DIR
  mkdir -p $WORKEXE_DIR
  mkdir -p "${WORKEXE_DIR}job"
  mkdir -p $OUTPUT_DIR
  mkdir -p $POSTPROCESS_DIR
  mkdir -p "${RESTART_DIR}run"
  mkdir -p "${RESTART_DIR}spinup"
  mkdir -p $DISTRIBUTOR_LOG_DIR
  mkdir -p $POINT_LOG_DIR

  if [[ ! -r $POINTLIST_PATH ]]; then
    echo "The pointlist is missing!"
    exit 1
  fi  
  
  cp ${POINTLIST_PATH} ${WORK_POINTLIST_PATH}
  #cp ${POINTLIST_PATH} ${WORK_DIR}pointlist.txt # do we want to keep an active pointlist as pointlist.txt?

  cp "${FDM_EXECECUTABLE_PATH}" "${WORKEXE_DIR}${FDM_EXECUTABLE}" #copies the model executable
  cp -r "${FDM_BASE}source" "$WORKEXE_DIR"       # copies the model sourcecode
  cp -r "${FDM_BASE}settings" "$WORKEXE_DIR"     # copies the settings files
  cp -r "${FDM_BASE}launch_job" "$WORKEXE_DIR"   # copies submission scripts for resubmission

  echo "Starting the model."
  bash "${FDM_BASE}launch_job/submit_job.sh" "$SETTINGS_DIR" "$SUBMISSION_ITERATION" "$WORK_POINTLIST_PATH"

fi
exit 0  
