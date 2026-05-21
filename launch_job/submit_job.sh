#!/bin/bash
# Non-interactive job submission — called by launch_job.sh (first run) and
# by the SLURM/offline job itself when resubmitting.
#
# Usage: submit_job.sh <settings_dir> <submission_iteration> <work_pointlist_path>

SETTINGS_DIR="$1"
SUBMISSION_ITERATION="$2"
WORK_POINTLIST_PATH="$3"

toml_get() {
    local file="$1" section="$2" key="$3"
    python3 -c "
import tomllib
with open('$file', 'rb') as f:
    d = tomllib.load(f)
print(d['$section']['$key'])
"
}

PROJECT_NAME=$(toml_get "${SETTINGS_DIR}run.toml" job project_name)
DOMAIN=$(toml_get "${SETTINGS_DIR}run.toml" job domain)
FORCING=$(toml_get "${SETTINGS_DIR}run.toml" job forcing)
RUN_TYPE=$(toml_get "${SETTINGS_DIR}run.toml" job run_type)
RELAUNCH=$(toml_get "${SETTINGS_DIR}run.toml" job relaunch)
FDM_EXECUTABLE=$(toml_get "${SETTINGS_DIR}run.toml" code-subdirectories fdm_executable)
DATA_DIR=$(toml_get "${SETTINGS_DIR}run.toml" directories data_dir)

WORK_DIR="${DATA_DIR}${PROJECT_NAME}/"
WORKEXE_DIR="${WORK_DIR}localcode/"
OUTPUT_DIR="${WORK_DIR}output/"
POINT_LOG_DIR="${WORK_DIR}logfiles/model_logfiles/"
DISTRIBUTOR_LOG_DIR="${WORK_DIR}logfiles/distributor_logfiles/"
DISTRLOG_FILE="${WORK_DIR}launch_${PROJECT_NAME}_${SUBMISSION_ITERATION}.log"

PREFIX="${DOMAIN}_${FORCING}"
NEXT_ITER=$(( SUBMISSION_ITERATION + 1 ))
NEXT_POINTLIST="${WORK_DIR}pointlist_${NEXT_ITER}.txt"

echo "${SUBMISSION_ITERATION}" > "${WORK_DIR}submission_iteration.txt"

N_POINTS=$(awk 'END{print NR}' "${WORK_POINTLIST_PATH}")

if [ "${RUN_TYPE}" == "ECMWF" ]; then

    MAX_NODES=$(toml_get "${SETTINGS_DIR}run.toml" ecmwf max_nodes)
    WALLTIME=$(toml_get "${SETTINGS_DIR}run.toml" ecmwf walltime)
    ACCOUNT=$(toml_get "${SETTINGS_DIR}run.toml" ecmwf account_no)
    MEM_PER_CPU=$(toml_get "${SETTINGS_DIR}run.toml" ecmwf memory_per_task)
    CORES_PER_NODE=128
    MAX_WORKERS=$(( MAX_NODES * CORES_PER_NODE - 1 ))  # -1 reserves a slot for the distributor
    NUM_WORKERS=$(( N_POINTS < MAX_WORKERS ? N_POINTS : MAX_WORKERS ))
    TOTAL_TASKS=$(( NUM_WORKERS + 1 ))                 # +1 for the distributor
    NUM_NODES=$(( (TOTAL_TASKS + CORES_PER_NODE - 1) / CORES_PER_NODE ))  # ceiling division

    cd "${WORK_DIR}"
    sbatch <<EOF
#!/bin/bash
#SBATCH --nodes=$NUM_NODES
#SBATCH --ntasks=$TOTAL_TASKS
#SBATCH --time=$WALLTIME
#SBATCH --account=$ACCOUNT
#SBATCH --mem-per-cpu=$MEM_PER_CPU
#SBATCH --output=$DISTRLOG_FILE

srun "${WORKEXE_DIR}${FDM_EXECUTABLE}" $WORK_POINTLIST_PATH $SETTINGS_DIR $POINT_LOG_DIR $DISTRIBUTOR_LOG_DIR

echo "\$(date +%c) Model run ended."

N_REMAINING=\$(python3 ${WORKEXE_DIR}launch_job/make_resubmit_pointlist.py \
    --output-dir $OUTPUT_DIR \
    --pointlist $WORK_POINTLIST_PATH \
    --next-pointlist $NEXT_POINTLIST \
    --prefix $PREFIX)

echo "\$(date +%c) \$N_REMAINING points remaining."

if [ \$N_REMAINING -gt 0 ] && [ \$N_REMAINING -lt $N_POINTS ] && [ "$RELAUNCH" == "yes" ]; then
    echo "\$(date +%c) Resubmitting for iteration $NEXT_ITER (\$N_REMAINING points remaining)."
    bash ${WORKEXE_DIR}launch_job/submit_job.sh $SETTINGS_DIR $NEXT_ITER $NEXT_POINTLIST
elif [ \$N_REMAINING -eq $N_POINTS ] && [ \$N_REMAINING -gt 0 ]; then
    echo "\$(date +%c) No progress made (0 points completed) — not resubmitting."
fi
EOF

    echo "$(date +%c) Job submitted to SLURM (iteration ${SUBMISSION_ITERATION}, ${N_POINTS} points, ${NUM_WORKERS} workers)."

elif [ "${RUN_TYPE}" == "offline" ]; then

    N_PROCS=$(toml_get "${SETTINGS_DIR}run.toml" offline n_procs)
    WORKERS=$(( N_PROCS - 1 ))
    WORKERS=$(( WORKERS < N_POINTS ? WORKERS : N_POINTS ))

    mpirun -n $(( WORKERS + 1 )) "${WORKEXE_DIR}${FDM_EXECUTABLE}" \
        $WORK_POINTLIST_PATH $SETTINGS_DIR $POINT_LOG_DIR $DISTRIBUTOR_LOG_DIR

    echo "$(date +%c) Model run ended."

    N_REMAINING=$(python3 "${WORKEXE_DIR}launch_job/make_resubmit_pointlist.py" \
        --output-dir "$OUTPUT_DIR" \
        --pointlist "$WORK_POINTLIST_PATH" \
        --next-pointlist "$NEXT_POINTLIST" \
        --prefix "$PREFIX")

    echo "$(date +%c) $N_REMAINING points remaining."

    if [ "$N_REMAINING" -gt 0 ] && [ "$N_REMAINING" -lt "$N_POINTS" ] && [ "$RELAUNCH" == "yes" ]; then
        echo "$(date +%c) Resubmitting for iteration $NEXT_ITER ($N_REMAINING points remaining)."
        bash "${WORKEXE_DIR}launch_job/submit_job.sh" "$SETTINGS_DIR" "$NEXT_ITER" "$NEXT_POINTLIST"
    elif [ "$N_REMAINING" -eq "$N_POINTS" ] && [ "$N_REMAINING" -gt 0 ]; then
        echo "$(date +%c) No progress made (0 points completed) — not resubmitting."
    fi

else
    echo "Run type -- ${RUN_TYPE} -- not recognized"
    exit 1
fi
