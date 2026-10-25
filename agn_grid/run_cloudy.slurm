#!/bin/bash
# Run CLOUDY on every model in a run list, one RM node per array task. Used for the timing
# pilot (the defaults below), production and every rerun.
#
# One CLOUDY process per core (CLOUDY is single-threaded), each with a walltime cap
# (RUN_CAP). Array task i of n takes every n-th line of the run list starting at line i,
# so slow corners spread across nodes; 3 tasks x 28 cores fit the 70 pilot runs in one wave.
# Each task writes its runs' wall times and exit status to <TIMING>_<task>.tsv
# (exit 124 = hit the cap). Afterwards, parse the outputs with:
#     python agn_grid/parse_cloudy_runs.py cloudy_runs/pilot_runs.txt \
#         --timing cloudy_runs/pilot_timing_*.tsv --csv cloudy_runs/pilot_runs.csv
#
# Generate the inputs first (make_incident_grid.py, then make_cloudy_inputs.py pilot --apply).
# Submit from the top of smock-recipes, which holds cloudy_runs/:  sbatch agn_grid/run_cloudy.slurm
# (one node for a short list:  sbatch --array=0 ...)
#
# The production grids (~21k runs) use the same script with overrides. Throttle the
# array to <= 10 nodes at a time (%10) so we never take over the 32-node RM partition:
#     RUN_LIST=cloudy_runs/production_runs.txt TIMING=cloudy_runs/production_timing \
#     sbatch --job-name=cloudy_prod --array=0-39%10 --time=20:00:00 \
#         --output=cloudy_runs/production_%A_%a.log agn_grid/run_cloudy.slurm
#
#SBATCH --job-name=cloudy_pilot
# The account and partition are the ones used on PSC's Vera: change them for your cluster.
#SBATCH --account=phy200026p
#SBATCH --partition=RM
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --array=0-2
#SBATCH --time=10:30:00
#SBATCH --output=cloudy_runs/pilot_%A_%a.log

set -euo pipefail

# Per-machine defaults, overridable from the environment.
CLOUDY_DIR=${CLOUDY_DIR:?set CLOUDY_DIR to your CLOUDY c23.01 folder (it holds source/cloudy.exe and data/)}
RUN_LIST=${RUN_LIST:-cloudy_runs/pilot_runs.txt}
TIMING=${TIMING:-cloudy_runs/pilot_timing}   # prefix; "_<task>.tsv" is appended
RUN_CAP=${RUN_CAP:-10h}          # per-run walltime cap (one wave, inside the job limit)
NPROC=${SLURM_NTASKS:-28}
TASK=${SLURM_ARRAY_TASK_ID:-0}
NTASK=${SLURM_ARRAY_TASK_COUNT:-1}

# CLOUDY finds its data files here; ./ lets it find linelist.dat in each model directory.
export CLOUDY_DATA_PATH="${CLOUDY_DIR}/data/:./"
export CLOUDY_EXE="${CLOUDY_DIR}/source/cloudy.exe"
TIMING_FILE=$(realpath "${TIMING}_${TASK}.tsv")
export RUN_CAP TIMING_FILE

echo -e "model_dir\tindex\twall_s\texit_code" > "$TIMING_FILE"

run_one() {
    local dir=$1 idx=$2 t0 rc
    t0=$(date +%s)
    cd "$dir"
    rc=0
    timeout "$RUN_CAP" "$CLOUDY_EXE" -r "$idx" || rc=$?
    echo -e "${dir}\t${idx}\t$(( $(date +%s) - t0 ))\t${rc}" >> "$TIMING_FILE"
}
export -f run_one

# This task's share of the run list (each line is "<model_dir> <index>").
MY_RUNS=$(awk -v t="$TASK" -v n="$NTASK" 'NF && (NR - 1) % n == t' "$RUN_LIST")
echo "$(date) task $TASK/$NTASK on $(hostname): $(echo "$MY_RUNS" | wc -l) runs on $NPROC cores, cap $RUN_CAP each"
echo "$MY_RUNS" | xargs -P "$NPROC" -L 1 bash -c 'run_one "$0" "$1"'
echo "$(date) done"
