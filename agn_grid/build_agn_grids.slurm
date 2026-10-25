#!/bin/bash
# Build every AGN grid with build_agn_grids.py --apply on one RM node: the dense core, the
# BLR core and the variation grids as three parallel processes (each builder process is
# serial). The dense core holds ~10 GB in memory while building and took ~1 h (21 Oct 2026).
#
# Run the dry run first and check each hole's stencil. Existing build directories stop the
# builder: move them aside (or pass --rebuild via BUILD_ARGS, which deletes them).
# Submit from the top of smock-recipes, which holds cloudy_runs/ (GRIDS picks a subset, e.g. GRIDS=dense_core):
#     sbatch agn_grid/build_agn_grids.slurm
#
#SBATCH --job-name=agn_grid_build
# The account and partition are the ones used on PSC's Vera: change them for your cluster.
#SBATCH --account=phy200026p
#SBATCH --partition=RM
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --time=04:00:00
#SBATCH --output=cloudy_runs/build_%j.log

set -euo pipefail

PYTHON=${PYTHON:-python}   # the active env's Python (smock's env, with syncretize)
# sbatch runs a spooled copy of this script, so the builder is found from the submit directory
BUILDER=${BUILDER:-agn_grid/build_agn_grids.py}
BUILD_ARGS=${BUILD_ARGS:-}
# substrings of grid names: the two cores, then every variation grid ("fiducial")
GRIDS=${GRIDS:-dense_core blr_core fiducial}

pids=()
for grids in $GRIDS; do
    "$PYTHON" "$BUILDER" --grids "$grids" --apply $BUILD_ARGS > "cloudy_runs/build_${SLURM_JOB_ID}_${grids}.log" 2>&1 &
    pids+=($!)
done
status=0
for pid in "${pids[@]}"; do
    wait "$pid" || status=1
done
exit $status
