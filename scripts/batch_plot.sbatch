#!/bin/bash
#SBATCH --job-name=bnf_batch_plot
#SBATCH --account=cli120
#SBATCH --partition=batch_all
#SBATCH --array=0-63%32
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=6G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/wolf2/arm/cli120/scratch/hengxiao80/bnf_batch_plot_%A_%a.out
#
# Render every frame of the comparison and/or single-panel plot series
# (see bnf_wrf_eval.batch_plot). Task k renders worklist[k::N] (strided --
# an even mix of plot types / case days / times per task). Restartable:
# frames already on disk are skipped, so resubmitting the whole array
# after a timeout or node failure just fills the gaps.
#
# --cpus-per-task=16 is NOT for compute (each task is serial): it caps packing at 8 tasks per
# 128-core node. With 2 cores, SLURM packed the whole array onto one node, whose GPFS link
# then starved every task reading 3-D fields -- e.g. the d3 reflectivity frames (job 1051284)
# took ~410 s each vs. ~16 s alone (2026-10-12). 2-D-only frames (OLR Tb) were unaffected.
#
# Pre-download the GOES/MRMS observation files FIRST so tasks don't all hit
# S3 (and race on the same missing path):
#   .venv/bin/python scripts/download_obs.py --hourly              # comparisons
#   .venv/bin/python scripts/download_obs.py --which goes,mrms --cadence-mrms 15
#
# Submit (from the main checkout):
#   WHICH=comparison sbatch scripts/batch_plot.sbatch
#   WHICH=single     sbatch scripts/batch_plot.sbatch
# Dry run (print one task's work list, render nothing):
#   .venv/bin/python -m bnf_wrf_eval.batch_plot --which comparison --list \
#     --task-index 0 --task-count 64

set -u

REPO=/gpfs/wolf2/arm/cli120/proj-shared/hengxiao80/bnf_wrf_eval
PY=$REPO/.venv/bin/python3
SRC=${SRC:-$REPO/src}

WHICH=${WHICH:-all}
OUTPUT_ROOT=${OUTPUT_ROOT:-$REPO/outputs}
EXTRA=${EXTRA:-}   # e.g. EXTRA="--start 2025-05-02T12:00 --end 2025-05-03T00:00"

module load intel/20.0.4 openmpi/4.1.4 netcdf-c/4.8.1 netcdf-fortran/4.5.3

export PYTHONPATH="$SRC${PYTHONPATH:+:$PYTHONPATH}"
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1

cd "$REPO"

ARRAY_COUNT=$(( ${SLURM_ARRAY_TASK_MAX:-0} - ${SLURM_ARRAY_TASK_MIN:-0} + 1 ))
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

echo "host=$(hostname) date=$(date -Is) task=${TASK_ID}/${ARRAY_COUNT} which=${WHICH}"

rc=1
for attempt in 1 2 3; do
    "$PY" -m bnf_wrf_eval.batch_plot \
        --which "$WHICH" \
        --output-root "$OUTPUT_ROOT" \
        --retry 1 \
        --task-index "$TASK_ID" \
        --task-count "$ARRAY_COUNT" \
        $EXTRA
    rc=$?
    [ "$rc" -eq 0 ] && break
    echo "task ${TASK_ID} attempt ${attempt} exited ${rc}; re-running (existing frames skipped)"
    sleep 5
done

echo "task ${TASK_ID} final rc=${rc} ($(date -Is))"
exit "$rc"
