#!/bin/bash
#SBATCH --job-name=bnf_run_compare
#SBATCH --account=cli120
#SBATCH --partition=batch_all
#SBATCH --array=0-15
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/wolf2/arm/cli120/scratch/hengxiao80/bnf_run_compare_%A_%a.out
#
# Multi-run comparison frames / time series (bnf_wrf_eval.run_compare).
# Task k renders worklist[k::N]; restartable (existing frames skipped).
# --cpus-per-task=16 is NOT for compute (each task is serial): it caps packing at 8 tasks per
# 128-core node. With 2 cores, SLURM packed all 64 tasks of four arrays onto one node, whose
# GPFS link then starved them -- a ~8 s frame took ~9 min (2026-10-12).
#   sbatch scripts/run_compare.sbatch            # default: panels + every pair field
#   WHICH=timeseries,pair_timeseries sbatch --array=0-0 scripts/run_compare.sbatch   # serial, ~25 min
# CONFIG= / OUTPUT_ROOT= / EXTRA= (extra run_compare flags) override the defaults.
set -u
REPO=/gpfs/wolf2/arm/cli120/proj-shared/hengxiao80/bnf_wrf_eval
PY=$REPO/.venv/bin/python3
WHICH=${WHICH:-panels_tb,panels_refl,pair_tb,pair_refl,pair_precip,pair_clwp,pair_rwp,pair_ciwp,pair_swp,pair_gwp}
CONFIG=${CONFIG:-$REPO/configs/da_20250520.json}
OUTPUT_ROOT=${OUTPUT_ROOT:-$REPO/outputs}
EXTRA=${EXTRA:-}   # e.g. EXTRA=--overwrite
module load intel/20.0.4 openmpi/4.1.4 netcdf-c/4.8.1 netcdf-fortran/4.5.3
export PYTHONPATH="${SRC:-$REPO/src}${PYTHONPATH:+:$PYTHONPATH}" MPLBACKEND=Agg OMP_NUM_THREADS=1
cd "$REPO"
echo "host=$(hostname) date=$(date -Is) which=${WHICH}"
rc=1
for attempt in 1 2 3; do
    "$PY" -m bnf_wrf_eval.run_compare --config "$CONFIG" --which "$WHICH" --output-root "$OUTPUT_ROOT" $EXTRA
    rc=$?; [ "$rc" -eq 0 ] && break
    echo "attempt ${attempt} exited ${rc}; re-running"; sleep 5
done
exit "$rc"
