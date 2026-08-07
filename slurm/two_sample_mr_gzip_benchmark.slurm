#!/bin/bash
# Bounded ordinary two-sample MR benchmark. This wrapper is intentionally
# disjoint from the production all-by-all runner and result tree.
#
# It discovers only moderate TSV/TXT/CSV.gz files under DATA_ROOT unless
# EXPOSURE_GZ and OUTCOME_GZ are supplied explicitly. It excludes tar/fat/huge
# paths and files outside the configured size window.

#SBATCH --job-name=fastmr-two-sample-gzip
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --account=sscm013902
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
module purge
module load languages/R/4.5.1
cd /user/work/fh6520/fastMR-test-goal
mkdir -p logs results/two_sample_mr_gzip_benchmark

export OUTPUT_DIR="${OUTPUT_DIR:-results/two_sample_mr_gzip_benchmark}"
export DATA_ROOT="${DATA_ROOT:-/user/work/fh6520}"
# 50 MiB to 750 MiB: large enough to be a real gzip I/O test, but bounded so
# two full decompressed tables remain comfortably below the 16 GiB allocation.
export MIN_BYTES="${MIN_BYTES:-52428800}"
export MAX_BYTES="${MAX_BYTES:-786432000}"

echo "hostname=$(hostname)"
echo "DATA_ROOT=${DATA_ROOT}"
echo "MIN_BYTES=${MIN_BYTES} MAX_BYTES=${MAX_BYTES}"
echo "EXPOSURE_GZ=${EXPOSURE_GZ:-<auto>}"
echo "OUTCOME_GZ=${OUTCOME_GZ:-<auto>}"

# /usr/bin/time provides the authoritative process peak RSS. Slurm sacct
# should be used after completion as the independent job-level check.
/usr/bin/time -v -o "${OUTPUT_DIR}/time_v.txt" \
  Rscript benchmarks/two_sample_mr_gzip_benchmark.R \
    --output-dir="${OUTPUT_DIR}" \
    ${EXPOSURE_GZ:+--exposure="${EXPOSURE_GZ}"} \
    ${OUTCOME_GZ:+--outcome="${OUTCOME_GZ}"} \
    --data-root="${DATA_ROOT}" \
    --min-bytes="${MIN_BYTES}" \
    --max-bytes="${MAX_BYTES}"

echo "benchmark output: ${OUTPUT_DIR}"
echo "post-run exact Slurm command:"
echo "sacct -X -j ${SLURM_JOB_ID} --format=JobID,State,Elapsed,TotalCPU,MaxRSS,MaxVMSize,ExitCode"
