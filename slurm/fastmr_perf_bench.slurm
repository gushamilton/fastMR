#!/bin/bash
# fastMR baseline-vs-candidate performance array for BluePebble.
#
# Each array task = one (SHA, replicate) cell.  Tasks are spread across nodes
# by Slurm; hostname/CPU are recorded per row so node effects can be modelled.
#
# Submit (do not run on the login node):
#   export BASELINE_SHA=9f573d4 CANDIDATE_SHA=<sha> REPLICATES=5
#   export SCENARIOS="grid_ivw,fast_mr_groups,sparse_ivw,steiger_filtering,clump_mock,compressed_pairwise"
#   sbatch --array=0-$((2*REPLICATES-1)) --export=ALL slurm/fastmr_perf_bench.slurm
# Even tasks run BASELINE_SHA, odd tasks CANDIDATE_SHA (interleaved).
# Results: ${BENCH_ROOT}/results/<run-id>/task_<id>.csv  (merge with rbind).

#SBATCH --job-name=fastmr-perf
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --account=sscm013902
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail
module purge
module load languages/R/4.5.1

REPO="${REPO:-/user/work/fh6520/fastMR}"           # git clone with both SHAs fetched
BENCH_ROOT="${BENCH_ROOT:-/user/work/fh6520/fastmr-bench}"
REPLICATES="${REPLICATES:-5}"
SCENARIOS="${SCENARIOS:-grid_ivw,grid_ivw_parquet,grid_mixed,fast_mr_groups,fast_mr_groups_mixed,sparse_ivw,sparse_ivw_mask,steiger_filtering,clump_mock,compressed_pairwise}"
SCALE="${SCALE:-full}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
RUN_ID="${RUN_ID:-${SLURM_ARRAY_JOB_ID:-manual}}"
TASK="${SLURM_ARRAY_TASK_ID:-0}"

# Even tasks run BASELINE_SHA, odd tasks CANDIDATE_SHA, so both SHAs are
# interleaved in time and spread over the same pool of nodes.
if (( TASK % 2 == 0 )); then SHA="${BASELINE_SHA}"; LABEL=baseline
else SHA="${CANDIDATE_SHA}"; LABEL=candidate; fi
REP=$((TASK / 2 + 1))
# Optional extra library (e.g. a CompreSSoR build) for compressed scenarios.
if [[ -n "${EXTRA_R_LIBS:-}" ]]; then export R_LIBS="${EXTRA_R_LIBS}${R_LIBS:+:${R_LIBS}}"; fi

mkdir -p "${BENCH_ROOT}/builds" "${BENCH_ROOT}/results/${RUN_ID}" logs
SRC="${BENCH_ROOT}/builds/src-${SHA}"
LIB="${BENCH_ROOT}/builds/lib-${SHA}"
# Build each SHA once; concurrent array tasks wait on the lock.
(
  flock 9
  if [[ ! -f "${LIB}/.built" ]]; then
    rm -rf "${LIB}"; mkdir -p "${LIB}"
    # Compute nodes may lack git: export ${SRC} on the login node beforehand
    # (git -C "$REPO" archive "$SHA" | tar -x -C "$SRC"); clone only as fallback.
    if [[ ! -d "${SRC}" ]]; then
      git -C "${REPO}" archive "${SHA}" | (mkdir -p "${SRC}" && tar -x -C "${SRC}")
    fi
    R CMD INSTALL --no-test-load --library="${LIB}" "${SRC}"
    touch "${LIB}/.built"
  fi
) 9>"${BENCH_ROOT}/builds/.lock-${SHA}"

# The benchmark script itself is taken from the CANDIDATE checkout's
# benchmarks/perf/ (or BENCH_SCRIPT) so both SHAs run identical harness code.
BENCH_SCRIPT="${BENCH_SCRIPT:-${BENCH_ROOT}/fastmr_perf_bench.R}"
OUT="${BENCH_ROOT}/results/${RUN_ID}/task_${TASK}.csv"
echo "hostname=$(hostname) sha=${SHA} label=${LABEL} rep=${REP} threads=${THREADS}"

IFS=',' read -ra SCN <<< "${SCENARIOS}"
for s in "${SCN[@]}"; do
  TMP="$(mktemp)"; TV="$(mktemp)"
  # One process per scenario so /usr/bin/time -v peak RSS is per scenario.
  /usr/bin/time -v -o "${TV}" Rscript "${BENCH_SCRIPT}" --scenario="${s}" --scale="${SCALE}" \
    --label="${LABEL}" --sha="${SHA}" --lib="${LIB}" --replicate="${REP}" \
    --threads="${THREADS}" --output="${TMP}" || echo "scenario ${s} failed" >&2
  RSS_KB="$(awk -F: '/Maximum resident set size/ {gsub(/ /,"",$2); print $2}' "${TV}")"
  if [[ -s "${TMP}" ]]; then
    Rscript -e 'a <- commandArgs(TRUE); x <- read.csv(a[1], stringsAsFactors = FALSE)
      x$process_peak_rss_mb <- as.numeric(a[2]) / 1024; x$slurm_job <- a[4]
      write.table(x, a[3], sep = ",", row.names = FALSE, append = file.exists(a[3]),
                  col.names = !file.exists(a[3]))' "${TMP}" "${RSS_KB:-NA}" "${OUT}" \
      "${SLURM_ARRAY_JOB_ID:-NA}_${TASK}"
  fi
  rm -f "${TMP}" "${TV}"
done
echo "results: ${OUT}"
echo "sacct -X -j ${SLURM_ARRAY_JOB_ID:-?}_${TASK} --format=JobID,NodeList,State,Elapsed,TotalCPU,MaxRSS"
