#!/bin/bash
#SBATCH --job-name=fmr-ns-bench
#SBATCH --partition=test
#SBATCH --time=00:58:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --account=sscm013902
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
module purge; module load languages/R/4.5.1
# Layout: $WORK/src-main and $WORK/src-new hold `git archive` exports of the
# two builds (BluePebble compute nodes have no git); results go under $WORK.
WORK=${WORK:-/user/work/fh6520/fastmr-native-speed}
cd "$WORK"; S=src-new/benchmarks/native-speed; mkdir -p logs equiv results
hostname; lscpu | grep -E "Model name|^CPU\(s\)|Thread|Core|Socket"; uptime; nproc
tag=${TAG:-$SLURM_JOB_ID}
# Whole node for quiet timings, but pin every run to 4 physical cores x 2
# hyperthreads (8 logical CPUs), the shape of an 8-CPU Slurm allocation.
CPUSET=$(for c in 0 1 2 3; do cat /sys/devices/system/cpu/cpu$c/topology/thread_siblings_list; done | paste -sd, -)
echo "pinned to CPUs: $CPUSET"
run() { taskset -c "$CPUSET" "$@"; }
mkdir -p results/$tag
# Interleave the two builds so node drift hits both equally.
for rep in 1 2; do
  for v in main new; do
    run Rscript $S/per-method.R lib-$v $v results/$tag/per-method-$v-$rep.csv > /dev/null
    run Rscript $S/batch.R lib-$v $v 8 results/$tag/batch8-$v-$rep.csv > /dev/null
  done
done
# Attribution: the new build with the hull path off (RNG overlap only).
FASTMR_HULL="0,1" run Rscript $S/batch.R lib-new new-nohull 8 results/$tag/batch8-new-nohull.csv > /dev/null
MODE_KS=3,10,30,100,300,1000 run Rscript $S/mode_paths.R src-new > results/$tag/mode-paths.txt 2>&1
echo done
