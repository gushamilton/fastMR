#!/bin/bash
#SBATCH --job-name=fmr-ns-final-equiv
#SBATCH --partition=test
#SBATCH --time=00:58:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --account=sscm013902
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
set -euo pipefail
module purge; module load languages/R/4.5.1
# Layout: $WORK/src-main and $WORK/src-new hold `git archive` exports of the
# two builds (BluePebble compute nodes have no git); results go under $WORK.
WORK=${WORK:-/user/work/fh6520/fastmr-native-speed}
cd "$WORK"; S=src-new/benchmarks/native-speed; mkdir -p logs equiv results
hostname; lscpu | grep "Model name"; uptime
for v in main new; do
  rm -rf lib-$v; mkdir -p lib-$v
  R CMD INSTALL --no-test-load --library=lib-$v src-$v > logs/build-$v.log 2>&1 && echo built $v
done
grep -c fopenmp logs/build-new.log || true
Rscript $S/equiv-run.R lib-main equiv/bp-main.rds
Rscript $S/equiv-run.R lib-new equiv/bp-new.rds
FASTMR_HULL="Inf,1" Rscript $S/equiv-run.R lib-new equiv/bp-new-hullInf-rec.rds
FASTMR_HULL="Inf,0" Rscript $S/equiv-run.R lib-new equiv/bp-new-hullInf-exact.rds
FASTMR_OVERLAP=0 Rscript $S/equiv-run.R lib-new equiv/bp-new-sync.rds
Rscript $S/equiv-compare.R equiv/bp-main.rds equiv/bp-new.rds equiv/bp-new-hullInf-rec.rds equiv/bp-new-hullInf-exact.rds equiv/bp-new-sync.rds
Rscript $S/audit.R src-new
