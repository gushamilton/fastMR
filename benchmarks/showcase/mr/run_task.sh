#!/bin/bash -l
#SBATCH -A sscm013902
#SBATCH -p short
#SBATCH -J mr_bench
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH --array=1-3
#SBATCH -o /user/work/fh6520/showcase/mr/logs/%x_%A_%a.out
# MR-compute suite: one replicate per array task. Arms are paired on the same node in randomised order.
# Env: LIB_FASTMR, LIB_TSMR, FASTMR_SHA, CLUMP_REF (bed prefix), CLUMP_PREF (pgen prefix), DRY=1, CAP (seconds)
module load languages/R/4.5.1 apps/plink2/2.00a68LM apps/plink1.9/1.90-b77
W=/user/work/fh6520/showcase/mr
REP=${SLURM_ARRAY_TASK_ID:-1}
DRY=${DRY:-0}; SUF=""; [ "$DRY" = 1 ] && SUF="_dry"
TD=$W/work/rep$REP$SUF; rm -rf $TD; mkdir -p $TD $W/results $W/logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 TMPDIR=$TD
export LIB_FASTMR=${LIB_FASTMR:-/user/work/fh6520/regress-audit/lib-fm-final}
export LIB_TSMR=${LIB_TSMR:-/user/work/fh6520/showcase/lib-tsmr}
export FASTMR_SHA=${FASTMR_SHA:-unknown}
export CLUMP_REF=${CLUMP_REF:-/user/work/fh6520/showcase/prep/ref/EUR_maf01}
export CLUMP_PREF=${CLUMP_PREF:-/user/work/fh6520/showcase/prep/ref/EUR_maf01}
export PLINK1=$(which plink) PLINK2=$(which plink2)
if [ "$SLURM_CPUS_PER_TASK" != "8" ]; then echo "SLURM_CPUS_PER_TASK != 8"; exit 2; fi
if [ ! -f "$CLUMP_PREF.pgen" ] && [ "$DRY" = 1 ]; then   # dry-run only: build a pgen from the bed once
  mkdir -p $W/dryref; $PLINK2 --bfile $CLUMP_REF --make-pgen --out $W/dryref/ref --threads 4 > $TD/mkpgen.log 2>&1; export CLUMP_PREF=$W/dryref/ref; fi
if [ ! -f "$LIB_TSMR/TwoSampleMR/DESCRIPTION" ]; then echo "NOTE: LIB_TSMR missing, using site TwoSampleMR"; export LIB_TSMR=""; fi
CPU=$(grep -m1 'model name' /proc/cpuinfo | sed 's/^[^:]*: //')
{ echo "host=$(hostname) cpu=$CPU nproc=$(env -u OMP_NUM_THREADS nproc) slurm_cpus=$SLURM_CPUS_PER_TASK slurm_job=$SLURM_JOB_ID"; R --version | head -1
  $PLINK1 --version; $PLINK2 --version; echo "LIB_FASTMR=$LIB_FASTMR FASTMR_SHA=$FASTMR_SHA LIB_TSMR=$LIB_TSMR CLUMP_REF=$CLUMP_REF CLUMP_PREF=$CLUMP_PREF"
  Rscript -e 'if (nzchar(Sys.getenv("LIB_TSMR"))) .libPaths(c(Sys.getenv("LIB_TSMR"), .libPaths())); for (p in c("TwoSampleMR","ieugwasr")) { d <- packageDescription(p); cat(p, d$Version, "sha=", if (!is.null(d$RemoteSha)) d$RemoteSha else "NA", "from", dirname(system.file(package=p)), "\n") }; .libPaths(c(Sys.getenv("LIB_FASTMR"), .libPaths())); d <- packageDescription("fastMR"); cat("fastMR", d$Version, "sha=", if (!is.null(d$RemoteSha)) d$RemoteSha else "NA", "\n"); print(sessionInfo()$BLAS)'
} > $W/results/meta_rep$REP$SUF.txt 2>&1
Rscript $W/driver.R --rep=$REP --out=$W/results --work=$TD --dry=$DRY --cap=${CAP:-1500}
