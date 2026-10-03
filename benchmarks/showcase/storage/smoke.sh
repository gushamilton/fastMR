#!/bin/bash
module load languages/R/4.5.1 bcftools/1.19-openblas-5yp2 htslib/1.19.1-jbbb
export STORAGE_ROOT=/user/work/fh6520/showcase/storage/smoke STORAGE_SCRIPTS=/user/work/fh6520/showcase/storage/scripts N_ROWS=300000
mkdir -p $STORAGE_ROOT/{logs,results,stores,work}
hostname
Rscript $STORAGE_SCRIPTS/prep.R
export SLURM_CPUS_PER_TASK=8 REP=0 RESULT_TAG=smoke LOOKUP_SIZES=25,1000,100000 OP_CAP=120
Rscript $STORAGE_SCRIPTS/run_rep.R
