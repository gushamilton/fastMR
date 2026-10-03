#!/bin/bash
#SBATCH --job-name=stor-cprw --account=sscm013902 --partition=short --constraint=Cascadelake
#SBATCH --cpus-per-task=8 --mem=64G --time=01:00:00 --array=1-3
#SBATCH --output=/user/work/fh6520/showcase/storage/logs/cprw-%A_%a.out
module load languages/R/4.5.1
export STORAGE_ROOT=/user/work/fh6520/showcase/storage STORAGE_SCRIPTS=/user/work/fh6520/showcase/storage/scripts
export LIB_CS=/user/work/fh6520/showcase/final/lib4 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
Rscript $STORAGE_SCRIPTS/cprwrite_rerun.R
