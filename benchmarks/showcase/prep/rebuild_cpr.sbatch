#!/bin/bash -l
#SBATCH -A sscm013902
#SBATCH -p test
#SBATCH -J cpr-rebuild
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 00:45:00
#SBATCH --array=1-4
#SBATCH -o /user/work/fh6520/showcase/prep/logs/cpr-rebuild-%A_%a.out
set -euo pipefail
module load languages/R/4.5.1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 LIB_CS=/user/work/fh6520/showcase/final/lib
export TMPDIR=/user/work/fh6520/showcase/prep/tmp_rebuild_$SLURM_ARRAY_TASK_ID; mkdir -p $TMPDIR
Rscript /user/work/fh6520/showcase/prep/rebuild_cpr.R
rm -rf $TMPDIR
