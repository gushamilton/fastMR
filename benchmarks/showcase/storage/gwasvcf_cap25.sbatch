#!/bin/bash
#SBATCH --job-name=stor-gvcf --account=sscm013902 --partition=short --constraint=Cascadelake
#SBATCH --cpus-per-task=8 --mem=64G --time=01:25:00
#SBATCH --output=/user/work/fh6520/showcase/storage/logs/gwasvcf-cap25-%j.out
module load languages/R/4.5.1 bcftools/1.19-openblas-5yp2 htslib/1.19.1-jbbb
export STORAGE_ROOT=/user/work/fh6520/showcase/storage STORAGE_SCRIPTS=/user/work/fh6520/showcase/storage/scripts LIB_CS=/user/work/fh6520/showcase/final/lib
Rscript $STORAGE_SCRIPTS/gwasvcf_cap25.R
