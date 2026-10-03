#!/bin/bash
#SBATCH --job-name=stor-diag
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --constraint=Cascadelake
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/user/work/fh6520/showcase/storage/logs/diag-%j.out
module load languages/R/4.5.1 bcftools/1.19-openblas-5yp2 htslib/1.19.1-jbbb
export STORAGE_ROOT=/user/work/fh6520/showcase/storage STORAGE_SCRIPTS=/user/work/fh6520/showcase/storage/scripts
hostname
Rscript $STORAGE_SCRIPTS/make_lists.R
Rscript $STORAGE_SCRIPTS/diag_cpr_read.R
