#!/bin/bash
#SBATCH -A sscm013902 -p short --constraint=Cascadelake -c 8 --mem=32G -t 00:20:00 -o /user/work/fh6520/showcase/storage/logs/chk-%j.out
module load languages/R/4.5.1 bcftools/1.19-openblas-5yp2 htslib/1.19.1-jbbb
export STORAGE_SCRIPTS=/user/work/fh6520/showcase/storage/scripts
Rscript $STORAGE_SCRIPTS/chk_harness.R
