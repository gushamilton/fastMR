#!/bin/bash -l
#SBATCH -A sscm013902 -p short --constraint=Cascadelake -J stor-v7prep -c 8 --mem=64G -t 01:00:00
#SBATCH -o /user/work/fh6520/showcase/storage/logs/prep-v7-%j.out
module load languages/R/4.5.1 bcftools/1.19-openblas-5yp2 htslib/1.19.1-jbbb
export STORAGE_ROOT=/user/work/fh6520/showcase/storage STORAGE_SCRIPTS=/user/work/fh6520/showcase/storage/scripts LIB_CS=/user/work/fh6520/showcase/final/lib6
mkdir -p $STORAGE_ROOT/{logs,results,stores,work}
Rscript $STORAGE_SCRIPTS/make_lists.R && PREP_FORMATS=cpr Rscript $STORAGE_SCRIPTS/prep.R
