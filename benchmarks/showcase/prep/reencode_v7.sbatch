#!/bin/bash -l
#SBATCH -A sscm013902 -p short --constraint=Cascadelake -J cpr-v7 -c 8 --mem=32G -t 01:30:00 --array=1-7
#SBATCH -o /user/work/fh6520/showcase/prep/logs/cpr-v7-%A_%a.out
set -euo pipefail
module load languages/R/4.5.1
export OMP_NUM_THREADS=1 LIB_CS=/user/work/fh6520/showcase/final/lib6
S=/user/work/fh6520/showcase/prep
# shared, idempotent layout: TSV/VCF and design tables are links to the originals
for p in "traits_v7 traits" "traits25_v7 traits25"; do set -- $p; mkdir -p $S/$1
  for x in tsv vcf traits.csv truth_theta.csv; do [ -e $S/$1/$x ] || ln -sfn $S/$2/$x $S/$1/$x; done; done
Rscript $S/reencode_v7.R
