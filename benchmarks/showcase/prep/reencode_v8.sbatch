#!/bin/bash -l
#SBATCH -A sscm013902 -p short --constraint=Cascadelake -J cpr-v8 -c 8 --mem=32G -t 01:30:00 --array=1-7
#SBATCH -o /user/work/fh6520/showcase/prep/logs/cpr-v8-%A_%a.out
# Re-encode the trait stores with lib7 (CompreSSoR 0.7.2 d158622) into traits_v8 / traits25_v8.
set -euo pipefail
module load languages/R/4.5.1
export OMP_NUM_THREADS=1 LIB_CS=/user/work/fh6520/showcase/final/lib7 STORE_TAG=v8
S=/user/work/fh6520/showcase/prep
for p in "traits_v8 traits" "traits25_v8 traits25"; do set -- $p; mkdir -p $S/$1
  for x in tsv vcf traits.csv truth_theta.csv; do [ -e $S/$1/$x ] || ln -sfn $S/$2/$x $S/$1/$x; done; done
Rscript $S/reencode_v7.R
