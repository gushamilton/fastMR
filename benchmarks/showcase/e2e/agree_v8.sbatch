#!/bin/bash -l
#SBATCH -A sscm013902 -p short --constraint=Cascadelake -J agree-v8 -c 8 --mem=32G -t 01:00:00
#SBATCH -o /user/work/fh6520/showcase/e2e/logs/agree-v8-%j.out
# Agreement (all replicates), the single-SNP Wald ratio check (rep 1), the v7 and v8 scaling models and the
# old-vs-new tables. STEPS selects a subset, e.g. --export=ALL,STEPS="wald compare".
module load languages/R/4.5.1
export OMP_NUM_THREADS=1 LIB_CS=/user/work/fh6520/showcase/final/lib7 LIB_FASTMR=/user/work/fh6520/showcase/final/lib7
E=/user/work/fh6520/showcase/e2e; P=/user/work/fh6520/showcase/prep; STEPS=${STEPS:-"agree wald scaling compare"}
# TAG (default v8) = which fastMR run; v9 uses lib8 (fix/clump-flag-preread).
TAG=${TAG:-v8}; export TAG NEW_TAG=$TAG
[ $TAG = v9 ] && export LIB_CS=/user/work/fh6520/showcase/final/lib8 LIB_FASTMR=/user/work/fh6520/showcase/final/lib8
has() { [[ " $STEPS " == *" $1 "* ]]; }
if has agree; then
  bash $E/assemble_v8.sh
  for r in 1 2 3; do
    for s in 1x1 10x10; do Rscript $E/agree.R $E/final10_$TAG/rep$r $s; done
    Rscript $E/agree.R $E/final25_$TAG/rep$r 25x25
  done
fi
if has wald; then
  TRAITS=$P/traits_v8 Rscript $E/wald_check.R $E/final10_$TAG/rep1 10x10
  TRAITS=$P/traits25_v8 Rscript $E/wald_check.R $E/final25_$TAG/rep1 25x25
fi
if has scaling; then
  SC_FIT=$E/final10 SC_OBS="$E/final25 $E/results_v7_traits50_v7" SC_TAG=v7 Rscript $E/scaling_model.R
  SC_FIT=$E/final10_$TAG SC_OBS="$E/final25_$TAG $E/results_${TAG}_traits50_v8" SC_TAG=$TAG Rscript $E/scaling_model.R
fi
has compare && Rscript $E/compare_v8.R
