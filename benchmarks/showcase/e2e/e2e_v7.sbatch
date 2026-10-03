#!/bin/bash -l
#SBATCH -A sscm013902 -p short --constraint=Cascadelake -J e2e-v7 -c 8 --mem=64G -t 03:00:00 --array=1-3%1
#SBATCH -o /user/work/fh6520/showcase/e2e/logs/e2e-v7-%A_%a.out
# Final fastMR arms (C1, C8) on lib6 + format-0.4.6 stores, at 1x1, 10x10, 25x25 and 50x50. Replicates run one at a time.
P=/user/work/fh6520/showcase/prep; E=/user/work/fh6520/showcase/e2e
[ -d $P/traits50_v7 ] || bash $E/make_traits50b.sh $P/traits25_v7 $P/traits50_v7
export LIB_CS=/user/work/fh6520/showcase/final/lib6 LIB_FASTMR=/user/work/fh6520/showcase/final/lib6 ARMS="C1 C8" CAP=1500
for spec in "1x1 10x10:traits_v7" "25x25:traits25_v7" "50x50:traits50_v7"; do
  SIZES=${spec%%:*} TRAITS=$P/${spec##*:} RESDIR=$E/results_v7_${spec##*:} bash $E/e2e.sbatch
done
