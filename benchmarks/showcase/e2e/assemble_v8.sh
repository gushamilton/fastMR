#!/bin/bash
# Assemble final10_v8 / final25_v8: the v7 TwoSampleMR and GWAS-VCF arm dirs (A, A8, B; same node type, code-independent)
# linked next to the new fastMR arm dirs (C1, C8 from results_v8_*), one dir per replicate, for agree.R / compare_v8.R.
set -euo pipefail
E=/user/work/fh6520/showcase/e2e
for spec in "final10:results_v8_traits_v8:1x1 10x10" "final25:results_v8_traits25_v8:25x25"; do
  IFS=: read old src sizes <<< "$spec"
  for r in 1 2 3; do
    d=$E/${old}_v8/rep$r; mkdir -p $d
    for s in $sizes; do
      for a in A A8 B; do [ -d $E/$old/rep$r/${s}_$a ] && ln -sfn $E/$old/rep$r/${s}_$a $d/${s}_$a; done
      for a in C1 C8; do ln -sfn $E/$src/rep$r/${s}_$a $d/${s}_$a; done
    done
    cp $E/$src/rep$r/cells.csv $d/cells_fastmr.csv
  done
done
ls -l $E/final10_v8/rep1 $E/final25_v8/rep1
