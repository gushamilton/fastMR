#!/bin/bash
# 50x50 timing point for the scaling model, without generating new GWAS: link the 50 traits of the
# 25x25 set under 50 exposure and 50 outcome names (exp26-50 = out01-25, out26-50 = exp01-25).
# Used only for fastMR timing; the pairs are not a meaningful MR design.
set -euo pipefail
S=${1:-/user/work/fh6520/showcase/prep/traits25}; D=${2:-/user/work/fh6520/showcase/prep/traits50}
rm -rf $D; mkdir -p $D/cpr
for i in $(seq 1 25); do
  e=$(printf "exp%02d" $i); o=$(printf "out%02d" $i); e2=$(printf "exp%02d" $((i+25))); o2=$(printf "out%02d" $((i+25)))
  ln -s $S/cpr/$e.cpr $D/cpr/$e.cpr; ln -s $S/cpr/$o.cpr $D/cpr/$o.cpr
  ln -s $S/cpr/$o.cpr $D/cpr/$e2.cpr; ln -s $S/cpr/$e.cpr $D/cpr/$o2.cpr
done
awk -F, 'NR==1 {print; next} {print; n[$1]=$3}
  END { for (i=1;i<=25;i++) { e=sprintf("exp%02d",i); o=sprintf("out%02d",i);
        printf "exp%02d,exposure,%s,\nout%02d,outcome,%s,\n", i+25, n[o], i+25, n[e] } }' $S/traits.csv > $D/traits.csv
ls $D/cpr | wc -l; wc -l < $D/traits.csv
