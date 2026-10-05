#!/bin/bash -l
#SBATCH -A sscm013902 -p test --constraint=Cascadelake -J prof-bx2 -c 8 --mem=16G -t 00:25:00
#SBATCH -o /user/work/fh6520/showcase/e2e/logs/prof-batch-extract2-%j.out
# Batch extract, 20 stores: v7 (lib6) vs d158622 default (one pass per store) vs d158622 with panel sharing on.
module load languages/R/4.5.1
grep -m1 'model name' /proc/cpuinfo | grep -q "Gold 6226R" || { echo wrong CPU; exit 3; }
F=/user/work/fh6520/showcase/final; P=/user/work/fh6520/showcase/prep; S=/user/work/fh6520/showcase/e2e/prof_batch_extract.R
for b in 6 8; do s=v8; [ $b = 6 ] && s=v7; echo "===== Rprof lib$b"; Rscript $S $F/lib$b $P/traits_$s/cpr 1000 1 >/dev/null; RPROF=1 Rscript $S $F/lib$b $P/traits_$s/cpr 1000 1; done
for rep in 1 2 3; do for K in 1000 100000; do for T in 1 8; do
  Rscript $S $F/lib6 $P/traits_v7/cpr $K $T
  Rscript $S $F/lib8 $P/traits_v8/cpr $K $T
  Rscript $S $F/lib8 $P/traits_v8/cpr $K $T share
done; done; done
