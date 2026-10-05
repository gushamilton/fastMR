#!/bin/bash -l
#SBATCH -A sscm013902 -p test --constraint=Cascadelake -J prof-bx -c 8 --mem=16G -t 00:25:00
#SBATCH -o /user/work/fh6520/showcase/e2e/logs/prof-batch-extract-%j.out
module load languages/R/4.5.1
grep -m1 'model name' /proc/cpuinfo | grep -q "Gold 6226R" || { echo wrong CPU; exit 3; }
F=/user/work/fh6520/showcase/final; P=/user/work/fh6520/showcase/prep
for rep in 1 2 3; do for K in 1000 100000; do for T in 1 8; do
  Rscript /user/work/fh6520/showcase/e2e/prof_batch_extract.R $F/lib6 $P/traits_v7/cpr $K $T
  Rscript /user/work/fh6520/showcase/e2e/prof_batch_extract.R $F/lib8 $P/traits_v8/cpr $K $T
  Rscript /user/work/fh6520/showcase/e2e/prof_batch_extract.R $F/lib8 $P/traits_v8/cpr $K $T noverify
done; done; done
