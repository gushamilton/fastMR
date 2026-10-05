#!/bin/bash -l
#SBATCH -A sscm013902 -p test --constraint=Cascadelake -J prof-clump -c 8 --mem=32G -t 00:20:00
#SBATCH -o /user/work/fh6520/showcase/e2e/logs/prof-clump-ab-%j.out
module load languages/R/4.5.1 apps/plink2/2.00a68LM
export OMP_NUM_THREADS=1
grep -m1 'model name' /proc/cpuinfo
F=/user/work/fh6520/showcase/final; P=/user/work/fh6520/showcase/prep
for b in v7 v8 v7 v8; do lib=$F/lib6; [ $b = v8 ] && lib=$F/lib7
  echo "===== $b"; Rscript /user/work/fh6520/showcase/e2e/prof_clump_ab.R $lib $P/traits25_$b; done
