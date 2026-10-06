#!/bin/bash -l
#SBATCH -A sscm013902 -p test --constraint=Cascadelake -J prof-verify -c 8 --mem=16G -t 00:15:00
#SBATCH -o /user/work/fh6520/showcase/e2e/logs/prof-verify-%j.out
module load languages/R/4.5.1
grep -m1 'model name' /proc/cpuinfo
for i in 1 2; do Rscript /user/work/fh6520/showcase/e2e/profile_verify.R /user/work/fh6520/showcase/final/lib7 /user/work/fh6520/showcase/prep/traits25_v8; done
