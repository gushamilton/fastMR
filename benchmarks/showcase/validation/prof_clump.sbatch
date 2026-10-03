#!/bin/bash -l
#SBATCH -A sscm013902 -p test -c 8 --mem=16G -t 00:20:00 -o /user/work/fh6520/showcase/e2e/logs/prof_clump_%j.out
module load languages/R/4.5.1 apps/plink2/2.00a68LM
Rscript /user/work/fh6520/showcase/e2e/prof_clump.R
