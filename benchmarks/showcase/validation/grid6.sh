#!/bin/bash
#SBATCH -A sscm013902 -p test -t 00:59:00 -c 4 --mem=24G -o grid6_%j.out
cd /user/work/fh6520/showcase/clumpdebug
module load languages/R/4.5.1 apps/plink2/2.00a6LM apps/plink1.9/1.90-b77 >/dev/null 2>&1
LIB=$1; B=$2; P=$3; shift 3
for mode in pfile bfile; do for r2 in 0.001 0.01 0.1; do for kb in 250 10000; do for seed in "$@"; do
 Rscript verify.R $LIB $B $P $r2 $kb 10 $seed $mode 2>&1 | grep "^ref="
done; done; done; done
