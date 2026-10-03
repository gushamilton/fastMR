#!/bin/bash
#SBATCH --job-name=sc-ref
#SBATCH --account=sscm013902
#SBATCH --partition=test
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/user/work/fh6520/showcase/prep/logs/ref-%j.out
set -euo pipefail
module load apps/plink2/2.00a68LM apps/plink1.9/1.90-b77 bcftools/1.19-openblas-5yp2 htslib/1.19.1-jbbb
P=/user/work/fh6520/showcase/prep; SRC=/user/work/fh6520/1kg_ref/grch38_superpops/EUR
cd $P/ref; T=$P/ref/tmp; mkdir -p $T
date
# autosomes only, biallelic ACGT SNPs, MAF>=0.01, IDs = chr:pos:ref:alt, dedup by ID
plink2 --bfile $SRC --chr 1-22 --snps-only just-acgt --max-alleles 2 --maf 0.01 \
  --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 8 \
  --rm-dup force-first --threads 16 --memory 56000 --make-pgen --out EUR_maf01
# bed from same pgen (A1=ALT, A2=REF)
plink2 --pfile EUR_maf01 --threads 16 --memory 56000 --make-bed --out EUR_maf01
rm -rf $T
wc -l EUR_maf01.bim EUR_maf01.fam EUR_maf01.pvar | tail -4
head -3 EUR_maf01.bim; grep -v '^##' EUR_maf01.pvar | head -3
ls -la
# sanity: plink1.9 can read the bed
plink --bfile EUR_maf01 --chr 22 --freq --out /dev/null 2>&1 | tail -3 || true
date
