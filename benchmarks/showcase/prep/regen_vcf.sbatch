#!/bin/bash
#SBATCH --job-name=sc-vcf
#SBATCH --account=sscm013902
#SBATCH --partition=test
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/user/work/fh6520/showcase/prep/logs/vcf-%A_%a.out
# Re-encode .vcf only: sbatch --array=1-7 --export=ALL,TASK_CHUNK=3,LIB_CS=...,LIB_FASTMR=... regen_vcf.sbatch
# Override LIB_CS / LIB_FASTMR / LIB_TSMR / SEED via --export.  Re-running skips finished traits.
set -euo pipefail
module load languages/R/4.5.1 apps/plink2/2.00a68LM bcftools/1.19-openblas-5yp2 htslib/1.19.1-jbbb
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NTHREADS=${SLURM_CPUS_PER_TASK}
cd /user/work/fh6520/showcase/prep
Rscript regen_vcf.R
