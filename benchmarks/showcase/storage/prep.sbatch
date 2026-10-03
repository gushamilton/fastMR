#!/bin/bash
#SBATCH --job-name=stor-prep
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --constraint=Cascadelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/user/work/fh6520/showcase/storage/logs/prep-%j.out
module load languages/R/4.5.1 bcftools/1.19-openblas-5yp2 htslib/1.19.1-jbbb
export STORAGE_ROOT=${STORAGE_ROOT:-/user/work/fh6520/showcase/storage}
export LIB_CS=${LIB_CS:-/user/work/fh6520/regress-audit/lib-cs-final}
export STORAGE_SCRIPTS=$STORAGE_ROOT/scripts
mkdir -p $STORAGE_ROOT/{logs,results,stores,work}
hostname; lscpu | grep 'Model name'
Rscript $STORAGE_SCRIPTS/prep.R
