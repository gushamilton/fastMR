#!/bin/bash
#SBATCH --job-name=stor-multi
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --constraint=Cascadelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --array=1-3
#SBATCH --output=/user/work/fh6520/showcase/storage/multi/logs/multi-%A_%a.out
# Multi-trait storage suite (N simulated GWAS per format). Override NTRAITS / KEY_SIZES / RESULT_TAG for a dry run.
module load languages/R/4.5.1 bcftools/1.19-openblas-5yp2 htslib/1.19.1-jbbb
export STORAGE_ROOT=/user/work/fh6520/showcase/storage STORAGE_SCRIPTS=/user/work/fh6520/showcase/storage/scripts
export MULTI_ROOT=$STORAGE_ROOT/multi MULTI_SCRIPTS=$STORAGE_ROOT/multi/scripts
export LIB_CS=${LIB_CS:-/user/work/fh6520/showcase/final/lib}
export REP=${REP:-$SLURM_ARRAY_TASK_ID}
hostname; lscpu | grep 'Model name'; echo "nproc=$(nproc) SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK LIB_CS=$LIB_CS"
Rscript $MULTI_SCRIPTS/make_keys.R && Rscript $MULTI_SCRIPTS/run_multi.R
