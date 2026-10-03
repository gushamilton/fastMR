#!/bin/bash -l
#SBATCH -A sscm013902
#SBATCH -p test
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 0:40:00
#SBATCH -o /user/work/fh6520/showcase/final/build2.log
module load languages/R/4.5.1
T=/user/work/fh6520/toolchains/rust-1.87.0
export PATH=$T/rustc/bin:$T/cargo/bin:$PATH LD_LIBRARY_PATH=$T/rustc/lib:$LD_LIBRARY_PATH CARGO_HOME=$(mktemp -d) CARGO_NET_OFFLINE=true
F=/user/work/fh6520/showcase/final; L=$F/lib2; rm -rf $L; mkdir -p $L
R CMD INSTALL --library=$L $F/src-cs && R CMD INSTALL --library=$L $F/src-fm2 || exit 1
Rscript -e ".libPaths(c('$L', .libPaths())); library(CompreSSoR); library(fastMR); cat('CompreSSoR', as.character(packageVersion('CompreSSoR')), 'fastMR', as.character(packageVersion('fastMR')), '\n')" && echo BUILD_OK
