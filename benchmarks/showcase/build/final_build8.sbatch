#!/bin/bash -l
#SBATCH -A sscm013902
#SBATCH -p test
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 0:50:00
#SBATCH -o /user/work/fh6520/showcase/final/build8.log
# lib8 = CompreSSoR d158622 + fastMR fix/clump-flag-preread (671c85f); lib7 (unfixed 543cae2) rebuilt for the A/B.
# Then the fastMR clumping guard tests against lib8.
module load languages/R/4.5.1
T=/user/work/fh6520/toolchains/rust-1.87.0
export PATH=$T/rustc/bin:$T/cargo/bin:$PATH LD_LIBRARY_PATH=$T/rustc/lib:$LD_LIBRARY_PATH CARGO_HOME=$(mktemp -d) CARGO_NET_OFFLINE=true
F=/user/work/fh6520/showcase/final
for v in 7 8; do L=$F/lib$v; rm -rf $L; mkdir -p $L
  R CMD INSTALL --library=$L $F/src-cs$v && R CMD INSTALL --library=$L $F/src-fm$v || exit 1
  Rscript -e ".libPaths(c('$L', .libPaths())); library(CompreSSoR); library(fastMR); cat('lib$v CompreSSoR', as.character(packageVersion('CompreSSoR')), 'fastMR', as.character(packageVersion('fastMR')), '\n')"
done
cd $F/src-fm8/tests/testthat && Rscript -e ".libPaths(c('$F/lib8', '/user/work/fh6520/r_packages', .libPaths())); library(testthat); library(fastMR);
  for (f in c('test-candidate-count-guard.R', list.files('.', pattern = '^test-.*clump.*[.]R$'))) { cat('==', f, '\n'); r <- as.data.frame(test_file(f, reporter = 'summary', stop_on_failure = FALSE)); print(colSums(r[c('failed','skipped','error','warning','passed')])) }" && echo BUILD_OK
