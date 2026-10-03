#!/bin/bash -l
#SBATCH -A sscm013902 -p short -c 8 --mem=24G -t 00:45:00 -o /user/work/fh6520/showcase/e2e/logs/clumpfix_%j.out
module load languages/R/4.5.1 apps/plink2/2.00a68LM
F=/user/work/fh6520/showcase/final; L=$F/lib-clumpfix
T=/user/work/fh6520/toolchains/rust-1.87.0
export PATH=$T/rustc/bin:$T/cargo/bin:$PATH LD_LIBRARY_PATH=$T/rustc/lib:$LD_LIBRARY_PATH CARGO_HOME=$(mktemp -d) CARGO_NET_OFFLINE=true
rm -rf $L; mkdir -p $L
R CMD INSTALL --library=$L $F/src-cs >/dev/null 2>&1 && R CMD INSTALL --library=$L $F/src-fm-cf > $F/build-clumpfix.log 2>&1 || { echo BUILD_FAIL; tail $F/build-clumpfix.log; exit 1; }
S=/user/work/fh6520/showcase/e2e/clumpfix_check.R
# alternate order so neither build always runs on a warm page cache
Rscript $S $F/lib old_warmup >/dev/null 2>&1
Rscript $S $F/lib-clumpfix new; Rscript $S $F/lib old; Rscript $S $F/lib-clumpfix new2
Rscript -e 'a <- readRDS("/user/work/fh6520/showcase/e2e/clumpfix_old.rds"); b <- readRDS("/user/work/fh6520/showcase/e2e/clumpfix_new.rds"); for (k in names(a)) cat(k, "identical:", identical(a[[k]], b[[k]]), "\n"); cat("small-cap == one-call:", identical(b$e10, b$e10_small_cap), "\n")'
