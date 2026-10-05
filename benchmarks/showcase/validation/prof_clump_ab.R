#!/usr/bin/env Rscript
# Where the 25x25 clump stage spends its time, per build: Rscript prof_clump_ab.R LIB TRAITS_DIR
# Same call as the e2e C8 arm (partition = "auto"), Rprof by.total top 30, plus the elapsed/CPU of the call.
a <- commandArgs(TRUE); .libPaths(c(a[1], .libPaths())); suppressMessages({library(CompreSSoR); library(fastMR); loadNamespace("arrow")})
REF <- "/user/work/fh6520/showcase/prep/ref/EUR_maf01"; PL2 <- Sys.which("plink2")
expo <- sprintf("exp%02d", 1:25); ef <- setNames(file.path(a[2], "cpr", paste0(expo, ".cpr")), expo)
cat("CompreSSoR", as.character(packageVersion("CompreSSoR")), "fastMR", as.character(packageVersion("fastMR")), "\n")
Rprof(pf <- tempfile(), interval = 0.02)
t <- system.time(cl <- fast_clump_compressed(ef, pvalue_threshold = 5e-8, candidate_source = "pvalue_flag", pvalue_order = "require_exact",
     partition = "auto", pfile = REF, plink2_bin = PL2, clump_kb = 10000, clump_r2 = 0.001, threads = 8, io_threads = 8))
Rprof(NULL); print(t); cat("instruments:", sum(lengths(cl$instruments)), " strategy:", paste(unlist(cl$diagnostics$auto), collapse = " "), "\n")
s <- summaryRprof(pf)$by.total; print(head(s[order(-s$total.time), ], 30))
