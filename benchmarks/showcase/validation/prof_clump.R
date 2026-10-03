.libPaths(c("/user/work/fh6520/showcase/final/lib", .libPaths()))
suppressMessages({library(CompreSSoR); library(fastMR); library(data.table)})
REF <- "/user/work/fh6520/showcase/prep/ref/EUR_maf01"; PL2 <- Sys.which("plink2")
ef <- c(exp01 = "/user/work/fh6520/showcase/prep/traits/cpr/exp01.cpr")
tm <- function(lab, e) { t <- system.time(v <- e)[["elapsed"]]; cat(sprintf("%-40s %7.2f s\n", lab, t)); invisible(v) }
cand <- tm("read_candidates (pvalue_flag)", read_candidates(ef[[1]], 5e-8, strategy = "pvalue_flag", threads = 8))
cat("candidates:", nrow(cand), "\n")
# bare plink2 cost: load the 9M-line pvar, extract a handful of SNPs, r2
ids <- head(fread(paste0(REF, ".pvar"), skip = "#CHROM", nrows = 50)$ID, 20)
f <- tempfile(); writeLines(ids, f)
tm("plink2 --extract 20 SNPs --r2-phased", system2(PL2, c("--pfile", REF, "--extract", f, "--r2-phased", "--threads", "8", "--out", tempfile()), stdout = FALSE, stderr = FALSE))
tm("plink2 --extract 20 SNPs --make-pgen", system2(PL2, c("--pfile", REF, "--extract", f, "--make-pgen", "--threads", "8", "--out", tempfile()), stdout = FALSE, stderr = FALSE))
Rprof(pf <- tempfile(), interval = 0.05)
cl <- tm("fast_clump_compressed (1 exposure)", fast_clump_compressed(ef, pvalue_threshold = 5e-8, candidate_source = "pvalue_flag",
     pvalue_order = "require_exact", partition = "graph", pfile = REF, plink2_bin = PL2, clump_kb = 10000, clump_r2 = 0.001, threads = 8, io_threads = 8))
Rprof(NULL); s <- summaryRprof(pf)$by.total; print(head(s[order(-s$total.time), ], 25))
