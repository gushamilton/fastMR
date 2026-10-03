# New graph clumping (one PLINK2 call) vs current build: identical instruments, wall time.
args <- commandArgs(TRUE); LIB <- args[1]; LAB <- args[2]
.libPaths(c(LIB, .libPaths()))
suppressMessages({library(CompreSSoR); library(fastMR)})
REF <- "/user/work/fh6520/showcase/prep/ref/EUR_maf01"; PL2 <- Sys.which("plink2")
TR <- "/user/work/fh6520/showcase/prep/traits/cpr_new"
run <- function(ids, ...) {
  ef <- setNames(file.path(TR, paste0(ids, ".cpr")), ids)
  t <- system.time(cl <- fast_clump_compressed(ef, pvalue_threshold = 5e-8, candidate_source = "pvalue_flag",
       pvalue_order = "require_exact", partition = "graph", pfile = REF, plink2_bin = PL2,
       clump_kb = 10000, clump_r2 = 0.001, threads = 8, io_threads = 8, ...))[["elapsed"]]
  list(t = t, inst = cl$instruments, calls = cl$diagnostics$graph_calls %||% NA, fb = cl$diagnostics$fallback %||% NA)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
out <- list(e1 = run("exp01"), e10 = run(sprintf("exp%02d", 1:10)),
            e10_small_cap = run(sprintf("exp%02d", 1:10), max_graph_pairs = 2e4))
for (k in names(out)) cat(sprintf("%s %-14s wall=%6.2fs graph_calls=%s fallback=%s n_inst=%d\n", LAB, k, out[[k]]$t,
                                  out[[k]]$calls, out[[k]]$fb, sum(lengths(out[[k]]$inst))))
saveRDS(lapply(out, `[[`, "inst"), sprintf("/user/work/fh6520/showcase/e2e/clumpfix_%s.rds", LAB))
