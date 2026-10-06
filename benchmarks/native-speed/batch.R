# Five-method batch wall time (seconds) at a given thread count.
#   Rscript batch.R <lib> <label> <threads> <out.csv>
args <- commandArgs(TRUE)
.libPaths(c(args[1], .libPaths())); suppressMessages(library(fastMR))
th <- as.integer(args[3])
# FASTMR_HULL="<max ratios>,<recurrence 0/1>" overrides the hull settings (new build only).
hull_env <- Sys.getenv("FASTMR_HULL")
if (nzchar(hull_env)) { h <- as.numeric(strsplit(hull_env, ",")[[1]]); fastMR:::fastmr_set_mode_hull_native(h[1], as.logical(h[2])) }
# FASTMR_OVERLAP=0/1 turns the double-buffered bootstrap batches off/on (new build only).
if (nzchar(Sys.getenv("FASTMR_OVERLAP"))) fastMR:::fastmr_set_bootstrap_overlap_native(Sys.getenv("FASTMR_OVERLAP") == "1")
reps <- as.integer(Sys.getenv("BATCH_REPS", "3"))
mk <- function(sizes, seed = 1) {
  set.seed(seed); g <- rep(seq_along(sizes), sizes); n <- length(g)
  bx <- rnorm(n, 0.1, 0.03) * sample(c(-1, 1), n, TRUE)
  data.frame(SNP = paste0("rs", seq_len(n)), id.exposure = paste0("e", g), id.outcome = "o",
             beta.exposure = bx, se.exposure = runif(n, 0.005, 0.02),
             beta.outcome = rnorm(length(sizes), 0, 0.3)[g] * bx + rnorm(n, 0, 0.01) + (runif(n) < 0.2) * rnorm(n, 0, 0.05),
             se.outcome = runif(n, 0.005, 0.02))
}
five <- c("wald_ratio", "egger", "weighted_median", "ivw", "weighted_mode")
set.seed(99)
work <- list(k3 = rep(3, 12000), k10 = rep(10, 6000), k30 = rep(30, 2400), k100 = rep(100, 800),
             mixed = sample(c(1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100), 8000, TRUE,
                            prob = c(10, 8, 8, 7, 6, 5, 5, 4, 3, 2, 2, 1, 0.5)))
out <- NULL
for (w in names(work)) {
  d <- mk(work[[w]])
  ts <- replicate(reps, system.time(fast_mr(d, methods = five, nboot = 1000, seed = 1, threads = th))[["elapsed"]])
  out <- rbind(out, data.frame(label = args[2], workload = w, pairs = length(work[[w]]), threads = th,
                               wall_s = min(ts), wall_s_median = median(ts)))
  print(tail(out, 1))
}
write.csv(out, args[4], row.names = FALSE)
