# Per-method native cost, microseconds per pair, 1 thread, nboot = 1000.
#   Rscript per-method.R <lib> <label> <out.csv>
args <- commandArgs(TRUE)
.libPaths(c(args[1], .libPaths())); suppressMessages(library(fastMR))
mk <- function(G, k, seed = 1) {
  set.seed(seed); n <- G * k
  bx <- rnorm(n, 0.1, 0.03) * sample(c(-1, 1), n, TRUE)
  grp <- rep(seq_len(G), each = k)
  data.frame(SNP = paste0("rs", rep(seq_len(k), G)), id.exposure = paste0("e", grp), id.outcome = "o",
             beta.exposure = bx, se.exposure = runif(n, 0.005, 0.02),
             beta.outcome = rnorm(G, 0, 0.3)[grp] * bx + rnorm(n, 0, 0.01) + (runif(n) < 0.2) * rnorm(n, 0, 0.05),
             se.outcome = runif(n, 0.005, 0.02))
}
meth <- list(wald_ratio = "wald_ratio", ivw = "ivw", egger = "egger", weighted_median = "weighted_median",
             weighted_mode = "weighted_mode", five = c("wald_ratio", "egger", "weighted_median", "ivw", "weighted_mode"))
out <- NULL
for (k in c(3, 10, 30, 100)) {
  G <- max(30, round(3000 / k))
  d <- mk(G, k)
  for (m in names(meth)) {
    ts <- replicate(3, { t <- system.time(fast_mr(d, methods = meth[[m]], nboot = 1000, seed = 1, threads = 1)); c(t[["elapsed"]], t[["user.self"]]) })
    out <- rbind(out, data.frame(label = args[2], k = k, G = G, method = m,
                                 us_per_pair_elapsed = 1e6 * min(ts[1, ]) / G, us_per_pair_user = 1e6 * min(ts[2, ]) / G))
  }
}
write.csv(out, args[3], row.names = FALSE)
print(out)
