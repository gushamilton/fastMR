# Synthetic batched-MVMR benchmark: E exposures x K outcomes x 2 covariates.
# Usage: Rscript synthetic.R <lib> <out.csv> [E] [K]
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[[1]], .libPaths()))
suppressPackageStartupMessages(library(fastMR))
out_csv <- args[[2]]
E <- if (length(args) >= 3) as.integer(args[[3]]) else 1000L
K <- if (length(args) >= 4) as.integer(args[[4]]) else 3000L
own <- 30L; cov_inst <- 200L; p <- 3L
set.seed(20261005)
U <- E * own + 2L * cov_inst
cat("panel", U, "x", K, "\n")
# Covariate effects at every SNP; exposure effects at its own instruments
# (strong) and at the covariate instruments (weak pleiotropy).
c1 <- rnorm(U, 0, 0.01); c2 <- rnorm(U, 0, 0.01)
cov_rows <- E * own + seq_len(2L * cov_inst)
c1[cov_rows[1:cov_inst]] <- rnorm(cov_inst, 0, 0.05)
c2[cov_rows[cov_inst + 1:cov_inst]] <- rnorm(cov_inst, 0, 0.05)
row_se <- runif(U, 0.005, 0.02)
scale <- runif(K, 0.7, 1.4)
S <- outer(row_se, scale)                     # proportional SEs (shared path applies)
gamma <- matrix(rnorm(2 * K, 0, 0.3), 2, K)
B <- cbind(c1, c2) %*% gamma + matrix(rnorm(U * K), U, K) * S
designs <- lapply(seq_len(E), function(e) {
  rows <- c((e - 1L) * own + seq_len(own), cov_rows)
  bx <- c(rnorm(own, 0, 0.05), rnorm(length(cov_rows), 0, 0.005))
  list(rows = rows, beta = cbind(exposure = bx, c1 = c1[rows], c2 = c2[rows]),
       se = cbind(rep(0.004, length(rows)), rep(0.003, length(rows)), rep(0.003, length(rows))))
})
names(designs) <- paste0("E", seq_len(E))
# Exposure-outcome causal effects for the first 30 instruments of each design.
for (e in seq_len(E)) {
  r <- designs[[e]]$rows[seq_len(own)]
  B[r, ] <- B[r, ] + outer(designs[[e]]$beta[seq_len(own), 1], rnorm(K, 0, 0.1))
}
gc()
timeit <- function(expr) { gc(); t <- system.time(force(expr))[["elapsed"]]; t }
rows <- list()
add <- function(...) { rows[[length(rows) + 1L]] <<- data.frame(...); print(tail(rows, 1)[[1]]) }
cpu <- tryCatch(sub(".*: *", "", grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1]), error = function(e) NA)
for (th in c(1L, 2L, 4L, 8L)) for (w in c("exact", "shared")) {
  res <- NULL
  t <- timeit(res <- suppressWarnings(fast_mvmr_ivw_batch(designs, B, S, threads = th, weights = w,
                                                          exposure_cor = diag(3), weak_f = -Inf)))
  add(case = "kernel_with_QA", weights = w, threads = th, E = E, K = K, p = p, seconds = t,
      pairs_per_second = E * K / t, shared_fits = sum(res$shared), cpu = cpu)
  if (th == 8L && w == "exact") ref8 <- res
  if (th == 8L && w == "shared") sh8 <- res
}
for (th in c(1L, 8L)) {
  d2 <- lapply(designs, function(d) d[c("rows", "beta")])
  t <- timeit(res <- fast_mvmr_ivw_batch(d2, B, S, threads = th, weak_f = -Inf))
  add(case = "kernel_no_diagnostics", weights = "exact", threads = th, E = E, K = K, p = p,
      seconds = t, pairs_per_second = E * K / t, shared_fits = 0, cpu = cpu)
}
cat("shared vs exact max rel diff b:", max(abs(sh8$b - ref8$b) / abs(ref8$b)), "\n")
# Baseline 1: fast_mr_multivariable() looped over outcomes, for a subset.
sub_e <- 3L
t <- timeit(for (e in seq_len(sub_e)) {
  d <- designs[[e]]
  for (k in seq_len(K)) fast_mr_multivariable(d$beta, B[d$rows, k], S[d$rows, k])
})
add(case = "fast_mr_multivariable_loop", weights = "exact", threads = 1L, E = sub_e, K = K, p = p,
    seconds = t, pairs_per_second = sub_e * K / t, shared_fits = 0, cpu = cpu)
# Agreement on that subset.
md <- 0
for (e in seq_len(sub_e)) {
  d <- designs[[e]]
  for (k in seq(1, K, length.out = 50)) {
    r <- fast_mr_multivariable(d$beta, B[d$rows, k], S[d$rows, k])
    md <- max(md, abs(r$b - ref8$b[e, k, ]) / abs(r$b), abs(r$se - ref8$se[e, k, ]) / r$se)
  }
}
cat("max relative difference vs fast_mr_multivariable:", md, "\n")
# Baseline 2: vectorised-over-outcomes R loop over exposures (the approach of
# F-platelet 06_mvmr.R wls()), for a subset.
wls_r <- function(X, rows) {
  By <- B[rows, , drop = FALSE]; W <- 1 / S[rows, , drop = FALSE]^2
  p <- ncol(X); XtWX <- array(0, c(p, p, K)); XtWy <- matrix(0, p, K)
  for (a in 1:p) { XtWy[a, ] <- colSums(W * X[, a] * By)
    for (b in a:p) { v <- colSums(W * X[, a] * X[, b]); XtWX[a, b, ] <- v; XtWX[b, a, ] <- v } }
  beta <- matrix(NA_real_, p, K)
  for (k in 1:K) beta[, k] <- solve(XtWX[, , k], XtWy[, k])
  res <- By - X %*% beta; colSums(W * res^2); beta
}
sub_e2 <- 50L
t <- timeit(for (e in seq_len(sub_e2)) wls_r(designs[[e]]$beta, designs[[e]]$rows))
add(case = "R_loop_vectorised_outcomes", weights = "exact", threads = 1L, E = sub_e2, K = K, p = p,
    seconds = t, pairs_per_second = sub_e2 * K / t, shared_fits = 0, cpu = cpu)
out <- do.call(rbind, rows)
out$extrapolated_full_seconds <- E * K / out$pairs_per_second
write.csv(out, out_csv, row.names = FALSE)
cat("max_rel_diff_vs_fast_mr_multivariable", md, "\n")
