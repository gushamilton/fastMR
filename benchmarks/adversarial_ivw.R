args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else getwd()
.libPaths(c(file.path(root, ".local", "Rlib"), .libPaths()))
suppressPackageStartupMessages(library(fastMR))
set.seed(20260808)
rows <- vector("list", 100L)
for (iteration in seq_len(100L)) {
  ne <- sample.int(8L, 1L) + 1L
  no <- sample.int(8L, 1L) + 1L
  nsnp <- sample.int(39L, 1L) + 1L
  x <- matrix(rnorm(ne * nsnp), nrow = ne)
  y <- matrix(rnorm(no * nsnp), nrow = no)
  sx <- matrix(runif(ne * nsnp, 0.01, 0.3), nrow = ne)
  sy <- matrix(runif(no * nsnp, 0.01, 0.3), nrow = no)
  if (iteration %% 4L == 0L) x[1L, 1L] <- 0
  if (iteration %% 4L == 0L) x[1L, 2L] <- -abs(x[1L, 2L])
  if (iteration %% 5L == 0L) x[2L, 3L] <- 1e-15
  if (iteration == 100L) x[,] <- 0
  direct <- fast_mr_grid(x, y, sx, sy,
                         methods = c("ivw", "ivw_fe", "ivw_mre"), nboot = 0, threads = 5)
  scalar <- fast_mr_grid(x, y, sx, sy,
                         methods = c("ivw", "ivw_fe", "ivw_mre", "uwr"), nboot = 0, threads = 1)
  scalar <- scalar[scalar$method_code %in% c("ivw", "ivw_fe", "ivw_mre"), ]
  numeric_columns <- c("b", "se", "pval")
  max_scalar_delta <- max(vapply(numeric_columns, function(name) {
    a <- direct[[name]]; b <- scalar[[name]]
    ok <- is.finite(a) & is.finite(b)
    if (any(ok)) max(abs(a[ok] - b[ok])) else 0
  }, numeric(1)))
  threaded <- fast_mr_grid(x, y, sx, sy,
                           methods = c("ivw", "ivw_fe", "ivw_mre"), nboot = 0, threads = 1)
  max_thread_delta <- max(vapply(numeric_columns, function(name) {
    a <- direct[[name]]; b <- threaded[[name]]
    ok <- is.finite(a) & is.finite(b)
    if (any(ok)) max(abs(a[ok] - b[ok])) else 0
  }, numeric(1)))
  rows[[iteration]] <- data.frame(iteration = iteration, exposures = ne, outcomes = no,
                                  snps = nsnp, max_scalar_delta = max_scalar_delta,
                                  max_thread_delta = max_thread_delta,
                                  scalar_gate = max_scalar_delta <= 1e-12,
                                  thread_gate = max_thread_delta == 0)
}
result <- do.call(rbind, rows)
dir.create(file.path(root, "outputs"), showWarnings = FALSE, recursive = TRUE)
write.csv(result, file.path(root, "outputs", "adversarial_ivw_threads.csv"), row.names = FALSE)
cat("iterations", nrow(result), "max scalar delta", max(result$max_scalar_delta),
    "max thread delta", max(result$max_thread_delta),
    "scalar failures", sum(!result$scalar_gate),
    "thread failures", sum(!result$thread_gate), "\n")
