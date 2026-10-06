# Rscript audit.R <fastMR package source root>
# Compiles audit.cpp against a private copy of src/fastmr.cpp (included verbatim).
pkg <- normalizePath(commandArgs(TRUE)[1])
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
work <- tempfile("fastmr-audit-"); dir.create(work)
file.copy(file.path(pkg, "src", "fastmr.cpp"), file.path(work, "fastmr_copy.h"))
file.copy(file.path(here, "audit.cpp"), file.path(work, "audit.cpp"))
Sys.setenv(PKG_CXXFLAGS = "-std=c++17 -O2")
Rcpp::sourceCpp(file.path(work, "audit.cpp"), rebuild = TRUE)
set.seed(1)
# 1. Recurrence kernel error vs the 2240 u bound, a = delta/h over [1e-5, 1].
a <- c(10^seq(-5, 0, length.out = 20000), runif(20000, 0, 1))
h <- 10^runif(length(a), -8, 8)
e <- recurrence_error_ulps(a * h, h)
cat(sprintf("recurrence: %d bandwidth ratios, max error %.1f u (bound 2240 u, guard uses 4096 u); 99.9%% quantile %.1f u\n",
            length(a), max(e), quantile(e, 0.999)))
# 2. Hull audit over adversarial draws.
tot <- c(certified = 0, fallbacks = 0, mismatches = 0); worst <- 0
for (rep in 1:400) {
  k <- sample(c(2:20, 24, 32, 48, 64, 100, 200), 1)
  draws <- 500
  scale <- 10^runif(1, -6, 6)
  type <- sample(c("normal", "outlier", "mirror", "discrete", "cluster"), 1)
  mu <- switch(type,
    normal = rnorm(k, 0.2, 0.1), outlier = c(rnorm(k - k %/% 3, 0.2, 0.05), rnorm(k %/% 3, 3, 1)),
    mirror = 0.37 + c(-1, 1)[(seq_len(k) %% 2) + 1] * rep(seq(0.1, 1, length.out = ceiling(k / 2)), each = 2)[seq_len(k)],
    discrete = sample(c(0.1, 0.2, 0.5), k, TRUE), cluster = rep(rnorm(2), length.out = k))
  se <- runif(k, 0.001, 0.1) * if (type %in% c("mirror", "discrete")) 10^runif(1, -9, -1) else 1
  if (type == "mirror") se <- rep(se[1], k)
  R <- scale * matrix(rnorm(k * draws, mu, se), k, draws)
  W <- matrix(if (runif(1) < 0.5) 1 / se^2 else 1, k, draws)
  r <- hull_audit(R, W, sample(c(0.25, 0.5, 1, 2), 1))
  tot <- tot + r[names(tot)]; worst <- max(worst, r[["worst_rel_discrepancy"]])
}
cat(sprintf("hull audit: %d draws certified, %d fell back, %d argmax mismatches vs FFT; worst |d_hull - d_fft| / max(y) = %.2e (u = %.2e; guard floor 1e-9)\n",
            tot[["certified"]], tot[["fallbacks"]], tot[["mismatches"]], worst, .Machine$double.eps / 2))
