# Rscript mode_paths.R <fastMR package source root>
# Compiles mode_paths.cpp against a private copy of src/fastmr.cpp (included verbatim).
pkg <- normalizePath(commandArgs(TRUE)[1])
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
work <- tempfile("fastmr-mode_paths-"); dir.create(work)
file.copy(file.path(pkg, "src", "fastmr.cpp"), file.path(work, "fastmr_copy.h"))
file.copy(file.path(here, "mode_paths.cpp"), file.path(work, "mode_paths.cpp"))
Sys.setenv(PKG_CXXFLAGS = "-std=c++17 -O2")
Rcpp::sourceCpp(file.path(work, "mode_paths.cpp"), rebuild = TRUE)
set.seed(3)
out <- NULL
ks <- if (nzchar(Sys.getenv("MODE_KS"))) as.numeric(strsplit(Sys.getenv("MODE_KS"), ",")[[1]]) else c(3, 5, 8, 10, 16, 24, 30, 40, 48, 64, 80, 100, 150, 200)
for (k in ks) for (sc in c("normal", "outlier")) {
  draws <- if (k <= 30) 6000 else if (k <= 200) 2000 else max(100, round(4e5 / k))
  mu <- rnorm(k, 0.2, 0.1)
  if (sc == "outlier") mu[seq_len(max(1, k %/% 3))] <- mu[seq_len(max(1, k %/% 3))] + 2
  se <- runif(k, 0.02, 0.1)
  R <- matrix(rnorm(k * draws, mu, se), k, draws)
  r <- mode_paths(R, 1 / se^2, rep(1, k), reps = 5)
  out <- rbind(out, data.frame(k = k, sc = sc, t(round(r, 1))))
}
options(width = 200); print(out, row.names = FALSE)
