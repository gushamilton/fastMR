.libPaths(c("/Users/fh6520/projects/fastMR/.local/Rlib", .libPaths())); suppressMessages(library(fastMR))
E <- 60L; S <- 60L; set.seed(1)
bx <- matrix(rnorm(E*S,.05,.02),E,S); by <- matrix(rnorm(E*S,0,.02),E,S)
sx <- matrix(runif(E*S,.005,.02),E,S); sy <- matrix(runif(E*S,.005,.02),E,S)
for (m in c("ivw","egger","weighted_median","simple_mode","weighted_mode","penalised_weighted_median","egger_bootstrap")) for (th in c(1L, 8L)) {
  t <- system.time(fast_mr_grid(bx,by,sx,sy,methods=m,nboot=100,seed=1,threads=th))[["elapsed"]]
  cat(sprintf("%-26s threads=%d %6.3fs  %.1f us/pair\n", m, th, t, 1e6*t/E^2))
}
