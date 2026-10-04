# Equivalence battery: run every case with one installed fastMR build and save
# each result together with the .Random.seed left behind.
#   Rscript equiv-run.R <lib> <out.rds> [quick]
# Compare two outputs with equiv-compare.R. Every case starts from a fixed
# .Random.seed, so the two builds see the same RNG stream.
args <- commandArgs(TRUE)
.libPaths(c(args[1], .libPaths()))
suppressMessages(library(fastMR))
quick <- length(args) >= 3 && args[3] == "quick"
has_hull <- exists("fastmr_set_mode_hull_native", asNamespace("fastMR"))
# FASTMR_HULL="<max ratios>,<recurrence 0/1>" overrides the hull settings.
hull_env <- Sys.getenv("FASTMR_HULL")
if (has_hull && nzchar(hull_env)) {
  h <- as.numeric(strsplit(hull_env, ",")[[1]])
  fastMR:::fastmr_set_mode_hull_native(h[1], as.logical(h[2]))
}

five <- c("wald_ratio", "egger", "weighted_median", "ivw", "weighted_mode")
all_methods <- c("ivw", "ivw_fe", "ivw_mre", "egger", "egger_bootstrap", "uwr", "sign",
                 "simple_median", "weighted_median", "penalised_weighted_median",
                 "simple_mode", "weighted_mode", "wald_ratio")

# Groups of the given sizes; scenario = clean, outliers, ties.
make_data <- function(sizes, scenario, seed) {
  set.seed(seed)
  g <- rep(seq_along(sizes), sizes)
  n <- length(g)
  bx <- sample(c(-1, 1), n, TRUE) * runif(n, 0.03, 0.12)
  sx <- runif(n, 0.004, 0.012)
  sy <- runif(n, 0.008, 0.02)
  by <- rnorm(length(sizes), 0, 0.3)[g] * bx + rnorm(n, 0, sy)
  if (scenario == "outliers") {
    hit <- runif(n) < 0.3
    by[hit] <- by[hit] + sample(c(-1, 1), sum(hit), TRUE) * runif(sum(hit), 0.05, 0.5)
  }
  if (scenario == "ties") {
    # Exactly repeated ratios and SEs within each group, plus mirrored ratios
    # with equal SEs (two density peaks that tie in exact arithmetic).
    for (j in seq_along(sizes)) {
      rows <- which(g == j)
      k <- length(rows)
      if (k >= 2) {
        half <- rows[seq_len(k %/% 2)]
        a <- seq(0.1, 1, length.out = length(half))
        r <- c(0.37 - a, 0.37 + a)
        rows2 <- rows[seq_along(r)]
        bx[rows2] <- 0.1
        by[rows2] <- 0.1 * r
        sx[rows2] <- 1e-9
        sy[rows2] <- 1e-9
        if (k %% 2 == 1) { bx[rows[k]] <- 0.1; by[rows[k]] <- 0.1 * 0.37; sx[rows[k]] <- 1e-9; sy[rows[k]] <- 1e-9 }
      }
    }
  }
  data.frame(SNP = paste0("rs", seq_len(n)), id.exposure = paste0("E", g), id.outcome = "O",
             beta.exposure = bx, se.exposure = sx, beta.outcome = by, se.outcome = sy,
             stringsAsFactors = FALSE)
}

cases <- list()
add <- function(name, f) cases[[name]] <<- f
ks <- c(2:12, 15, 16, 17, 20, 24, 30, 32, 40, 50, 63, 64, 65, 80, 100, 128, 200, 300, 500)
if (quick) ks <- c(2, 3, 5, 10, 30, 64, 65, 100, 500)
threads_set <- c(1L, 2L, 4L, 8L)
for (scenario in c("clean", "outliers", "ties")) {
  # One group per k (k = 2..500) plus a mixed-size batch, five-method set.
  for (k in ks) {
    G <- if (k <= 10) 6 else if (k <= 64) 3 else 1
    nb <- if (k >= 200) 200 else 500
    for (seed in list(NULL, 11L)) for (th in threads_set) {
      add(sprintf("%s/k%d/G%d/seed%s/t%d", scenario, k, G, format(seed), th), local({
        k <- k; G <- G; seed <- seed; th <- th; nb <- nb; scenario <- scenario
        function() fast_mr(make_data(rep(k, G), scenario, 100 + k), methods = five,
                           nboot = nb, seed = seed, threads = th)
      }))
    }
  }
  mixed <- c(1, 2, 3, 4, 5, 3, 7, 10, 2, 30, 3, 64, 65, 12, 100, 3, 5, 8, 1, 200)
  for (m in list(five, all_methods, c("simple_mode", "weighted_mode"), "weighted_mode"))
    for (seed in list(NULL, 7L)) for (th in threads_set) for (phi in c(1, 0.5)) {
      add(sprintf("%s/mixed/%s/seed%s/t%d/phi%g", scenario, paste(m, collapse = "+"), format(seed), th, phi), local({
        m <- m; seed <- seed; th <- th; phi <- phi; scenario <- scenario
        function() fast_mr(make_data(mixed, scenario, 5), methods = m, nboot = 300,
                           seed = seed, threads = th, phi = phi)
      }))
    }
}
# Many batches (double-buffered pipeline with tiny buffers) and the forced
# parallel paths on small inputs.
for (bd in c(2000, 50000, 2^23)) for (ws in c(0, 1)) for (seed in list(NULL, 3L)) for (th in threads_set) {
  add(sprintf("batches/bd%g/ws%g/seed%s/t%d", bd, ws, format(seed), th), local({
    bd <- bd; ws <- ws; seed <- seed; th <- th
    function() {
      old <- options(fastMR.bootstrap_batch_draws = bd); on.exit(options(old), add = TRUE)
      prev <- fastMR:::fastmr_set_work_scale_native(ws); on.exit(fastMR:::fastmr_set_work_scale_native(prev), add = TRUE)
      fast_mr(make_data(rep(c(3, 10, 30, 2, 5), 40), "outliers", 9), methods = all_methods,
              nboot = 100, seed = seed, threads = th)
    }
  }))
}
# Every normal.kind / RNG kind on the batched path.
for (kind in list(c("Mersenne-Twister", "Inversion"), c("Mersenne-Twister", "Box-Muller"),
                  c("Mersenne-Twister", "Kinderman-Ramage"), c("Mersenne-Twister", "Ahrens-Dieter"),
                  c("L'Ecuyer-CMRG", "Inversion"), c("Knuth-TAOCP-2002", "Inversion"),
                  c("Wichmann-Hill", "Box-Muller")))
  for (seed in list(NULL, 5L)) for (th in threads_set) {
    add(sprintf("rngkind/%s/seed%s/t%d", paste(kind, collapse = "+"), format(seed), th), local({
      kind <- kind; seed <- seed; th <- th
      function() {
        old <- RNGkind(); on.exit(RNGkind(old[1], old[2], old[3]), add = TRUE)
        suppressWarnings(RNGkind(kind[1], kind[2]))
        set.seed(42)
        prev <- fastMR:::fastmr_set_work_scale_native(0); on.exit(fastMR:::fastmr_set_work_scale_native(prev), add = TRUE)
        r <- fast_mr(make_data(rep(c(3, 8, 20), 10), "clean", 4), methods = five,
                     nboot = 100, seed = seed, threads = th)
        list(r, .Random.seed)
      }
    }))
  }
# Seed sweep: 40 random datasets and seeds (seeded and unseeded) per thread count.
for (s in 1:40) for (th in threads_set) {
  add(sprintf("seeds/s%d/t%d", s, th), local({
    s <- s; th <- th
    function() {
      set.seed(1000 + s)
      sizes <- sample(c(2:12, 20, 30, 64, 65, 100), sample(5:25, 1), TRUE)
      scenario <- sample(c("clean", "outliers", "ties"), 1)
      seed <- if (s %% 2 == 0) s else NULL
      set.seed(s * 7919)
      fast_mr(make_data(sizes, scenario, s), methods = c(five, "simple_mode", "egger_bootstrap"),
              nboot = 150, seed = seed, threads = th)
    }
  }))
}
# Shared grid (OpenMP pair loop) and single-pair paths.
for (k in c(3, 10, 64, 100)) for (th in threads_set) for (seed in list(NULL, 2L)) {
  add(sprintf("grid/k%d/seed%s/t%d", k, format(seed), th), local({
    k <- k; th <- th; seed <- seed
    function() {
      set.seed(k)
      E <- 3; O <- 2
      bx <- matrix(sample(c(-1, 1), E * k, TRUE) * runif(E * k, 0.03, 0.12), E, k)
      by <- matrix(rnorm(O * k, 0, 0.02), O, k)
      fast_mr_grid(bx, by, matrix(runif(E * k, 0.004, 0.012), E, k), matrix(runif(O * k, 0.008, 0.02), O, k),
                   methods = c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode"),
                   nboot = 200, seed = seed, threads = th)
    }
  }))
}
for (k in c(3, 4, 10, 64, 65, 500)) for (scenario in c("clean", "outliers", "ties")) {
  add(sprintf("single/%s/k%d", scenario, k), local({
    k <- k; scenario <- scenario
    function() fast_mr(make_data(k, scenario, k), methods = all_methods, nboot = 1000, threads = 1)
  }))
}

if (has_hull) fastMR:::fastmr_mode_path_counts_native(TRUE)
out <- vector("list", length(cases)); names(out) <- names(cases)
t0 <- Sys.time()
for (nm in names(cases)) {
  set.seed(20261004)
  value <- cases[[nm]]()
  out[[nm]] <- list(value = value, state = get(".Random.seed", envir = .GlobalEnv))
}
counts <- fastMR:::fastmr_mode_path_counts_native(TRUE)
saveRDS(list(results = out, counts = counts, info = list(lib = args[1], R = R.version.string,
             omp = system.file(package = "fastMR"), seconds = as.numeric(Sys.time() - t0, units = "secs"))), args[2])
cat(length(out), "cases;", format(Sys.time() - t0), "; mode path counts:",
    paste(names(counts), counts, sep = "=", collapse = " "), "\n")
