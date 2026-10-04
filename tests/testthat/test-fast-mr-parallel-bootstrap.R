# Threaded bootstrap groups (fastmr_run_bootstrap_groups) must be identical()
# and representation-identical to the per-group reference for every thread count,
# seeded and unseeded, and must leave R's RNG in exactly the same state.

boot_fixture <- function(seed = 21L) {
  set.seed(seed)
  sizes <- c(1L, 2L, 3L, 4L, 1L, 3L, 30L, 2L, 12L, 3L, 7L, 50L)
  g <- rep(seq_along(sizes), sizes)
  n <- length(g)
  d <- data.frame(
    SNP = unlist(lapply(sizes, function(k) paste0("rs", sample(1000, k)))),
    id.exposure = paste0("E", (g - 1L) %% 4L),
    id.outcome = paste0("O", (g - 1L) %/% 4L),
    beta.exposure = rnorm(n, 0.2, 0.05),
    beta.outcome = rnorm(n, 0.1, 0.05),
    se.exposure = runif(n, 0.01, 0.05),
    se.outcome = runif(n, 0.01, 0.05),
    stringsAsFactors = FALSE
  )
  d$mr_keep <- TRUE
  # A zero exposure beta drops a ratio but keeps the row for Egger/IVW.
  d$beta.exposure[which(g == 7L)[3]] <- 0
  d
}

boot_method_sets <- list(
  c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode"),
  c("simple_median", "weighted_median"),
  c("weighted_median", "penalised_weighted_median", "egger_bootstrap"),
  c("egger_bootstrap", "penalised_weighted_median", "ivw"),
  c("simple_mode"),
  c("weighted_mode", "sign", "uwr", "wald_ratio", "ivw_fe", "ivw_mre"),
  c("ivw", "ivw_fe", "ivw_mre", "egger", "egger_bootstrap", "uwr", "sign",
    "simple_median", "weighted_median", "penalised_weighted_median",
    "simple_mode", "weighted_mode", "wald_ratio")
)

# Run f() from a fixed .Random.seed; return the value and the RNG state after.
with_rng <- function(f, state_seed = 99L) {
  set.seed(state_seed)
  value <- f()
  list(value = value, state = get(".Random.seed", envir = .GlobalEnv))
}

expect_boot_identical <- function(d, threads = c(1L, 2L, 4L), ...) {
  for (seed in list(NULL, 7)) {
    ref <- with_rng(function() fast_mr_reference(d, seed = seed, ...))
    for (k in threads) {
      new <- with_rng(function() fast_mr(d, seed = seed, threads = k, ...))
      expect_identical(new$value, ref$value)
      expect_same_representation(new$value, ref$value)
      expect_identical(new$state, ref$state)
    }
  }
}

test_that("threaded bootstrap methods are identical to the reference", {
  d <- boot_fixture()
  for (m in boot_method_sets) {
    expect_boot_identical(d, methods = m, nboot = 40)
  }
  expect_boot_identical(d, methods = c("simple_mode", "weighted_mode"), nboot = 25, phi = 0.7)
  expect_boot_identical(d, methods = c("penalised_weighted_median"), nboot = 25, penk = 3)
})

test_that("bootstrap batching does not change results or RNG state", {
  d <- boot_fixture()
  m <- boot_method_sets[[7]]
  ref <- with_rng(function() fast_mr_reference(d, methods = m, nboot = 30))
  ref_seeded <- fast_mr_reference(d, methods = m, nboot = 30, seed = 4)
  old <- options(fastMR.bootstrap_batch_draws = 500)
  on.exit(options(old))
  for (k in c(1L, 3L)) {
    new <- with_rng(function() fast_mr(d, methods = m, nboot = 30, threads = k))
    expect_identical(new, ref)
    expect_identical(fast_mr(d, methods = m, nboot = 30, seed = 4, threads = k), ref_seeded)
  }
})

test_that("mr_keep, duplicates, NA ids and tiny groups match under threads", {
  d <- boot_fixture()
  d <- rbind(d, d[c(1, 5, 6, 20), ], d[20:40, ])
  d$beta.exposure[nrow(d) - 3L] <- d$beta.exposure[nrow(d) - 3L] * 2
  set.seed(5)
  d$mr_keep <- runif(nrow(d)) > 0.25
  d$mr_keep[1:3] <- FALSE
  d$mr_keep[4] <- NA
  d$se.exposure[2] <- 0
  d$id.exposure[30:33] <- NA
  d$exposure <- paste("exposure", d$id.exposure)
  d <- d[sample(nrow(d)), ]
  rownames(d) <- NULL
  expect_boot_identical(d, methods = boot_method_sets[[1]], nboot = 30)
  expect_boot_identical(d, methods = boot_method_sets[[3]], nboot = 30)
  # A group whose rows are all mr_keep = FALSE.
  e <- boot_fixture()
  e$mr_keep[e$id.exposure == "E1" & e$id.outcome == "O0"] <- FALSE
  expect_boot_identical(e, methods = boot_method_sets[[1]], nboot = 20)
  # Single group, and only groups too small to bootstrap.
  one <- e[e$id.exposure == "E2" & e$id.outcome == "O1", ]
  expect_boot_identical(one, methods = boot_method_sets[[7]], nboot = 20)
  expect_boot_identical(e[1:3, ], methods = boot_method_sets[[1]], nboot = 20)
})

test_that("unseeded runs without .Random.seed create one as before", {
  d <- boot_fixture()
  if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  fast_mr(d[1:3, ], methods = "weighted_median", nboot = 10, threads = 2)
  expect_true(exists(".Random.seed", envir = .GlobalEnv))
  rm(".Random.seed", envir = .GlobalEnv)
  fast_mr(d, methods = "weighted_median", nboot = 10, seed = 1, threads = 2)
  expect_false(exists(".Random.seed", envir = .GlobalEnv))
})

test_that("streamed and buffered batches mix without changing results or RNG state", {
  # Budgets that force every group to stream serially (0), a mix of streamed
  # single-group batches and buffered multi-group batches (600), and one
  # buffered batch (Inf), for each thread count, seeded and unseeded.
  d <- boot_fixture()
  for (budget in c(0, 600, Inf)) {
    old <- options(fastMR.bootstrap_batch_draws = budget)
    expect_boot_identical(d, threads = c(1L, 2L, 5L), methods = boot_method_sets[[1]], nboot = 20)
    expect_boot_identical(d, threads = c(1L, 3L), methods = boot_method_sets[[7]], nboot = 15)
    options(old)
  }
})

test_that("double-buffered batches keep the serial draws, results and RNG state", {
  # Tiny budgets give many buffered batches, so drawing batch b + 1 overlaps
  # computing batch b throughout; every RNG kind takes the same pipeline.
  d <- boot_fixture()
  d <- rbind(d, transform(d, id.outcome = paste0(id.outcome, "b")))
  old_kind <- RNGkind()
  on.exit(RNGkind(old_kind[1], old_kind[2], old_kind[3]), add = TRUE)
  for (kind in c("Inversion", "Box-Muller")) {
    RNGkind(normal.kind = kind)
    for (budget in c(150, 2000)) {
      old <- options(fastMR.bootstrap_batch_draws = budget)
      expect_boot_identical(d, threads = c(1L, 2L, 4L, 8L), methods = boot_method_sets[[1]], nboot = 12)
      previous <- fastMR:::fastmr_set_bootstrap_overlap_native(FALSE)
      expect_boot_identical(d, threads = c(2L, 8L), methods = boot_method_sets[[1]], nboot = 12)
      fastMR:::fastmr_set_bootstrap_overlap_native(previous)
      options(old)
    }
  }
})
