# fast_mr() groups rows once and batches non-RNG methods natively; both paths
# must be identical() to the original per-group implementation.

groups_fixture <- function(seed = 11L) {
  set.seed(seed)
  sizes <- c(1L, 2L, 1L, 3L, 7L, 12L, 2L, 30L, 5L, 1L)
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
  d
}

expect_same_as_reference <- function(d, ...) {
  new <- fast_mr(d, ...)
  old <- fast_mr_reference(d, ...)
  expect_identical(new, old)
  invisible(new)
}

method_sets <- list(
  "ivw",
  c("ivw", "ivw_fe", "ivw_mre", "uwr", "sign", "wald_ratio"),
  c("ivw", "egger", "simple_median", "weighted_median"),
  c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode"),
  c("egger_bootstrap", "penalised_weighted_median", "ivw")
)

test_that("grouped fast_mr is identical to the reference for all method sets", {
  d <- groups_fixture()
  for (m in method_sets) {
    for (threads in c(1L, 4L)) {
      expect_same_as_reference(d, methods = m, nboot = 0, threads = threads)
      expect_same_as_reference(d, methods = m, nboot = 50, seed = 7, threads = threads)
    }
  }
  expect_same_as_reference(d, methods = c("ivw", "simple_mode", "weighted_mode"),
                           nboot = 40, seed = 3, phi = 0.8)
  expect_same_as_reference(d, methods = c("ivw", "penalised_weighted_median"),
                           nboot = 40, seed = 3, penk = 5)
})

test_that("duplicate SNPs and mr_keep FALSE rows are handled identically", {
  d <- groups_fixture()
  d <- rbind(d, d[c(1, 5, 6, 20), ], d[20:40, ])
  d$beta.exposure[nrow(d) - 3L] <- d$beta.exposure[nrow(d) - 3L] * 2
  set.seed(5)
  d$mr_keep <- runif(nrow(d)) > 0.25
  d$mr_keep[1:3] <- FALSE
  d$mr_keep[4] <- NA
  # a non-kept row may carry unusable values
  d$se.exposure[2] <- 0
  d <- d[sample(nrow(d)), ]
  rownames(d) <- NULL
  for (m in method_sets[c(1, 2, 4)]) {
    expect_same_as_reference(d, methods = m, nboot = 0)
    expect_same_as_reference(d, methods = m, nboot = 30, seed = 1)
  }
  d$mr_keep <- FALSE
  expect_same_as_reference(d, methods = c("ivw", "egger"), nboot = 0)
})

test_that("NA, empty, factor and missing ids give the same groups", {
  d <- groups_fixture()
  d$id.exposure[1:5] <- NA
  d$id.outcome[3:8] <- ""
  d$id.outcome[20:25] <- NA
  expect_same_as_reference(d, methods = c("ivw", "egger"), nboot = 0)
  expect_same_as_reference(d, methods = c("ivw", "weighted_median"), nboot = 20, seed = 2)
  f <- groups_fixture()
  f$id.exposure <- factor(f$id.exposure, levels = c("E3", "E2", "E1", "E0"))
  f$id.outcome <- factor(f$id.outcome)
  expect_same_as_reference(f, methods = c("ivw", "egger"), nboot = 0)
  expect_same_as_reference(f, methods = c("ivw", "simple_mode"), nboot = 20, seed = 2)
  g <- groups_fixture()
  g$id.exposure <- NULL
  expect_same_as_reference(g, methods = "ivw", nboot = 0)
  g$id.outcome <- NULL
  expect_same_as_reference(g, methods = c("ivw", "egger"), nboot = 0)
  h <- groups_fixture()
  h$id.exposure <- c("a\rb", "a")[1 + (seq_len(nrow(h)) %% 2)]
  h$id.outcome <- c("c", "b\rc")[1 + (seq_len(nrow(h)) %% 3 == 0)]
  expect_same_as_reference(h, methods = "ivw", nboot = 0)
})

test_that("custom labels, single group, tiny and empty inputs are identical", {
  d <- groups_fixture()
  d$exposure <- paste("exposure", d$id.exposure)
  d$outcome <- paste("outcome", d$id.outcome)
  expect_same_as_reference(d, methods = c("ivw", "egger", "weighted_median"), nboot = 0)
  one <- d[d$id.exposure == "E0" & d$id.outcome == "O0", ]
  expect_same_as_reference(one, methods = c("ivw", "egger"), nboot = 0)
  expect_same_as_reference(one[1, ], methods = c("ivw", "wald_ratio"), nboot = 0)
  expect_same_as_reference(d[0, ], methods = "ivw", nboot = 0)
})

test_that("output files and row names are preserved", {
  d <- groups_fixture()
  out <- fast_mr(d, methods = c("ivw", "egger"), nboot = 0)
  expect_identical(rownames(out), as.character(seq_len(nrow(out))))
  expect_identical(attributes(out)$row.names, attributes(fast_mr_reference(d, methods = c("ivw", "egger"), nboot = 0))$row.names)
})

test_that("grouped fast_mr serializes byte-identically to the reference", {
  d <- groups_fixture()
  for (m in method_sets) {
    new <- fast_mr(d, methods = m, nboot = 0)
    old <- fast_mr_reference(d, methods = m, nboot = 0)
    expect_identical(serialize(new, NULL), serialize(old, NULL))
  }
})
