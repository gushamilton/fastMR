# Vectorised implementations must match the pre-change code exactly.

make_steiger_data <- function(G = 12L, per = 5L, seed = 7L) {
  set.seed(seed)
  n <- G * per
  grp <- rep(seq_len(G), each = per)
  d <- data.frame(
    SNP = paste0("rs", seq_len(n)),
    id.exposure = paste0("E", (grp %% 4L) + 1L),
    id.outcome = paste0("O", (grp %/% 4L) + 1L),
    exposure = paste0("exp", (grp %% 4L) + 1L),
    outcome = paste0("out", (grp %/% 4L) + 1L),
    beta.exposure = rnorm(n, 0.1, 0.05), se.exposure = runif(n, 0.01, 0.05),
    beta.outcome = rnorm(n, 0.05, 0.05), se.outcome = runif(n, 0.01, 0.05),
    eaf.exposure = runif(n, 0.1, 0.9), eaf.outcome = runif(n, 0.1, 0.9),
    samplesize.exposure = 5000 + (grp %% 4L) * 1000,
    samplesize.outcome = 8000,
    pval.exposure = runif(n, 1e-12, 0.05), pval.outcome = runif(n),
    mr_keep = TRUE, stringsAsFactors = FALSE)
  d <- d[order(d$id.exposure, d$id.outcome), ]
  rownames(d) <- NULL
  d[sample(nrow(d)), ]  # interleave groups, keep non-trivial row names
}

expect_same_steiger <- function(d) {
  suppressWarnings(
    expect_identical(fast_mr_steiger_filtering(d), old_steiger_filtering(d)))
}

with_old_groups <- function(code) {
  testthat::local_mocked_bindings(
    fastmr_diagnostic_groups = old_diagnostic_groups, .package = "fastMR")
  code
}

test_that("steiger filtering matches old code: mixed units per group", {
  d <- make_steiger_data()
  u <- c("SD", "log odds", "SD units", NA, "mg/dL")
  d$units.exposure <- u[(as.integer(factor(d$id.exposure)) %% 5L) + 1L]
  d$units.outcome <- "SD"
  d$units.outcome[d$id.outcome == "O2"] <- "log odds"
  d$prevalence.exposure <- 0.1
  d$ncase.exposure <- 1000; d$ncontrol.exposure <- 4000
  d$prevalence.outcome <- 0.05
  d$ncase.outcome <- 2000; d$ncontrol.outcome <- 6000
  expect_same_steiger(d)
  d$prevalence.outcome <- NULL
  d$ncase.outcome <- NULL
  expect_same_steiger(d)
})

test_that("steiger filtering matches old code: all input shapes", {
  d <- make_steiger_data()
  expect_same_steiger(d)                                    # no units columns
  d$units.exposure <- "SD"; d$units.outcome <- "SD"
  expect_same_steiger(d)
  d$units.exposure <- "log odds"
  d$prevalence.exposure <- 0.1; d$ncase.exposure <- 10; d$ncontrol.exposure <- 90
  expect_same_steiger(d)
  d$prevalence.exposure <- NULL
  expect_same_steiger(d)                                    # no prevalence
  d$units.exposure <- NULL; d$units.outcome <- NULL
  r <- make_steiger_data()
  r$rsq.exposure <- runif(nrow(r), -0.1, 0.1)
  r$rsq.outcome <- as.character(runif(nrow(r)))
  r$rsq.outcome[3] <- NA
  expect_same_steiger(r)                                    # supplied rsq
  m <- make_steiger_data()
  expect_same_steiger(m[setdiff(names(m), c("eaf.exposure", "eaf.outcome"))])
  expect_same_steiger(m[setdiff(names(m), c("samplesize.outcome", "pval.exposure", "units.exposure"))])
})

test_that("steiger filtering matches old code: invalid rows, NA/factor ids, singletons", {
  d <- make_steiger_data()
  d$beta.exposure[2] <- NA; d$se.outcome[5] <- 0; d$samplesize.exposure[7] <- 2
  d$pval.exposure[4] <- 0; d$eaf.outcome[9] <- 1.5
  expect_same_steiger(d)
  n <- make_steiger_data()
  n$id.exposure[c(1, 10)] <- NA; n$id.outcome[3] <- NA
  expect_error(fast_mr_steiger_filtering(n), "unique labels")
  n$exposure <- n$id.exposure; n$exposure[is.na(n$exposure)] <- "x"
  expect_identical(tryCatch(fast_mr_steiger_filtering(n), error = conditionMessage),
                   tryCatch(old_steiger_filtering(n), error = conditionMessage))
  f <- make_steiger_data()
  f$id.exposure <- factor(f$id.exposure); f$id.outcome <- factor(f$id.outcome)
  f$exposure <- factor(f$exposure)
  expect_same_steiger(f)
  s <- make_steiger_data(G = 6L, per = 1L)
  expect_same_steiger(s)                                    # single-row groups
  expect_same_steiger(make_steiger_data(G = 1L, per = 6L))  # one group
  expect_same_steiger(d[0, ])
  e <- make_steiger_data(); e$id.exposure <- NULL
  expect_identical(tryCatch(fast_mr_steiger_filtering(e), error = conditionMessage),
                   tryCatch(old_steiger_filtering(e), error = conditionMessage))
  e$exposure <- "exp"
  expect_same_steiger(e)
})

test_that("steiger filtering keeps the unique label and units error", {
  d <- make_steiger_data()
  d$units.exposure <- "SD"
  d$units.exposure[d$id.exposure == "E1"][2] <- "log odds"
  expect_error(fast_mr_steiger_filtering(d), "unique labels and units")
  expect_error(old_steiger_filtering(d), "unique labels and units")
  d2 <- make_steiger_data()
  d2$outcome[1] <- "other"
  expect_error(fast_mr_steiger_filtering(d2), "unique labels and units")
  d3 <- make_steiger_data(); d3$exposure <- NULL
  expect_error(fast_mr_steiger_filtering(d3), "unique labels and units")
  expect_error(old_steiger_filtering(d3), "unique labels and units")
})

test_that("diagnostic groups match old grouping", {
  d <- make_steiger_data()
  d$exposure[c(3, 8)] <- NA
  d$id.outcome[c(2, 6)] <- NA
  d$id.exposure[11] <- ""
  expect_identical(fastmr_diagnostic_groups(d), old_diagnostic_groups(d))
  f <- make_steiger_data(); f$id.exposure <- factor(f$id.exposure)
  expect_identical(fastmr_diagnostic_groups(f), old_diagnostic_groups(f))
  g <- make_steiger_data(); g$exposure <- NULL; g$outcome <- NULL
  expect_identical(fastmr_diagnostic_groups(g), old_diagnostic_groups(g))
  h <- g; h$id.exposure <- NULL
  expect_identical(fastmr_diagnostic_groups(h), old_diagnostic_groups(h))
  expect_identical(fastmr_diagnostic_groups(d[0, ]), old_diagnostic_groups(d[0, ]))
})

test_that("diagnostics using groups match old grouping", {
  d <- make_steiger_data(G = 8L, per = 6L)
  d$r.exposure <- runif(nrow(d), 0.01, 0.1)
  d$r.outcome <- runif(nrow(d), 0.01, 0.1)
  fns <- list(
    het = function(x) fast_mr_heterogeneity(x),
    pleio = function(x) fast_mr_pleiotropy_test(x),
    single = function(x) fast_mr_singlesnp(x),
    loo = function(x) fast_mr_leaveoneout(x),
    dir = function(x) fast_mr_directionality_test(x))
  for (f in fns) {
    new <- f(d)
    old <- with_old_groups(f(d))
    expect_identical(new, old)
  }
})
