# The batched diagnostics must reproduce the former per-group implementations
# (helper-old-diagnostics.R) exactly: values, types, column and row order,
# attributes, row.names form, and the error raised first.

batched_base <- function(sizes, seed = 1L, labels = TRUE, samplesize = "samplesize") {
  set.seed(seed)
  G <- length(sizes)
  n <- sum(sizes)
  group <- rep(seq_len(G), sizes)
  d <- data.frame(
    SNP = paste0("rs", unlist(lapply(sizes, seq_len))),
    beta.exposure = rnorm(n, 0.1, 0.05),
    se.exposure = runif(n, 0.01, 0.03),
    beta.outcome = rnorm(n, 0.03, 0.04),
    se.outcome = runif(n, 0.01, 0.03),
    id.exposure = paste0("E", (group - 1L) %/% 2L + 1L),
    id.outcome = paste0("O", (group - 1L) %% 2L + 1L),
    mr_keep = TRUE,
    stringsAsFactors = FALSE
  )
  if (labels) {
    d$exposure <- paste0("Exposure ", (group - 1L) %/% 2L + 1L)
    d$outcome <- paste0("Outcome ", (group - 1L) %% 2L + 1L)
  }
  if (!is.null(samplesize)) d[[samplesize]] <- 1000 + group
  d
}

batched_cases <- function() {
  cases <- list()
  cases$many <- batched_base(rep(c(5L, 8L, 12L, 4L), 6L))
  cases$one <- batched_base(15L)
  cases$tiny <- batched_base(c(1L, 2L, 3L, 1L, 2L, 3L, 4L), seed = 2L)
  cases$single_row <- batched_base(1L, seed = 3L)
  dup <- batched_base(c(6L, 5L, 4L), seed = 4L)
  dup <- rbind(dup, dup[c(2, 8, 8), ], dup[1, ])
  dup$beta.outcome[nrow(dup) - 1L] <- 0.5
  cases$duplicates <- dup
  bad <- batched_base(c(6L, 7L, 5L), seed = 5L)
  bad$mr_keep[c(2, 9, 10)] <- FALSE
  bad$beta.exposure[2] <- NA
  bad$se.outcome[9] <- Inf
  bad$se.exposure[10] <- 0
  bad$mr_keep[14] <- NA
  bad$beta.exposure[16] <- 0
  # A dropped first row with a SNP that is kept later, and an unkept group.
  bad$mr_keep[1] <- FALSE
  bad$SNP[1] <- bad$SNP[3]
  cases$nonfinite <- bad
  allbad <- batched_base(c(4L, 5L), seed = 6L)
  allbad$mr_keep[1:4] <- FALSE
  cases$all_dropped <- allbad
  na_snp <- batched_base(c(5L, 6L, 4L), seed = 7L)
  na_snp$mr_keep[c(1, 8, 13)] <- FALSE
  na_snp$SNP[c(1, 8)] <- NA
  na_snp$SNP[13] <- ""
  cases$na_snp <- na_snp
  na_blank <- na_snp
  na_blank$id.exposure <- NULL
  na_blank$id.outcome <- NULL
  cases$na_snp_no_ids <- na_blank
  first_dup <- batched_base(c(4L, 4L), seed = 8L)
  first_dup$SNP[2] <- first_dup$SNP[1]
  first_dup$mr_keep[2] <- FALSE
  first_dup$SNP[3] <- NA
  first_dup$mr_keep[3] <- FALSE
  cases$na_after_first <- first_dup
  cases$no_labels <- batched_base(c(5L, 3L, 6L), seed = 9L, labels = FALSE)
  cases$no_samplesize <- batched_base(c(5L, 3L, 1L), seed = 10L, samplesize = NULL)
  cases$samplesize_outcome <- batched_base(c(4L, 6L), seed = 11L,
                                           samplesize = "samplesize.outcome")
  cases$sample_size <- batched_base(c(4L, 6L), seed = 12L, samplesize = "sample_size")
  odd_size <- batched_base(c(4L, 5L, 3L), seed = 13L)
  odd_size$samplesize <- as.character(odd_size$samplesize)
  odd_size$samplesize[c(1, 6)] <- c("n/a", NA)
  cases$character_samplesize <- odd_size
  no_ids <- batched_base(c(7L), seed = 14L, labels = FALSE)
  no_ids$id.exposure <- NULL
  no_ids$id.outcome <- NULL
  cases$no_ids <- no_ids
  only_outcome <- batched_base(c(4L, 5L, 6L), seed = 15L)
  only_outcome$id.exposure <- NULL
  cases$only_outcome_id <- only_outcome
  na_ids <- batched_base(c(4L, 5L, 6L, 3L), seed = 16L)
  na_ids$id.exposure[c(1, 2)] <- NA
  na_ids$exposure[c(5, 10)] <- NA
  na_ids$outcome[11] <- NA
  cases$na_ids <- na_ids
  fac <- batched_base(c(5L, 4L, 6L), seed = 17L)
  for (column in c("SNP", "id.exposure", "id.outcome", "exposure", "outcome")) {
    fac[[column]] <- factor(fac[[column]])
  }
  fac$samplesize <- factor(fac$samplesize)
  fac$beta.outcome <- factor(format(fac$beta.outcome, digits = 6))
  fac$mr_keep <- factor(ifelse(seq_len(nrow(fac)) %% 5L == 0L, "FALSE", "TRUE"))
  cases$factors <- fac
  shuffled <- batched_base(rep(c(6L, 3L, 9L), 4L), seed = 18L)
  set.seed(19)
  cases$interleaved <- shuffled[sample.int(nrow(shuffled)), ]
  no_keep <- batched_base(c(5L, 6L), seed = 20L)
  no_keep$mr_keep <- NULL
  cases$no_mr_keep <- no_keep
  cases$large_pair <- batched_base(c(120L, 3L), seed = 21L)
  cases$empty <- batched_base(c(3L, 2L), seed = 22L)[0, ]
  cases
}

batched_with_steiger <- function(d) {
  set.seed(23)
  n <- nrow(d)
  d$r.exposure <- runif(n, 0.01, 0.1)
  d$r.outcome <- runif(n, 0.001, 0.05)
  d$pval.exposure <- 10^-runif(n, 5, 20)
  d$pval.outcome <- runif(n)
  d$samplesize.exposure <- 5e4 + seq_len(n)
  d$samplesize.outcome <- rep(2e4, n)
  if (n) {
    d$r.exposure[seq(1, n, by = 3)] <- NA
    d$r.outcome[seq(2, n, by = 4)] <- NA
    d$samplesize.exposure[seq(1, n, by = 5)] <- NA
  }
  d
}

batched_outcome <- function(expr) {
  tryCatch(expr, error = function(e) structure(list(message = conditionMessage(e)),
                                               class = "batched_error"))
}

expect_batched_same <- function(new, old, label) {
  if (inherits(old, "batched_error") || inherits(new, "batched_error")) {
    expect_identical(new, old, label = label)
  } else {
    expect_same_representation(new, old)
  }
}

batched_call <- function(new, old, ...) {
  extra <- list(...)
  list(new = function(d, threads) do.call(new, c(list(d), extra, list(threads = threads))),
       old = function(d, threads = 1) do.call(old, c(list(d), extra, list(threads = threads))))
}

batched_steiger_call <- function(drop = NULL) {
  prepare <- function(d) {
    if (!is.data.frame(d)) return(d)
    d <- batched_with_steiger(d)
    d[drop] <- NULL
    d
  }
  list(new = function(d, threads) fast_mr_directionality_test(prepare(d)),
       old = function(d, threads = 1) old_fast_mr_directionality_test(prepare(d)))
}

batched_calls <- list(
  heterogeneity = batched_call(fast_mr_heterogeneity, old_fast_mr_heterogeneity),
  heterogeneity_all = batched_call(
    fast_mr_heterogeneity, old_fast_mr_heterogeneity,
    methods = c("uwr", "ivw_mre", "egger", "ivw", "ivw_fe")),
  pleiotropy = batched_call(fast_mr_pleiotropy_test, old_fast_mr_pleiotropy_test),
  singlesnp = batched_call(fast_mr_singlesnp, old_fast_mr_singlesnp),
  singlesnp_methods = batched_call(
    fast_mr_singlesnp, old_fast_mr_singlesnp,
    all_method = c("weighted_median", "mr_ivw_fe", "simple_mode")),
  loo_ivw = batched_call(fast_mr_leaveoneout, old_fast_mr_leaveoneout),
  loo_ivw_fe = batched_call(fast_mr_leaveoneout, old_fast_mr_leaveoneout, method = "ivw_fe"),
  loo_ivw_mre = batched_call(fast_mr_leaveoneout, old_fast_mr_leaveoneout, method = "ivw_mre"),
  loo_uwr = batched_call(fast_mr_leaveoneout, old_fast_mr_leaveoneout, method = "uwr"),
  loo_egger = batched_call(fast_mr_leaveoneout, old_fast_mr_leaveoneout,
                           method = "mr_egger_regression"),
  directionality = batched_steiger_call(),
  directionality_pn = batched_steiger_call("r.exposure")
)

# The former Egger leave-one-out read the sample size of an empty frame when a
# pair had a single SNP (and no other rows), which errored; it is now NA.
batched_old_crash <- function(old) {
  inherits(old, "batched_error") && grepl("subscript out of bounds", old$message)
}

test_that("batched diagnostics are identical to the per-group implementations", {
  cases <- batched_cases()
  for (case in names(cases)) {
    d <- cases[[case]]
    for (call in names(batched_calls)) {
      label <- paste(case, call)
      old <- batched_outcome(batched_calls[[call]]$old(d))
      new <- batched_outcome(batched_calls[[call]]$new(d, 1))
      if (call == "loo_egger" && batched_old_crash(old)) {
        # Compare with the former output without the sample-size column.
        sizes <- intersect(names(d), c("samplesize.outcome", "samplesize", "sample_size"))
        old <- batched_calls[[call]]$old(d[setdiff(names(d), sizes)])
        expect_true(anyNA(new$samplesize))
        old$samplesize <- new$samplesize
      }
      expect_batched_same(new, old, label)
      threaded <- batched_outcome(batched_calls[[call]]$new(d, 4))
      expect_batched_same(threaded, new, paste(label, "threads"))
    }
  }
})

test_that("Egger leave-one-out of single-SNP pairs reports NA sample sizes", {
  d <- batched_base(c(1L, 4L, 1L), seed = 30L)
  expect_error(old_fast_mr_leaveoneout(d, method = "egger"), "subscript out of bounds")
  new <- fast_mr_leaveoneout(d, method = "egger")
  d$samplesize <- NULL
  reference <- old_fast_mr_leaveoneout(d, method = "egger")
  expect_identical(new$samplesize, c(NA, 1001, rep(1002, 5), NA, 1003))
  reference$samplesize <- new$samplesize
  expect_same_representation(new, reference)
})

test_that("batched diagnostics raise the error the per-group loop raised first", {
  d <- batched_base(c(4L, 5L, 6L), seed = 31L)
  value_late <- d; value_late$beta.exposure[12] <- NA
  snp_early <- value_late; snp_early$SNP[2] <- ""
  value_early <- snp_early; value_early$se.outcome[3] <- -1
  missing_column <- d; missing_column$se.outcome <- NULL
  not_numeric <- d; not_numeric$beta.outcome <- "x"
  inputs <- list(clean = d, value_late = value_late, snp_early = snp_early,
                 value_early = value_early, missing = missing_column,
                 not_numeric = not_numeric, list = as.list(d), empty = d[0, ])
  for (input in names(inputs)) {
    for (threads in list(1, 0, "x")) {
      for (call in names(batched_calls)) {
        old <- batched_outcome(batched_calls[[call]]$old(inputs[[input]], threads))
        new <- batched_outcome(batched_calls[[call]]$new(inputs[[input]], threads))
        expect_batched_same(new, old, paste(input, call, threads))
      }
    }
  }
})

test_that("RNG-free fast_mr() groups are identical for every thread count", {
  d <- batched_base(rep(c(3L, 7L, 12L, 25L, 2L), 8L), seed = 32L)
  methods <- c("ivw", "ivw_fe", "ivw_mre", "egger", "uwr", "sign", "simple_median",
               "weighted_median", "penalised_weighted_median", "simple_mode",
               "weighted_mode", "wald_ratio")
  for (method in c(as.list(methods), list(methods))) {
    serial <- fast_mr(d, methods = method, nboot = 0, threads = 1)
    for (threads in c(2, 8)) {
      expect_same_representation(fast_mr(d, methods = method, nboot = 0, threads = threads),
                                 serial)
    }
  }
  # RNG-free methods with nboot > 0 also take the batched path.
  expect_same_representation(fast_mr(d, methods = c("ivw", "egger"), nboot = 100, threads = 8),
                             fast_mr(d, methods = c("ivw", "egger"), nboot = 100, threads = 1))
})

test_that("drop-one native fits equal fits of the reduced groups", {
  d <- batched_base(c(6L, 4L, 1L), seed = 33L)
  offsets <- c(0L, 6L, 10L, 11L)
  args <- list(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome)
  methods <- c("egger", "ivw", "weighted_median", "sign")
  jobs <- list(c(1L, 0L), c(1L, 5L), c(2L, -1L), c(2L, 2L), c(3L, 0L), c(NA, -1L))
  native <- do.call(fastmr_run_groups_drop_native, c(list(offsets), args, list(
    vapply(jobs, `[`, integer(1), 1L), vapply(jobs, `[`, integer(1), 2L), methods,
    threads = 3L)))
  for (j in seq_along(jobs)) {
    g <- jobs[[j]][[1L]]
    rows <- if (is.na(g)) integer() else seq.int(offsets[g] + 1L, offsets[g + 1L])
    if (!is.na(g) && jobs[[j]][[2L]] >= 0L) rows <- rows[-(jobs[[j]][[2L]] + 1L)]
    one <- fastmr_run_groups_native(c(0L, length(rows)), d$beta.exposure[rows],
                                    d$beta.outcome[rows], d$se.exposure[rows],
                                    d$se.outcome[rows], methods)
    k <- (j - 1L) * length(methods) + seq_along(methods)
    expect_identical(lapply(native, `[`, k), one)
  }
  expect_error(fastmr_run_groups_drop_native(offsets, args[[1]], args[[2]], args[[3]], args[[4]],
                                             4L, 0L, "ivw"), "job_group")
  expect_error(fastmr_run_groups_drop_native(offsets, args[[1]], args[[2]], args[[3]], args[[4]],
                                             3L, 1L, "ivw"), "job_drop")
})

test_that("grouped native sums and means equal sum() and mean()", {
  set.seed(34)
  x <- c(rnorm(50) * 10^runif(50, -10, 10), NA, NaN, 1e308, 1e308, -Inf, 3)
  x <- x[sample.int(length(x))]
  offsets <- c(0L, 0L, 1L, 7L, 20L, 33L, 56L)
  groups <- lapply(seq_len(length(offsets) - 1L),
                   function(i) x[seq_len(offsets[i + 1L] - offsets[i]) + offsets[i]])
  expect_identical(fastmr_group_sum_native(offsets, x, TRUE),
                   vapply(groups, sum, numeric(1), na.rm = TRUE))
  expect_identical(fastmr_group_sum_native(offsets, x, FALSE),
                   vapply(groups, sum, numeric(1)))
  expect_identical(fastmr_group_mean_native(offsets, x),
                   vapply(groups, mean, numeric(1), na.rm = TRUE))
  big <- c(1e308, 1.5e308, 1e308, NA)
  expect_identical(fastmr_group_mean_native(c(0L, 4L), big), mean(big, na.rm = TRUE))
  n <- 5e4 + seq_len(1000) / 3
  expect_identical(fastmr_group_mean_native(c(0L, 1000L), n), mean(n))
})
