# The vectorised pair assembly behind fast_mr_compressed()'s pairwise and
# Steiger paths must reproduce the former per-pair loop exactly
# (helper-old-compressed-pairwise.R holds verbatim copies of it).

# In-memory stand-ins for decoded stores: exposures hold (most of) their
# instruments, outcomes (most of) the instrument union, with some missing
# keys, non-finite betas and non-positive standard errors.
pairwise_assembly_study <- function(n_exposure = 6L, n_outcome = 5L, seed = 1L,
                                    problems = TRUE) {
  set.seed(seed)
  pool <- sprintf("1:%d:A:G", seq_len(60L) * 100L)
  store <- function(keys) {
    n <- length(keys)
    data.frame(
      chromosome = "1",
      beta = stats::rnorm(n, 0, 0.1),
      standard_error = stats::runif(n, 0.01, 0.05),
      effect_allele_frequency = stats::runif(n, 0.05, 0.95),
      p_value = stats::runif(n),
      variant_key = keys,
      stringsAsFactors = FALSE
    )
  }
  instruments <- lapply(seq_len(n_exposure), function(e) {
    sample(pool, sample(2:8, 1L))
  })
  names(instruments) <- sprintf("exp%02d", seq_len(n_exposure))
  exposure_data <- lapply(instruments, function(keys) {
    present <- if (problems) keys[stats::runif(length(keys)) > 0.15] else keys
    store(sample(present))
  })
  union_keys <- unique(unlist(instruments, use.names = FALSE))
  outcome_data <- lapply(seq_len(n_outcome), function(o) {
    present <- if (problems) union_keys[stats::runif(length(union_keys)) > 0.1] else union_keys
    store(sample(present))
  })
  names(outcome_data) <- sprintf("out%02d", seq_len(n_outcome))
  if (problems) {
    exposure_data[[2L]]$beta[[1L]] <- NA_real_
    exposure_data[[3L]]$standard_error[[1L]] <- 0
    outcome_data[[1L]]$beta[1:3] <- Inf
    outcome_data[[2L]]$standard_error[2:4] <- -1
    # One exposure with no instrument found at all.
    exposure_data[[n_exposure]] <- exposure_data[[n_exposure]][0, , drop = FALSE]
  }
  list(exposure_data = exposure_data, outcome_data = outcome_data,
       instruments = instruments)
}

capture_run <- function(expr) {
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(expr, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) structure(conditionMessage(e), class = "captured_error")
  )
  list(value = value, warnings = warnings)
}

test_that("pairwise assembly is identical to the per-pair loop", {
  withr::local_options(fastMR.warning_pairs = Inf)
  for (seed in 1:6) {
    for (problems in c(TRUE, FALSE)) {
      study <- pairwise_assembly_study(seed = seed, problems = problems)
      for (strict in c(FALSE, TRUE)) {
        for (minimum_snps in c(1L, 3L, 6L)) {
          run <- function(f) capture_run(f(
            study$exposure_data, study$outcome_data, study$instruments,
            minimum_snps, strict
          ))
          old <- run(old_compressed_pairwise_data)
          new <- run(fastMR:::fastmr_compressed_pairwise_data)
          if (!inherits(new$value, "captured_error")) new$value$index <- NULL
          expect_identical(new, old, info = paste(seed, problems, strict, minimum_snps))
        }
      }
    }
  }
})

test_that("strict mode reports the first failing pair the loop met", {
  study <- pairwise_assembly_study(seed = 3L, problems = FALSE)
  instruments <- study$instruments
  scenarios <- list(
    exposure = function(s) {
      s$exposure_data$exp03 <- s$exposure_data$exp03[-1L, ]
      s
    },
    outcome_before_exposure = function(s) {
      s$outcome_data$out02 <- s$outcome_data$out02[
        !s$outcome_data$out02$variant_key %in% s$instruments$exp01, ]
      s$exposure_data$exp04 <- s$exposure_data$exp04[-1L, ]
      s
    },
    invalid = function(s) {
      hit <- match(s$instruments$exp02[[1L]], s$outcome_data$out03$variant_key)
      s$outcome_data$out03$standard_error[[hit]] <- NA_real_
      s
    }
  )
  for (name in names(scenarios)) {
    s <- scenarios[[name]](study)
    for (minimum_snps in c(1L, 5L)) {
      run <- function(f) capture_run(f(
        s$exposure_data, s$outcome_data, s$instruments, minimum_snps, TRUE
      ))
      old <- run(old_compressed_pairwise_data)
      new <- run(fastMR:::fastmr_compressed_pairwise_data)
      expect_s3_class(old$value, "captured_error")
      expect_identical(new, old, info = paste(name, minimum_snps))
    }
  }
})

test_that("fast_mr on the assembled table matches for every method set and thread count", {
  withr::local_options(fastMR.warning_pairs = Inf)
  study <- pairwise_assembly_study(n_exposure = 5L, n_outcome = 4L, seed = 11L)
  old <- suppressWarnings(old_compressed_pairwise_data(
    study$exposure_data, study$outcome_data, study$instruments, 1L, FALSE
  ))
  new <- suppressWarnings(fastMR:::fastmr_compressed_pairwise_data(
    study$exposure_data, study$outcome_data, study$instruments, 1L, FALSE
  ))
  method_sets <- list(
    "ivw",
    c("wald_ratio", "egger", "weighted_median", "ivw", "weighted_mode"),
    fastMR:::fastmr_method_registry()$code
  )
  for (methods in method_sets) {
    for (threads in c(1L, 8L)) {
      fit <- function(data) fast_mr(data, methods = methods, nboot = 40,
                                    seed = 7, threads = threads)
      expect_identical(fit(new$data), fit(old$data),
                       info = paste(paste(methods, collapse = ","), threads))
    }
  }
})

test_that("compressed Steiger assembly is identical to the per-pair loop", {
  for (seed in 1:4) {
    study <- pairwise_assembly_study(seed = seed)
    options <- list(
      samplesize_exposure = stats::setNames(
        seq(1e5, by = 1e4, length.out = length(study$exposure_data)),
        names(study$exposure_data)
      ),
      samplesize_outcome = stats::setNames(
        rep(5e4, length(study$outcome_data)), names(study$outcome_data)
      ),
      binary = if (seed %% 2L) NULL else data.frame(
        id = "out02", ncase = 3000, ncontrol = 27000, prevalence = 0.1,
        stringsAsFactors = FALSE
      )
    )
    if (!is.null(options$binary)) options$samplesize_outcome[["out02"]] <- NA_real_
    for (minimum_snps in c(1L, 4L, 50L)) {
      run <- function(f, ...) f(study$exposure_data, study$outcome_data,
                                study$instruments, minimum_snps, options, ...)
      expected <- run(old_compressed_steiger)
      expect_identical(run(fastMR:::fastmr_compressed_steiger), expected,
                       info = paste(seed, minimum_snps))
      index <- fastMR:::fastmr_compressed_pair_index(
        study$exposure_data, study$outcome_data, study$instruments
      )
      expect_identical(run(fastMR:::fastmr_compressed_steiger, index = index),
                       expected, info = paste(seed, minimum_snps, "index"))
    }
  }
})

test_that("a store missing a gathered column is an error, not a misaligned table", {
  study <- pairwise_assembly_study(seed = 2L, problems = FALSE)
  broken <- study
  broken$outcome_data$out03$standard_error <- NULL
  expect_error(
    fastMR:::fastmr_compressed_pairwise_data(
      broken$exposure_data, broken$outcome_data, broken$instruments, 1L, FALSE
    ),
    "outcome store\\(s\\) lack a complete 'standard_error' column: out03"
  )
  broken <- study
  broken$exposure_data$exp02$effect_allele_frequency <- NULL
  options <- list(
    samplesize_exposure = stats::setNames(rep(1e5, 6), names(study$exposure_data)),
    samplesize_outcome = stats::setNames(rep(5e4, 5), names(study$outcome_data)),
    binary = NULL
  )
  expect_error(
    fastMR:::fastmr_compressed_steiger(
      broken$exposure_data, broken$outcome_data, broken$instruments, 1L, options
    ),
    "exposure store\\(s\\) lack a complete 'effect_allele_frequency' column: exp02"
  )
})

test_that("long omission warnings are truncated with the total count", {
  study <- pairwise_assembly_study(n_exposure = 12L, n_outcome = 10L, seed = 5L)
  run <- function() capture_run(fastMR:::fastmr_compressed_pairwise_data(
    study$exposure_data, study$outcome_data, study$instruments, 3L, FALSE
  ))
  full <- withr::with_options(list(fastMR.warning_pairs = Inf), run())
  full$value$index <- NULL
  old <- capture_run(old_compressed_pairwise_data(
    study$exposure_data, study$outcome_data, study$instruments, 3L, FALSE
  ))
  expect_identical(full, old)
  short <- withr::with_options(list(fastMR.warning_pairs = 3), run())
  short$value$index <- NULL
  expect_identical(short$value, full$value)
  expect_length(short$warnings, length(full$warnings))
  for (i in seq_along(full$warnings)) {
    entries <- strsplit(sub("^[^:]*: ", "", full$warnings[[i]]), "; ")[[1L]]
    prefix <- sub(": .*$", ": ", full$warnings[[i]])
    expected <- if (length(entries) > 3L) {
      paste0(prefix, paste(entries[1:3], collapse = "; "), "; ... and ",
             length(entries) - 3L, " more (", length(entries),
             " in total; per-pair counts are in ",
             "attr(result, \"compressed_input\")$counts)")
    } else {
      full$warnings[[i]]
    }
    expect_identical(short$warnings[[i]], expected)
  }
  none <- withr::with_options(list(fastMR.warning_pairs = 0), run())
  expect_true(all(grepl(": \\.\\.\\. and [0-9]+ more \\(", none$warnings)))
  # The default limit (20) keeps short listings verbatim and truncates long ones.
  default <- run()
  long <- lengths(strsplit(full$warnings, "; ")) > 20L
  expect_true(any(long))
  expect_identical(default$warnings[!long], full$warnings[!long])
  expect_true(all(grepl("more \\([0-9]+ in total", default$warnings[long])))
  expect_s3_class(
    withr::with_options(list(fastMR.warning_pairs = -1), run())$value,
    "captured_error"
  )
})

test_that("pairwise assembly scales linearly to 50k pairs", {
  skip_on_cran()
  # 250 exposures x 200 outcomes = 50,000 pairs.  On BluePebble (Cascade
  # Lake) the per-pair loop took 58 s with a 239 MB peak R heap here; the
  # vectorised assembly took 0.09 s with a 64 MB peak.
  set.seed(42)
  n_exposure <- 250L
  n_outcome <- 200L
  pool <- sprintf("2:%d:C:T", seq_len(2000L) * 10L)
  instruments <- lapply(seq_len(n_exposure), function(e) sample(pool, 3L))
  names(instruments) <- sprintf("e%03d", seq_len(n_exposure))
  union_keys <- unique(unlist(instruments, use.names = FALSE))
  store <- function(keys) data.frame(
    beta = stats::rnorm(length(keys)), standard_error = 0.1,
    variant_key = keys, stringsAsFactors = FALSE
  )
  exposure_data <- lapply(instruments, store)
  outcome_data <- lapply(seq_len(n_outcome), function(o) store(union_keys))
  names(outcome_data) <- sprintf("o%03d", seq_len(n_outcome))
  gc(reset = TRUE)
  before <- sum(gc()[, 2L])
  elapsed <- system.time(
    out <- fastMR:::fastmr_compressed_pairwise_data(
      exposure_data, outcome_data, instruments, 1L, TRUE
    )
  )[["elapsed"]]
  usage <- gc()
  peak <- sum(usage[, ncol(usage)]) - before
  expect_identical(nrow(out$counts), 50000L)
  expect_identical(nrow(out$data), 150000L)
  expect_lt(elapsed, 10)
  # Peak R heap growth in MB (the long table itself is ~10 MB).
  expect_lt(peak, 150)
})

test_that("fast_mr_compressed pairwise path matches the per-pair loop end to end", {
  skip_if_compressor_unavailable()
  withr::local_options(fastMR.warning_pairs = Inf)
  stores <- vapply(c(1, 1.3, 0.7, -0.4, 0.2), function(multiplier) {
    path <- tempfile("fastmr-pairwise-assembly-")
    CompreSSoR::compress_sumstats(compressor_canonical_fixture(multiplier), path,
                                  overwrite = TRUE)
    path
  }, character(1))
  identity <- compressor_canonical_fixture()
  keys <- CompreSSoR::compressor_variant_key(
    identity$chromosome, identity$base_pair_location,
    identity$other_allele, identity$effect_allele
  )
  exposures <- setNames(stores[1:2], c("exposure_a", "exposure_b"))
  outcomes <- setNames(stores[3:5], c("outcome_a", "outcome_b", "outcome_c"))
  absent <- "2:200000000:A:C"
  instruments <- list(
    exposure_a = c(keys[c(2L, 7L, 14L, 25L, 40L, 61L)], absent),
    exposure_b = keys[c(3L, 8L, 19L)]
  )
  n_exposure <- c(exposure_a = 50000, exposure_b = 80000)
  n_outcome <- c(outcome_a = 120000, outcome_b = 30000, outcome_c = 60000)
  methods <- c("wald_ratio", "egger", "weighted_median", "ivw", "weighted_mode")
  # The old pipeline: the same reads, the verbatim per-pair loop and Steiger.
  columns <- c("chromosome", "base_pair_location", "effect_allele", "other_allele",
               "beta", "standard_error", "effect_allele_frequency", "p_value")
  union_keys <- unique(unlist(instruments, use.names = FALSE))
  all_data <- fastMR:::fastmr_io_map(
    c(unname(exposures), unname(outcomes)),
    c(unname(instruments), rep(list(union_keys), length(outcomes))),
    columns, 1L
  )
  exposure_data <- setNames(all_data[1:2], names(exposures))
  outcome_data <- setNames(all_data[3:5], names(outcomes))
  options <- list(samplesize_exposure = n_exposure, samplesize_outcome = n_outcome,
                  binary = NULL)
  for (threads in c(1L, 8L)) {
    for (minimum_snps in c(1L, 4L)) {
      expected <- capture_run({
        old <- old_compressed_pairwise_data(exposure_data, outcome_data,
                                            instruments, minimum_snps, FALSE)
        list(
          result = fast_mr(old$data, methods = methods, nboot = 50, seed = 3,
                           threads = threads),
          counts = old$counts,
          steiger = old_compressed_steiger(exposure_data, outcome_data, instruments,
                                           minimum_snps, options)
        )
      })
      got <- capture_run({
        result <- fast_mr_compressed(
          exposures, outcomes, instruments, methods = methods, nboot = 50,
          seed = 3, threads = threads, minimum_snps = minimum_snps,
          strict = FALSE, steiger = TRUE, samplesize_exposure = n_exposure,
          samplesize_outcome = n_outcome
        )
        meta <- attr(result, "compressed_input")
        expect_identical(meta$estimator_path, "pairwise")
        steiger <- attr(result, "steiger")
        attr(result, "compressed_input") <- NULL
        attr(result, "steiger") <- NULL
        list(result = result, counts = meta$counts, steiger = steiger)
      })
      expect_identical(got, expected, info = paste(threads, minimum_snps))
    }
  }
})
