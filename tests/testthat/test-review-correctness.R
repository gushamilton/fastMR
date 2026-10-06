# Regression tests for the correctness review (findings 4, 6, 7, 8, 10, 11, 15).

# --- #4: seed + group - 1 must stay a valid R integer --------------------------

test_that("a seed whose per-pair stream would overflow is rejected before any work", {
  set.seed(1)
  d <- data.frame(SNP = paste0("rs", rep(1:4, 3)), id.exposure = rep(paste0("e", 1:3), each = 4),
                  id.outcome = "o", beta.exposure = rnorm(12, 0.1, 0.02), se.exposure = 0.01,
                  beta.outcome = rnorm(12, 0.02, 0.01), se.outcome = 0.01)
  before <- .Random.seed
  expect_error(fast_mr(d, methods = "weighted_median", nboot = 10, seed = .Machine$integer.max - 1L),
               "exceeds the largest R integer")
  expect_identical(.Random.seed, before)
  # the largest seed that fits: groups use max - 2, max - 1, max
  ok <- fast_mr(d, methods = "weighted_median", nboot = 10, seed = .Machine$integer.max - 2L)
  expect_identical(nrow(ok), 3L)
  expect_identical(ok$se[3], fast_mr(d[9:12, ], methods = "weighted_median", nboot = 10,
                                     seed = .Machine$integer.max)$se)
  # seeds outside the integer range fail validation in every entry point
  expect_error(fast_mr(d, methods = "weighted_median", nboot = 10, seed = 2^31), "seed must lie in")
  expect_error(fast_mr(d, methods = "weighted_median", nboot = 10, seed = -2^31), "seed must lie in")
  # RNG-free runs do not use the seed streams
  expect_no_error(fast_mr(d, methods = "ivw", nboot = 10, seed = .Machine$integer.max))
})

# --- #15: grid loop index (64-bit); mapping unchanged ---------------------------

test_that("threaded grid maps pair indices to exposure/outcome as the serial path does", {
  set.seed(15)
  k <- 6
  eb <- matrix(rnorm(5 * k, 0.1, 0.05), 5); es <- matrix(runif(5 * k, 0.01, 0.02), 5)
  ob <- matrix(rnorm(3 * k, 0.05, 0.05), 3); os <- matrix(runif(3 * k, 0.01, 0.02), 3)
  m <- c("ivw", "egger", "weighted_median")
  one <- fast_mr_grid(eb, ob, es, os, methods = m, nboot = 50, seed = 3, threads = 1)
  two <- fast_mr_grid(eb, ob, es, os, methods = m, nboot = 50, seed = 3, threads = 2)
  expect_identical(one, two)
  direct <- fast_mr(data.frame(SNP = paste0("s", 1:k), beta.exposure = eb[4, ], se.exposure = es[4, ],
                               beta.outcome = ob[2, ], se.outcome = os[2, ]), methods = "ivw", nboot = 0)
  row <- one[one$id.exposure == "4" & one$id.outcome == "2" & one$method == direct$method, ]
  expect_equal(row$b, direct$b, tolerance = 1e-12)
})

# --- #7 / #8: clumping order, duplicates, absent candidates, build mismatch ------

review_clump_mocks <- function(present, edges, env = parent.frame()) {
  # edges: list of SNP pairs in LD; every other pair is not in LD.
  in_ld <- function(a, b) any(vapply(edges, function(e) setequal(e, c(a, b)), logical(1)))
  testthat::local_mocked_bindings(
    fastmr_clump_reference_ids = function(snps, ...) intersect(snps, present),
    fastmr_clump_run_graph = function(snps, ...) {
      g <- if (length(snps) > 1L) t(utils::combn(length(snps), 2L)) else matrix(integer(), 0, 2)
      ok <- if (nrow(g)) mapply(function(i, j) in_ld(snps[i], snps[j]), g[, 1], g[, 2]) else logical()
      list(lead = g[ok, 1], target = g[ok, 2])
    },
    fastmr_clump_run_frontier = function(leads, targets, ...) {
      g <- expand.grid(lead = unique(leads), target = unique(targets), stringsAsFactors = FALSE)
      g <- g[g$lead != g$target & g$lead %in% present & g$target %in% present, , drop = FALSE]
      if (nrow(g)) g <- g[mapply(in_ld, g$lead, g$target), , drop = FALSE]
      rownames(g) <- NULL
      g
    },
    fastmr_clump_plink_version = function(...) "PLINK v2.0.0-mock",
    .package = "fastMR", .env = env)
}

test_that("p = 0 ties are led by the largest |z|, not the SNP-ID string", {
  skip_on_os("windows")
  z <- c(40, 90, 60, 39)
  dat <- data.frame(SNP = c("1:100500:A:C", "1:99000:A:C", "1:100200:A:C", "1:100900:A:C"),
                    id.exposure = "prot", beta.exposure = z * 0.01, se.exposure = 0.01,
                    pval.exposure = 2 * stats::pnorm(-z), chr_name = "1",
                    chrom_start = c(100500, 99000, 100200, 100900), stringsAsFactors = FALSE)
  expect_true(all(dat$pval.exposure == 0))
  snps <- dat$SNP
  review_clump_mocks(present = snps, edges = utils::combn(snps, 2L, simplify = FALSE))
  g <- fast_clump_data_graph(dat, clump_kb = 10000, clump_r2 = 0.001, bfile = "mock", plink2_bin = "/bin/true")
  expect_identical(g$instruments$prot, "1:99000:A:C")   # z = 90
  lr <- fast_clump_data_lead_rows(dat, clump_kb = 10000, clump_r2 = 0.001, bfile = "mock",
                                  plink2_bin = "/bin/true")
  expect_identical(lr$data, g$data)
  bt <- fast_clump_data_batched(dat, clump_kb = 10000, clump_r2 = 0.001, bfile = "mock",
                                plink2_bin = "/bin/true")
  expect_identical(bt$data, g$data)
  # per-exposure hands PLINK the greedy rank: the strongest variant is rank 1
  ranked <- NULL
  testthat::local_mocked_bindings(
    fastmr_clump_make_subset = function(snps, ...) list(
      reference_args = c("--pfile", "sub"),
      pvar = data.frame(id = dat$SNP, chr = "1", pos = dat$chrom_start, stringsAsFactors = FALSE)),
    fastmr_clump_run_clump = function(snps, ...) { ranked <<- snps; list(ids = snps[1L], pvar = NULL) },
    fastmr_clump_run_lead_graph = function(leads, snps, ...) {
      list(lead = rep(match(leads, snps), length(snps) - 1L),
           target = setdiff(seq_along(snps), match(leads, snps)), invalid_r2 = 0, invalid_example = "")
    },
    .package = "fastMR")
  pe <- fast_clump_data_per_exposure(dat, clump_kb = 10000, clump_r2 = 0.001, bfile = "mock",
                                     plink2_bin = "/bin/true")
  expect_identical(ranked, dat$SNP[order(-z)])
  expect_identical(pe$data, g$data)
  # the rest of the tie (equal p, equal |z|) still falls back to SNP order
  expect_identical(fastMR:::fastmr_clump_order(c(0, 0, 0), c("b", "a", "c"), tiebreak = c(-5, -5, -9)),
                   c(3L, 2L, 1L))
})

test_that("without p ties the clumping order ignores |z|", {
  set.seed(7)
  p <- runif(50, 0, 1e-6)
  snp <- sprintf("rs%02d", sample(50))
  tb <- -abs(rnorm(50, 0, 30))
  expect_identical(fastMR:::fastmr_clump_order(p, snp, tiebreak = tb), order(p, snp, method = "radix"))
  rank <- sample(50)
  expect_identical(fastMR:::fastmr_clump_order(p, snp, rank = rank, tiebreak = tb),
                   order(rank, p, snp, method = "radix"))
  dat <- data.frame(SNP = snp, id.exposure = "e", pval.exposure = p, beta.exposure = tb, se.exposure = 1)
  expect_null(fastMR:::fastmr_clump_tiebreak(dat[c("SNP", "pval.exposure")], "pval.exposure"))
  expect_equal(fastMR:::fastmr_clump_tiebreak(dat, "pval.exposure"), -abs(tb))
})

test_that("a repeated (exposure, SNP) is ordered by its best p; every row is returned", {
  skip_on_os("windows")
  dat <- data.frame(SNP = c("rsA", "rsB", "rsA"), id.exposure = "e", pval.exposure = c(1e-8, 1e-10, 1e-30),
                    chr_name = "1", chrom_start = c(100, 200, 100), outcome = c("o1", "o1", "o2"))
  review_clump_mocks(present = c("rsA", "rsB"), edges = list(c("rsA", "rsB")))
  g <- fast_clump_data_graph(dat, bfile = "x", plink2_bin = "/bin/true")
  expect_identical(g$instruments$e, c("rsA", "rsA"))  # one entry per returned row, as before
  expect_identical(g$data, dat[c(1L, 3L), ])
  lr <- fast_clump_data_lead_rows(dat, bfile = "x", plink2_bin = "/bin/true")
  expect_identical(lr$data, g$data)
  ld <- matrix(c(1, 0.9, 0.9, 1), 2, dimnames = list(c("rsA", "rsB"), c("rsA", "rsB")))
  local <- fast_clump_data(dat, ld_matrix = ld)
  expect_identical(local, dat[c(1L, 3L), ])
  keep <- fastMR:::fastmr_clump_dedup(transform(dat, pval.exposure = c(NA, 1e-10, 1e-9)), "pval.exposure")
  expect_identical(keep, c(FALSE, TRUE, TRUE))
})

test_that("candidates absent from the LD reference are counted, warned about, optionally dropped", {
  skip_on_os("windows")
  dat <- data.frame(SNP = c("rs1", "rs2", "rs3", "rs4"), id.exposure = c("e", "e", "e", "f"),
                    pval.exposure = c(1e-10, 1e-9, 1e-8, 1e-8), chr_name = "1",
                    chrom_start = c(100, 200, 300, 400))
  review_clump_mocks(present = c("rs1", "rs2"), edges = list(c("rs1", "rs2")))
  expect_warning(keep <- fast_clump_data_graph(dat, bfile = "x", plink2_bin = "/bin/true"),
                 "2 of 4 eligible candidate rows \\(2 variants\\) are absent .* kept unclumped")
  expect_identical(keep$instruments, list(e = c("rs1", "rs3"), f = "rs4"))
  expect_identical(keep$diagnostics$absent_from_reference, 2L)
  expect_warning(drop <- fast_clump_data_graph(dat, bfile = "x", plink2_bin = "/bin/true", absent = "drop"),
                 "were dropped")
  expect_identical(drop$instruments, list(e = "rs1"))
  expect_error(fast_clump_data_graph(dat[3:4, ], bfile = "x", plink2_bin = "/bin/true"),
               "none of the 2 eligible candidate rows \\(2 variants\\) is in the LD reference")
  expect_error(fast_clump_data_graph(dat, bfile = "x", plink2_bin = "/bin/true", absent = "nope"))
  # per-exposure (and auto) agree, and pass `absent` on when delegating
  testthat::local_mocked_bindings(
    fastmr_clump_make_subset = function(snps, ...) list(
      reference_args = c("--pfile", "sub"),
      pvar = data.frame(id = c("rs1", "rs2"), chr = "1", pos = c(100, 200), stringsAsFactors = FALSE)),
    fastmr_clump_run_clump = function(snps, ...) list(ids = snps[1L], pvar = NULL),
    fastmr_clump_run_lead_graph = function(leads, snps, ...) {
      list(lead = 1L, target = 2L, invalid_r2 = 0, invalid_example = "")
    },
    .package = "fastMR")
  for (mode in c("keep", "drop")) {
    expect_warning(pe <- fast_clump_data_per_exposure(dat, bfile = "x", plink2_bin = "/bin/true",
                                                      absent = mode), "absent from the LD reference")
    expect_null(pe$diagnostics$delegated)
    expect_identical(pe$data, if (mode == "keep") keep$data else drop$data)
    expect_warning(au <- fast_clump_data_auto(dat, bfile = "x", plink2_bin = "/bin/true", absent = mode),
                   "absent from the LD reference")
    expect_identical(au$data, pe$data)
  }
  expect_warning(dg <- fast_clump_data_per_exposure(dat, clump_r2 = 0, bfile = "x", plink2_bin = "/bin/true",
                                                    absent = "drop"), "were dropped")
  expect_match(dg$diagnostics$delegated, "clump_r2")
  expect_identical(dg$data, drop$data)
})

test_that("a genome-build (position) mismatch with the reference is warned about", {
  skip_on_os("windows")
  dat <- data.frame(SNP = paste0("rs", 1:6), id.exposure = "e", pval.exposure = 10^-(20:15),
                    chr_name = "1", chrom_start = (1:6) * 1e6)
  review_clump_mocks(present = dat$SNP, edges = list())
  testthat::local_mocked_bindings(
    fastmr_clump_auto_plan = function(...) list(strategy = "per_exposure", reason = "forced"),
    fastmr_clump_make_subset = function(snps, ...) list(
      reference_args = "x",
      pvar = data.frame(id = dat$SNP, chr = "1", pos = dat$chrom_start + 150000, stringsAsFactors = FALSE)),
    .package = "fastMR")
  expect_warning(res <- fast_clump_data_auto(dat, bfile = "ref", plink2_bin = "/bin/true"),
                 "6 of 6 candidate variants found in the LD reference have a different position")
  expect_identical(res$diagnostics$delegated, "reference_position_mismatch")
  expect_identical(nrow(res$data), 6L)
})

# --- #10: identical exposure alleles --------------------------------------------

test_that("an exposure with identical alleles is removed (two outcome alleles) and never flipped", {
  e <- data.frame(SNP = c("rs1", "rs2", "rs3", "rs4"), beta.exposure = 0.1, se.exposure = 0.01,
                  effect_allele.exposure = c("A", "A", "a", "C"), other_allele.exposure = c("A", "A", "A", "T"),
                  eaf.exposure = 0.3, id.exposure = "e", exposure = "e", pval.exposure = 1e-10)
  o <- data.frame(SNP = c("rs1", "rs2", "rs3", "rs4"), beta.outcome = 0.2, se.outcome = 0.02,
                  effect_allele.outcome = c("G", "A", "A", "T"), other_allele.outcome = c("A", NA, "G", "C"),
                  eaf.outcome = 0.3, id.outcome = "o", outcome = "o", pval.outcome = 1e-3)
  for (action in 1:3) {
    h <- fast_harmonise_data(e, o, action = action)
    h <- h[order(h$SNP), ]
    # rs2 has only the outcome effect allele: TwoSampleMR keeps it (subject to
    # the ambiguity rules), so only its outcome effect stays unflipped.
    expect_identical(h$remove, c(TRUE, FALSE, TRUE, FALSE))
    expect_false(any(h$mr_keep[c(1, 3)]))
    expect_true(h$mr_keep[4])
    expect_identical(h$beta.outcome, c(0.2, 0.2, 0.2, -0.2))   # rs4 is a genuine switch
    expect_identical(h$eaf.outcome[1:3], c(0.3, 0.3, 0.3))
  }
  skip_if_not_installed("TwoSampleMR")
  for (action in 1:3) {
    t <- suppressMessages(TwoSampleMR::harmonise_data(e, o, action = action))
    t <- t[order(t$SNP), ]
    f <- fast_harmonise_data(e, o, action = action)
    f <- f[order(f$SNP), ]
    expect_identical(f$remove, t$remove)
    expect_identical(f$mr_keep, t$mr_keep)
    expect_equal(f$beta.outcome, t$beta.outcome)
    expect_equal(f$eaf.outcome, t$eaf.outcome)
  }
})

# --- #11: compressed pipeline robustness ----------------------------------------

test_that("the sparse-IVW memory estimate includes the outcome x union matrices", {
  b <- fastMR:::fastmr_sparse_ivw_dispatch_bytes(2940, 2940, 50000)
  expect_equal(b, 40 * 2940^2 + 24 * 2940 * 50000 + 8 * 2940^2)
  expect_gt(b, 3.5e9)
  expect_true(fastMR:::fastmr_sparse_ivw_fits(2940, 2940, 50000))
  withr::local_options(fastMR.sparse_ivw_max_memory_mb = 1024)
  expect_false(fastMR:::fastmr_sparse_ivw_fits(2940, 2940, 50000))
  expect_true(fastMR:::fastmr_sparse_ivw_fits(10, 10, 100))
})

review_compressed_stores <- function() {
  vapply(c(1, 1.3, 0.7), function(m) {
    p <- tempfile("fm-review-")
    CompreSSoR::compress_sumstats(compressor_canonical_fixture(m), p, overwrite = TRUE)
    p
  }, character(1))
}

review_keys <- function() {
  id <- compressor_canonical_fixture()
  CompreSSoR::compressor_variant_key(id$chromosome, id$base_pair_location, id$other_allele, id$effect_allele)
}

test_that("an empty instrument set is dropped with strict = FALSE and an error with strict = TRUE", {
  skip_if_compressor_unavailable()
  skip_on_os("windows")
  stores <- review_compressed_stores()
  keys <- review_keys()
  ex <- stats::setNames(stores[1:2], c("ea", "eb"))
  out <- c(oa = stores[[3]])
  expect_warning(r <- fast_mr_compressed(ex, out, list(ea = keys[1:5], eb = character()),
                                         methods = c("ivw", "egger"), strict = FALSE),
                 "dropping exposure\\(s\\) with an empty instrument set: eb")
  ref <- fast_mr_compressed(ex[1], out, list(ea = keys[1:5]), methods = c("ivw", "egger"))
  attr(r, "compressed_input") <- NULL
  attr(ref, "compressed_input") <- NULL
  expect_identical(r, ref)
  expect_error(fast_mr_compressed(ex, out, list(ea = keys[1:5], eb = character()), methods = "ivw"),
               "instrument set is empty for exposure\\(s\\): eb")
  expect_error(fast_mr_compressed(ex, out, list(ea = character(), eb = character()), methods = "ivw",
                                  strict = FALSE),
               "no exposure has any instruments")
})

test_that("estimator = 'pairwise' is overridden by the shared-grid path; memory cap selects pairwise", {
  skip_if_compressor_unavailable()
  skip_on_os("windows")
  stores <- review_compressed_stores()
  keys <- review_keys()
  ex <- stats::setNames(stores[1:2], c("ea", "eb"))
  out <- c(oa = stores[[3]])
  shared <- fast_mr_compressed(ex, out, keys[1:10], methods = "ivw", estimator = "pairwise")
  expect_identical(attr(shared, "compressed_input")$estimator_path, "shared_instrument_grid")
  sets <- list(ea = keys[1:10], eb = keys[5:20])
  sparse <- fast_mr_compressed(ex, out, sets, methods = "ivw")
  expect_identical(attr(sparse, "compressed_input")$estimator_path, "sparse_ivw")
  withr::local_options(fastMR.sparse_ivw_max_memory_mb = 1e-6)
  capped <- fast_mr_compressed(ex, out, sets, methods = "ivw")
  expect_identical(attr(capped, "compressed_input")$estimator_path, "pairwise")
  expect_equal(capped$b, sparse$b, tolerance = 1e-12)
})

# #6 (row-id/key guard, fallback short read) is tested in
# test-candidate-count-guard.R, next to its store fixture.
