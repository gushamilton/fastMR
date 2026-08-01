test_that("local harmonisation aligns strands and marks ambiguous palindromes", {
  exposure <- data.frame(
    SNP = paste0("rs", 1:4), id.exposure = "E", beta.exposure = c(.1, .2, .3, .4),
    se.exposure = .05, effect_allele.exposure = c("A", "A", "C", "A"),
    other_allele.exposure = c("G", "T", "G", "C"), eaf.exposure = c(.2, .2, .5, .3)
  )
  outcome <- data.frame(
    SNP = paste0("rs", 1:4), id.outcome = "O", beta.outcome = c(.3, .4, .5, .6),
    se.outcome = .06, effect_allele.outcome = c("T", "T", "G", "G"),
    other_allele.outcome = c("C", "A", "C", "T"), eaf.outcome = c(.2, .8, .5, .7)
  )
  result <- fast_harmonise_data(exposure, outcome, action = 2)
  result <- result[match(exposure$SNP, result$SNP), ]
  expect_equal(result$beta.outcome, c(.3, -.4, -.5, -.6), tolerance = 0)
  expect_equal(result$mr_keep, c(TRUE, TRUE, FALSE, TRUE))
  expect_true(all(result$palindromic[2:3]))
  expect_true(result$ambiguous[[3]])
})

test_that("harmonisation supports all native allele-information cases", {
  cases <- list(
    `2-2` = data.frame(effect_allele.exposure = "A", other_allele.exposure = "G",
                       effect_allele.outcome = "G", other_allele.outcome = "A"),
    `2-1` = data.frame(effect_allele.exposure = "A", other_allele.exposure = "G",
                       effect_allele.outcome = "G", other_allele.outcome = NA_character_),
    `1-2` = data.frame(effect_allele.exposure = "G", other_allele.exposure = NA_character_,
                       effect_allele.outcome = "A", other_allele.outcome = "G"),
    `1-1` = data.frame(effect_allele.exposure = "A", other_allele.exposure = NA_character_,
                       effect_allele.outcome = "A", other_allele.outcome = NA_character_)
  )
  for (case in cases) {
    exposure <- data.frame(SNP = "rs_case", id.exposure = "E", beta.exposure = .2,
      se.exposure = .05, eaf.exposure = .2, case[1, 1:2, drop = FALSE])
    outcome <- data.frame(SNP = "rs_case", id.outcome = "O", beta.outcome = .1,
      se.outcome = .04, eaf.outcome = .8, case[1, 3:4, drop = FALSE])
    result <- fast_harmonise_data(exposure, outcome, action = 1)
    expect_equal(nrow(result), 1L)
    expect_false(result$remove)
    expect_true(result$mr_keep)
  }

  exposure <- data.frame(SNP = "rs_indel", id.exposure = "E", beta.exposure = .2,
    se.exposure = .05, eaf.exposure = .2, effect_allele.exposure = "A",
    other_allele.exposure = "AT")
  outcome <- data.frame(SNP = "rs_indel", id.outcome = "O", beta.outcome = .1,
    se.outcome = .04, eaf.outcome = .8, effect_allele.outcome = "I",
    other_allele.outcome = "D")
  result <- fast_harmonise_data(exposure, outcome, action = 1)
  expect_false(result$remove)
  expect_equal(result$effect_allele.outcome, "A")
  expect_equal(result$other_allele.outcome, "AT")
  expect_equal(result$beta.outcome, -.1)

  outcome$se.outcome <- NA_real_
  result <- fast_harmonise_data(exposure, outcome, action = 1)
  expect_false(result$mr_keep)
  expect_true("samplesize.outcome" %in% names(result))
})

test_that("action vector applies one native policy per outcome", {
  exposure <- data.frame(SNP = paste0("rs_action", 1:3), id.exposure = "E",
    beta.exposure = .2, se.exposure = .05, eaf.exposure = c(.2, .5, .2),
    effect_allele.exposure = c("A", "A", "A"), other_allele.exposure = c("G", "T", "G"))
  outcome <- data.frame(SNP = paste0("rs_action", 1:3),
    id.outcome = c("O1", "O2", "O3"), beta.outcome = .1, se.outcome = .04,
    eaf.outcome = c(.2, .5, .2), effect_allele.outcome = c("A", "T", "A"),
    other_allele.outcome = c("G", "A", "C"))
  result <- fast_harmonise_data(exposure, outcome, action = c(1, 2, 3))
  result <- result[match(outcome$SNP, result$SNP), ]
  expect_equal(result$mr_keep, c(TRUE, FALSE, FALSE))
})

test_that("local LD-matrix clumping retains the best independent index SNPs", {
  dat <- data.frame(SNP = paste0("rs", 1:4), id.exposure = "E",
                    pval.exposure = c(1e-8, 1e-6, 1e-5, 1e-4),
                    chr_name = 1, chrom_start = c(100, 200, 500000, 600000))
  ld <- diag(4)
  ld[1, 2] <- ld[2, 1] <- .9
  rownames(ld) <- colnames(ld) <- dat$SNP
  result <- fast_clump_data(dat, clump_kb = 10, clump_r2 = .5, ld_matrix = ld)
  expect_equal(result$SNP, c("rs1", "rs3", "rs4"))
})
