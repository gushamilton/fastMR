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
