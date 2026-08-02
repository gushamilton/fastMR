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

test_that("action vectors follow native SNP-sorted outcome order", {
  exposure <- data.frame(SNP = c("rs1", "rs2"), id.exposure = "E",
    beta.exposure = .2, se.exposure = .05, eaf.exposure = .2,
    effect_allele.exposure = c("A", "A"), other_allele.exposure = c("T", "T"))
  outcome <- data.frame(SNP = c("rs2", "rs1"), id.outcome = c("O2", "O1"),
    beta.outcome = .1, se.outcome = .04, eaf.outcome = .8,
    effect_allele.outcome = c("T", "T"), other_allele.outcome = c("A", "A"))
  result <- fast_harmonise_data(exposure, outcome, action = c(1, 3))
  result <- result[match(exposure$SNP, result$SNP), ]
  expect_equal(result$mr_keep, c(TRUE, FALSE))
})

test_that("harmonisation handles native empty alleles, equal alleles, and factors", {
  exposure <- data.frame(SNP = "rs_factor", id.exposure = "E",
    beta.exposure = factor("0.2"), se.exposure = factor("0.05"), eaf.exposure = .2,
    effect_allele.exposure = "A", other_allele.exposure = "G")
  outcome <- data.frame(SNP = "rs_factor", id.outcome = "O",
    beta.outcome = factor("0.1"), se.outcome = factor("0.04"), eaf.outcome = .2,
    effect_allele.outcome = "A", other_allele.outcome = "")
  result <- fast_harmonise_data(exposure, outcome, action = 1)
  expect_equal(result$beta.exposure, .2)
  expect_equal(result$beta.outcome, .1)
  expect_false(result$remove)
  expect_true(result$mr_keep)

  exposure$effect_allele.exposure <- "A"
  exposure$other_allele.exposure <- "A"
  outcome$effect_allele.outcome <- "A"
  result <- fast_harmonise_data(exposure, outcome, action = 1)
  expect_false(result$remove)
})

test_that("action length is validated against all merged outcomes", {
  exposure <- data.frame(SNP = c("rs_valid", "rs_missing"), id.exposure = "E",
    beta.exposure = c(.2, .3), se.exposure = c(.05, .06), eaf.exposure = c(.2, .3),
    effect_allele.exposure = c("A", NA), other_allele.exposure = c("G", NA))
  outcome <- data.frame(SNP = c("rs_valid", "rs_missing"),
    id.outcome = c("O1", "O2"), beta.outcome = c(.1, .2), se.outcome = c(.04, .05),
    eaf.outcome = c(.2, .3), effect_allele.outcome = c("A", NA),
    other_allele.outcome = c("G", NA))
  result <- fast_harmonise_data(exposure, outcome, action = c(1, 2))
  expect_equal(result$SNP, "rs_valid")
  expect_error(fast_harmonise_data(exposure, outcome, action = c(1, 2, 3)),
               "one value per unique id.outcome")
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

test_that("local clumping converts factor p-values and handles short groups", {
  dat <- data.frame(SNP = paste0("rs", 1:3), id.exposure = "E",
                    pval.exposure = factor(c("1e-8", "1e-2", "1e-6")))
  ld <- diag(3)
  rownames(ld) <- colnames(ld) <- dat$SNP
  result <- fast_clump_data(dat, clump_p1 = 1e-5, ld_matrix = ld)
  expect_equal(result$SNP, c("rs1", "rs3"))
  expect_equal(fast_clump_data(dat[1, , drop = FALSE], ld_matrix = ld)$SNP, "rs1")
  expect_error(fast_clump_data(dat[c(1, 1), ], ld_matrix = ld),
               "unique within each id.exposure")
  dat$SNP[1] <- NA_character_
  expect_error(fast_clump_data(dat, ld_matrix = ld), "non-missing")
})

test_that("PLINK clumping passes an explicit clump field and preserves groups", {
  skip_on_os("windows")
  plink <- tempfile("fastMR_plink_stub_")
  writeLines(c(
    "#!/bin/sh",
    "input=''",
    "out=''",
    "bfile=''",
    "field=''",
    "while [ \"$#\" -gt 0 ]; do",
    "  case \"$1\" in",
    "    --bfile) bfile=\"$2\"; shift 2;;",
    "    --clump) input=\"$2\"; shift 2;;",
    "    --clump-field) field=\"$2\"; shift 2;;",
    "    --out) out=\"$2\"; shift 2;;",
    "    *) shift;;",
    "  esac",
    "done",
    "if [ \"$bfile\" != \"/tmp/panel with spaces\" ] || [ \"$field\" != \"P\" ]; then exit 42; fi",
    "printf 'CHR SNP BP P NSIG S05 S01 S001 S0001\n' > \"${out}.clumped\"",
    "awk 'NR == 2 { print \"1\", $1, 100, $2, 1, 1, 1, 1, 1 }' \"$input\" >> \"${out}.clumped\""
  ), plink)
  Sys.chmod(plink, "0755")
  dat <- data.frame(SNP = c("rs1", "rs2", "rs3"),
    id.exposure = c("E1", "E1", "E2"), pval.exposure = c(1e-8, 1e-6, 1e-7))
  result <- fast_clump_data(dat, bfile = "/tmp/panel with spaces", plink_bin = plink)
  expect_equal(result$SNP, c("rs1", "rs3"))
})

test_that("PLINK failures include a useful status", {
  skip_on_os("windows")
  plink <- tempfile("fastMR_plink_failure_")
  writeLines(c("#!/bin/sh", "echo synthetic PLINK failure >&2", "exit 17"), plink)
  Sys.chmod(plink, "0755")
  dat <- data.frame(SNP = c("rs1", "rs2"), id.exposure = "E",
                    pval.exposure = c(1e-8, 1e-6))
  expect_error(fast_clump_data(dat, bfile = "panel", plink_bin = plink),
               "exit status 17.*synthetic PLINK failure")
})
