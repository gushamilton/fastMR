# Reference MAF floor (fastMR#19): candidates that are monomorphic in, or
# absent from, the LD reference have undefined LD with every lead.  They were
# never pruned and survived as spurious "independent" instruments.

plink2_for_tests <- function() {
  path <- Sys.getenv("FASTMR_PLINK2", unset = "")
  if (!nzchar(path)) path <- Sys.which("plink2")
  if (!nzchar(path) || !file.exists(path)) return("")
  unname(path)
}

skip_if_no_plink2 <- function() {
  skip_on_os("windows")
  skip_if(!nzchar(plink2_for_tests()), "PLINK2 not available (set FASTMR_PLINK2)")
  skip_if(!nzchar(Sys.which("zstd")) && !nzchar(Sys.which("zstdcat")), "zstd not available")
}

# 40-sample synthetic reference.  rs2 duplicates rs1 (r2 = 1), rs3 is
# monomorphic and sits between them, rs4 is exactly uncorrelated with rs1,
# and rs6 is on chromosome 2.  rs5 is deliberately absent.
synthetic_reference <- function() {
  plink2 <- plink2_for_tests()
  dir <- tempfile("fastMR_ref_maf_")
  dir.create(dir)
  n <- 40L
  gt <- function(dosage) c("0/0", "0/1", "1/1")[dosage + 1L]
  g1 <- rep(c(0L, 2L), n / 2L)
  g4 <- rep(c(0L, 0L, 2L, 2L), n / 4L)
  rows <- list(
    c("1", "1000", "rs1", gt(g1)),
    c("1", "2000", "rs2", gt(g1)),
    c("1", "3000", "rs3", gt(rep(0L, n))),
    c("1", "4000", "rs4", gt(g4)),
    c("2", "1000", "rs6", gt(rep(c(0L, 1L), n / 2L)))
  )
  vcf <- file.path(dir, "ref.vcf")
  samples <- sprintf("S%02d", seq_len(n))
  lines <- c(
    "##fileformat=VCFv4.2", "##contig=<ID=1>", "##contig=<ID=2>",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples),
          collapse = "\t"),
    vapply(rows, function(r) paste(c(r[1:3], "A", "G", ".", ".", ".", "GT", r[-(1:3)]),
                                   collapse = "\t"), character(1))
  )
  writeLines(lines, vcf)
  out <- system2(plink2, c("--vcf", vcf, "--make-bed", "--out", file.path(dir, "ref")),
                 stdout = TRUE, stderr = TRUE)
  if (!file.exists(file.path(dir, "ref.bed"))) stop("could not build synthetic reference: ", paste(out, collapse = "\n"))
  file.path(dir, "ref")
}

reference_candidates <- function() {
  data.frame(
    SNP = c("rs1", "rs2", "rs3", "rs4", "rs5", "rs3", "rs1", "rs6"),
    id.exposure = c("E1", "E1", "E1", "E1", "E1", "E2", "E2", "E2"),
    pval.exposure = c(1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-12, 1e-9, 1e-8),
    chr_name = c(1, 1, 1, 1, 1, 1, 1, 2),
    chrom_start = c(1000, 2000, 3000, 4000, 5000, 3000, 1000, 1000)
  )
}

test_that("legacy clumping (min_ref_maf = 0) keeps monomorphic and absent candidates", {
  skip_if_no_plink2()
  ref <- synthetic_reference()
  dat <- reference_candidates()
  expect_message(
    old <- fast_clump_data_batched(dat, bfile = ref, plink2_bin = plink2_for_tests(),
                                   clump_r2 = 0.1, min_ref_maf = 0),
    "floor disabled"
  )
  # This is the result of the pre-#19 implementation: rs3 (monomorphic) and
  # rs5 (absent) cannot be pruned, and the monomorphic E2 lead rs3 cannot
  # prune rs1.
  expect_equal(old$instruments$E1, c("rs1", "rs3", "rs4", "rs5"))
  expect_equal(old$instruments$E2, c("rs3", "rs1", "rs6"))
  summary <- attr(old, "reference_maf")
  expect_identical(summary, old$diagnostics$reference_maf)
  expect_equal(summary$min_ref_maf, 0)
  expect_equal(summary$absent_rows, 1)
  expect_equal(summary$below_floor_rows, 0)
  expect_equal(summary$kept_rows, nrow(dat))
  expect_false(summary$absent_dropped)
})

test_that("the default reference MAF floor drops monomorphic and absent candidates", {
  skip_if_no_plink2()
  ref <- synthetic_reference()
  dat <- reference_candidates()
  expect_message(
    new <- fast_clump_data_batched(dat, bfile = ref, plink2_bin = plink2_for_tests(),
                                   clump_r2 = 0.1),
    "kept 5 of 8 candidate rows \\(1 absent from reference, 2 below floor dropped\\)"
  )
  expect_equal(new$instruments$E1, c("rs1", "rs4"))
  expect_equal(new$instruments$E2, c("rs1", "rs6"))
  expect_false(any(new$data$SNP %in% c("rs3", "rs5")))
  summary <- attr(new, "reference_maf")
  expect_equal(summary$min_ref_maf, 0.01)
  expect_equal(summary$candidate_rows, 8)
  expect_equal(summary$absent_rows, 1)
  expect_equal(summary$below_floor_rows, 2)
  expect_equal(summary$kept_rows, 5)
  expect_equal(summary$absent_variants, 1)
  expect_equal(summary$below_floor_variants, 1)
  expect_true(summary$absent_dropped)

  # A floor above every reference MAF drops everything except MAF 0.5 variants.
  strict <- suppressMessages(fast_clump_data_batched(
    dat, bfile = ref, plink2_bin = plink2_for_tests(), clump_r2 = 0.1, min_ref_maf = 0.3
  ))
  expect_equal(strict$instruments$E1, c("rs1", "rs4"))
  expect_equal(strict$instruments$E2, "rs1")
  expect_equal(attr(strict, "reference_maf")$below_floor_rows, 3)
})

test_that("lead-row and chromosome-partitioned clumpers apply the same floor", {
  skip_if_no_plink2()
  ref <- synthetic_reference()
  dat <- reference_candidates()
  plink2 <- plink2_for_tests()
  for (fun in list(fast_clump_data_lead_rows, fast_clump_data_batched_chromosomal)) {
    old <- suppressMessages(fun(dat, bfile = ref, plink2_bin = plink2, clump_r2 = 0.1, min_ref_maf = 0))
    expect_equal(old$instruments$E1, c("rs1", "rs3", "rs4", "rs5"))
    expect_equal(old$instruments$E2, c("rs3", "rs1", "rs6"))
    new <- suppressMessages(fun(dat, bfile = ref, plink2_bin = plink2, clump_r2 = 0.1))
    expect_equal(new$instruments$E1, c("rs1", "rs4"))
    expect_equal(new$instruments$E2, c("rs1", "rs6"))
    summary <- attr(new, "reference_maf")
    expect_identical(summary, new$diagnostics$reference_maf)
    expect_equal(c(summary$absent_rows, summary$below_floor_rows, summary$kept_rows), c(1, 2, 5))
  }
})

test_that("min_ref_maf is validated", {
  dat <- reference_candidates()
  expect_error(fast_clump_data_batched(dat, bfile = "panel", plink2_bin = "plink2", min_ref_maf = 0.5),
               "min_ref_maf")
  expect_error(fast_clump_data(dat, ld_matrix = diag(1), min_ref_maf = -1), "min_ref_maf")
})

test_that("PLINK 1.9 clumping drops candidates that are rare or absent in the reference", {
  skip_on_os("windows")
  plink <- tempfile("fastMR_plink_freq_stub_")
  writeLines(c(
    "#!/bin/sh",
    "case \" $* \" in *\" --freq \"*)",
    "  fx=''; fo=''; prev=''",
    "  for a in \"$@\"; do case \"$prev\" in --extract) fx=\"$a\";; --out) fo=\"$a\";; esac; prev=\"$a\"; done",
    "  printf ' CHR SNP A1 A2 MAF NCHROBS\\n   1 rs1 G A 0.4 100\\n   1 rs2 G A 0 100\\n   1 rs3 G A 0.005 100\\n' > \"${fo}.frq\"",
    "  exit 0;;",
    "esac",
    "input=''; out=''",
    "while [ \"$#\" -gt 0 ]; do case \"$1\" in --clump) input=\"$2\"; shift 2;; --out) out=\"$2\"; shift 2;; *) shift;; esac; done",
    "printf 'CHR SNP BP P NSIG S05 S01 S001 S0001\\n' > \"${out}.clumped\"",
    "awk 'NR > 1 { print \"1\", $1, 100, $2, 1, 1, 1, 1, 1 }' \"$input\" >> \"${out}.clumped\""
  ), plink)
  Sys.chmod(plink, "0755")
  dat <- data.frame(SNP = c("rs1", "rs2", "rs3", "rs4"), id.exposure = "E",
                    pval.exposure = c(1e-8, 1e-7, 1e-6, 1e-5))
  expect_message(new <- fast_clump_data(dat, bfile = "panel", plink_bin = plink),
                 "kept 1 of 4")
  expect_equal(new$SNP, "rs1")
  expect_equal(attr(new, "reference_maf")$absent_rows, 1)
  expect_equal(attr(new, "reference_maf")$below_floor_rows, 2)
  old <- suppressMessages(fast_clump_data(dat, bfile = "panel", plink_bin = plink, min_ref_maf = 0))
  expect_equal(old$SNP, dat$SNP)
  expect_equal(attr(old, "reference_maf")$absent_rows, 1)
  lenient <- suppressMessages(fast_clump_data(dat, bfile = "panel", plink_bin = plink, min_ref_maf = 0.001))
  expect_equal(attr(lenient, "reference_maf")$below_floor_rows, 1)
})

test_that("LD-matrix clumping never treats undefined LD as independence", {
  dat <- data.frame(SNP = paste0("rs", 1:3), id.exposure = "E",
                    pval.exposure = c(1e-8, 1e-7, 1e-6),
                    chr_name = 1, chrom_start = c(100, 200, 300))
  ld <- diag(3)
  ld[1, 2] <- ld[2, 1] <- NaN  # e.g. rs2 monomorphic in the panel
  rownames(ld) <- colnames(ld) <- dat$SNP
  expect_message(new <- fast_clump_data(dat, clump_kb = 10, clump_r2 = 0.5, ld_matrix = ld),
                 "1 dropped for undefined LD")
  expect_equal(new$SNP, c("rs1", "rs3"))
  expect_equal(attr(new, "reference_maf")$undefined_ld_rows, 1)
  old <- fast_clump_data(dat, clump_kb = 10, clump_r2 = 0.5, ld_matrix = ld, min_ref_maf = 0)
  expect_equal(old$SNP, c("rs1", "rs2", "rs3"))
})
