test_that("batched PLINK2 clumping preserves exposure-specific greedy decisions", {
  skip_on_os("windows")
  plink2 <- tempfile("fastMR_plink2_stub_")
  writeLines(c(
    "#!/bin/sh",
    "out=''",
    "while [ \"$#\" -gt 0 ]; do",
    "  case \"$1\" in",
    "    --out) out=\"$2\"; shift 2;;",
    "    *) shift;;",
    "  esac",
    "done",
    "printf '1\\t100\\trs1\\t1\\t200\\trs2\\t0.90\\n1\\t200\\trs2\\t1\\t100\\trs1\\t0.90\\n' > \"${out}.vcor\"",
    "zstd -q -f \"${out}.vcor\" -o \"${out}.vcor.zst\"",
    "rm -f \"${out}.vcor\""
  ), plink2)
  Sys.chmod(plink2, "0755")
  dat <- data.frame(
    SNP = c("rs1", "rs2", "rs3", "rs2", "rs3"),
    id.exposure = c("E1", "E1", "E1", "E2", "E2"),
    pval.exposure = c(1e-8, 1e-7, 1e-6, 1e-8, 2e-8),
    chr_name = 1, chrom_start = c(100, 200, 500, 200, 500)
  )
  result <- fast_clump_data_batched(dat, bfile = "panel", plink2_bin = plink2)
  expect_equal(result$instruments$E1, c("rs1", "rs3"))
  expect_equal(result$instruments$E2, c("rs2", "rs3"))
  expect_true(isTRUE(result$diagnostics$exact))
  expect_gt(result$diagnostics$plink_calls, 0)
})

test_that("chromosome-partitioned clumping aggregates exact results", {
  skip_on_os("windows")
  plink2 <- tempfile("fastMR_plink2_chr_stub_")
  writeLines(c(
    "#!/bin/sh",
    "out=''",
    "while [ \"$#\" -gt 0 ]; do",
    "  case \"$1\" in",
    "    --out) out=\"$2\"; shift 2;;",
    "    *) shift;;",
    "  esac",
    "done",
    "printf '1\\t100\\trs1\\t1\\t200\\trs2\\t0.90\\n1\\t200\\trs2\\t1\\t100\\trs1\\t0.90\\n' > \"${out}.vcor\"",
    "zstd -q -f \"${out}.vcor\" -o \"${out}.vcor.zst\"",
    "rm -f \"${out}.vcor\""
  ), plink2)
  Sys.chmod(plink2, "0755")
  dat <- data.frame(
    SNP = c("rs1", "rs2", "rs3", "rs2", "rs3"),
    id.exposure = c("E1", "E1", "E1", "E2", "E2"),
    pval.exposure = c(1e-8, 1e-7, 1e-6, 1e-8, 2e-8),
    chr_name = c(1, 1, 1, 2, 2),
    chrom_start = c(100, 200, 500, 200, 500)
  )
  result <- fast_clump_data_batched_chromosomal(
    dat, bfile = "panel", plink2_bin = plink2
  )
  expect_equal(result$instruments$E1, c("rs1", "rs3"))
  expect_equal(result$instruments$E2, c("rs2", "rs3"))
  expect_identical(result$diagnostics$partition, "chromosome")
  expect_true(isTRUE(result$diagnostics$exact))
})

test_that("lead-row clumping reuses a row for exposures with the same lead", {
  skip_on_os("windows")
  plink2 <- tempfile("fastMR_plink2_lead_row_stub_")
  writeLines(c(
    "#!/bin/sh",
    "out=''",
    "while [ \"$#\" -gt 0 ]; do",
    "  case \"$1\" in",
    "    --out) out=\"$2\"; shift 2;;",
    "    *) shift;;",
    "  esac",
    "done",
    "printf '1\\t100\\trs1\\t1\\t200\\trs2\\t0.90\\n' > \"${out}.vcor\"",
    "zstd -q -f \"${out}.vcor\" -o \"${out}.vcor.zst\"",
    "rm -f \"${out}.vcor\""
  ), plink2)
  Sys.chmod(plink2, "0755")
  dat <- data.frame(
    SNP = c("rs1", "rs2", "rs1", "rs2"),
    id.exposure = c("E1", "E1", "E2", "E2"),
    pval.exposure = c(1e-8, 1e-7, 2e-8, 2e-7),
    chr_name = 1, chrom_start = c(100, 200, 100, 200)
  )
  result <- fast_clump_data_lead_rows(dat, bfile = "panel", plink2_bin = plink2)
  expect_equal(result$instruments$E1, "rs1")
  expect_equal(result$instruments$E2, "rs1")
  expect_equal(result$diagnostics$plink_calls, 1)
  expect_equal(result$diagnostics$unique_leads, 1)
  expect_identical(result$diagnostics$strategy, "lead_row")
  expect_true(isTRUE(result$diagnostics$exact))
})

test_that("lead-row follow-up queries keep the lead in --extract", {
  skip_on_os("windows")
  plink2 <- tempfile("fastMR_plink2_lead_extract_stub_")
  writeLines(c(
    "#!/bin/sh",
    "out=''",
    "lead_file=''",
    "extract_file=''",
    "while [ \"$#\" -gt 0 ]; do",
    "  case \"$1\" in",
    "    --out) out=\"$2\"; shift 2;;",
    "    --ld-snp-list) lead_file=\"$2\"; shift 2;;",
    "    --extract) extract_file=\"$2\"; shift 2;;",
    "    *) shift;;",
    "  esac",
    "done",
    "lead=$(head -n 1 \"$lead_file\")",
    "if [ \"$lead\" = rs1 ]; then",
    "  grep -Fxq rs1 \"$extract_file\" || exit 9",
    "fi",
    "printf '' > \"${out}.vcor\"",
    "if [ \"$lead\" = rs1 ] && grep -Fxq rs2 \"$extract_file\"; then",
    "  printf '1\\t100\\trs1\\t1\\t200\\trs2\\t0.90\\n' > \"${out}.vcor\"",
    "fi",
    "zstd -q -f \"${out}.vcor\" -o \"${out}.vcor.zst\""
  ), plink2)
  Sys.chmod(plink2, "0755")
  dat <- data.frame(
    SNP = c("rs1", "rs2", "rs0", "rs1", "rs3"),
    id.exposure = c("E1", "E1", "E2", "E2", "E2"),
    pval.exposure = c(1e-8, 1e-7, 1e-9, 2e-8, 3e-8),
    chr_name = 1, chrom_start = c(100, 200, 50, 100, 300)
  )
  result <- fast_clump_data_lead_rows(dat, bfile = "panel", plink2_bin = plink2)
  expect_equal(result$instruments$E1, "rs1")
  expect_equal(result$instruments$E2, c("rs0", "rs1", "rs3"))
  expect_true(isTRUE(result$diagnostics$exact))
})

test_that("batched clumping enforces bounded work", {
  skip_on_os("windows")
  plink2 <- tempfile("fastMR_plink2_limit_stub_")
  writeLines(c("#!/bin/sh", "exit 0"), plink2)
  Sys.chmod(plink2, "0755")
  dat <- data.frame(SNP = paste0("rs", 1:3), id.exposure = "E",
                    pval.exposure = c(1e-8, 2e-8, 3e-8))
  expect_error(
    fast_clump_data_batched(dat, bfile = "panel", plink2_bin = plink2,
                            max_pair_requests = 1),
    "pair-request limit"
  )
})

test_that("compressed candidate extraction records reconstructed p-value provenance", {
  skip_if_compressor_unavailable()
  skip_on_os("windows")
  input <- data.frame(
    chromosome = "1", base_pair_location = 1:3,
    reference_allele = "A", alternate_allele = c("C", "G", "T"),
    effect_allele = c("C", "G", "T"), other_allele = "A",
    beta = c(1, .1, .2), standard_error = .1,
    effect_allele_frequency = .2
  )
  store <- tempfile("fastMR-clump-store-")
  CompreSSoR::compress_sumstats(input, store, overwrite = TRUE)
  plink2 <- tempfile("fastMR_plink2_empty_stub-")
  writeLines(c(
    "#!/bin/sh", "out=''",
    "while [ \"$#\" -gt 0 ]; do case \"$1\" in --out) out=\"$2\"; shift 2;; *) shift;; esac; done",
    "printf '' > \"${out}.vcor\"",
    "zstd -q -f \"${out}.vcor\" -o \"${out}.vcor.zst\""
  ), plink2)
  Sys.chmod(plink2, "0755")
  result <- fast_clump_compressed(
    c(exposure = store), candidate_source = "full", pvalue_threshold = 1,
    bfile = "panel", plink2_bin = plink2
  )
  expect_false(result$diagnostics$compressed_input$pvalue_order_exact)
  expect_identical(result$diagnostics$compressed_input$pvalue_order, "reconstructed")
})
