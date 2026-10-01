# Exact LD-graph clumping against a mocked LD oracle (no PLINK needed).
#
# LD model: r2 = 1 - |dbp| / 10000 for same-chromosome pairs, so r2 is exactly
# representable at the threshold (0.5 at 5 kb) and SNPs at 1 kb spacing sit
# exactly on the window edge for clump_kb = 5.

graph_ld_reference <- function() {
  chr <- rep(c("1", "2", "X"), each = 40)
  bp <- rep(seq(1000, 40000, by = 1000), 3)
  data.frame(SNP = sprintf("%s:%d:A:C", chr, bp), chr = chr, bp = bp,
             stringsAsFactors = FALSE)
}

graph_make_dat <- function(E, seed, ref = graph_ld_reference()) {
  set.seed(seed)
  pool <- c(1e-8, 1e-8, 1e-7, 1e-6, 1e-6, 1e-5)  # many exact ties in p
  rows <- lapply(seq_len(E), function(e) {
    ii <- sort(sample.int(nrow(ref), sample(5:60, 1L)))
    data.frame(SNP = ref$SNP[ii], id.exposure = sprintf("e%03d", e),
               pval.exposure = sample(pool, length(ii), TRUE),
               chr_name = ref$chr[ii], chrom_start = ref$bp[ii],
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

graph_pairs <- function(a, b, ref, kb, r2) {
  # all (a_i, b_j) pairs on the same chromosome within the window with r2 >= r2
  ia <- match(a, ref$SNP); ib <- match(b, ref$SNP)
  g <- expand.grid(i = seq_along(ia), j = seq_along(ib))
  g <- g[ia[g$i] != ib[g$j] & ref$chr[ia[g$i]] == ref$chr[ib[g$j]], , drop = FALSE]
  d <- abs(ref$bp[ia[g$i]] - ref$bp[ib[g$j]])
  ok <- d <= kb * 1000 & (1 - d / 10000) >= r2
  data.frame(lead = a[g$i[ok]], target = b[g$j[ok]], stringsAsFactors = FALSE)
}

with_ld_oracle <- function(ref, code) {
  calls <- new.env()
  calls$graph <- 0L; calls$frontier <- 0L
  testthat::local_mocked_bindings(
    fastmr_clump_run_graph = function(snps, reference_args, plink2_bin, clump_kb,
                                      clump_r2, threads, workdir, tag, ...) {
      calls$graph <- calls$graph + 1L
      u <- unique(snps)
      p <- graph_pairs(u, u, ref, clump_kb, clump_r2)
      p[as.integer(factor(p$lead, u)) < as.integer(factor(p$target, u)), , drop = FALSE]  # once per pair
    },
    fastmr_clump_run_frontier = function(leads, targets, reference_args, plink2_bin,
                                         clump_kb, clump_r2, threads, workdir, round) {
      calls$frontier <- calls$frontier + 1L
      graph_pairs(unique(leads), unique(targets), ref, clump_kb, clump_r2)
    },
    fastmr_clump_plink_version = function(plink2_bin) "PLINK v2.0.0-mock",
    .package = "fastMR", .env = parent.frame()
  )
  calls
}

graph_run_all <- function(dat, kb, r2, p1 = 1, ...) {
  common <- list(clump_kb = kb, clump_r2 = r2, clump_p1 = p1,
                 bfile = "mock", plink2_bin = "/bin/true")
  list(
    graph = do.call(fast_clump_data_graph, c(list(dat), common, list(...))),
    global = do.call(fast_clump_data_batched, c(list(dat), common)),
    lead_row = do.call(fast_clump_data_lead_rows, c(list(dat), common)),
    chromosome = do.call(fast_clump_data_batched_chromosomal, c(list(dat), common))
  )
}

test_that("graph clumping is identical to global and lead_row for E = 1, 10, 100", {
  skip_on_os("windows")
  ref <- graph_ld_reference()
  for (E in c(1L, 10L, 100L)) {
    dat <- graph_make_dat(E, seed = E, ref)
    with_ld_oracle(ref)
    for (cfg in list(c(5, 0.5), c(8, 0.2), c(2, 0.5))) {
      res <- graph_run_all(dat, kb = cfg[1], r2 = cfg[2], p1 = 1e-6)
      expect_identical(res$graph$instruments, res$global$instruments)
      expect_identical(res$graph$instruments, res$lead_row$instruments)
      expect_identical(res$graph$data, res$global$data)
      expect_identical(res$graph$data, res$lead_row$data)
      expect_identical(res$graph$instruments, res$chromosome$instruments)
      expect_true(isTRUE(res$graph$diagnostics$exact))
      expect_false(res$graph$diagnostics$fallback)
    }
  }
})

test_that("graph clumping makes one LD call per chromosome", {
  skip_on_os("windows")
  ref <- graph_ld_reference()
  dat <- graph_make_dat(50L, 3L, ref)
  calls <- with_ld_oracle(ref)
  res <- fast_clump_data_graph(dat, clump_kb = 5, clump_r2 = 0.5, bfile = "mock",
                               plink2_bin = "/bin/true")
  expect_identical(calls$graph, length(unique(dat$chr_name)))
  expect_identical(calls$frontier, 0L)
  expect_identical(res$diagnostics$plink_calls, calls$graph)
  expect_identical(res$diagnostics$partition, "graph")
  expect_identical(res$diagnostics$ld_provenance$plink2_version, "PLINK v2.0.0-mock")
  expect_true(any(grepl("--ld-window-r2 0.5", res$diagnostics$ld_provenance$flags, fixed = TRUE)))
  expect_true(any(grepl("--ld-window ", res$diagnostics$ld_provenance$flags, fixed = TRUE)))
})

test_that("p ties, window edges and r2 exactly at threshold match PLINK semantics", {
  skip_on_os("windows")
  ref <- graph_ld_reference()
  # chain 1:1000 .. 1:7000: the SNP 5 kb away has r2 == 0.5 (== threshold,
  # inclusive) exactly at the 5 kb window edge; 6 kb is outside.
  s <- ref$SNP[ref$chr == "1"][1:7]
  dat <- data.frame(SNP = s, id.exposure = "E", pval.exposure = 1e-8,
                    chr_name = "1", chrom_start = seq(1000, 7000, by = 1000))
  with_ld_oracle(ref)
  res <- graph_run_all(dat, kb = 5, r2 = 0.5)
  expect_identical(res$graph$instruments, res$global$instruments)
  expect_identical(res$graph$instruments, res$lead_row$instruments)
  # all tied: SNP order decides; first SNP kills 2..6 (<= 5 kb), not 7 (6 kb).
  expect_identical(res$graph$instruments$E, c(s[1], s[7]))
  # a second exposure whose order differs through p only
  dat2 <- rbind(dat, transform(dat, id.exposure = "F", pval.exposure = c(1e-7, 1e-9, 1e-8, 1e-8, 1e-6, 1e-8, 1e-8)))
  res2 <- graph_run_all(dat2, kb = 5, r2 = 0.5)
  expect_identical(res2$graph$data, res2$lead_row$data)
  expect_identical(res2$graph$instruments, res2$global$instruments)
})

test_that("multi-chromosome and ineligible rows are handled", {
  skip_on_os("windows")
  ref <- graph_ld_reference()
  dat <- graph_make_dat(20L, 11L, ref)
  dat$pval.exposure[c(1, 5, 9)] <- NA
  dat <- rbind(dat, dat[3:6, ])  # duplicated exposure/SNP rows are restored
  with_ld_oracle(ref)
  res <- graph_run_all(dat, kb = 5, r2 = 0.3, p1 = 1e-6)
  expect_identical(res$graph$data, res$global$data)
  expect_identical(res$graph$data, res$lead_row$data)
})

test_that("graph clumping falls back to lead_row above max_graph_pairs", {
  skip_on_os("windows")
  ref <- graph_ld_reference()
  dat <- graph_make_dat(10L, 4L, ref)
  calls <- with_ld_oracle(ref)
  res <- fast_clump_data_graph(dat, clump_kb = 5, clump_r2 = 0.5, bfile = "mock",
                               plink2_bin = "/bin/true", max_graph_pairs = 1)
  lead <- fast_clump_data_lead_rows(dat, clump_kb = 5, clump_r2 = 0.5, bfile = "mock",
                                    plink2_bin = "/bin/true")
  expect_identical(res$data, lead$data)
  expect_identical(res$instruments, lead$instruments)
  expect_true(res$diagnostics$fallback)
  expect_equal(length(res$diagnostics$fallbacks), 3L)
  expect_match(res$diagnostics$fallbacks[[1]]$reason, "max_graph_pairs")
  expect_identical(calls$graph, 0L)
  expect_gt(calls$frontier, 0L)
  # mixed: only the chromosome above the cap falls back
  mixed <- fast_clump_data_graph(dat, clump_kb = 5, clump_r2 = 0.5, bfile = "mock",
                                 plink2_bin = "/bin/true", max_graph_pairs = 400)
  expect_identical(mixed$data, lead$data)
  # missing positions -> fallback with a stated reason
  nopos <- dat[setdiff(names(dat), c("chr_name", "chrom_start"))]
  res3 <- fast_clump_data_graph(nopos, clump_kb = 5, clump_r2 = 0.5, bfile = "mock",
                                plink2_bin = "/bin/true")
  expect_identical(res3$diagnostics$fallbacks[[1]]$reason, "missing_positions")
})

test_that("C++ graph pass and vcor parser behave", {
  keep <- fastMR:::.fastmr_graph_clump(4L, c(0L, 1L), c(1L, 2L), c(0L, 1L, 2L, 3L, 2L, 0L), c(0L, 4L, 6L))
  # exposure 1: vertices 0,1,2,3 with edges 0-1, 1-2: keep 0, kill 1, keep 2, keep 3
  # exposure 2: vertices 2,0 (no edge between them)
  expect_identical(keep, c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE))
  ids <- fastMR:::.fastmr_vcor_ids(c("#CHROM_A\tPOS_A\tID_A\tCHROM_B\tPOS_B\tID_B\tUNPHASED_R2",
                                     "1\t100\trs1\t1\t200\trs2\t0.9", "1\t100\trs1\t1\t300\trs3"))
  expect_identical(ids$lead, c("rs1", "rs1"))
  expect_identical(ids$target, c("rs2", "rs3"))
})

test_that("graph partition runs through the PLINK2 argument surface", {
  skip_on_os("windows")
  plink2 <- tempfile("fastMR_plink2_graph_stub_")
  log <- tempfile("graph_args_")
  writeLines(c(
    "#!/bin/sh", "out=''",
    sprintf("echo \"$@\" >> %s", log),
    "while [ \"$#\" -gt 0 ]; do case \"$1\" in --out) out=\"$2\"; shift 2;; *) shift;; esac; done",
    "[ -n \"$out\" ] || exit 0",
    "printf '1\\t1000\\tA\\t1\\t2000\\tB\\t0.9\\n' > \"${out}.vcor\"",
    "zstd -q -f \"${out}.vcor\" -o \"${out}.vcor.zst\"; rm -f \"${out}.vcor\""
  ), plink2)
  Sys.chmod(plink2, "0755")
  skip_if(!nzchar(Sys.which("zstd")) && !nzchar(Sys.which("zstdcat")), "zstd unavailable")
  dat <- data.frame(SNP = c("A", "B", "C"), id.exposure = "E", pval.exposure = c(1e-8, 1e-7, 1e-6),
                    chr_name = "1", chrom_start = c(1000, 2000, 90000))
  res <- fast_clump_data_graph(dat, clump_kb = 500, clump_r2 = 0.01, bfile = "panel", plink2_bin = plink2)
  expect_identical(res$instruments$E, c("A", "C"))
  args <- readLines(log)
  expect_length(args, 2L)  # --version + one all-pairs call
  call <- args[grepl("--r2-unphased", args)]
  expect_match(call, "--r2-unphased zs --ld-window-kb 500 --ld-window 1000000000 --ld-window-r2 0.01", fixed = TRUE)
  expect_match(call, "--extract", fixed = TRUE)
})

test_that("compressed candidate reading is identical with and without read_candidates", {
  skip_if_compressor_unavailable()
  skip_on_os("windows")
  skip_if_not(exists("read_candidates", asNamespace("CompreSSoR")), "read_candidates unavailable")
  set.seed(2)
  V <- 3000L
  input <- data.frame(
    chromosome = "1", base_pair_location = seq.int(100001L, length.out = V),
    reference_allele = "A", alternate_allele = "C", effect_allele = "C", other_allele = "A",
    beta = rnorm(V, 0, 0.05), standard_error = runif(V, 0.01, 0.03),
    effect_allele_frequency = runif(V, 0.05, 0.95)
  )
  input$beta[c(5, 900, 2000)] <- c(0.9, -0.8, 0.7)
  paths <- c(a = tempfile("fm-a-"), b = tempfile("fm-b-"))
  CompreSSoR::compress_sumstats(input, paths[["a"]], overwrite = TRUE)
  input$beta <- rev(input$beta)
  CompreSSoR::compress_sumstats(input, paths[["b"]], overwrite = TRUE)
  legacy <- function(path, thr) {
    n <- CompreSSoR::open_compressor(path)$manifest$n_rows
    x <- CompreSSoR::read_sumstats(path, variants = seq.int(0L, n - 1L),
      columns = c("chromosome", "base_pair_location", "effect_allele", "other_allele", "p_value"))
    key <- CompreSSoR::compressor_variant_key(x$chromosome, x$base_pair_location, x$other_allele, x$effect_allele)
    keep <- is.finite(x$p_value) & x$p_value <= thr
    data.frame(SNP = key[keep], pval.exposure = x$p_value[keep])
  }
  thrs <- c(1e-300, 1e-3, 0.05, 1)
  cand <- function(thr) fastMR:::fastmr_compressed_candidate_data(paths, names(paths), thr, "full", "reconstructed", 1L)
  new <- lapply(thrs, cand)
  testthat::local_mocked_bindings(fastmr_have_compressor_fn = function(name) FALSE, .package = "fastMR")
  old <- lapply(thrs, cand)
  for (k in seq_along(thrs)) {
    expect_identical(new[[k]]$data, old[[k]]$data)
    for (nm in names(paths)) {
      ref <- legacy(paths[[nm]], thrs[[k]])
      got <- new[[k]]$data[new[[k]]$data$id.exposure == nm, c("SNP", "pval.exposure")]
      rownames(got) <- NULL
      expect_identical(got, ref)
    }
  }
})

test_that("require_exact uses the CompreSSoR rank domain when available", {
  skip_if_compressor_unavailable()
  skip_on_os("windows")
  skip_if_not(exists("read_pvalue_order", asNamespace("CompreSSoR")), "read_pvalue_order unavailable")
  set.seed(3)
  V <- 400L
  input <- data.frame(
    chromosome = "1", base_pair_location = seq.int(100001L, length.out = V),
    reference_allele = "A", alternate_allele = "C", effect_allele = "C", other_allele = "A",
    beta = c(rnorm(V - 5, 0, 0.02), rep(0.4, 5)), standard_error = 0.05,
    effect_allele_frequency = 0.3
  )
  plain <- tempfile("fm-plain-"); exact <- tempfile("fm-exact-")
  CompreSSoR::compress_sumstats(input, plain, overwrite = TRUE)
  expect_error(
    fast_clump_compressed(c(x = plain), pvalue_order = "require_exact", pvalue_threshold = 1e-3,
                          bfile = "panel", plink2_bin = "/bin/true", candidate_source = "full"),
    "exact p-value ordering domain")
  ok <- tryCatch({
    CompreSSoR::compress_sumstats(input, exact, overwrite = TRUE, pvalue_order = TRUE,
                                  pvalue_order_threshold = 0.01)
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(ok, "cannot write an exact-order store with this CompreSSoR")
  cand <- fastMR:::fastmr_compressed_candidate_data(c(x = exact), "x", 1e-3, "full", "require_exact", 1L)
  expect_true(cand$exact)
  expect_true("pvalue_rank" %in% names(cand$data))
  expect_false(anyNA(cand$data$pvalue_rank))
})

test_that(".fastmr_vcor_ids is GC-safe under gctorture", {
  skip_on_cran()
  n <- 25
  lines <- c("#CHROM_A\tPOS_A\tID_A\tCHROM_B\tPOS_B\tID_B\tUNPHASED_R2",
             sprintf("22\t%d\tidA_%d_x\t22\t%d\tidB_%d_y\t0.5", 1:n, 1:n, 1:n, 1:n))
  gctorture(TRUE)
  r <- tryCatch(fastMR:::.fastmr_vcor_ids(lines), finally = gctorture(FALSE))
  expect_identical(r$lead, sprintf("idA_%d_x", 1:n))
  expect_identical(r$target, sprintf("idB_%d_y", 1:n))
})

test_that("graph run treats 'No variants remaining' as an empty graph", {
  skip_on_os("windows")
  plink2 <- tempfile("fastMR_plink2_novar_")
  writeLines(c("#!/bin/sh",
               "case \"$*\" in *--version*) echo 'PLINK v2.0.0-stub'; exit 0;; esac",
               "echo 'Error: No variants remaining after --extract.'", "exit 3"), plink2)
  Sys.chmod(plink2, "0755")
  wd <- tempfile("wd"); dir.create(wd)
  ld <- fastMR:::fastmr_clump_run_graph(c("A", "B"), c("--bfile", "panel"), plink2, 500, 0.01, 1L, wd, "1")
  expect_identical(nrow(ld), 0L)
  expect_named(ld, c("lead", "target"))
  # other failures still abort
  writeLines(c("#!/bin/sh", "echo 'Error: boom'", "exit 3"), plink2)
  expect_error(fastMR:::fastmr_clump_run_graph(c("A", "B"), c("--bfile", "panel"), plink2, 500, 0.01, 1L, wd, "1"),
               "failed")
})

test_that("graph clumping skips PLINK when no candidate pairs exist", {
  skip_on_os("windows")
  plink2 <- tempfile("fastMR_plink2_never_")
  log <- tempfile("never_args_")
  writeLines(c("#!/bin/sh", sprintf("echo \"$@\" >> %s", log),
               "case \"$*\" in *--version*) echo 'PLINK v2.0.0-stub'; exit 0;; esac",
               "echo 'Error: No variants remaining after --extract.'", "exit 3"), plink2)
  Sys.chmod(plink2, "0755")
  dat <- data.frame(SNP = c("A", "B"), id.exposure = "E", pval.exposure = c(1e-8, 1e-7),
                    chr_name = "1", chrom_start = c(1000, 900000))
  res <- fast_clump_data_graph(dat, clump_kb = 1, clump_r2 = 0.01, bfile = "panel", plink2_bin = plink2)
  expect_identical(res$instruments$E, c("A", "B"))
  expect_false(any(grepl("--r2-unphased", readLines(log))))
})
