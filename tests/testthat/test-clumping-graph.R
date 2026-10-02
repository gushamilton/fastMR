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
      p <- p[as.integer(factor(p$lead, u)) < as.integer(factor(p$target, u)), , drop = FALSE]  # once per pair
      list(lead = match(p$lead, snps), target = match(p$target, snps))
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

test_that("graph clumping makes one LD call for all chromosomes", {
  skip_on_os("windows")
  ref <- graph_ld_reference()
  dat <- graph_make_dat(50L, 3L, ref)
  calls <- with_ld_oracle(ref)
  res <- fast_clump_data_graph(dat, clump_kb = 5, clump_r2 = 0.5, bfile = "mock",
                               plink2_bin = "/bin/true")
  expect_gt(length(unique(dat$chr_name)), 1L)
  expect_identical(calls$graph, 1L)
  expect_identical(calls$frontier, 0L)
  expect_identical(res$diagnostics$plink_calls, calls$graph)
  expect_identical(res$diagnostics$partition, "graph")
  expect_identical(res$diagnostics$ld_provenance$plink2_version, "PLINK v2.0.0-mock")
  expect_true(any(grepl("--ld-window-r2 0.5", res$diagnostics$ld_provenance$flags, fixed = TRUE)))
  expect_true(any(grepl("--ld-window ", res$diagnostics$ld_provenance$flags, fixed = TRUE)))
})

test_that("max_graph_pairs splits chromosomes into several LD calls with identical results", {
  skip_on_os("windows")
  ref <- graph_ld_reference()
  dat <- graph_make_dat(30L, 5L, ref)
  est <- vapply(split(dat$chrom_start, dat$chr_name), function(bp)
    fastMR:::fastmr_graph_pair_estimate(unique(bp), 5000), numeric(1))
  calls <- with_ld_oracle(ref)
  one <- fast_clump_data_graph(dat, clump_kb = 5, clump_r2 = 0.5, bfile = "mock", plink2_bin = "/bin/true")
  calls$graph <- 0L
  split_res <- fast_clump_data_graph(dat, clump_kb = 5, clump_r2 = 0.5, bfile = "mock",
                                     plink2_bin = "/bin/true", max_graph_pairs = max(est))
  lead <- fast_clump_data_lead_rows(dat, clump_kb = 5, clump_r2 = 0.5, bfile = "mock", plink2_bin = "/bin/true")
  expect_gt(calls$graph, 1L)
  expect_lte(calls$graph, length(est))
  expect_false(split_res$diagnostics$fallback)
  expect_identical(split_res$instruments, one$instruments)
  expect_identical(split_res$data, lead$data)
  expect_identical(split_res$diagnostics$graph_edges, one$diagnostics$graph_edges)
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
    "printf '#ID_A\\tID_B\\tUNPHASED_R2\\nA\\tB\\t0.9\\n' > \"${out}.vcor\""
  ), plink2)
  Sys.chmod(plink2, "0755")
  dat <- data.frame(SNP = c("A", "B", "C"), id.exposure = "E", pval.exposure = c(1e-8, 1e-7, 1e-6),
                    chr_name = "1", chrom_start = c(1000, 2000, 90000))
  res <- fast_clump_data_graph(dat, clump_kb = 500, clump_r2 = 0.01, bfile = "panel", plink2_bin = plink2)
  expect_identical(res$instruments$E, c("A", "C"))
  args <- readLines(log)
  expect_length(args, 2L)  # --version + one all-pairs call
  call <- args[grepl("--r2-phased", args)]
  expect_match(call, "--r2-phased cols=id --ld-window-kb 500 --ld-window 1000000000 --ld-window-r2 0.01", fixed = TRUE)
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
  # CompreSSoR >= the exact-cis-order change writes the domain by default; opt out here.
  CompreSSoR::compress_sumstats(input, plain, overwrite = TRUE, pvalue_order = FALSE)
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

test_that("graph run treats 'No variants remaining' as an empty graph", {
  skip_on_os("windows")
  plink2 <- tempfile("fastMR_plink2_novar_")
  writeLines(c("#!/bin/sh",
               "case \"$*\" in *--version*) echo 'PLINK v2.0.0-stub'; exit 0;; esac",
               "echo 'Error: No variants remaining after --extract.'", "exit 3"), plink2)
  Sys.chmod(plink2, "0755")
  wd <- tempfile("wd"); dir.create(wd)
  ld <- fastMR:::fastmr_clump_run_graph(c("A", "B"), c("--bfile", "panel"), plink2, 500, 0.01, 1L, wd, "1")
  expect_identical(ld, list(lead = integer(), target = integer()))
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
  expect_false(any(grepl("--r2-phased", readLines(log))))
})

write_vcor <- function(lines, final_newline = TRUE, eol = "\n") {
  f <- tempfile(fileext = ".vcor")
  writeBin(charToRaw(paste0(paste(lines, collapse = eol), if (final_newline) eol else "")), f)
  f
}

test_that("native vcor reader maps IDs, honours header columns, final line and CRLF", {
  ids <- c("rs1", "rs2", "rs3", "rs4")
  f <- write_vcor(c("#CHROM_A\tPOS_A\tID_A\tCHROM_B\tPOS_B\tID_B\tUNPHASED_R2",
                    "1\t100\trs1\t1\t200\trs2\t0.9", "", "1\t100\trs3\t1\t300\trs4\t0.5"),
                  final_newline = FALSE)
  r <- fastMR:::.fastmr_vcor_read(f, ids)
  expect_identical(r$lead, c(1L, 3L)); expect_identical(r$target, c(2L, 4L))
  f2 <- write_vcor(c("#ID_A\tID_B\tUNPHASED_R2", "rs1\trs2\t0.9", "rs3\trs4"), eol = "\r\n")
  r2 <- fastMR:::.fastmr_vcor_read(f2, ids)
  expect_identical(r2$lead, c(1L, 3L)); expect_identical(r2$target, c(2L, 4L))
  # ID_B as the last column with CRLF
  f3 <- write_vcor(c("#ID_A\tID_B", "rs1\trs2"), eol = "\r\n")
  expect_identical(fastMR:::.fastmr_vcor_read(f3, ids)$target, 2L)
})

test_that("native vcor reader errors on unknown IDs, short lines, bad header; empty is empty", {
  ids <- c("a", "b")
  expect_error(fastMR:::.fastmr_vcor_read(write_vcor(c("#ID_A\tID_B\tR", "a\tzzz\t1")), ids), "not among")
  expect_error(fastMR:::.fastmr_vcor_read(write_vcor(c("#ID_A\tID_B\tR", "a\tb\t1", "a")), ids), "malformed")
  expect_error(fastMR:::.fastmr_vcor_read(write_vcor(c("A\tB\t1", "a\tb\t1")), ids), "unrecognised")
  expect_error(fastMR:::.fastmr_vcor_read(write_vcor(c("#ID_B\tID_A\tR", "a\tb\t1")), ids), "unrecognised")
  expect_error(fastMR:::.fastmr_vcor_read(tempfile(), ids), "cannot open")
  empty <- list(lead = integer(), target = integer(), invalid_r2 = 0, invalid_example = "")
  r <- fastMR:::.fastmr_vcor_read(write_vcor("#ID_A\tID_B\tUNPHASED_R2"), ids)
  expect_identical(r, empty)
  f0 <- tempfile(); file.create(f0)
  expect_identical(fastMR:::.fastmr_vcor_read(f0, ids), empty)
})

test_that("native vcor reader on a >1e6-line file equals a pure-R reference", {
  skip_on_cran()
  set.seed(1)
  nv <- 5000L; n <- 1200003L
  ids <- sprintf("rs%d", seq_len(nv))
  a <- sample.int(nv, n, TRUE); b <- sample.int(nv, n, TRUE)
  f <- write_vcor(c("#CHROM_A\tPOS_A\tID_A\tCHROM_B\tPOS_B\tID_B\tUNPHASED_R2",
                    sprintf("1\t%d\t%s\t1\t%d\t%s\t0.5", a, ids[a], b, ids[b])))
  ref <- strsplit(readLines(f)[-1L], "\t", fixed = TRUE)
  ref_a <- match(vapply(ref[seq_len(1000)], `[`, "", 3L), ids)
  r <- fastMR:::.fastmr_vcor_read(f, ids)
  expect_identical(r$lead, a); expect_identical(r$target, b)
  expect_identical(r$lead[1:1000], ref_a)
  tab <- utils::read.delim(f, comment.char = "", header = TRUE, colClasses = "character")
  expect_identical(r$lead, match(tab$ID_A, ids)); expect_identical(r$target, match(tab$ID_B, ids))
})

test_that("native vcor reader is GC-safe under gctorture", {
  skip_on_cran()
  n <- 25
  ids <- c(sprintf("idA_%d_x", 1:n), sprintf("idB_%d_y", 1:n))
  f <- write_vcor(c("#CHROM_A\tPOS_A\tID_A\tCHROM_B\tPOS_B\tID_B\tUNPHASED_R2",
                    sprintf("22\t%d\tidA_%d_x\t22\t%d\tidB_%d_y\t0.5", 1:n, 1:n, 1:n, 1:n)))
  gctorture(TRUE)
  r <- tryCatch(fastMR:::.fastmr_vcor_read(f, ids), finally = gctorture(FALSE))
  expect_identical(r$lead, 1:n); expect_identical(r$target, as.integer(n + 1:n))
})

test_that("graph run parses plain PLINK output into vertex ids", {
  skip_on_os("windows")
  n <- 25
  tab <- c("#ID_A\tID_B\tUNPHASED_R2", sprintf("v%d\tw%d\t0.5", 1:n, 1:n))
  plink2 <- tempfile("fastMR_plink2_cols_")
  src <- tempfile(); writeLines(tab, src)
  writeLines(c("#!/bin/sh", "out=''",
               "while [ \"$#\" -gt 0 ]; do case \"$1\" in --out) out=\"$2\"; shift 2;; *) shift;; esac; done",
               sprintf("cp %s \"${out}.vcor\"", src)), plink2)
  Sys.chmod(plink2, "0755")
  wd <- tempfile("wd"); dir.create(wd)
  snps <- c(sprintf("v%d", 1:n), sprintf("w%d", 1:n))
  ld <- fastMR:::fastmr_clump_run_graph(snps, c("--bfile", "p"), plink2, 500, 0.01, 1L, wd, "1")
  expect_identical(ld$lead, 1:n)
  expect_identical(ld$target, as.integer(n + 1:n))
  expect_length(list.files(wd, pattern = "vcor"), 0L)  # temp table removed
  expect_error(fastMR:::fastmr_clump_run_graph(c("v1", "w1"), c("--bfile", "p"), plink2, 500, 0.01, 1L, wd, "1"),
               "not among")
})

test_that("graph LD parse refuses an unrecognised .vcor header", {
  skip_on_os("windows")
  plink2 <- tempfile("fastMR_plink2_badheader_")
  writeLines(c(
    "#!/bin/sh", "out=''",
    "while [ \"$#\" -gt 0 ]; do case \"$1\" in --out) out=\"$2\"; shift 2;; *) shift;; esac; done",
    "[ -n \"$out\" ] || exit 0",
    "printf 'A\\tB\\t0.9\\n' > \"${out}.vcor\""
  ), plink2)
  Sys.chmod(plink2, "0755")
  dat <- data.frame(SNP = c("A", "B"), id.exposure = "E", pval.exposure = c(1e-8, 1e-7),
                    chr_name = "1", chrom_start = c(1000, 2000))
  expect_error(fast_clump_data_graph(dat, clump_kb = 500, clump_r2 = 0.01, bfile = "panel", plink2_bin = plink2))
})

test_that("one-pass compressed candidates equal the three-pass fallback (3-exposure batch, exact order)", {
  skip_if_compressor_unavailable()
  skip_on_os("windows")
  skip_if_not(fastMR:::fastmr_have_one_pass_candidates(), "one-pass candidates unavailable")
  V <- 600L
  mk <- function(seed) {
    set.seed(seed)
    z <- rnorm(V, 0, 1.3); z[sample.int(V, 40L)] <- sample(c(-1, 1), 40L, TRUE) * runif(40L, 3, 12)
    data.frame(
      chromosome = rep(c("1", "2"), each = V / 2L),
      base_pair_location = rep(seq.int(100001L, length.out = V / 2L), 2L),
      reference_allele = "A", alternate_allele = "C", effect_allele = "C", other_allele = "A",
      beta = z * 0.05, standard_error = 0.05, effect_allele_frequency = 0.3)
  }
  paths <- vapply(1:3, function(i) {
    p <- tempfile("fm-onepass-")
    CompreSSoR::compress_sumstats(mk(i), p, overwrite = TRUE, pvalue_order = TRUE,
                                  pvalue_order_threshold = 0.01)
    p
  }, character(1))
  labels <- c("a", "b", "c")
  for (po in c("reconstructed", "require_exact")) {
    new <- fastMR:::fastmr_compressed_candidate_data(paths, labels, 1e-3, "full", po, 1L)
    testthat::local_mocked_bindings(fastmr_have_one_pass_candidates = function() FALSE,
                                    .package = "fastMR")
    old <- fastMR:::fastmr_compressed_candidate_data(paths, labels, 1e-3, "full", po, 1L)
    expect_gt(nrow(new$data), 0L)
    expect_identical(new, old)
  }
})

# Regression: PLINK --clump (1.9 and 2) uses haplotype-frequency ("phased", EM
# when phase is unknown) r2.  The graph and lead-row LD queries used
# --r2-unphased (squared dosage correlation), which disagrees for rare
# variants: on 1000G EUR, rs55942980 -> rs10157022 has unphased r2 0.00192 but
# phased r2 0.000548, so at clump_r2 = 0.001 fastMR dropped a SNP PLINK keeps.
# Fixture (20 unphased samples): A-B unphased 0.396 / phased 0.103 (PLINK
# keeps B); C-D unphased 0.116 / phased 0.306 (PLINK drops D).
test_that("graph and lead-row clumping use the same r2 as PLINK --clump", {
  skip_on_os("windows")
  plink2 <- Sys.getenv("FASTMR_PLINK2", Sys.which("plink2"))
  skip_if(!nzchar(plink2), "plink2 not available")
  g <- list(A = c(0,0,0,1,0,1,0,0,0,0,0,1,0,1,1,1,2,0,0,0),
            B = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0),
            C = c(1,1,1,1,1,0,1,2,1,2,1,1,0,1,0,0,1,1,1,1),
            D = c(1,0,1,0,0,0,0,1,1,0,1,1,0,1,0,0,0,0,1,0))
  chr <- c(A = "1", B = "1", C = "2", D = "2")
  bp <- c(A = 1000, B = 2000, C = 1000, D = 2000)
  dir <- tempfile("fastMR_r2fx_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  vcf <- file.path(dir, "fx.vcf")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=1,length=100000>", "##contig=<ID=2,length=100000>",
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
               paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                       sprintf("s%02d", 1:20)), collapse = "\t"),
               vapply(names(g), function(s) paste(c(chr[[s]], bp[[s]], s, "A", "G", ".", ".", ".", "GT",
                                                    c("0/0", "0/1", "1/1")[g[[s]] + 1]), collapse = "\t"), "")),
             vcf)
  pref <- file.path(dir, "fx")
  st <- system2(plink2, c("--vcf", shQuote(vcf), "--make-pgen", "--out", shQuote(pref)), stdout = FALSE, stderr = FALSE)
  skip_if(st != 0L, "plink2 could not build the fixture")
  dat <- data.frame(SNP = names(g), id.exposure = "E", pval.exposure = c(1e-10, 1e-9, 1e-10, 1e-9),
                    chr_name = unname(chr), chrom_start = unname(bp), stringsAsFactors = FALSE)
  # PLINK2's own --clump defines the expected set.
  pfile <- file.path(dir, "p.txt")
  write.table(data.frame(SNP = dat$SNP, P = dat$pval.exposure), pfile, row.names = FALSE, quote = FALSE)
  system2(plink2, c("--pfile", shQuote(pref), "--clump", shQuote(pfile), "--clump-p1", "1", "--clump-p2", "1",
                    "--clump-r2", "0.2", "--clump-kb", "250", "--out", shQuote(file.path(dir, "c"))),
          stdout = FALSE, stderr = FALSE)
  clumps <- utils::read.table(file.path(dir, "c.clumps"), header = FALSE, comment.char = "#")
  expect_setequal(clumps[[3]], c("A", "B", "C"))
  res <- fast_clump_data_graph(dat, clump_kb = 250, clump_r2 = 0.2, pfile = pref, plink2_bin = plink2)
  expect_setequal(res$instruments$E, c("A", "B", "C"))
  skip_if(!nzchar(Sys.which("zstdcat")) && !nzchar(Sys.which("zstd")), "zstd not available")
  lr <- fast_clump_data_lead_rows(dat, clump_kb = 250, clump_r2 = 0.2, pfile = pref, plink2_bin = plink2)
  expect_setequal(lr$instruments$E, c("A", "B", "C"))
})

test_that("impossible r2 > 1 from PLINK2 (2.00a6.8 --r2-phased bug) is detected and warned", {
  ids <- c("rs74865827", "rs139603618", "rs3")
  f <- write_vcor(c("#ID_A\tID_B\tPHASED_R2", "rs74865827\trs139603618\t96.1295", "rs74865827\trs3\t0.5",
                    "rs139603618\trs3\t1"))
  r <- fastMR:::.fastmr_vcor_read(f, ids)
  expect_identical(r$lead, c(1L, 1L, 2L))
  expect_identical(r$invalid_r2, 1)
  expect_identical(r$invalid_example, "rs74865827 rs139603618 r2=96.1295")
  expect_warning(fastMR:::fastmr_clump_warn_invalid_r2(r$invalid_r2, r$invalid_example), "r2 > 1")
  expect_silent(fastMR:::fastmr_clump_warn_invalid_r2(0, ""))
  ok <- fastMR:::.fastmr_vcor_read(write_vcor(c("#ID_A\tID_B\tUNPHASED_R2", "rs74865827\trs3\t1.0000001")), ids)
  expect_identical(ok$invalid_r2, 0)
})
