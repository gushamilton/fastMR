# Per-exposure --clump partition and auto dispatch, against a mocked LD oracle
# like test-clumping-graph.R's (r2 = 1 - |dbp| / 10000 on one chromosome,
# --ld-window-r2 threshold applied with PLINK2's (1 - 2^-44) tolerance)
# plus a PLINK2 --clump emulator that applies PLINK2's own argument parsing
# (bp_radius = (int)(kb * 1000 * (1 + eps) - 1), edge iff r2 > r2 * (1 + eps))
# and greedy order by P.  A real-PLINK2 fixture checks the boundaries at the end.

pe_eps <- 2^-44

pe_ref <- function() {
  chr <- rep(c("1", "2", "X"), each = 40)
  bp <- rep(seq(1000, 40000, by = 1000), 3)
  data.frame(SNP = sprintf("%s:%d:A:C", chr, bp), chr = chr, bp = bp, stringsAsFactors = FALSE)
}

pe_make_dat <- function(E, seed, ref = pe_ref(), absent = 0L) {
  set.seed(seed)
  pool <- c(1e-8, 1e-8, 1e-7, 1e-6, 1e-6, 1e-5)
  rows <- lapply(seq_len(E), function(e) {
    ii <- sort(sample.int(nrow(ref), sample(5:60, 1L)))
    d <- data.frame(SNP = ref$SNP[ii], id.exposure = sprintf("e%03d", e),
                    pval.exposure = sample(pool, length(ii), TRUE),
                    chr_name = ref$chr[ii], chrom_start = ref$bp[ii], stringsAsFactors = FALSE)
    if (absent > 0L) {
      # absent from the reference, but placed inside LD windows of present SNPs
      a <- sample(1:39, absent)
      d <- rbind(d, data.frame(SNP = sprintf("1:%d:G:T", a * 1000 + 500), id.exposure = sprintf("e%03d", e),
                               pval.exposure = sample(pool, absent, TRUE), chr_name = "1",
                               chrom_start = a * 1000 + 500, stringsAsFactors = FALSE))
    }
    d
  })
  do.call(rbind, rows)
}

pe_pairs <- function(a, b, ref, kb, r2) {
  ia <- match(a, ref$SNP); ib <- match(b, ref$SNP)
  ok_a <- !is.na(ia); ok_b <- !is.na(ib)
  a <- a[ok_a]; ia <- ia[ok_a]; b <- b[ok_b]; ib <- ib[ok_b]
  if (!length(a) || !length(b)) return(data.frame(lead = character(), target = character()))
  g <- expand.grid(i = seq_along(ia), j = seq_along(ib))
  g <- g[ia[g$i] != ib[g$j] & ref$chr[ia[g$i]] == ref$chr[ib[g$j]], , drop = FALSE]
  d <- abs(ref$bp[ia[g$i]] - ref$bp[ib[g$j]])
  ok <- d <= kb * 1000 & (1 - d / 10000) >= r2 * (1 - pe_eps)   # PLINK2 --ld-window-r2 parsing
  data.frame(lead = a[g$i[ok]], target = b[g$j[ok]], stringsAsFactors = FALSE)
}

# Emulates `plink2 --clump` on the reference `ref` (exact r2 model).
pe_emulate_clump <- function(snps, ref, kb_arg, r2_arg) {
  radius <- trunc(as.numeric(kb_arg) * 1000 * (1 + pe_eps) - 1)
  cut <- as.numeric(r2_arg) * (1 + pe_eps)
  s <- snps[snps %in% ref$SNP]          # PLINK2 ignores IDs missing from the dataset
  i <- match(s, ref$SNP)
  dead <- rep(FALSE, length(s)); lead <- rep(FALSE, length(s))
  for (k in seq_along(s)) {             # `snps` is in P (greedy) order
    if (dead[k]) next
    lead[k] <- TRUE
    d <- abs(ref$bp[i] - ref$bp[i[k]])
    dead <- dead | (ref$chr[i] == ref$chr[i[k]] & d <= radius & (1 - d / 10000) > cut)
  }
  s[lead]
}

with_pe_oracle <- function(ref, code) {
  calls <- new.env()
  calls$graph <- 0L; calls$frontier <- 0L; calls$subset <- 0L; calls$clump <- 0L; calls$cert <- 0L
  calls$extract <- logical()
  testthat::local_mocked_bindings(
    fastmr_clump_run_graph = function(snps, reference_args, plink2_bin, clump_kb,
                                      clump_r2, threads, workdir, tag, ...) {
      calls$graph <- calls$graph + 1L
      u <- unique(snps)
      p <- pe_pairs(u, u, ref, clump_kb, clump_r2)
      p <- p[as.integer(factor(p$lead, u)) < as.integer(factor(p$target, u)), , drop = FALSE]
      list(lead = match(p$lead, snps), target = match(p$target, snps))
    },
    fastmr_clump_run_frontier = function(leads, targets, reference_args, plink2_bin,
                                         clump_kb, clump_r2, threads, workdir, round) {
      calls$frontier <- calls$frontier + 1L
      pe_pairs(unique(leads), unique(targets), ref, clump_kb, clump_r2)
    },
    fastmr_clump_make_subset = function(snps, reference_args, plink2_bin, threads, stem) {
      calls$subset <- calls$subset + 1L
      hit <- ref[ref$SNP %in% snps, , drop = FALSE]
      if (!nrow(hit)) return(NULL)
      list(reference_args = c("--pfile", stem),
           pvar = data.frame(id = hit$SNP, chr = hit$chr, pos = hit$bp, stringsAsFactors = FALSE))
    },
    fastmr_clump_run_clump = function(snps, reference_args, plink2_bin, kb_arg, r2_arg, stem,
                                      extract = FALSE) {
      calls$clump <- calls$clump + 1L
      calls$extract <- c(calls$extract, extract)
      hit <- ref[ref$SNP %in% snps, , drop = FALSE]
      list(ids = pe_emulate_clump(snps, ref, kb_arg, r2_arg),
           pvar = if (extract) data.frame(id = hit$SNP, chr = hit$chr, pos = hit$bp,
                                          stringsAsFactors = FALSE) else NULL)
    },
    fastmr_clump_run_lead_graph = function(leads, snps, reference_args, plink2_bin, clump_kb,
                                           clump_r2, threads, workdir) {
      calls$cert <- calls$cert + 1L
      p <- pe_pairs(unique(leads), unique(snps), ref, clump_kb, clump_r2)
      list(lead = match(p$lead, snps), target = match(p$target, snps), invalid_r2 = 0, invalid_example = "")
    },
    fastmr_clump_plink_version = function(plink2_bin) "PLINK v2.0.0-mock",
    .package = "fastMR", .env = parent.frame()
  )
  calls
}

pe_common <- function(kb, r2, p1 = 1) list(clump_kb = kb, clump_r2 = r2, clump_p1 = p1,
                                           bfile = "mock", plink2_bin = "/bin/true")

test_that("--clump argument translation reproduces the graph's inclusive window and r2", {
  for (kb in c(0.0004, 0.001, 0.0015, 1, 4.9999, 5, 250, 250.00001, 10000, 1e6)) {
    a <- fastMR:::fastmr_clump_per_exposure_args(kb, 0.1)
    plink_radius <- trunc(as.numeric(a$kb_arg) * 1000 * (1 + pe_eps) - 1)
    graph_fmt <- as.numeric(format(kb, trim = TRUE, scientific = FALSE))
    graph_radius <- min(floor(graph_fmt * 1000 * (1 + pe_eps)), floor(kb * 1000))
    expect_identical(plink_radius, graph_radius)
    expect_gte(as.numeric(a$kb_arg), 0.001)   # PLINK2 rejects --clump-kb < 0.001
  }
  for (r2 in c(1e-6, 0.001, 0.01, 0.1, 0.2, 0.25, 0.5, 0.99, 1)) {
    a <- fastMR:::fastmr_clump_per_exposure_args(10000, r2)
    t <- r2 * (1 - pe_eps)                         # --ld-window-r2: edge iff r2 >= t
    cut <- as.numeric(a$r2_arg) * (1 + pe_eps)     # --clump: edge iff r2 > cut
    expect_lt(cut, t)
    expect_gt(cut, t * (1 - 2e-14))
    expect_true(r2 > cut)                          # r2 exactly at the threshold is in LD
    expect_lt(as.numeric(a$r2_arg), 1 - pe_eps)    # PLINK2 rejects --clump-r2 >= 1 - eps
  }
  expect_null(fastMR:::fastmr_clump_per_exposure_args(10000, 0))
  expect_identical(fastMR:::fastmr_clump_per_exposure_args(10000, 0.001)$radius_bp, 1e7)
})

test_that("per-exposure clumping is identical to graph and lead_row for E = 1, 10, 100", {
  skip_on_os("windows")
  ref <- pe_ref()
  for (E in c(1L, 10L, 100L)) {
    dat <- pe_make_dat(E, seed = E, ref)
    calls <- with_pe_oracle(ref)
    for (cfg in list(c(5, 0.5), c(8, 0.2), c(2, 0.5), c(4.9999, 0.6), c(10, 0.01))) {
      common <- pe_common(cfg[1], cfg[2], p1 = 1e-6)
      g <- do.call(fast_clump_data_graph, c(list(dat), common))
      lr <- do.call(fast_clump_data_lead_rows, c(list(dat), common))
      for (sub in c("auto", "always", "never")) {
        pe <- do.call(fast_clump_data_per_exposure, c(list(dat), common, list(subset = sub)))
        expect_identical(pe$diagnostics$partition, "per_exposure")
        expect_null(pe$diagnostics$delegated)
        expect_identical(pe$instruments, g$instruments)
        expect_identical(pe$data, g$data)
        expect_identical(pe$data, lr$data)
        expect_identical(pe$diagnostics$subset, sub != "never")
      }
      au <- do.call(fast_clump_data_auto, c(list(dat), common))
      expect_identical(au$data, g$data)
      expect_true(au$diagnostics$auto$strategy %in% c("graph", "per_exposure"))
      expect_identical(au$diagnostics$partition, au$diagnostics$auto$used)
    }
  }
})

test_that("window edge, r2 exactly at threshold and p ties match the graph", {
  skip_on_os("windows")
  ref <- pe_ref()
  s <- ref$SNP[ref$chr == "1"][1:7]
  dat <- data.frame(SNP = s, id.exposure = "E", pval.exposure = 1e-8,
                    chr_name = "1", chrom_start = seq(1000, 7000, by = 1000), stringsAsFactors = FALSE)
  with_pe_oracle(ref)
  for (cfg in list(c(5, 0.5), c(5, 0.1), c(4.999, 0.1), c(5, 0.5000001), c(6, 0.4), c(6, 0.4000001))) {
    common <- pe_common(cfg[1], cfg[2])
    g <- do.call(fast_clump_data_graph, c(list(dat), common))
    for (sub in c("always", "never")) {
      pe <- do.call(fast_clump_data_per_exposure, c(list(dat), common, list(subset = sub)))
      expect_identical(pe$instruments, g$instruments)
    }
  }
  # 5 kb window, r2 == 0.5 at exactly 5 kb: inclusive on both, so SNP 1 kills 2..6.
  pe <- fast_clump_data_per_exposure(dat, clump_kb = 5, clump_r2 = 0.5, bfile = "mock", plink2_bin = "/bin/true")
  expect_identical(pe$instruments$E, c(s[1], s[7]))
  # ties broken by SNP ID, a second exposure ordered by p, and the same SNPs in both
  dat2 <- rbind(dat, transform(dat, id.exposure = "F", pval.exposure = c(1e-7, 1e-9, 1e-8, 1e-8, 1e-6, 1e-8, 1e-8)))
  g2 <- fast_clump_data_graph(dat2, clump_kb = 5, clump_r2 = 0.5, bfile = "mock", plink2_bin = "/bin/true")
  pe2 <- fast_clump_data_per_exposure(dat2, clump_kb = 5, clump_r2 = 0.5, bfile = "mock", plink2_bin = "/bin/true")
  expect_identical(pe2$data, g2$data)
  # an exact rank column overrides p
  dat3 <- transform(dat2, pvalue_rank = c(7:1, 1:7))
  g3 <- fast_clump_data_graph(dat3, clump_kb = 5, clump_r2 = 0.5, bfile = "mock", plink2_bin = "/bin/true")
  pe3 <- fast_clump_data_per_exposure(dat3, clump_kb = 5, clump_r2 = 0.5, bfile = "mock", plink2_bin = "/bin/true")
  expect_identical(pe3$data, g3$data)
  expect_identical(pe3$instruments$E, c(s[1], s[7]))   # rank 1 = s[7] kills s[2..6]; data order kept
})

test_that("candidates absent from the reference are retained, duplicates and ineligible rows handled", {
  skip_on_os("windows")
  ref <- pe_ref()
  dat <- pe_make_dat(12L, 7L, ref, absent = 3L)
  dat$pval.exposure[c(2, 9, 30)] <- NA
  dat <- rbind(dat, dat[4:8, ])
  calls <- with_pe_oracle(ref)
  common <- pe_common(6, 0.3, p1 = 1e-6)
  g <- do.call(fast_clump_data_graph, c(list(dat), common))
  lr <- do.call(fast_clump_data_lead_rows, c(list(dat), common))
  for (sub in c("always", "never")) {
    pe <- do.call(fast_clump_data_per_exposure, c(list(dat), common, list(subset = sub)))
    expect_identical(pe$data, g$data)
    expect_identical(pe$data, lr$data)
    expect_gt(pe$diagnostics$absent_from_reference, 0L)
  }
  # nothing in the reference: every eligible candidate is retained, no --clump runs
  none <- data.frame(SNP = c("9:1:A:C", "9:2:A:C"), id.exposure = c("a", "b"), pval.exposure = 1e-9,
                     chr_name = "9", chrom_start = c(1, 2))
  calls$clump <- 0L
  pe <- fast_clump_data_per_exposure(none, clump_kb = 5, clump_r2 = 0.5, bfile = "mock", plink2_bin = "/bin/true")
  expect_identical(pe$data, none)
  expect_identical(calls$clump, 0L)
  # no eligible rows
  empty <- fast_clump_data_per_exposure(transform(none, pval.exposure = 0.5), clump_kb = 5, clump_r2 = 0.5,
                                        clump_p1 = 1e-3, bfile = "mock", plink2_bin = "/bin/true")
  expect_identical(nrow(empty$data), 0L)
})

test_that("per-exposure delegates to the graph when it cannot be exact", {
  skip_on_os("windows")
  ref <- pe_ref()
  dat <- pe_make_dat(5L, 2L, ref)
  with_pe_oracle(ref)
  g <- fast_clump_data_graph(dat, clump_kb = 5, clump_r2 = 0.5, bfile = "mock", plink2_bin = "/bin/true")
  moved <- dat
  moved$chrom_start[moved$SNP == moved$SNP[1]] <- moved$chrom_start[1] + 1
  r <- fast_clump_data_per_exposure(moved, clump_kb = 5, clump_r2 = 0.5, bfile = "mock", plink2_bin = "/bin/true",
                                    subset = "always")
  expect_identical(r$diagnostics$delegated, "reference_position_mismatch")
  expect_identical(r$diagnostics$partition, "graph")
  relabel <- dat
  relabel$chr_name[relabel$chr_name == "2"] <- "1"
  r <- fast_clump_data_per_exposure(relabel, clump_kb = 5, clump_r2 = 0.5, bfile = "mock", plink2_bin = "/bin/true",
                                    subset = "always")
  expect_identical(r$diagnostics$delegated, "reference_chromosome_mismatch")
  r <- fast_clump_data_per_exposure(dat, clump_kb = 5, clump_r2 = 0, bfile = "mock", plink2_bin = "/bin/true")
  expect_match(r$diagnostics$delegated, "clump_r2")
  nopos <- dat[setdiff(names(dat), "chrom_start")]
  r <- fast_clump_data_per_exposure(nopos, clump_kb = 5, clump_r2 = 0.5, bfile = "mock", plink2_bin = "/bin/true")
  expect_identical(r$diagnostics$delegated, "missing_positions")
  # chr label renamed consistently (e.g. "chr1" in the data): still exact via the mapping
  ren <- dat
  ren$chr_name <- paste0("chr", ren$chr_name)
  r <- fast_clump_data_per_exposure(ren, clump_kb = 5, clump_r2 = 0.5, bfile = "mock", plink2_bin = "/bin/true",
                                    subset = "always")
  expect_null(r$diagnostics$delegated)
  expect_identical(r$instruments, g$instruments)
})

test_that("a --clump result that disagrees with the graph's LD is caught by the certificate", {
  skip_on_os("windows")
  ref <- pe_ref()
  dat <- pe_make_dat(10L, 21L, ref)
  calls <- with_pe_oracle(ref)
  g <- fast_clump_data_graph(dat, clump_kb = 5, clump_r2 = 0.5, bfile = "mock", plink2_bin = "/bin/true")
  ok <- fast_clump_data_per_exposure(dat, clump_kb = 5, clump_r2 = 0.5, bfile = "mock", plink2_bin = "/bin/true",
                                     subset = "always")
  expect_true(ok$diagnostics$certificate$verified)
  expect_gt(ok$diagnostics$certificate$leads, 0L)
  expect_gt(calls$cert, 0L)
  # --clump that sees one spurious LD pair (as PLINK2 2.00a6.8 does for some |D'| = 1 pairs)
  testthat::local_mocked_bindings(
    fastmr_clump_run_clump = function(snps, reference_args, plink2_bin, kb_arg, r2_arg, stem, extract = FALSE) {
      ids <- pe_emulate_clump(snps, ref, kb_arg, r2_arg)
      if (length(ids) > 1L) ids <- ids[-length(ids)]   # the last lead "clumped" by a bogus pair
      hit <- ref[ref$SNP %in% snps, , drop = FALSE]
      list(ids = ids, pvar = if (extract) data.frame(id = hit$SNP, chr = hit$chr, pos = hit$bp) else NULL)
    }, .package = "fastMR")
  bad <- fast_clump_data_per_exposure(dat, clump_kb = 5, clump_r2 = 0.5, bfile = "mock", plink2_bin = "/bin/true",
                                      subset = "always")
  expect_match(bad$diagnostics$delegated, "certificate_mismatch")
  expect_identical(bad$diagnostics$partition, "graph")
  expect_identical(bad$data, g$data)
})

test_that("auto dispatch follows the cost model and records it", {
  skip_on_os("windows")
  ref <- pe_ref()
  with_pe_oracle(ref)
  one <- pe_make_dat(1L, 9L, ref)
  a1 <- fast_clump_data_auto(one, clump_kb = 5, clump_r2 = 0.5, bfile = "mock", plink2_bin = "/bin/true")
  expect_identical(a1$diagnostics$auto$strategy, "per_exposure")   # E = 1: subset + one --clump
  expect_true(a1$diagnostics$subset)
  plan <- fastMR:::fastmr_clump_auto_plan
  big <- data.frame(SNP = paste0("s", 1:2000), id.exposure = rep(sprintf("e%02d", 1:10), 200),
                    pval.exposure = 1e-9, chr_name = "1", chrom_start = seq_len(2000) * 10)
  expect_identical(plan(big, 10000, 1, 8L)$strategy, "per_exposure")   # ~2e6 pairs
  sparse <- transform(big, chrom_start = seq_len(2000) * 1e7, id.exposure = sprintf("e%03d", rep(1:100, 20)))
  expect_identical(plan(sparse, 10, 1, 1L)$strategy, "graph")          # no pairs at all
  expect_identical(plan(sparse, 10, 1, 1L)$estimated_pairs, 0)
  # w = max(1, min(threads, E) / 3.5); costs follow fastmr_clump_auto_model
  p <- plan(big, 10000, 1, 8L)
  m <- fastMR:::fastmr_clump_auto_model
  expect_equal(p$effective_workers, 8 / 3.5)
  expect_equal(p$cost_per_exposure, m$per_exposure_fixed + m$per_exposure_each * 10 / (8 / 3.5))
  expect_equal(p$cost_graph, m$graph_fixed + m$graph_per_pair * p$estimated_pairs + m$graph_per_row * 2000)
  # a per-exposure failure falls back to the graph with the error recorded
  testthat::local_mocked_bindings(fastmr_clump_run_clump = function(...) stop("boom"), .package = "fastMR")
  a2 <- fast_clump_data_auto(one, clump_kb = 5, clump_r2 = 0.5, bfile = "mock", plink2_bin = "/bin/true")
  expect_identical(a2$diagnostics$auto$strategy, "graph")
  expect_match(a2$diagnostics$auto$per_exposure_error, "boom")
  expect_identical(a2$data, a1$data)
})

test_that("per-exposure runs through the PLINK2 argument surface", {
  skip_on_os("windows")
  plink2 <- tempfile("fastMR_plink2_pe_stub_")
  log <- tempfile("pe_args_")
  # --make-pgen/--make-just-pvar write a .pvar of A, B, C; --clump keeps A and C.
  writeLines(c(
    "#!/bin/sh", "out=''", "clump=0", "pvar=0",
    sprintf("echo \"$@\" >> %s", log),
    "case \"$*\" in *--version*) echo 'PLINK v2.0.0-stub'; exit 0;; esac",
    "cert=0",
    "while [ \"$#\" -gt 0 ]; do case \"$1\" in --out) out=\"$2\"; shift 2;; --clump) clump=1; shift 2;; --ld-snp-list) cert=1; shift 2;; --make-pgen|--make-just-pvar) pvar=1; shift;; *) shift;; esac; done",
    "if [ $cert = 1 ]; then printf '#ID_A\\tID_B\\tPHASED_R2\\nA\\tB\\t0.9\\n' > \"${out}.vcor\"; exit 0; fi",
    "if [ $pvar = 1 ]; then printf '##fileformat\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\n1\\t1000\\tA\\tC\\tG\\n1\\t2000\\tB\\tC\\tG\\n1\\t90000\\tC\\tC\\tG\\n' > \"${out}.pvar\"; fi",
    "if [ $clump = 1 ]; then printf '#CHROM\\tPOS\\tID\\tP\\tTOTAL\\n1\\t1000\\tA\\t0.25\\t1\\n1\\t90000\\tC\\t0.75\\t0\\n' > \"${out}.clumps\"; fi"
  ), plink2)
  Sys.chmod(plink2, "0755")
  dat <- data.frame(SNP = c("A", "B", "C", "D"), id.exposure = "E", pval.exposure = c(1e-8, 1e-7, 1e-6, 1e-9),
                    chr_name = c("1", "1", "1", "5"), chrom_start = c(1000, 2000, 90000, 5))
  for (sub in c("never", "always")) {
    unlink(log)
    res <- fast_clump_data_per_exposure(dat, clump_kb = 500, clump_r2 = 0.01, bfile = "panel",
                                        plink2_bin = plink2, subset = sub)
    expect_identical(res$instruments$E, c("A", "C", "D"))   # D is absent from the reference
    args <- readLines(log)
    expect_length(args[grepl("--ld-snp-list", args)], 1L)   # the certificate query
    expect_match(args[grepl("--ld-snp-list", args)], "--r2-phased cols=id --ld-snp-list", fixed = TRUE)
    call <- args[grepl("--clump ", args)]
    expect_length(call, 1L)
    expect_match(call, "--clump-id-field SNP --clump-p-field P --clump-p1 1 --clump-p2 1", fixed = TRUE)
    expect_match(call, "--clump-r2 0.0099999999999987616 --clump-kb 500.00099999999998 --threads 1", fixed = TRUE)
    expect_identical(any(grepl("--make-pgen", args)), sub == "always")
  }
})

# Real PLINK2: phased haplotypes with r2 exactly 0.25 (A-B) and perfect LD
# (C-D) exactly 5 kb apart; the graph partition is the reference.
test_that("per-exposure and graph agree with real PLINK2 at the r2 and window boundaries", {
  skip_on_os("windows")
  plink2 <- Sys.getenv("FASTMR_PLINK2", Sys.which("plink2"))
  skip_if(!nzchar(plink2), "plink2 not available")
  hapA <- c(1, 1, 1, 1, 0, 0, 0, 0); hapB <- c(1, 1, 1, 0, 0, 0, 0, 1)   # pA = pB = .5, pAB = 3/8: r2 = 0.25
  hapC <- c(1, 0, 1, 0, 0, 1, 1, 0); hapD <- hapC                       # r2 = 1
  rep5 <- function(h) rep(h, 5)
  haps <- list(A = rep5(hapA), B = rep5(hapB), C = rep5(hapC), D = rep5(hapD), F = rep5(hapC))
  bp <- c(A = 10000, B = 11000, C = 20000, D = 25000, F = 30001)
  dir <- tempfile("fastMR_pefx_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  gt <- function(h) paste0(h[c(TRUE, FALSE)], "|", h[c(FALSE, TRUE)])
  vcf <- file.path(dir, "fx.vcf")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=1,length=100000>",
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
               paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                       sprintf("s%02d", 1:20)), collapse = "\t"),
               vapply(names(haps), function(s) paste(c("1", bp[[s]], s, "A", "G", ".", ".", ".", "GT",
                                                       gt(haps[[s]])), collapse = "\t"), "")), vcf)
  pref <- file.path(dir, "fx")
  st <- system2(plink2, c("--vcf", shQuote(vcf), "--make-pgen", "--out", shQuote(pref)), stdout = FALSE, stderr = FALSE)
  skip_if(st != 0L, "plink2 could not build the fixture")
  dat <- data.frame(SNP = names(haps), id.exposure = "E", pval.exposure = c(1e-10, 1e-9, 1e-10, 1e-9, 1e-8),
                    chr_name = "1", chrom_start = unname(bp), stringsAsFactors = FALSE)
  dat <- rbind(dat, transform(dat, id.exposure = "G", pval.exposure = rev(pval.exposure)))
  for (cfg in list(c(5, 0.25), c(5, 0.2500001), c(4.999, 0.25), c(10.001, 0.9), c(10, 1), c(250, 0.001))) {
    g <- fast_clump_data_graph(dat, clump_kb = cfg[1], clump_r2 = cfg[2], pfile = pref, plink2_bin = plink2)
    for (sub in c("never", "always")) {
      pe <- suppressWarnings(fast_clump_data_per_exposure(dat, clump_kb = cfg[1], clump_r2 = cfg[2], pfile = pref,
                                                          plink2_bin = plink2, subset = sub))
      expect_null(pe$diagnostics$delegated)
      expect_identical(pe$data, g$data, info = paste(cfg, collapse = "/"))
    }
  }
  # the boundaries are really exercised: r2 == 0.25 is in LD, 5 kb is inside a 5 kb window
  g <- fast_clump_data_graph(dat, clump_kb = 5, clump_r2 = 0.25, pfile = pref, plink2_bin = plink2)
  expect_identical(g$instruments$E, c("A", "C", "F"))
  g <- fast_clump_data_graph(dat, clump_kb = 5, clump_r2 = 0.2500001, pfile = pref, plink2_bin = plink2)
  expect_identical(g$instruments$E, c("A", "B", "C", "F"))
  g <- fast_clump_data_graph(dat, clump_kb = 4.999, clump_r2 = 0.25, pfile = pref, plink2_bin = plink2)
  expect_identical(g$instruments$E, c("A", "C", "D", "F"))
})
