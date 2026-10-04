# Exact LD-graph multi-exposure clumping.
#
# One PLINK2 all-pairs call over the union of candidate SNPs (split into a few
# calls only when the estimated pair count exceeds max_graph_pairs) replaces
# the per-lead (or per-frontier-round) process launches.  The pair set
# is the symmetric closure of the pairs PLINK2 reports at r2 >= clump_r2 within
# clump_kb, which is exactly the set of pairs the lead-row and frontier
# strategies query; each exposure is then clumped greedily, in the same
# (p, SNP) order, against that graph in C++.
#
# The LD statistic is `--r2-phased` (haplotype-frequency r2, EM-estimated when
# phase is unknown), which is what both `plink --clump` (1.9) and
# `plink2 --clump` use by default.  `--r2-unphased` (squared dosage
# correlation) is a different estimator: for rare variants it can sit on the
# other side of a low clump_r2 threshold and change the retained set.

# Internal LD oracle for the graph strategy (mockable in tests/benchmarks).
# Returns a list of 1-based integer vertex ids `lead` and `target` (positions
# in `snps`; each reported pair once or twice; the caller symmetrises).  PLINK2
# writes an UNCOMPRESSED .vcor (about 30 bytes/edge, e.g. ~200 MB for 6M edges)
# into `workdir`; it is parsed in C++ and deleted.  No zstd is needed.
fastmr_clump_run_graph <- function(snps, reference_args, plink2_bin, clump_kb,
                                   clump_r2, threads, workdir, tag,
                                   ld_window_variants = 1e9) {
  stem <- file.path(workdir, paste0("graph_", tag))
  extract_file <- paste0(stem, ".extract.txt")
  writeLines(as.character(snps), extract_file)
  args <- c(reference_args, "--extract", fastmr_clump_quote(extract_file),
            "--r2-phased", "cols=id", "--ld-window-kb", format(clump_kb, trim = TRUE, scientific = FALSE),
            "--ld-window", format(ld_window_variants, trim = TRUE, scientific = FALSE),
            "--ld-window-r2", format(clump_r2, trim = TRUE),
            "--threads", as.integer(threads), "--out", fastmr_clump_quote(stem))
  output <- tryCatch(suppressWarnings(system2(plink2_bin, args, stdout = TRUE, stderr = TRUE)),
                     error = function(e) structure(character(), status = 1L, error = conditionMessage(e)))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  path <- paste0(stem, ".vcor")
  if (status != 0L && any(grepl("No variants remaining after", output, fixed = TRUE))) {
    # None of these candidates is in the reference: no LD edges (lead-row
    # semantics keep such SNPs, so the graph must not abort).
    unlink(extract_file)
    return(list(lead = integer(), target = integer()))
  }
  if (status != 0L || !file.exists(path)) {
    detail <- attr(output, "error")
    if (is.null(detail)) detail <- paste(utils::tail(output, 8L), collapse = " | ")
    stop("PLINK2 all-pairs LD graph query failed (", tag, "): ", detail, call. = FALSE)
  }
  on.exit(unlink(c(path, extract_file)), add = TRUE)
  # Parsed natively (buffered fread, hash map ID -> vertex): an ID outside
  # `snps` is an error.  Returns 1-based integer vertex ids.
  ld <- .fastmr_vcor_read(path, as.character(snps))
  fastmr_clump_warn_invalid_r2(ld$invalid_r2, ld$invalid_example)
  ld
}

# Internal (mockable): the subset of `snps` present in the reference, from one
# PLINK2 --extract --write-snplist query.  character() when none is present.
fastmr_clump_reference_ids <- function(snps, reference_args, plink2_bin, threads, stem) {
  extract_file <- paste0(stem, ".extract.txt")
  writeLines(unique(as.character(snps)), extract_file)
  on.exit(unlink(c(extract_file, paste0(stem, c(".snplist", ".log")))), add = TRUE)
  output <- fastmr_clump_system2(plink2_bin, c(
    reference_args, "--extract", fastmr_clump_quote(extract_file), "--write-snplist", "allow-dups",
    "--threads", as.integer(threads), "--out", fastmr_clump_quote(stem)))
  if (fastmr_clump_no_variants(output)) return(character())
  path <- paste0(stem, ".snplist")
  if (attr(output, "status") != 0L || !file.exists(path)) {
    fastmr_clump_plink_failed(output, "reference membership query (--write-snplist)")
  }
  unique(readLines(path))
}

# PLINK2 2.00a6.8 (Jan 2025) --r2-phased (and its --clump) reports impossible
# values (r2 ~ 96, D' ~ -220) for some pairs whose true |D'| is 1; 2.00a6
# (Oct 2024) and PLINK 1.9 give the correct small r2.  Such a pair passes any
# --ld-window-r2 filter, so the lower-ranked SNP is clumped away exactly as
# PLINK2's own --clump does, but not as PLINK 1.9 does.
fastmr_clump_warn_invalid_r2 <- function(n, example) {
  if (!length(n) || !isTRUE(n > 0)) return(invisible(FALSE))
  warning(sprintf(paste0("PLINK2 reported %.0f LD pair(s) with r2 > 1 (e.g. %s); this is a PLINK2 ",
                         "--r2-phased bug (seen in 2.00a6.8). Those pairs were treated as in LD, as ",
                         "PLINK2 --clump would; use a PLINK2 build without the bug to match PLINK 1.9."),
                  n, example), call. = FALSE)
  invisible(TRUE)
}

fastmr_clump_plink_version <- function(plink2_bin) {
  out <- tryCatch(suppressWarnings(system2(plink2_bin, "--version", stdout = TRUE, stderr = TRUE)),
                  error = function(e) character())
  if (!length(out) || !is.null(attr(out, "status"))) return(NA_character_)
  out[[1L]]
}

# Number of unordered pairs within `window` bp, from sorted positions.
fastmr_graph_pair_estimate <- function(bp, window) {
  bp <- sort(bp)
  if (length(bp) < 2L) return(0)
  sum(as.numeric(findInterval(bp + window, bp) - seq_along(bp)))
}

#' Exact LD-graph multi-exposure clumping
#'
#' Clump every exposure with a single PLINK2 all-pairs LD call over the union
#' of candidate SNPs (`--r2-phased` with `--ld-window-kb` equal to `clump_kb`,
#' `--ld-window-r2` equal to `clump_r2` and a very large `--ld-window` variant
#' count), then run each exposure's greedy clump in C++ against the
#' symmetrised graph.  Retained instruments are identical to
#' [fast_clump_data_lead_rows()] and [fast_clump_data_batched()].
#'
#' The expected number of pairs is estimated from sorted positions before any
#' PLINK call.  Chromosomes share a call while their summed estimate stays
#' within `max_graph_pairs`.  A chromosome whose own estimate exceeds it (or
#' any data set lacking complete positions) is clumped with
#' [fast_clump_data_lead_rows()] instead, and the reason is recorded in
#' `diagnostics$fallbacks`.
#'
#' Within an exposure, candidates are taken in order of p (exact rank first
#' when a `pvalue_rank` column is present); equal p, including p = 0 from
#' underflow at |z| above about 38, is broken by larger |z| (from
#' `beta.exposure / se.exposure` when present) and then by SNP ID.  When an
#' (exposure, SNP) pair appears in several rows, the smallest p orders it, and
#' every row of a retained pair is returned.
#'
#' Candidate membership in the reference is checked with one PLINK2
#' `--write-snplist` query.  Candidates absent from the reference have no LD
#' edges; see `absent`.  If no eligible candidate is in the reference, the call
#' stops (a SNP-ID scheme mismatch is the usual cause).
#'
#' @inheritParams fast_clump_data_lead_rows
#' @param max_graph_pairs Cap on the estimated number of candidate pairs
#'   within the window per PLINK2 call; a chromosome above it falls back to
#'   lead-row mode.
#' @param absent What to do with eligible candidates absent from the LD
#'   reference: `"keep"` (default; retained unclumped, the historical
#'   behaviour) or `"drop"` (removed, as TwoSampleMR local clumping does).
#'   Either way a warning reports how many there are.
#' @return A list with `data`, named `instruments`, and `diagnostics`
#'   (including `ld_provenance` and `absent_from_reference`).
#' @export
fast_clump_data_graph <- function(
    dat, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1,
    bfile = NULL, pfile = NULL, plink2_bin = NULL, threads = 1L,
    max_graph_pairs = 5e7, max_pair_requests = 2e8, max_target_variants = 2e6,
    max_rounds = 10000L, workdir = NULL, reference_manifest = NULL,
    absent = c("keep", "drop")) {
  if (!is.data.frame(dat) || !"SNP" %in% names(dat)) stop("dat must contain SNP", call. = FALSE)
  absent <- match.arg(absent)
  clump_kb <- fastmr_clump_number(clump_kb, "clump_kb", 0)
  clump_r2 <- fastmr_clump_number(clump_r2, "clump_r2", 0, 1)
  clump_p1 <- fastmr_clump_number(clump_p1, "clump_p1", 0, 1)
  threads <- as.integer(fastmr_clump_number(threads, "threads", 1))
  max_graph_pairs <- fastmr_clump_number(max_graph_pairs, "max_graph_pairs", 0)
  if (is.null(plink2_bin)) plink2_bin <- Sys.which("plink2")
  if (!nzchar(plink2_bin)) stop("PLINK2 executable not found; provide plink2_bin", call. = FALSE)
  reference_args <- fastmr_clump_reference_args(bfile, pfile)
  if (!is.null(reference_manifest) &&
      (length(reference_manifest) != 1L || !is.character(reference_manifest) ||
       is.na(reference_manifest) || !file.exists(reference_manifest))) {
    stop("reference_manifest must be one existing file when supplied", call. = FALSE)
  }
  reference_md5 <- if (is.null(reference_manifest)) NULL else unname(tools::md5sum(reference_manifest))
  pcol <- if ("pval.exposure" %in% names(dat)) "pval.exposure" else if ("pval.outcome" %in% names(dat)) "pval.outcome" else NULL
  if (is.null(pcol)) {
    dat$pval.exposure <- 0.99
    pcol <- "pval.exposure"
  }
  if (!"id.exposure" %in% names(dat)) dat$id.exposure <- "exposure"
  if (anyNA(dat$SNP) || any(!nzchar(trimws(as.character(dat$SNP)))) || anyNA(dat$id.exposure)) {
    stop("SNP and id.exposure must be non-missing and non-empty", call. = FALSE)
  }
  original <- dat
  dat <- dat[fastmr_clump_dedup(dat, pcol), , drop = FALSE]
  p <- suppressWarnings(as.numeric(as.character(dat[[pcol]])))
  tiebreak <- fastmr_clump_tiebreak(dat, pcol)
  snp <- as.character(dat$SNP)
  expo <- as.character(dat$id.exposure)
  rank <- dat[["pvalue_rank"]]
  position <- fastmr_clump_position(dat)
  eligible <- is.finite(p) & p <= clump_p1
  window <- clump_kb * 1000
  workdir_owned <- is.null(workdir)
  if (workdir_owned) workdir <- tempfile("fastMR_graph_")
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  if (workdir_owned) on.exit(unlink(workdir, recursive = TRUE, force = TRUE), add = TRUE)

  retained_key <- character()
  fallbacks <- list()
  chromosomes <- list()
  graph_calls <- 0L
  fallback_calls <- 0L
  fallback_pairs <- 0
  fallback_rounds <- 0L
  total_edges <- 0
  lead_row_fallback <- function(rows, reason, label) {
    part <- original[paste(as.character(original$id.exposure), as.character(original$SNP), sep = "\r") %in%
                       paste(expo[rows], snp[rows], sep = "\r"), , drop = FALSE]
    ans <- fast_clump_data_lead_rows(
      part, clump_kb = clump_kb, clump_r2 = clump_r2, clump_p1 = clump_p1,
      bfile = bfile, pfile = pfile, plink2_bin = plink2_bin, threads = threads,
      max_pair_requests = max_pair_requests, max_target_variants = max_target_variants,
      max_rounds = max_rounds, workdir = file.path(workdir, paste0("fb_", length(fallbacks) + 1L)))
    fallback_calls <<- fallback_calls + ans$diagnostics$plink_calls
    fallback_pairs <<- fallback_pairs + ans$diagnostics$logical_pairs
    fallback_rounds <<- fallback_rounds + ans$diagnostics$rounds
    fallbacks[[length(fallbacks) + 1L]] <<- list(scope = label, reason = reason, rows = length(rows),
                                                  strategy = "lead_row")
    paste(as.character(ans$data$id.exposure), as.character(ans$data$SNP), sep = "\r")
  }

  elig_rows <- which(eligible)
  absent_rows <- integer()
  if (length(elig_rows)) {
    present <- fastmr_clump_reference_ids(unique(snp[elig_rows]), reference_args, plink2_bin,
                                          threads, file.path(workdir, "reference_ids"))
    absent_rows <- elig_rows[!snp[elig_rows] %in% present]
    fastmr_clump_check_absent(length(absent_rows), length(elig_rows),
                              length(unique(snp[absent_rows])), absent)
  }
  positions_ok <- !anyNA(position$chr[elig_rows]) && all(is.finite(position$bp[elig_rows]))
  if (length(elig_rows) && !positions_ok) {
    retained_key <- lead_row_fallback(elig_rows, "missing_positions", "all")
  } else if (length(elig_rows)) {
    # Plan every chromosome first, then query LD for as many chromosomes per
    # PLINK2 call as fit under max_graph_pairs: each call re-parses the whole
    # reference .pvar (~1 s for 9M variants), which dominated one-call-per-
    # chromosome runs.  PLINK2 --r2 only pairs variants on the same
    # chromosome, so the union query returns exactly the per-chromosome edges,
    # and the per-call cap keeps the parsed edge count bounded as before.
    plan <- list()
    for (cc in unique(position$chr[elig_rows])) {
      rows <- elig_rows[position$chr[elig_rows] == cc]
      usnp <- unique(snp[rows])
      ubp <- position$bp[rows][match(usnp, snp[rows])]
      est <- fastmr_graph_pair_estimate(ubp, window)
      info <- list(snps = length(usnp), rows = length(rows), estimated_pairs = est)
      if (est > max_graph_pairs) {
        reason <- sprintf("estimated_pairs %.0f > max_graph_pairs %.0f", est, max_graph_pairs)
        info$fallback <- reason
        retained_key <- c(retained_key, lead_row_fallback(rows, reason, paste0("chr", cc)))
        chromosomes[[cc]] <- info
        next
      }
      plan[[length(plan) + 1L]] <- list(cc = cc, rows = rows, usnp = usnp, ubp = ubp, est = est,
                                         info = info, a = integer(), b = integer())
    }
    # Greedy batches in chromosome order; a batch never repeats a SNP ID, so
    # every reported ID maps to exactly one chromosome.
    batch <- integer(length(plan)); nb <- 0L; load <- Inf; seen <- character()
    for (i in seq_along(plan)) {
      if (plan[[i]]$est <= 0) next
      if (load + plan[[i]]$est > max_graph_pairs || any(plan[[i]]$usnp %in% seen)) {
        nb <- nb + 1L; load <- 0; seen <- character()
      }
      batch[i] <- nb; load <- load + plan[[i]]$est; seen <- c(seen, plan[[i]]$usnp)
    }
    for (k in seq_len(nb)) {
      members <- which(batch == k)
      sizes <- vapply(plan[members], function(x) length(x$usnp), integer(1))
      offsets <- c(0L, cumsum(sizes))
      tag <- if (length(members) == 1L) gsub("[^A-Za-z0-9_.-]", "_", plan[[members]]$cc) else paste0("batch", k)
      graph_calls <- graph_calls + 1L
      ld <- fastmr_clump_run_graph(unlist(lapply(plan[members], `[[`, "usnp"), use.names = FALSE),
                                   reference_args, plink2_bin, clump_kb, clump_r2, threads, workdir, tag)
      ca <- findInterval(ld$lead - 1L, offsets, rightmost.closed = FALSE)
      cb <- findInterval(ld$target - 1L, offsets, rightmost.closed = FALSE)
      same <- ca == cb
      for (j in seq_along(members)) {
        sel <- same & ca == j
        plan[[members[j]]]$a <- ld$lead[sel] - offsets[j]
        plan[[members[j]]]$b <- ld$target[sel] - offsets[j]
      }
    }
    for (x in plan) {
      # est == 0: no candidate pair within the window, so the graph is empty.
      rows <- x$rows; usnp <- x$usnp; ubp <- x$ubp; info <- x$info
      a <- x$a
      b <- x$b
      ok <- a != b
      a <- a[ok]; b <- b[ok]
      near <- abs(ubp[a] - ubp[b]) <= window
      a <- a[near]; b <- b[near]
      info$edges <- length(a)
      total_edges <- total_edges + length(a)
      vtx <- match(snp[rows], usnp)
      e_id <- match(expo[rows], unique(expo[rows]))
      ord <- fastmr_clump_order(p[rows], snp[rows], rank[rows], tiebreak[rows], group = e_id)
      rows <- rows[ord]; vtx <- vtx[ord]; e_id <- e_id[ord]
      starts <- c(0L, cumsum(tabulate(e_id)))
      keep <- .fastmr_graph_clump(length(usnp), a - 1L, b - 1L, vtx - 1L, as.integer(starts))
      retained_key <- c(retained_key, paste(expo[rows[keep]], snp[rows[keep]], sep = "\r"))
      chromosomes[[x$cc]] <- info
    }
  }
  if (identical(absent, "drop") && length(absent_rows)) {
    retained_key <- setdiff(retained_key, paste(expo[absent_rows], snp[absent_rows], sep = "\r"))
  }
  original_key <- paste(as.character(original$id.exposure), as.character(original$SNP), sep = "\r")
  result <- original[original_key %in% retained_key, , drop = FALSE]
  instruments <- lapply(split(result$SNP, result$id.exposure, drop = TRUE), as.character)
  ld_provenance <- list(
    plink2_version = fastmr_clump_plink_version(plink2_bin),
    mode = "all_pairs_graph",
    flags = c("--r2-phased cols=id", paste("--ld-window-kb", format(clump_kb, trim = TRUE, scientific = FALSE)),
              "--ld-window 1000000000", paste("--ld-window-r2", format(clump_r2, trim = TRUE))),
    reference = paste(gsub("'", "", reference_args), collapse = " "),
    reference_manifest_md5 = reference_md5,
    clump_kb = clump_kb, clump_r2 = clump_r2, clump_p1 = clump_p1
  )
  list(data = result, instruments = instruments,
       diagnostics = list(
         exposures = length(unique(as.character(original$id.exposure))),
         candidate_rows = nrow(dat), retained = length(retained_key),
         rounds = fallback_rounds,
         plink_calls = graph_calls + fallback_calls + as.integer(length(elig_rows) > 0L),
         graph_calls = graph_calls, fallback_plink_calls = fallback_calls,
         graph_edges = total_edges, logical_pairs = total_edges + fallback_pairs,
         exact = TRUE, fallback = length(fallbacks) > 0L, fallbacks = fallbacks,
         partition = "graph", chromosomes = chromosomes,
         absent_from_reference = length(absent_rows), absent = absent,
         max_graph_pairs = max_graph_pairs,
         reference_manifest_md5 = reference_md5, ld_provenance = ld_provenance))
}
