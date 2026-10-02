# Exact LD-graph multi-exposure clumping.
#
# One PLINK2 all-pairs call per chromosome over the union of candidate SNPs
# replaces the per-lead (or per-frontier-round) process launches.  The pair set
# is the symmetric closure of the pairs PLINK2 reports at r2 >= clump_r2 within
# clump_kb, which is exactly the set of pairs the lead-row and frontier
# strategies query; each exposure is then clumped greedily, in the same
# (p, SNP) order, against that graph in C++.

# Internal LD oracle for the graph strategy (mockable in tests/benchmarks).
# Returns a data frame with columns `lead` and `target` (each reported pair
# once or twice; the caller symmetrises).
fastmr_clump_run_graph <- function(snps, reference_args, plink2_bin, clump_kb,
                                   clump_r2, threads, workdir, tag,
                                   ld_window_variants = 1e9) {
  stem <- file.path(workdir, paste0("graph_", tag))
  extract_file <- paste0(stem, ".extract.txt")
  writeLines(as.character(snps), extract_file)
  args <- c(reference_args, "--extract", fastmr_clump_quote(extract_file),
            "--r2-unphased", "zs", "cols=id", "--ld-window-kb", format(clump_kb, trim = TRUE, scientific = FALSE),
            "--ld-window", format(ld_window_variants, trim = TRUE, scientific = FALSE),
            "--ld-window-r2", format(clump_r2, trim = TRUE),
            "--threads", as.integer(threads), "--out", fastmr_clump_quote(stem))
  output <- tryCatch(suppressWarnings(system2(plink2_bin, args, stdout = TRUE, stderr = TRUE)),
                     error = function(e) structure(character(), status = 1L, error = conditionMessage(e)))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  path <- paste0(stem, ".vcor.zst")
  if (status != 0L && any(grepl("No variants remaining after", output, fixed = TRUE))) {
    # None of these candidates is in the reference: no LD edges (lead-row
    # semantics keep such SNPs, so the graph must not abort).
    unlink(extract_file)
    return(data.frame(lead = character(), target = character(), stringsAsFactors = FALSE))
  }
  if (status != 0L || !file.exists(path)) {
    detail <- attr(output, "error")
    if (is.null(detail)) detail <- paste(utils::tail(output, 8L), collapse = " | ")
    stop("PLINK2 all-pairs LD graph query failed (", tag, "): ", detail, call. = FALSE)
  }
  on.exit(unlink(c(path, extract_file)), add = TRUE)
  zstdcat <- Sys.which("zstdcat")
  if (nzchar(zstdcat)) cmd_args <- shQuote(path) else {
    zstdcat <- Sys.which("zstd")
    if (!nzchar(zstdcat)) stop("PLINK produced a .vcor.zst file but neither zstdcat nor zstd is available", call. = FALSE)
    cmd_args <- c("-dc", shQuote(path))
  }
  # Stream the decompressed table in blocks so that only a block of lines (not
  # the whole ~200 B/edge line vector) is ever resident; each block is reduced
  # to its two ID columns, which share interned CHARSXPs.
  # "rb": a text-mode pipe silently drops an unterminated final line.
  con <- pipe(paste(shQuote(zstdcat), paste(cmd_args, collapse = " ")), "rb")
  closed <- FALSE
  on.exit(if (!closed) try(close(con), silent = TRUE), add = TRUE)
  fa <- NA_integer_; fb <- NA_integer_  # set from the header
  first <- TRUE
  leads <- list(); targets <- list()
  repeat {
    block <- readLines(con, n = 1000000L, warn = FALSE)
    if (!length(block)) break
    if (first) {
      first <- FALSE
      if (startsWith(block[[1L]], "#")) {
        header <- strsplit(sub("^#", "", block[[1L]]), "\t", fixed = TRUE)[[1L]]
        ia <- match("ID_A", header); ib <- match("ID_B", header)
        if (!is.na(ia) && !is.na(ib) && ia < ib) { fa <- ia; fb <- ib }
      }
      if (is.na(fa)) {
        stop("unrecognised PLINK2 .vcor header (need ID_A and ID_B): ",
             substr(block[[1L]], 1L, 200L), call. = FALSE)
      }
    }
    ids <- .fastmr_vcor_ids(block, fa, fb)
    data_lines <- sum(nzchar(block) & !startsWith(block, "#"))
    if (length(ids$lead) != data_lines) {
      stop("malformed PLINK2 .vcor line(s): ", data_lines - length(ids$lead),
           " line(s) have fewer than ", fb, " fields", call. = FALSE)
    }
    leads[[length(leads) + 1L]] <- ids$lead
    targets[[length(targets) + 1L]] <- ids$target
  }
  status <- close(con)
  closed <- TRUE
  if (!is.null(status) && !identical(as.integer(status), 0L)) {
    stop("could not decompress PLINK LD output", call. = FALSE)
  }
  data.frame(lead = as.character(unlist(leads, use.names = FALSE)),
             target = as.character(unlist(targets, use.names = FALSE)),
             stringsAsFactors = FALSE)
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
#' Clump every exposure with a single PLINK2 all-pairs LD call per chromosome
#' over the union of candidate SNPs (`--r2-unphased` with `--ld-window-kb`
#' equal to `clump_kb`, `--ld-window-r2` equal to `clump_r2` and a very large
#' `--ld-window` variant count), then run each exposure's greedy clump in C++
#' against the symmetrised graph.  Retained instruments are identical to
#' [fast_clump_data_lead_rows()] and [fast_clump_data_batched()].
#'
#' The expected number of pairs is estimated from sorted positions before any
#' PLINK call.  A chromosome whose estimate exceeds `max_graph_pairs` (or any
#' data set lacking complete positions) is clumped with
#' [fast_clump_data_lead_rows()] instead, and the reason is recorded in
#' `diagnostics$fallbacks`.
#'
#' @inheritParams fast_clump_data_lead_rows
#' @param max_graph_pairs Per-chromosome cap on the estimated number of
#'   candidate pairs within the window before falling back to lead-row mode.
#' @return A list with `data`, named `instruments`, and `diagnostics`
#'   (including `ld_provenance`).
#' @export
fast_clump_data_graph <- function(
    dat, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1,
    bfile = NULL, pfile = NULL, plink2_bin = NULL, threads = 1L,
    max_graph_pairs = 5e7, max_pair_requests = 2e8, max_target_variants = 2e6,
    max_rounds = 10000L, workdir = NULL, reference_manifest = NULL) {
  if (!is.data.frame(dat) || !"SNP" %in% names(dat)) stop("dat must contain SNP", call. = FALSE)
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
  dedup <- !duplicated(paste(as.character(dat$id.exposure), as.character(dat$SNP), sep = "\r"))
  dat <- dat[dedup, , drop = FALSE]
  p <- suppressWarnings(as.numeric(as.character(dat[[pcol]])))
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
  positions_ok <- !anyNA(position$chr[elig_rows]) && all(is.finite(position$bp[elig_rows]))
  if (length(elig_rows) && !positions_ok) {
    retained_key <- lead_row_fallback(elig_rows, "missing_positions", "all")
  } else if (length(elig_rows)) {
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
      tag <- gsub("[^A-Za-z0-9_.-]", "_", cc)
      ld <- if (est > 0) {
        graph_calls <- graph_calls + 1L
        fastmr_clump_run_graph(usnp, reference_args, plink2_bin, clump_kb, clump_r2,
                               threads, workdir, tag)
      } else {
        # No candidate pair within the window: the graph is empty.
        data.frame(lead = character(), target = character(), stringsAsFactors = FALSE)
      }
      a <- match(as.character(ld$lead), usnp)
      b <- match(as.character(ld$target), usnp)
      ok <- !is.na(a) & !is.na(b) & a != b
      a <- a[ok]; b <- b[ok]
      near <- abs(ubp[a] - ubp[b]) <= window
      a <- a[near]; b <- b[near]
      info$edges <- length(a)
      total_edges <- total_edges + length(a)
      vtx <- match(snp[rows], usnp)
      e_id <- match(expo[rows], unique(expo[rows]))
      ord <- if (is.null(rank)) order(e_id, p[rows], snp[rows], method = "radix")
             else order(e_id, rank[rows], p[rows], snp[rows], method = "radix")
      rows <- rows[ord]; vtx <- vtx[ord]; e_id <- e_id[ord]
      starts <- c(0L, cumsum(tabulate(e_id)))
      keep <- .fastmr_graph_clump(length(usnp), a - 1L, b - 1L, vtx - 1L, as.integer(starts))
      retained_key <- c(retained_key, paste(expo[rows[keep]], snp[rows[keep]], sep = "\r"))
      chromosomes[[cc]] <- info
    }
  }
  original_key <- paste(as.character(original$id.exposure), as.character(original$SNP), sep = "\r")
  result <- original[original_key %in% retained_key, , drop = FALSE]
  instruments <- lapply(split(result$SNP, result$id.exposure, drop = TRUE), as.character)
  ld_provenance <- list(
    plink2_version = fastmr_clump_plink_version(plink2_bin),
    mode = "all_pairs_graph",
    flags = c("--r2-unphased zs cols=id", paste("--ld-window-kb", format(clump_kb, trim = TRUE, scientific = FALSE)),
              "--ld-window 1000000000", paste("--ld-window-r2", format(clump_r2, trim = TRUE))),
    reference = paste(gsub("'", "", reference_args), collapse = " "),
    reference_manifest_md5 = reference_md5,
    clump_kb = clump_kb, clump_r2 = clump_r2, clump_p1 = clump_p1
  )
  list(data = result, instruments = instruments,
       diagnostics = list(
         exposures = length(unique(as.character(original$id.exposure))),
         candidate_rows = nrow(dat), retained = length(retained_key),
         rounds = fallback_rounds, plink_calls = graph_calls + fallback_calls,
         graph_calls = graph_calls, fallback_plink_calls = fallback_calls,
         graph_edges = total_edges, logical_pairs = total_edges + fallback_pairs,
         exact = TRUE, fallback = length(fallbacks) > 0L, fallbacks = fallbacks,
         partition = "graph", chromosomes = chromosomes,
         max_graph_pairs = max_graph_pairs,
         reference_manifest_md5 = reference_md5, ld_provenance = ld_provenance))
}
