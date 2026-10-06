# Exact per-exposure clumping on a candidate-only reference, and the size-based
# dispatch between it and the all-pairs graph.
#
# At permissive settings (r2 0.001, 10 Mb) the all-pairs candidate graph is
# nearly complete (tens of millions of edges for a few exposures), while the
# greedy clump only ever needs the pairs that involve a lead.  PLINK2's own
# --clump evaluates exactly those, so for many candidate pairs it is cheaper to
# extract the union of candidates once (`--extract --make-pgen`) and run one
# single-threaded `--clump` per exposure on that subset, in parallel.
#
# Exactness against the graph partition rests on four details:
# * P is a unique exact rank (the graph's (rank, p, SNP) order), never the
#   p-value itself: real p-values have many ties, which PLINK2 breaks by file
#   order rather than by SNP ID.
# * The window and r2 arguments are translated so that PLINK2 --clump's strict
#   comparisons (r2 > r2_thresh * (1 + eps), |dbp| <= kb * 1000 * (1 + eps) - 1)
#   reproduce the graph's inclusive ones (--ld-window-r2: r2 >= r2 * (1 - eps);
#   --ld-window-kb plus |dbp| <= clump_kb * 1000), see
#   fastmr_clump_per_exposure_args().
# * Candidates absent from the reference are retained (they have no LD edges
#   in the graph either), or dropped with absent = "drop".
# * Candidate positions must equal the reference positions (the graph windows
#   on the data positions, --clump on the reference); otherwise, or when the
#   reference has duplicate candidate IDs, the call is delegated to the graph.
#
# --clump and the graph's --r2-phased are separate PLINK2 code paths, so the
# per-exposure leads are certified against the graph's own LD statistic with
# one lead-restricted --r2-phased query (see fast_clump_data_per_exposure()).

# PLINK2's kSmallEpsilon (2^-44), used when it parses --clump-kb/--clump-r2 and
# --ld-window-kb/--ld-window-r2.
fastmr_plink2_small_epsilon <- 2^-44

# --clump-kb and --clump-r2 values that make PLINK2 --clump decide "in LD
# within the window" exactly as the graph partition does for the same
# clump_kb/clump_r2.  Returns NULL when no exact translation exists
# (clump_r2 <= 0: --clump would drop r2 == 0 pairs that the graph keeps).
fastmr_clump_per_exposure_args <- function(clump_kb, clump_r2) {
  eps <- fastmr_plink2_small_epsilon
  # The graph passes format(x) (7 significant digits) to PLINK2 and then
  # filters |dbp| <= clump_kb * 1000 on the (verified equal) data positions.
  kb_graph <- as.numeric(format(clump_kb, trim = TRUE, scientific = FALSE))
  r2_graph <- as.numeric(format(clump_r2, trim = TRUE))
  if (!is.finite(r2_graph) || r2_graph <= 0) return(NULL)
  radius <- min(floor(kb_graph * 1000 * (1 + eps)), floor(clump_kb * 1000), 2147483646)
  # --clump: bp_radius = (int)(kb * 1000 * (1 + eps) - 1); passing
  # (radius + 1) / 1000 gives exactly `radius` (the eps bump dominates the
  # parse error of a 17-significant-digit decimal).
  kb_arg <- sprintf("%.17g", (radius + 1) / 1000)
  # graph edge: r2 >= t, t = r2_graph * (1 - eps).  --clump edge:
  # r2 > arg * (1 + eps).  Aim the --clump cut a relative 1e-14 below t: PLINK2
  # parses decimals to within a few ulps (~1e-15), so the cut stays strictly
  # below t and only r2 within [t * (1 - 1e-14), t) -- far below the precision
  # of an EM r2 estimate -- could be classified differently.
  t <- r2_graph * (1 - eps)
  r2_arg <- sprintf("%.17g", t * (1 - 1e-14) / (1 + eps))
  list(kb_arg = kb_arg, r2_arg = r2_arg, radius_bp = radius, r2_cut = t)
}

# Reads ID/CHROM/POS from a .pvar (header line starting with #CHROM).
fastmr_clump_read_pvar <- function(path) {
  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)
  skip <- 0L
  header <- NULL
  repeat {
    line <- readLines(con, n = 1L)
    if (!length(line)) break
    if (startsWith(line, "##")) { skip <- skip + 1L; next }
    header <- line
    break
  }
  if (is.null(header) || !startsWith(header, "#CHROM")) {
    stop("unrecognised .pvar header in ", path, call. = FALSE)
  }
  cols <- strsplit(sub("^#", "", header), "\t", fixed = TRUE)[[1L]]
  need <- match(c("CHROM", "POS", "ID"), cols)
  if (anyNA(need)) stop("a .pvar lacks CHROM/POS/ID columns: ", path, call. = FALSE)
  classes <- rep("NULL", length(cols))
  classes[need] <- c("character", "numeric", "character")
  x <- utils::read.table(path, sep = "\t", header = FALSE, skip = skip + 1L,
                         colClasses = classes, comment.char = "", quote = "",
                         na.strings = character(), col.names = cols)
  data.frame(id = x$ID, chr = x$CHROM, pos = x$POS, stringsAsFactors = FALSE)
}

fastmr_clump_plink_failed <- function(output, what) {
  detail <- attr(output, "error")
  if (is.null(detail)) detail <- paste(utils::tail(output, 8L), collapse = " | ")
  stop("PLINK2 ", what, " failed: ", detail, call. = FALSE)
}

fastmr_clump_system2 <- function(plink2_bin, args) {
  output <- tryCatch(suppressWarnings(system2(plink2_bin, args, stdout = TRUE, stderr = TRUE)),
                     error = function(e) structure(character(), status = 1L, error = conditionMessage(e)))
  status <- attr(output, "status")
  attr(output, "status") <- if (is.null(status)) 0L else as.integer(status)
  output
}

fastmr_clump_no_variants <- function(output) {
  attr(output, "status") != 0L && any(grepl("No variants remaining after", output, fixed = TRUE))
}

# Internal (mockable): extract the candidate union into a small .pgen.
# Returns list(reference_args, pvar) or NULL when no candidate is in the
# reference.
fastmr_clump_make_subset <- function(snps, reference_args, plink2_bin, threads, stem) {
  extract_file <- paste0(stem, ".extract.txt")
  writeLines(as.character(snps), extract_file)
  on.exit(unlink(extract_file), add = TRUE)
  output <- fastmr_clump_system2(plink2_bin, c(
    reference_args, "--extract", fastmr_clump_quote(extract_file), "--make-pgen",
    "--threads", as.integer(threads), "--out", fastmr_clump_quote(stem)))
  if (fastmr_clump_no_variants(output)) return(NULL)
  if (attr(output, "status") != 0L || !file.exists(paste0(stem, ".pvar"))) {
    fastmr_clump_plink_failed(output, "candidate reference subset")
  }
  list(reference_args = c("--pfile", fastmr_clump_quote(stem)),
       pvar = fastmr_clump_read_pvar(paste0(stem, ".pvar")))
}

# Internal (mockable): one PLINK2 --clump.  `snps` are in greedy order (best
# first); P is their exact rank.  With `extract = TRUE` (no subset) the call
# also restricts the reference to `snps` and writes their .pvar so positions
# can be checked.  Returns list(ids = index variant IDs, pvar = data frame or
# NULL); ids = character() and pvar with zero rows when none is in the
# reference.
fastmr_clump_run_clump <- function(snps, reference_args, plink2_bin, kb_arg, r2_arg,
                                   stem, extract = FALSE) {
  clump_file <- paste0(stem, ".p.txt")
  n <- length(snps)
  writeLines(c("SNP\tP", paste0(snps, "\t", sprintf("%.17g", seq_len(n) / (n + 1)))), clump_file)
  args <- reference_args
  if (extract) {
    extract_file <- paste0(stem, ".extract.txt")
    writeLines(as.character(snps), extract_file)
    args <- c(args, "--extract", fastmr_clump_quote(extract_file), "--make-just-pvar")
  }
  # On the small candidate subset, cap PLINK2's workspace: by default every
  # one of the concurrent processes reserves half of the machine's RAM.
  if (!extract) args <- c(args, "--memory", "2048")
  args <- c(args, "--clump", fastmr_clump_quote(clump_file),
            "--clump-id-field", "SNP", "--clump-p-field", "P",
            "--clump-p1", "1", "--clump-p2", "1",
            "--clump-r2", r2_arg, "--clump-kb", kb_arg,
            "--threads", "1", "--out", fastmr_clump_quote(stem))
  output <- fastmr_clump_system2(plink2_bin, args)
  empty <- data.frame(id = character(), chr = character(), pos = numeric(), stringsAsFactors = FALSE)
  if (extract && fastmr_clump_no_variants(output)) return(list(ids = character(), pvar = empty))
  clumps <- paste0(stem, ".clumps")
  if (attr(output, "status") != 0L || !file.exists(clumps)) {
    fastmr_clump_plink_failed(output, "--clump")
  }
  lines <- readLines(clumps)
  lines <- lines[nzchar(lines)]
  if (!length(lines) || !startsWith(lines[[1L]], "#")) stop("unrecognised PLINK2 .clumps file", call. = FALSE)
  cols <- strsplit(sub("^#", "", lines[[1L]]), "\t", fixed = TRUE)[[1L]]
  id_col <- match("ID", cols)
  if (is.na(id_col)) stop("PLINK2 .clumps file has no ID column", call. = FALSE)
  fields <- strsplit(lines[-1L], "\t", fixed = TRUE)
  ids <- vapply(fields, function(f) if (length(f) >= id_col) f[[id_col]] else NA_character_, character(1))
  if (anyNA(ids)) stop("malformed PLINK2 .clumps file", call. = FALSE)
  pvar <- if (extract) fastmr_clump_read_pvar(paste0(stem, ".pvar")) else NULL
  list(ids = ids, pvar = pvar)
}

# Internal (mockable): the graph partition's LD query (same --r2-phased flags
# as fastmr_clump_run_graph()) restricted to pairs with one variant in
# `leads`.  Returns 1-based vertex ids into `snps` plus the r2 > 1 count.
fastmr_clump_run_lead_graph <- function(leads, snps, reference_args, plink2_bin, clump_kb,
                                        clump_r2, threads, workdir) {
  stem <- file.path(workdir, "cert")
  lead_file <- paste0(stem, ".leads.txt")
  writeLines(as.character(leads), lead_file)
  on.exit(unlink(c(lead_file, paste0(stem, ".vcor"))), add = TRUE)
  output <- fastmr_clump_system2(plink2_bin, c(
    reference_args, "--r2-phased", "cols=id", "--ld-snp-list", fastmr_clump_quote(lead_file),
    "--ld-window-kb", format(clump_kb, trim = TRUE, scientific = FALSE),
    "--ld-window", "1000000000", "--ld-window-r2", format(clump_r2, trim = TRUE),
    "--threads", as.integer(threads), "--out", fastmr_clump_quote(stem)))
  path <- paste0(stem, ".vcor")
  if (attr(output, "status") != 0L || !file.exists(path)) {
    fastmr_clump_plink_failed(output, "lead LD certificate query")
  }
  .fastmr_vcor_read(path, as.character(snps))
}

# Warns about a position or chromosome-label mismatch between the candidates
# and the reference (typically a different genome build): the per-exposure
# call is then delegated to the graph partition, which windows on the data
# positions, so the mismatch would otherwise go unreported.
fastmr_clump_warn_reference_mismatch <- function(reason, pvar, snp, chr, bp) {
  at <- match(snp, pvar$id)
  hit <- !is.na(at)
  if (identical(reason, "reference_position_mismatch")) {
    u <- !duplicated(snp) & hit
    n_bad <- sum(pvar$pos[at[u]] != bp[u])
    warning(sprintf(paste0("%d of %d candidate variants found in the LD reference have a different ",
                           "position there (different genome build?); per-exposure clumping was ",
                           "delegated to the graph partition, which windows on the data positions"),
                    as.integer(n_bad), as.integer(sum(u))), call. = FALSE)
  } else if (identical(reason, "reference_chromosome_mismatch")) {
    warning(paste0("candidate chromosome labels do not map one-to-one onto the LD reference's ",
                   "(different genome build or chromosome naming?); per-exposure clumping was ",
                   "delegated to the graph partition, which uses the data chromosomes"), call. = FALSE)
  }
  invisible(NULL)
}

# Checks that the reference describes the candidates exactly as the data does:
# unique IDs, equal positions and a one-to-one chromosome-label mapping.
# Returns NULL when consistent, else a reason string.
fastmr_clump_reference_mismatch <- function(pvar, snp, chr, bp) {
  if (anyDuplicated(pvar$id)) return("duplicate_reference_ids")
  at <- match(snp, pvar$id)
  hit <- !is.na(at)
  if (!any(hit)) return(NULL)
  if (any(pvar$pos[at[hit]] != bp[hit])) return("reference_position_mismatch")
  pairs <- unique(data.frame(d = chr[hit], r = pvar$chr[at[hit]], stringsAsFactors = FALSE))
  if (anyDuplicated(pairs$d) || anyDuplicated(pairs$r)) return("reference_chromosome_mismatch")
  NULL
}

#' Exact per-exposure clumping on a candidate-only reference
#'
#' Extracts the union of candidate SNPs from the reference once
#' (`--extract --make-pgen`), then runs one single-threaded PLINK2 `--clump`
#' per exposure on that subset, `threads` processes at a time.  The P column
#' given to PLINK2 is each candidate's exact rank in the greedy order used by
#' [fast_clump_data_graph()] (exact CompreSSoR rank when `pvalue_rank` is
#' present, then p, then SNP ID), and the `--clump-kb`/`--clump-r2` arguments
#' are translated so that PLINK2's window and r2 comparisons match the graph
#' partition's inclusive ones.  Candidates absent from the reference are
#' handled as `absent` says (retained by default, with a warning).  The leads are then certified with one `--r2-phased
#' --ld-snp-list` query (the graph partition's own LD statistic and flags,
#' restricted to pairs involving a lead): the greedy pass over that
#' lead-incident graph must keep exactly the `--clump` leads, which proves the
#' result equals [fast_clump_data_graph()] (and [fast_clump_data_lead_rows()]).
#' This matters because `--clump` computes r2 in a separate code path: PLINK2
#' 2.00a6.8's `--clump` mis-estimates some rare-variant pairs with |D'| = 1.
#' Pairs reported with r2 > 1 are warned about as in the graph partition.
#'
#' The call is delegated to [fast_clump_data_graph()] (reason in
#' `diagnostics$delegated`) when the certificate fails, when candidate
#' positions are incomplete or differ from the reference, when the reference
#' has duplicate candidate IDs, or when `clump_r2 <= 0`.  A position or
#' chromosome-label mismatch with the reference (for example a different
#' genome build) is also reported with a warning, since the graph partition
#' then windows on the data positions.
#'
#' @inheritParams fast_clump_data_graph
#' @param subset `"auto"` or `"always"` (default behaviour: extract the
#'   candidate subset), or `"never"` (clump and certify on the full reference
#'   with `--extract`; slower, since the reference is parsed twice).
#' @param ... Forwarded to [fast_clump_data_graph()] when the call is
#'   delegated (for example `max_graph_pairs`).
#' @return A list with `data`, named `instruments`, and `diagnostics`.
#' @export
fast_clump_data_per_exposure <- function(
    dat, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1,
    bfile = NULL, pfile = NULL, plink2_bin = NULL, threads = 1L,
    subset = c("auto", "always", "never"), workdir = NULL,
    reference_manifest = NULL, absent = c("keep", "drop"), ...) {
  if (!is.data.frame(dat) || !"SNP" %in% names(dat)) stop("dat must contain SNP", call. = FALSE)
  absent <- match.arg(absent)
  clump_kb <- fastmr_clump_number(clump_kb, "clump_kb", 0)
  clump_r2 <- fastmr_clump_number(clump_r2, "clump_r2", 0, 1)
  clump_p1 <- fastmr_clump_number(clump_p1, "clump_p1", 0, 1)
  threads <- as.integer(fastmr_clump_number(threads, "threads", 1))
  subset <- match.arg(subset)
  if (is.null(plink2_bin)) plink2_bin <- Sys.which("plink2")
  if (!nzchar(plink2_bin)) stop("PLINK2 executable not found; provide plink2_bin", call. = FALSE)
  reference_args <- fastmr_clump_reference_args(bfile, pfile)
  if (!is.null(reference_manifest) &&
      (length(reference_manifest) != 1L || !is.character(reference_manifest) ||
       is.na(reference_manifest) || !file.exists(reference_manifest))) {
    stop("reference_manifest must be one existing file when supplied", call. = FALSE)
  }
  reference_md5 <- if (is.null(reference_manifest)) NULL else unname(tools::md5sum(reference_manifest))
  delegate <- function(reason) {
    ans <- fast_clump_data_graph(original_input, clump_kb = clump_kb, clump_r2 = clump_r2,
                                 clump_p1 = clump_p1, bfile = bfile, pfile = pfile,
                                 plink2_bin = plink2_bin, threads = threads, workdir = workdir,
                                 reference_manifest = reference_manifest, absent = absent, ...)
    ans$diagnostics$delegated <- reason
    ans$diagnostics$partition_requested <- "per_exposure"
    ans
  }
  original_input <- dat
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
  elig <- which(is.finite(p) & p <= clump_p1)

  args <- fastmr_clump_per_exposure_args(clump_kb, clump_r2)
  if (is.null(args)) return(delegate("clump_r2 <= 0 has no exact --clump translation"))
  if (length(elig)) {
    if (anyNA(position$chr[elig]) || !all(is.finite(position$bp[elig]))) {
      return(delegate("missing_positions"))
    }
    first <- match(snp[elig], snp[elig])
    if (any(position$chr[elig] != position$chr[elig][first]) ||
        any(position$bp[elig] != position$bp[elig][first])) {
      return(delegate("inconsistent_candidate_positions"))
    }
  }

  workdir_owned <- is.null(workdir)
  if (workdir_owned) workdir <- tempfile("fastMR_clump_pe_")
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  if (workdir_owned) on.exit(unlink(workdir, recursive = TRUE, force = TRUE), add = TRUE)

  e_ids <- unique(expo[elig])
  # Always subset by default: without a subset the certificate query has to
  # parse the full reference a second time (~1 s for 9M variants), which costs
  # more than the subset itself even for one exposure.
  use_subset <- subset != "never"
  plink_calls <- 0L
  subset_variants <- NA_integer_
  retained_key <- character()
  n_absent <- 0L
  absent_snps <- character()
  certificate <- list(leads = 0L, lead_edges = 0L, verified = TRUE)
  if (length(elig)) {
    usnp <- unique(snp[elig])
    clump_reference <- reference_args
    present <- NULL
    if (use_subset) {
      plink_calls <- plink_calls + 1L
      sub <- fastmr_clump_make_subset(usnp, reference_args, plink2_bin, threads,
                                      file.path(workdir, "candidates"))
      if (is.null(sub)) {
        present <- character()
        subset_variants <- 0L
      } else {
        bad <- fastmr_clump_reference_mismatch(sub$pvar, snp[elig], position$chr[elig], position$bp[elig])
        if (!is.null(bad)) {
          fastmr_clump_warn_reference_mismatch(bad, sub$pvar, snp[elig], position$chr[elig], position$bp[elig])
          return(delegate(bad))
        }
        present <- sub$pvar$id
        subset_variants <- nrow(sub$pvar)
        clump_reference <- sub$reference_args
      }
    }
    ord <- fastmr_clump_order(p[elig], snp[elig], rank[elig], tiebreak[elig],
                              group = match(expo[elig], e_ids))
    rows <- elig[ord]
    groups <- split(rows, factor(expo[rows], levels = e_ids))
    # Reference membership once per row (not once per exposure against the
    # whole subset: E x |subset| hashing dominated at E = 1000).
    present_row <- if (is.null(present)) NULL else {
      x <- logical(length(snp)); x[elig] <- snp[elig] %in% present; x
    }
    jobs <- lapply(seq_along(groups), function(k) {
      r <- groups[[k]]
      if (!is.null(present_row)) r <- r[present_row[r]]
      list(k = k, snps = snp[r])
    })
    run_job <- function(job) {
      if (!length(job$snps)) return(list(ids = character(), pvar = NULL))
      tryCatch(fastmr_clump_run_clump(job$snps, clump_reference, plink2_bin, args$kb_arg, args$r2_arg,
                                      file.path(workdir, paste0("e", job$k)), extract = !use_subset),
               error = function(e) e)
    }
    todo <- which(vapply(jobs, function(j) length(j$snps) > 0L, logical(1)))
    workers <- min(threads, max(1L, length(todo)))
    out <- vector("list", length(jobs))
    if (length(todo)) {
      res <- if (.Platform$OS.type != "windows" && workers > 1L) {
        parallel::mclapply(jobs[todo], run_job, mc.cores = workers, mc.preschedule = TRUE)
      } else {
        lapply(jobs[todo], run_job)
      }
      out[todo] <- res
    }
    plink_calls <- plink_calls + length(todo)
    failed <- vapply(out[todo], function(x) inherits(x, "error") || !is.list(x) || is.null(x$ids), logical(1))
    if (any(failed)) {
      bad <- out[todo][[which(failed)[1L]]]
      stop(if (inherits(bad, "error")) conditionMessage(bad) else "a per-exposure --clump worker failed",
           call. = FALSE)
    }
    kept <- vector("list", length(groups))
    ordered_present <- vector("list", length(groups))
    for (k in seq_along(groups)) {
      r <- groups[[k]]
      s <- snp[r]
      in_ref <- if (use_subset) present_row[r] else {
        pv <- out[[k]]$pvar
        if (!is.null(pv) && nrow(pv)) {
          bad <- fastmr_clump_reference_mismatch(pv, s, position$chr[r], position$bp[r])
          if (!is.null(bad)) {
            fastmr_clump_warn_reference_mismatch(bad, pv, s, position$chr[r], position$bp[r])
            return(delegate(bad))
          }
        }
        if (is.null(pv)) rep(FALSE, length(s)) else s %in% pv$id
      }
      unknown <- setdiff(out[[k]]$ids, s[in_ref])
      if (length(unknown)) stop("PLINK2 --clump reported an index variant that was not a candidate: ",
                                unknown[[1L]], call. = FALSE)
      n_absent <- n_absent + sum(!in_ref)
      absent_snps <- c(absent_snps, s[!in_ref])
      keep <- (if (identical(absent, "drop")) FALSE else !in_ref) | s %in% out[[k]]$ids
      kept[[k]] <- paste(e_ids[[k]], s[keep], sep = "\r")
      ordered_present[[k]] <- s[in_ref]
    }
    retained_key <- unlist(kept, use.names = FALSE)
    fastmr_clump_check_absent(n_absent, length(elig), length(unique(absent_snps)), absent)
    # Certificate: the graph partition's own LD statistic (--r2-phased with
    # the graph's flags) for every pair involving a --clump lead.  The greedy
    # pass only ever consults edges to kept SNPs, so if the C++ greedy over
    # this lead-incident graph keeps exactly the --clump leads in every
    # exposure, the result equals the all-pairs graph result.  --clump and
    # --r2-phased are separate PLINK2 code paths (2.00a6.8's --clump
    # mis-estimates r2 for some |D'| = 1 rare-variant pairs that --r2-phased
    # gets right), so this is checked rather than assumed.
    leads <- unique(unlist(lapply(out, function(x) x$ids), use.names = FALSE))
    if (length(leads)) {
      vertices <- unique(unlist(ordered_present, use.names = FALSE))
      cert_reference <- if (use_subset) clump_reference else
        c(clump_reference, "--extract", fastmr_clump_quote(file.path(workdir, "cert.extract.txt")))
      if (!use_subset) writeLines(vertices, file.path(workdir, "cert.extract.txt"))
      plink_calls <- plink_calls + 1L
      ld <- fastmr_clump_run_lead_graph(leads, vertices, cert_reference, plink2_bin, clump_kb, clump_r2,
                                        threads, workdir)
      vbp <- position$bp[elig][match(vertices, snp[elig])]
      a <- ld$lead; b <- ld$target
      ok <- a != b & abs(vbp[a] - vbp[b]) <= clump_kb * 1000
      a <- a[ok]; b <- b[ok]
      vtx <- match(unlist(ordered_present, use.names = FALSE), vertices)
      starts <- c(0L, cumsum(lengths(ordered_present)))
      gkeep <- .fastmr_graph_clump(length(vertices), a - 1L, b - 1L, vtx - 1L, as.integer(starts))
      gk <- split(gkeep, rep(seq_along(ordered_present), lengths(ordered_present)))
      for (k in which(lengths(ordered_present) > 0L)) {
        g_leads <- ordered_present[[k]][gk[[as.character(k)]]]
        if (!setequal(g_leads, out[[k]]$ids)) {
          return(delegate(sprintf("per_exposure_certificate_mismatch (exposure %s)", e_ids[[k]])))
        }
      }
      certificate <- list(leads = length(leads), lead_edges = length(a), verified = TRUE)
      fastmr_clump_warn_invalid_r2(ld$invalid_r2, ld$invalid_example)
    }
  }
  version <- fastmr_clump_plink_version(plink2_bin)
  original_key <- paste(as.character(original$id.exposure), as.character(original$SNP), sep = "\r")
  result <- original[original_key %in% retained_key, , drop = FALSE]
  instruments <- lapply(split(result$SNP, result$id.exposure, drop = TRUE), as.character)
  ld_provenance <- list(
    plink2_version = version,
    mode = if (use_subset) "per_exposure_clump_candidate_subset" else "per_exposure_clump",
    flags = c("--clump (P = exact greedy rank)", "--clump-p1 1", "--clump-p2 1",
              paste("--clump-r2", args$r2_arg), paste("--clump-kb", args$kb_arg)),
    window_bp_inclusive = args$radius_bp, r2_inclusive_cut = args$r2_cut,
    reference = paste(gsub("'", "", reference_args), collapse = " "),
    reference_manifest_md5 = reference_md5,
    clump_kb = clump_kb, clump_r2 = clump_r2, clump_p1 = clump_p1
  )
  list(data = result, instruments = instruments,
       diagnostics = list(
         exposures = length(unique(as.character(original$id.exposure))),
         candidate_rows = nrow(dat), retained = length(retained_key),
         rounds = 0L, plink_calls = plink_calls, subset = use_subset,
         subset_variants = subset_variants, absent_from_reference = n_absent, absent = absent,
         certificate = certificate,
         workers = if (length(elig)) min(threads, max(1L, length(e_ids))) else 0L,
         exact = TRUE, fallback = FALSE, fallbacks = list(),
         partition = "per_exposure",
         reference_manifest_md5 = reference_md5, ld_provenance = ld_provenance))
}

# Measured cost model (BluePebble Cascade Lake nodes, 9M-variant 1000G EUR
# reference, clump-suite E = 1..1000 and the 10 x 10 showcase; seconds).
fastmr_clump_auto_model <- list(
  graph_fixed = 1.1,          # PLINK2 start + reference .pvar parse
  graph_per_pair = 0.75e-6,   # all-pairs LD + .vcor parse + R, per estimated pair
  graph_per_row = 2e-6,       # per eligible candidate row (R bookkeeping)
  per_exposure_fixed = 1.0,   # candidate subset + lead certificate query
  per_exposure_each = 0.0155, # one single-threaded --clump on the subset
  worker_divisor = 3.5        # concurrent --clump processes scale ~ threads / 3.5
)

# Chooses "graph" or "per_exposure" from the estimated candidate pair count.
fastmr_clump_auto_plan <- function(dat, clump_kb, clump_p1, threads,
                                   model = fastmr_clump_auto_model) {
  pcol <- if ("pval.exposure" %in% names(dat)) "pval.exposure" else if ("pval.outcome" %in% names(dat)) "pval.outcome" else NULL
  # Estimate only: skip as.character() on numeric p (0.4 s per 300k rows).
  p <- if (is.null(pcol)) rep(0.99, nrow(dat)) else if (is.numeric(dat[[pcol]])) dat[[pcol]]
       else suppressWarnings(as.numeric(as.character(dat[[pcol]])))
  elig <- which(is.finite(p) & p <= clump_p1)
  expo <- if ("id.exposure" %in% names(dat)) as.character(dat$id.exposure) else rep("exposure", nrow(dat))
  E <- length(unique(expo[elig]))
  # Estimate only: avoid fastmr_clump_position()'s as.character() round trip.
  position <- list(
    chr = if ("chr_name" %in% names(dat)) as.character(dat$chr_name) else rep(NA_character_, nrow(dat)),
    bp = if (!"chrom_start" %in% names(dat)) rep(NA_real_, nrow(dat))
         else if (is.numeric(dat$chrom_start)) as.numeric(dat$chrom_start)
         else suppressWarnings(as.numeric(as.character(dat$chrom_start))))
  plan <- list(exposures = E, threads = threads, estimated_pairs = NA_real_,
               effective_workers = NA_real_, cost_graph = NA_real_,
               cost_per_exposure = NA_real_, strategy = "graph", reason = "")
  if (!length(elig)) {
    plan$reason <- "no eligible candidates"
    return(plan)
  }
  if (anyNA(position$chr[elig]) || !all(is.finite(position$bp[elig]))) {
    plan$reason <- "missing positions (graph falls back to lead_row)"
    return(plan)
  }
  snp <- as.character(dat$SNP[elig])
  u <- !duplicated(snp)
  est <- sum(vapply(split(position$bp[elig][u], position$chr[elig][u]),
                    fastmr_graph_pair_estimate, numeric(1), window = clump_kb * 1000))
  # Concurrent single-threaded PLINK2 processes scaled at about threads / 3.5
  # on 8-thread (4-core, hyper-threaded) allocations.
  w <- max(1, min(threads, E) / model$worker_divisor)
  plan$estimated_pairs <- est
  plan$effective_workers <- w
  plan$cost_graph <- model$graph_fixed + model$graph_per_pair * est + model$graph_per_row * length(elig)
  plan$cost_per_exposure <- model$per_exposure_fixed + model$per_exposure_each * E / w
  if (plan$cost_per_exposure < plan$cost_graph) {
    plan$strategy <- "per_exposure"
    plan$reason <- sprintf("estimated pairs %.3g: per-exposure cost %.2f s < graph %.2f s",
                           est, plan$cost_per_exposure, plan$cost_graph)
  } else {
    plan$reason <- sprintf("estimated pairs %.3g: graph cost %.2f s <= per-exposure %.2f s",
                           est, plan$cost_graph, plan$cost_per_exposure)
  }
  plan
}

#' Exact multi-exposure clumping with size-based strategy dispatch
#'
#' Chooses between [fast_clump_data_graph()] (one all-pairs LD call over the
#' candidate union) and [fast_clump_data_per_exposure()] (one PLINK2 `--clump`
#' per exposure on a candidate-only reference) from a cost model measured on
#' a 9M-variant reference: graph about `1.1 s + 0.75 us x P + 2 us x rows`,
#' with `P` the estimated candidate pairs within `clump_kb` (from sorted
#' positions) and `rows` the eligible candidate rows; per-exposure about
#' `1.0 s + 15.5 ms x E / w`, with `w = max(1, min(threads, E) / 3.5)`
#' concurrent workers.  The cheaper prediction wins.  Both strategies return
#' identical instruments; the choice and the model inputs
#' are recorded in `diagnostics$auto`.  If the per-exposure run fails, the
#' graph strategy is used and the error is recorded in
#' `diagnostics$auto$per_exposure_error`.
#'
#' @inheritParams fast_clump_data_graph
#' @param ... Forwarded to the chosen strategy (for example
#'   `max_graph_pairs`).
#' @return A list with `data`, named `instruments`, and `diagnostics`.
#' @export
fast_clump_data_auto <- function(
    dat, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1,
    bfile = NULL, pfile = NULL, plink2_bin = NULL, threads = 1L,
    workdir = NULL, reference_manifest = NULL, ...) {
  if (!is.data.frame(dat) || !"SNP" %in% names(dat)) stop("dat must contain SNP", call. = FALSE)
  clump_kb <- fastmr_clump_number(clump_kb, "clump_kb", 0)
  clump_p1 <- fastmr_clump_number(clump_p1, "clump_p1", 0, 1)
  threads <- as.integer(fastmr_clump_number(threads, "threads", 1))
  plan <- fastmr_clump_auto_plan(dat, clump_kb, clump_p1, threads)
  res <- NULL
  if (identical(plan$strategy, "per_exposure")) {
    res <- tryCatch(fast_clump_data_per_exposure(
      dat, clump_kb = clump_kb, clump_r2 = clump_r2, clump_p1 = clump_p1, bfile = bfile,
      pfile = pfile, plink2_bin = plink2_bin, threads = threads, workdir = workdir,
      reference_manifest = reference_manifest, ...), error = function(e) e)
    if (inherits(res, "error")) {
      plan$per_exposure_error <- conditionMessage(res)
      plan$strategy <- "graph"
      res <- NULL
    }
  }
  if (is.null(res)) {
    graph <- function(..., subset) fast_clump_data_graph(...)   # drop a per-exposure-only argument
    res <- graph(dat, clump_kb = clump_kb, clump_r2 = clump_r2, clump_p1 = clump_p1, bfile = bfile,
                 pfile = pfile, plink2_bin = plink2_bin, threads = threads, workdir = workdir,
                 reference_manifest = reference_manifest, ...)
  }
  plan$used <- res$diagnostics$partition
  res$diagnostics$auto <- plan
  res$diagnostics$partition_requested <- "auto"
  res
}
