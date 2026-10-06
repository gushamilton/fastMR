# Batched, exact greedy LD clumping for many exposure-specific candidate sets.
#
# The implementation deliberately keeps the exposure state in R but asks
# PLINK2 for LD rows for all current leads in one call per frontier round.  A
# pair is never requested twice: negative (tested below threshold) and
# positive (reported by PLINK) decisions are cached separately.  This is
# equivalent to running PLINK clumping independently for every exposure, but
# shares genotype/LD work whenever exposures have the same lead or target.

fastmr_clump_number <- function(x, name, lower = -Inf, upper = Inf) {
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || !is.finite(x) ||
      x < lower || x > upper) {
    stop(name, " must be one finite value in [", lower, ", ", upper, "]", call. = FALSE)
  }
  as.numeric(x)
}

fastmr_clump_default <- function(x, value) if (is.null(x)) value else x

fastmr_clump_position <- function(dat) {
  chr <- if ("chr_name" %in% names(dat)) as.character(dat$chr_name) else rep(NA_character_, nrow(dat))
  bp <- if ("chrom_start" %in% names(dat)) suppressWarnings(as.numeric(as.character(dat$chrom_start))) else rep(NA_real_, nrow(dat))
  list(chr = chr, bp = bp)
}

# Greedy order: exact CompreSSoR rank when the optional `pvalue_rank` column
# is present, then p, then the tie-break (-|z|, see fastmr_clump_tiebreak())
# when available, then SNP with C-locale ties.  `group`, when given, is the
# leading sort key (exposure).  Without ties in p the order is (p, SNP).
fastmr_clump_order <- function(p, snp, rank = NULL, tiebreak = NULL, group = NULL) {
  keys <- c(if (!is.null(group)) list(group), if (!is.null(rank)) list(rank), list(p),
            if (!is.null(tiebreak)) list(tiebreak), list(snp))
  do.call(order, c(keys, list(method = "radix")))
}

# Tie-break for equal p: -|z| (strongest first; NA, sorted last, when
# unknown).  p underflows to 0 for |z| above ~38 (CompreSSoR reconstructed p
# and many cis-pQTL files), and (p, SNP) would then lead with the
# lexicographically first SNP.  |z| comes from the internal `.fastmr_abs_z`
# column (compressed candidates), else beta/se on the p-value's side.
# NULL when unavailable, which keeps the (p, SNP) order.
fastmr_clump_tiebreak <- function(dat, pcol) {
  num <- function(x) if (is.numeric(x)) x else suppressWarnings(as.numeric(as.character(x)))
  if (".fastmr_abs_z" %in% names(dat)) {
    z <- num(dat[[".fastmr_abs_z"]])
  } else {
    side <- if (identical(pcol, "pval.outcome")) "outcome" else "exposure"
    b <- paste0("beta.", side)
    s <- paste0("se.", side)
    if (!all(c(b, s) %in% names(dat))) return(NULL)
    z <- num(dat[[b]]) / num(dat[[s]])
  }
  z <- -abs(z)
  z[!is.finite(z)] <- NA_real_
  z
}

# Rows kept when an (id.exposure, SNP) pair repeats: the one with the
# smallest p (the first such row on ties; rows without a usable p last).
# Logical over the rows of `dat`; all TRUE when nothing repeats.
fastmr_clump_dedup <- function(dat, pcol) {
  key <- paste(as.character(dat$id.exposure), as.character(dat$SNP), sep = "\r")
  keep <- !duplicated(key)
  if (all(keep)) return(keep)
  p <- suppressWarnings(as.numeric(as.character(dat[[pcol]])))
  p[!is.finite(p)] <- NA_real_
  code <- match(key, key)
  o <- order(code, p, seq_along(key), method = "radix")
  keep <- logical(length(key))
  keep[o[!duplicated(code[o])]] <- TRUE
  keep
}

# Candidates absent from the LD reference have no LD edges: kept unclumped
# (absent = "keep", the historical behaviour) or dropped (absent = "drop", as
# TwoSampleMR / ieugwasr local clumping does).  Errors when every eligible
# candidate is absent (almost always a SNP-ID scheme mismatch), otherwise
# warns with the count.
fastmr_clump_check_absent <- function(n_absent, n_rows, absent_variants, absent) {
  if (!n_rows || !n_absent) return(invisible(NULL))
  hint <- paste0("check that SNP IDs follow the reference's ID scheme ",
                 "(rsID vs chr:pos:ref:alt, allele order, 'chr' prefix)")
  if (n_absent == n_rows) {
    stop("none of the ", n_rows, " eligible candidate rows (", absent_variants,
         " variants) is in the LD reference; ", hint, call. = FALSE)
  }
  warning(sprintf(paste0("%d of %d eligible candidate rows (%d variants) are absent from the LD ",
                         "reference and were %s; %s"),
                  as.integer(n_absent), as.integer(n_rows), as.integer(absent_variants),
                  if (identical(absent, "drop")) "dropped (absent = \"drop\")"
                  else "kept unclumped (absent = \"keep\"; use absent = \"drop\" to drop them)",
                  hint), call. = FALSE)
  invisible(NULL)
}

fastmr_clump_pair_key <- function(a, b) {
  ifelse(a < b, paste0(a, "\r", b), paste0(b, "\r", a))
}

fastmr_clump_read_vcor <- function(path, zstdcat = NULL) {
  if (!file.exists(path)) return(data.frame(lead = character(), target = character(), stringsAsFactors = FALSE))
  if (is.null(zstdcat)) zstdcat <- Sys.which("zstdcat")
  if (!nzchar(zstdcat)) {
    zstdcat <- Sys.which("zstd")
    if (!nzchar(zstdcat)) stop("PLINK produced a .vcor.zst file but neither zstdcat nor zstd is available", call. = FALSE)
    args <- c("-dc", shQuote(path))
  } else {
    args <- shQuote(path)
  }
  lines <- tryCatch(system2(zstdcat, args, stdout = TRUE, stderr = TRUE),
                    error = function(e) stop("could not read PLINK LD output: ", conditionMessage(e), call. = FALSE))
  status <- attr(lines, "status")
  if (!is.null(status) && status != 0L) stop("could not decompress PLINK LD output", call. = FALSE)
  lines <- lines[nzchar(lines) & !grepl("^#", lines)]
  if (!length(lines)) return(data.frame(lead = character(), target = character(), stringsAsFactors = FALSE))
  fields <- strsplit(lines, "\t", fixed = TRUE)
  fields <- fields[vapply(fields, length, integer(1)) >= 6L]
  if (!length(fields)) return(data.frame(lead = character(), target = character(), stringsAsFactors = FALSE))
  r2 <- suppressWarnings(as.numeric(vapply(fields, function(f) if (length(f) >= 7L) f[[7L]] else NA_character_, character(1))))
  bad <- which(r2 > 1 + 1e-6)
  if (length(bad)) {
    fastmr_clump_warn_invalid_r2(length(bad), paste(fields[[bad[1L]]][c(3L, 6L, 7L)], collapse = " "))
  }
  data.frame(lead = vapply(fields, `[[`, character(1), 3L),
             target = vapply(fields, `[[`, character(1), 6L),
             stringsAsFactors = FALSE)
}

fastmr_clump_quote <- function(x) {
  shQuote(as.character(x), type = if (.Platform$OS.type == "windows") "cmd" else "sh")
}

fastmr_clump_reference_args <- function(bfile = NULL, pfile = NULL) {
  if (!is.null(bfile) && !is.null(pfile)) stop("supply only one of bfile or pfile", call. = FALSE)
  if (is.null(bfile) && is.null(pfile)) stop("supply a PLINK bfile or pfile reference", call. = FALSE)
  if (!is.null(bfile)) c("--bfile", fastmr_clump_quote(bfile)) else c("--pfile", fastmr_clump_quote(pfile))
}

fastmr_clump_run_frontier <- function(leads, targets, reference_args, plink2_bin,
                                      clump_kb, clump_r2, threads, workdir, round) {
  stem <- file.path(workdir, sprintf("frontier_%06d", round))
  lead_file <- paste0(stem, ".leads.txt")
  target_file <- paste0(stem, ".targets.txt")
  writeLines(unique(as.character(leads)), lead_file)
  writeLines(unique(as.character(targets)), target_file)
  args <- c(reference_args, "--extract", fastmr_clump_quote(target_file),
            "--ld-snp-list", fastmr_clump_quote(lead_file),
            "--r2-phased", "zs", "--ld-window-kb", format(clump_kb, trim = TRUE),
            "--ld-window-r2", format(clump_r2, trim = TRUE), "--threads", as.integer(threads),
            "--out", fastmr_clump_quote(stem))
  output <- tryCatch(suppressWarnings(system2(plink2_bin, args, stdout = TRUE, stderr = TRUE)),
                     error = function(e) structure(character(), status = 1L, error = conditionMessage(e)))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  path <- paste0(stem, ".vcor.zst")
  if (status != 0L || !file.exists(path)) {
    detail <- attr(output, "error")
    if (is.null(detail)) detail <- paste(utils::tail(output, 8L), collapse = " | ")
    stop("batched PLINK2 LD query failed in round ", round, ": ", detail, call. = FALSE)
  }
  fastmr_clump_read_vcor(path)
}

#' Batched multi-exposure PLINK2 LD clumping
#'
#' Clump each exposure independently while sharing PLINK2 LD queries across
#' all exposures.  The function is an opt-in replacement for repeatedly
#' calling [fast_clump_data()] with a PLINK reference.  It retains the exact
#' greedy index-SNP decisions of the per-exposure workflow and returns the
#' original rows for retained SNPs.
#'
#' @param dat Data frame containing `SNP`, `id.exposure`, and a p-value column.
#' @param clump_kb Maximum index/target distance in kilobases.
#' @param clump_r2 Minimum LD r-squared for removing a target.
#' @param clump_p1 Maximum p-value for an index SNP.
#' @param bfile PLINK binary reference prefix, or use `pfile`.
#' @param pfile PLINK2 pgen reference prefix, or use `bfile`.
#' @param plink2_bin PLINK2 executable.  Defaults to `plink2` on `PATH`.
#' @param threads Threads passed to each PLINK2 query.
#' @param max_pair_requests Safety limit on the number of logical LD pairs.
#' @param max_target_variants Safety limit on a single frontier target union.
#' @param max_rounds Safety limit on frontier rounds.
#' @param on_limit Either `"error"` (default) or `"fallback"`; fallback uses
#'   the existing per-exposure PLINK workflow and requires `bfile`.
#' @param workdir Optional directory for temporary frontier files.
#' @param reference_manifest Optional reference-panel manifest whose MD5 is
#'   recorded in diagnostics.
#' @return A list with `data`, `instruments`, and `diagnostics`.
#' @export
fast_clump_data_batched <- function(
    dat, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1,
    bfile = NULL, pfile = NULL, plink2_bin = NULL, threads = 1L,
    max_pair_requests = 2e8, max_target_variants = 2e6, max_rounds = 10000L,
    on_limit = c("error", "fallback"), workdir = NULL,
    reference_manifest = NULL) {
  if (!is.data.frame(dat) || !"SNP" %in% names(dat)) stop("dat must contain SNP", call. = FALSE)
  clump_kb <- fastmr_clump_number(clump_kb, "clump_kb", 0)
  clump_r2 <- fastmr_clump_number(clump_r2, "clump_r2", 0, 1)
  clump_p1 <- fastmr_clump_number(clump_p1, "clump_p1", 0, 1)
  threads <- as.integer(fastmr_clump_number(threads, "threads", 1))
  max_pair_requests <- fastmr_clump_number(max_pair_requests, "max_pair_requests", 1)
  max_target_variants <- as.integer(fastmr_clump_number(max_target_variants, "max_target_variants", 1))
  max_rounds <- as.integer(fastmr_clump_number(max_rounds, "max_rounds", 1))
  on_limit <- match.arg(on_limit)
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
  position <- fastmr_clump_position(dat)
  exposure_ids <- unique(as.character(dat$id.exposure))
  states <- lapply(exposure_ids, function(id) {
    ii <- which(as.character(dat$id.exposure) == id)
    ii <- ii[is.finite(p[ii]) & p[ii] <= clump_p1]
    ii <- ii[fastmr_clump_order(p[ii], as.character(dat$SNP[ii]), dat[["pvalue_rank"]][ii],
                                tiebreak[ii])]
    list(index = ii, dead = rep(FALSE, length(ii)))
  })
  names(states) <- exposure_ids
  retained <- logical(nrow(dat))
  pair_tested <- new.env(hash = TRUE, parent = emptyenv())
  pair_positive <- new.env(hash = TRUE, parent = emptyenv())
  n_pairs <- 0
  n_rounds <- 0L
  n_calls <- 0L
  workdir_owned <- is.null(workdir)
  if (workdir_owned) workdir <- tempfile("fastMR_batched_")
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  cleanup <- if (workdir_owned) on.exit(unlink(workdir, recursive = TRUE, force = TRUE), add = TRUE) else NULL
  current_live <- function(state) {
    if (!length(state$index)) return(NA_integer_)
    hit <- which(!state$dead)
    if (!length(hit)) NA_integer_ else state$index[hit[1L]]
  }
  limit <- function(message) {
    if (on_limit == "error") stop(message, call. = FALSE)
    if (is.null(bfile)) stop(message, "; fallback requires bfile", call. = FALSE)
    warning(message, "; falling back to per-exposure PLINK clumping", call. = FALSE)
    result <- fast_clump_data(original, clump_kb = clump_kb, clump_r2 = clump_r2,
                              clump_p1 = clump_p1, bfile = bfile)
    instruments <- lapply(split(result$SNP, result$id.exposure, drop = TRUE), as.character)
    return(list(data = result, instruments = instruments,
                diagnostics = list(exposures = length(exposure_ids), rounds = n_rounds,
                                   plink_calls = n_calls, logical_pairs = n_pairs,
                                   exact = TRUE, fallback = TRUE,
                                   reference_manifest_md5 = reference_md5)))
  }
  repeat {
    leads <- vapply(states, current_live, integer(1))
    if (all(is.na(leads))) break
    n_rounds <- n_rounds + 1L
    if (n_rounds > max_rounds) return(limit("batched clumping exceeded max_rounds"))
    lead_groups <- split(seq_along(leads)[!is.na(leads)], dat$SNP[leads[!is.na(leads)]])
    lead_names <- names(lead_groups)
    target_ids <- character()
    for (lead in lead_names) {
      target <- character()
      for (exposure in lead_groups[[lead]]) {
        state <- states[[exposure]]
        live <- state$index[!state$dead]
        if (!length(live)) next
        same <- is.na(position$chr[leads[exposure]]) | is.na(position$chr[live]) |
          position$chr[leads[exposure]] == position$chr[live]
        close <- is.na(position$bp[leads[exposure]]) | is.na(position$bp[live]) |
          abs(position$bp[leads[exposure]] - position$bp[live]) <= clump_kb * 1000
        target <- c(target, as.character(dat$SNP[live[same & close]]))
      }
      target_ids <- c(target_ids, target)
    }
    target_ids <- unique(target_ids)
    if (length(target_ids) > max_target_variants) return(limit("batched clumping target union exceeded max_target_variants"))
    # Each lead is compared with the union of all relevant live targets.  The
    # cache means only previously unseen logical pairs count towards the cap.
    pair_keys <- unlist(lapply(lead_names, function(lead) fastmr_clump_pair_key(lead, target_ids)), use.names = FALSE)
    new_pairs <- pair_keys[!vapply(pair_keys, exists, logical(1), envir = pair_tested, inherits = FALSE)]
    if (n_pairs + length(new_pairs) > max_pair_requests) return(limit("batched clumping pair-request limit exceeded"))
    n_pairs <- n_pairs + length(new_pairs)
    if (length(new_pairs)) for (key in new_pairs) assign(key, TRUE, envir = pair_tested)
    ld <- fastmr_clump_run_frontier(lead_names, target_ids, reference_args, plink2_bin,
                                    clump_kb, clump_r2, threads, workdir, n_rounds)
    n_calls <- n_calls + 1L
    if (nrow(ld)) {
      for (j in seq_len(nrow(ld))) assign(fastmr_clump_pair_key(ld$lead[j], ld$target[j]), TRUE, envir = pair_positive)
    }
    for (exposure in seq_along(states)) {
      lead_index <- leads[exposure]
      if (is.na(lead_index)) next
      lead <- as.character(dat$SNP[lead_index])
      retained[lead_index] <- TRUE
      state <- states[[exposure]]
      state$dead[match(lead_index, state$index)] <- TRUE
      live <- which(!state$dead)
      if (length(live)) {
        candidate <- state$index[live]
        same <- is.na(position$chr[lead_index]) | is.na(position$chr[candidate]) |
          position$chr[lead_index] == position$chr[candidate]
        close <- is.na(position$bp[lead_index]) | is.na(position$bp[candidate]) |
          abs(position$bp[lead_index] - position$bp[candidate]) <= clump_kb * 1000
        candidate <- candidate[same & close]
        if (length(candidate)) {
          keys <- fastmr_clump_pair_key(lead, as.character(dat$SNP[candidate]))
          blocked <- vapply(keys, exists, logical(1), envir = pair_positive, inherits = FALSE)
          state$dead[match(candidate, state$index)] <- blocked
        }
      }
      states[[exposure]] <- state
    }
  }
  retained_key <- paste(as.character(dat$id.exposure[retained]), as.character(dat$SNP[retained]), sep = "\r")
  original_key <- paste(as.character(original$id.exposure), as.character(original$SNP), sep = "\r")
  result <- original[original_key %in% retained_key, , drop = FALSE]
  instruments <- lapply(split(result$SNP, result$id.exposure, drop = TRUE), as.character)
  list(data = result, instruments = instruments,
       diagnostics = list(exposures = length(exposure_ids), candidate_rows = nrow(dat),
                          retained = nrow(dat[retained, , drop = FALSE]), rounds = n_rounds,
                          plink_calls = n_calls, logical_pairs = n_pairs,
                          positive_pairs = length(ls(pair_positive)), exact = TRUE,
                          fallback = FALSE, reference_manifest_md5 = reference_md5))
}

#' Exact lead-row LD clumping with a shared pair cache
#'
#' This variant asks PLINK for one LD row whenever a SNP first becomes the
#' current lead of one or more exposures.  Exposures sharing that lead share
#' the row, and the symmetric pair cache prevents an LD decision from being
#' requested again when the same pair is encountered later.  Unlike the global
#' frontier implementation, each query contains only targets relevant to that
#' lead group; this is the direct implementation of the lead-row strategy.
#'
#' @param dat Data frame containing `SNP`, `id.exposure`, and a p-value column.
#' @param clump_kb Maximum index/target distance in kilobases.
#' @param clump_r2 Minimum LD r-squared for removing a target.
#' @param clump_p1 Maximum p-value for an index SNP.
#' @param bfile PLINK binary reference prefix, or use `pfile`.
#' @param pfile PLINK2 pgen reference prefix, or use `bfile`.
#' @param plink2_bin PLINK2 executable. Defaults to `plink2` on `PATH`.
#' @param threads Threads passed to each PLINK query.
#' @param max_pair_requests Safety limit on logical LD pairs.
#' @param max_target_variants Safety limit for one lead-row target set.
#' @param max_rounds Safety limit on frontier rounds.
#' @param on_limit Either `"error"` (default) or `"fallback"`.
#' @param workdir Optional directory for query files.
#' @param reference_manifest Optional reference-panel manifest whose MD5 is
#'   recorded in diagnostics.
#' @return A list with `data`, named `instruments`, and `diagnostics`.
#' @export
fast_clump_data_lead_rows <- function(
    dat, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1,
    bfile = NULL, pfile = NULL, plink2_bin = NULL, threads = 1L,
    max_pair_requests = 2e8, max_target_variants = 2e6, max_rounds = 10000L,
    on_limit = c("error", "fallback"), workdir = NULL,
    reference_manifest = NULL) {
  if (!is.data.frame(dat) || !"SNP" %in% names(dat)) stop("dat must contain SNP", call. = FALSE)
  clump_kb <- fastmr_clump_number(clump_kb, "clump_kb", 0)
  clump_r2 <- fastmr_clump_number(clump_r2, "clump_r2", 0, 1)
  clump_p1 <- fastmr_clump_number(clump_p1, "clump_p1", 0, 1)
  threads <- as.integer(fastmr_clump_number(threads, "threads", 1))
  max_pair_requests <- fastmr_clump_number(max_pair_requests, "max_pair_requests", 1)
  max_target_variants <- as.integer(fastmr_clump_number(max_target_variants, "max_target_variants", 1))
  max_rounds <- as.integer(fastmr_clump_number(max_rounds, "max_rounds", 1))
  on_limit <- match.arg(on_limit)
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
  position <- fastmr_clump_position(dat)
  exposure_ids <- unique(as.character(dat$id.exposure))
  states <- lapply(exposure_ids, function(id) {
    ii <- which(as.character(dat$id.exposure) == id)
    ii <- ii[is.finite(p[ii]) & p[ii] <= clump_p1]
    ii <- ii[fastmr_clump_order(p[ii], as.character(dat$SNP[ii]), dat[["pvalue_rank"]][ii],
                                tiebreak[ii])]
    list(index = ii, dead = rep(FALSE, length(ii)))
  })
  names(states) <- exposure_ids
  retained <- logical(nrow(dat))
  pair_tested <- new.env(hash = TRUE, parent = emptyenv())
  pair_positive <- new.env(hash = TRUE, parent = emptyenv())
  n_pairs <- 0
  n_rounds <- 0L
  n_calls <- 0L
  n_unique_leads <- 0L
  workdir_owned <- is.null(workdir)
  if (workdir_owned) workdir <- tempfile("fastMR_lead_rows_")
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  if (workdir_owned) on.exit(unlink(workdir, recursive = TRUE, force = TRUE), add = TRUE)
  current_live <- function(state) {
    if (!length(state$index)) return(NA_integer_)
    hit <- which(!state$dead)
    if (!length(hit)) NA_integer_ else state$index[hit[1L]]
  }
  limit <- function(message) {
    if (on_limit == "error") stop(message, call. = FALSE)
    if (is.null(bfile)) stop(message, "; fallback requires bfile", call. = FALSE)
    warning(message, "; falling back to per-exposure PLINK clumping", call. = FALSE)
    result <- fast_clump_data(original, clump_kb = clump_kb, clump_r2 = clump_r2,
                              clump_p1 = clump_p1, bfile = bfile)
    instruments <- lapply(split(result$SNP, result$id.exposure, drop = TRUE), as.character)
    list(data = result, instruments = instruments,
         diagnostics = list(exposures = length(exposure_ids), rounds = n_rounds,
                            plink_calls = n_calls, logical_pairs = n_pairs,
                            unique_leads = n_unique_leads, exact = TRUE,
                            fallback = TRUE, strategy = "lead_row",
                            reference_manifest_md5 = reference_md5))
  }
  repeat {
    leads <- vapply(states, current_live, integer(1))
    if (all(is.na(leads))) break
    n_rounds <- n_rounds + 1L
    if (n_rounds > max_rounds) return(limit("lead-row clumping exceeded max_rounds"))
    live_exposures <- which(!is.na(leads))
    lead_groups <- split(live_exposures, as.character(dat$SNP[leads[live_exposures]]))
    for (lead in names(lead_groups)) {
      group <- lead_groups[[lead]]
      target_parts <- lapply(group, function(exposure) {
        state <- states[[exposure]]
        live <- state$index[!state$dead]
        if (!length(live)) return(character())
        lead_index <- leads[exposure]
        same <- is.na(position$chr[lead_index]) | is.na(position$chr[live]) |
          position$chr[lead_index] == position$chr[live]
        close <- is.na(position$bp[lead_index]) | is.na(position$bp[live]) |
          abs(position$bp[lead_index] - position$bp[live]) <= clump_kb * 1000
        as.character(dat$SNP[live[same & close]])
      })
      targets <- unique(unlist(target_parts, use.names = FALSE))
      if (!length(targets)) next
      keys <- fastmr_clump_pair_key(lead, targets)
      missing <- !vapply(keys, exists, logical(1), envir = pair_tested, inherits = FALSE)
      targets <- targets[missing]
      keys <- keys[missing]
      if (!length(targets)) next
      # PLINK applies --extract before --ld-snp-list.  The current lead may
      # already have its self-pair in the tested cache, but it must still be
      # present in the extraction set or PLINK will silently drop the row.
      targets <- unique(c(lead, targets))
      if (length(targets) > max_target_variants) {
        return(limit("lead-row target set exceeded max_target_variants"))
      }
      if (n_pairs + length(targets) > max_pair_requests) {
        return(limit("lead-row pair-request limit exceeded"))
      }
      for (key in keys) assign(key, TRUE, envir = pair_tested)
      n_pairs <- n_pairs + length(keys)
      n_unique_leads <- n_unique_leads + 1L
      n_calls <- n_calls + 1L
      ld <- fastmr_clump_run_frontier(lead, targets, reference_args, plink2_bin,
                                      clump_kb, clump_r2, threads, workdir, n_calls)
      if (nrow(ld)) {
        for (j in seq_len(nrow(ld))) {
          assign(fastmr_clump_pair_key(ld$lead[j], ld$target[j]), TRUE,
                 envir = pair_positive)
        }
      }
    }
    for (exposure in seq_along(states)) {
      lead_index <- leads[exposure]
      if (is.na(lead_index)) next
      lead <- as.character(dat$SNP[lead_index])
      retained[lead_index] <- TRUE
      state <- states[[exposure]]
      state$dead[match(lead_index, state$index)] <- TRUE
      live <- which(!state$dead)
      if (length(live)) {
        candidate <- state$index[live]
        same <- is.na(position$chr[lead_index]) | is.na(position$chr[candidate]) |
          position$chr[lead_index] == position$chr[candidate]
        close <- is.na(position$bp[lead_index]) | is.na(position$bp[candidate]) |
          abs(position$bp[lead_index] - position$bp[candidate]) <= clump_kb * 1000
        candidate <- candidate[same & close]
        if (length(candidate)) {
          keys <- fastmr_clump_pair_key(lead, as.character(dat$SNP[candidate]))
          blocked <- vapply(keys, exists, logical(1), envir = pair_positive, inherits = FALSE)
          state$dead[match(candidate, state$index)] <- blocked
        }
      }
      states[[exposure]] <- state
    }
  }
  retained_key <- paste(as.character(dat$id.exposure[retained]), as.character(dat$SNP[retained]), sep = "\r")
  original_key <- paste(as.character(original$id.exposure), as.character(original$SNP), sep = "\r")
  result <- original[original_key %in% retained_key, , drop = FALSE]
  instruments <- lapply(split(result$SNP, result$id.exposure, drop = TRUE), as.character)
  list(data = result, instruments = instruments,
       diagnostics = list(exposures = length(exposure_ids), candidate_rows = nrow(dat),
                          retained = nrow(dat[retained, , drop = FALSE]), rounds = n_rounds,
                          plink_calls = n_calls, logical_pairs = n_pairs,
                          positive_pairs = length(ls(pair_positive)), unique_leads = n_unique_leads,
                          exact = TRUE, fallback = FALSE, strategy = "lead_row",
                          reference_manifest_md5 = reference_md5))
}

#' Chromosome-partitioned batched PLINK2 LD clumping
#'
#' This is the scalable form of [fast_clump_data_batched()].  Variants on
#' different chromosomes cannot be in LD, so each chromosome can be clumped
#' independently.  The partition keeps the PLINK target union local to one
#' chromosome while preserving the exact exposure-specific greedy decisions.
#'
#' @param dat Data frame containing `SNP`, `id.exposure`, p-values,
#'   `chr_name`, and `chrom_start`.
#' @param ... Arguments forwarded to [fast_clump_data_batched()].
#' @return A list with `data`, named `instruments`, and aggregated diagnostics.
#' @export
fast_clump_data_batched_chromosomal <- function(dat, ...) {
  if (!is.data.frame(dat) || !all(c("chr_name", "chrom_start") %in% names(dat))) {
    stop("chromosome-partitioned clumping requires chr_name and chrom_start", call. = FALSE)
  }
  chr <- as.character(dat[["chr_name"]])
  bp <- suppressWarnings(as.numeric(as.character(dat[["chrom_start"]])))
  if (anyNA(chr) || any(!nzchar(trimws(chr))) || any(!is.finite(bp))) {
    stop("chr_name and chrom_start must be complete for chromosome-partitioned clumping", call. = FALSE)
  }
  dots <- list(...)
  workdir <- dots$workdir
  workdir_owned <- is.null(workdir)
  if (workdir_owned) workdir <- tempfile("fastMR_chromosomal_")
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  if (workdir_owned) on.exit(unlink(workdir, recursive = TRUE, force = TRUE), add = TRUE)
  chromosomes <- unique(chr)
  pieces <- lapply(seq_along(chromosomes), function(k) {
    cc <- chromosomes[[k]]
    idx <- which(chr == cc)
    part <- dat[idx, , drop = FALSE]
    part$.fastmr_row_id <- idx
    part_dots <- dots
    part_dots$workdir <- file.path(workdir, paste0("chr_", gsub("[^A-Za-z0-9_.-]", "_", cc)))
    ans <- do.call(fast_clump_data_batched, c(list(dat = part), part_dots))
    ans$chromosome <- cc
    ans
  })
  names(pieces) <- chromosomes
  retained <- lapply(pieces, `[[`, "data")
  retained <- retained[vapply(retained, nrow, integer(1)) > 0L]
  if (length(retained)) {
    combined <- do.call(rbind, retained)
    combined <- combined[order(combined$.fastmr_row_id), , drop = FALSE]
    combined$.fastmr_row_id <- NULL
    rownames(combined) <- NULL
  } else {
    combined <- dat[FALSE, , drop = FALSE]
  }
  instruments <- lapply(split(combined$SNP, combined$id.exposure, drop = TRUE), as.character)
  diagnostics <- lapply(pieces, `[[`, "diagnostics")
  sum_diag <- function(name, default = 0) {
    vals <- vapply(diagnostics, function(x) if (is.null(x[[name]])) default else x[[name]], numeric(1))
    sum(vals)
  }
  reference_manifest <- dots$reference_manifest
  list(
    data = combined,
    instruments = instruments,
    diagnostics = list(
      exposures = length(unique(as.character(if ("id.exposure" %in% names(dat)) dat$id.exposure else "exposure"))),
      candidate_rows = nrow(dat), retained = nrow(combined),
      rounds = sum_diag("rounds"), plink_calls = sum_diag("plink_calls"),
      logical_pairs = sum_diag("logical_pairs"), positive_pairs = sum_diag("positive_pairs"),
      exact = all(vapply(diagnostics, function(x) isTRUE(x$exact), logical(1))),
      fallback = any(vapply(diagnostics, function(x) isTRUE(x$fallback), logical(1))),
      partition = "chromosome", chromosomes = diagnostics,
      reference_manifest_md5 = if (!is.null(reference_manifest) && file.exists(reference_manifest))
        unname(tools::md5sum(reference_manifest)) else NULL
    )
  )
}

fastmr_have_compressor_fn <- function(name) {
  exists(name, envir = asNamespace("CompreSSoR"), inherits = FALSE)
}

# CompreSSoR with one-pass candidate extraction (values, key, p and exact rank
# for only the candidate rows; p bit-identical to read_sumstats()).
fastmr_have_one_pass_candidates <- function() {
  fastmr_have_compressor_fn("compressor_capabilities") &&
    fastmr_have_compressor_fn("read_candidates_batch") &&
    "candidates_one_pass" %in% CompreSSoR::compressor_capabilities()
}

# CompreSSoR whose batch candidate reader accepts strategy = "pvalue_flag".
fastmr_have_flag_candidates_batch <- function() {
  fastmr_have_compressor_fn("read_candidates_batch") &&
    "pvalue_flag" %in% eval(formals(CompreSSoR::read_candidates_batch)$strategy)
}

# Indirection over CompreSSoR::read_candidates_batch() (lets tests substitute a
# faulty batch reader).
fastmr_read_candidates_batch <- function(...) CompreSSoR::read_candidates_batch(...)

# CompreSSoR whose read_candidates_batch() is known to return the per-store
# result for every batch composition and to stop rather than drop rows.  The
# full-store batch path has no manifest count to check against, so it is used
# only with such a build; otherwise the per-store reader runs.
fastmr_have_checked_candidates_batch <- function() {
  fastmr_have_compressor_fn("compressor_capabilities") &&
    "candidates_batch_rows_checked" %in% CompreSSoR::compressor_capabilities()
}

# Number of rows flagged in a store's p-value flag domain: the manifest's
# recorded count, or (when the manifest lacks it) the flag stream itself.
fastmr_store_flag_count <- function(store, io_threads = 1L) {
  domain <- fastmr_clump_default(store$manifest$domains, list())$pvalue_flag
  n <- suppressWarnings(as.numeric(domain$hit_rows))
  if (length(n) == 1L && is.finite(n) && n >= 0) return(n)
  length(CompreSSoR::read_pvalue_flag(store, as = "row_ids", threads = io_threads))
}

# Zero-based row ids of a store's p-value flag domain.
fastmr_store_flag_rows <- function(store, io_threads = 1L) {
  as.integer(CompreSSoR::read_pvalue_flag(store, as = "row_ids", threads = io_threads))
}

# NULL when every key's position field (canonical keys are
# chromosome:position:ref:alt) equals the row's base_pair_location column;
# otherwise a short description.  Catches identity decoded against the wrong
# variant panel.
fastmr_batch_key_problem <- function(x, label) {
  if (!nrow(x)) return(NULL)
  if (!all(c("key", "base_pair_location") %in% names(x))) {
    return(paste0("store '", label, "' returned no key/base_pair_location columns"))
  }
  parts <- strsplit(as.character(x[["key"]]), ":", fixed = TRUE)
  key_pos <- suppressWarnings(as.numeric(vapply(parts, function(f) if (length(f) >= 2L) f[[2L]] else NA_character_,
                                                character(1))))
  col_pos <- suppressWarnings(as.numeric(x[["base_pair_location"]]))
  bad <- is.na(key_pos) | is.na(col_pos) | key_pos != col_pos
  if (any(bad)) {
    return(paste0("store '", label, "' returned ", sum(bad),
                  " key(s) whose position does not match base_pair_location"))
  }
  NULL
}

# NULL when a read_candidates_batch(strategy = "pvalue_flag") result has, for
# every store, exactly the store's flagged rows (read at the flag's own
# threshold, so before the user threshold is applied) with keys consistent
# with their positions; otherwise a short description of the first mismatch.
# `flag_rows`, when given, holds each store's flagged row ids
# (fastmr_store_flag_rows()); the returned row ids must equal them as a set.
fastmr_flag_batch_problem <- function(got, stores, labels, io_threads = 1L, flag_rows = NULL) {
  if (!is.list(got) || length(got) != length(stores)) {
    return(paste0("returned ", if (is.list(got)) length(got) else 0L,
                  " tables for ", length(stores), " stores"))
  }
  if (!is.null(names(got)) && !identical(names(got), labels)) {
    return("store names out of order")
  }
  for (i in seq_along(stores)) {
    x <- got[[i]]
    if (!is.data.frame(x)) return(paste0("store '", labels[[i]], "' returned no table"))
    expected <- if (is.null(flag_rows)) fastmr_store_flag_count(stores[[i]], io_threads) else
      length(flag_rows[[i]])
    if (nrow(x) != expected || ("row" %in% names(x) && anyDuplicated(x[["row"]]))) {
      return(paste0("store '", labels[[i]], "' returned ", nrow(x), " of ", expected,
                    " flagged rows"))
    }
    if (!is.null(flag_rows)) {
      if (!"row" %in% names(x)) return(paste0("store '", labels[[i]], "' returned no row ids"))
      if (!identical(sort(as.integer(x[["row"]])), sort(as.integer(flag_rows[[i]])))) {
        return(paste0("store '", labels[[i]], "' returned row ids that are not its flagged rows"))
      }
    }
    problem <- fastmr_batch_key_problem(x, labels[[i]])
    if (!is.null(problem)) return(problem)
  }
  NULL
}

# NULL when a full-store read_candidates_batch() result is internally
# consistent (one table per store, unique row ids, keys matching positions).
fastmr_full_batch_problem <- function(got, labels) {
  if (!is.list(got) || length(got) != length(labels)) {
    return(paste0("returned ", if (is.list(got)) length(got) else 0L,
                  " tables for ", length(labels), " stores"))
  }
  for (i in seq_along(got)) {
    x <- got[[i]]
    if (!is.data.frame(x)) return(paste0("store '", labels[[i]], "' returned no table"))
    if (!"row" %in% names(x) || anyDuplicated(x[["row"]])) {
      return(paste0("store '", labels[[i]], "' returned missing or duplicated row ids"))
    }
    problem <- fastmr_batch_key_problem(x, labels[[i]])
    if (!is.null(problem)) return(problem)
  }
  NULL
}

# |z| of the kept candidate rows: the clumpers' tie-break for equal p
# (reconstructed p underflows to 0 above |z| ~ 38).  Internal column
# `.fastmr_abs_z`, removed again by fast_clump_compressed().
fastmr_candidate_abs_z <- function(x, keep) {
  if (!"z" %in% names(x)) return(rep(NA_real_, sum(keep)))
  abs(suppressWarnings(as.numeric(x[["z"]][keep])))
}

# Candidate table from a read_candidates_batch() result: native row order,
# p <= pvalue_threshold, and the exact rank when requested.
fastmr_candidates_from_batch <- function(got, labels, pvalue_threshold, exact_order) {
  lapply(seq_along(got), function(i) {
    x <- got[[i]]
    x <- x[order(x$row), , drop = FALSE]
    p <- suppressWarnings(as.numeric(x[["p_value"]]))
    keep <- is.finite(p) & p <= pvalue_threshold
    out <- data.frame(SNP = x[["key"]][keep], id.exposure = rep(labels[[i]], sum(keep)),
                      pval.exposure = p[keep],
                      chr_name = as.character(x[["chromosome"]][keep]),
                      chrom_start = as.numeric(x[["base_pair_location"]][keep]),
                      .fastmr_abs_z = fastmr_candidate_abs_z(x, keep),
                      stringsAsFactors = FALSE, check.names = FALSE)
    if (exact_order) {
      rank <- x[["exact_rank"]][keep]
      if (anyNA(rank) || any(rank <= 0L)) {
        stop("store '", labels[[i]], "' has candidates without an exact p-value rank; cannot guarantee exact clumping order",
             call. = FALSE)
      }
      out$pvalue_rank <- as.integer(rank)
    }
    out
  })
}

fastmr_compressed_candidate_data <- function(paths, labels, pvalue_threshold,
                                              candidate_source, pvalue_order,
                                              io_threads) {
  fastmr_require_compressor()
  stores <- lapply(paths, CompreSSoR::open_compressor)
  names(stores) <- labels
  flag_thresholds <- NULL
  if (identical(candidate_source, "pvalue_flag")) {
    if (!"read_pvalue_flag" %in% getNamespaceExports("CompreSSoR")) {
      stop("candidate_source='pvalue_flag' requires the current CompreSSoR read_pvalue_flag() API; see CompreSSoR#45", call. = FALSE)
    }
    flag_thresholds <- vapply(stores, function(store) {
      domains <- fastmr_clump_default(store$manifest$domains, list())
      domain <- domains$pvalue_flag
      if (is.null(domain) || !isTRUE(fastmr_clump_default(domain$enabled, TRUE))) {
        stop("a requested store has no p-value flag domain; use candidate_source='full' or rebuild it with pvalue_flag=TRUE", call. = FALSE)
      }
      threshold <- as.numeric(domain$threshold)
      if (length(threshold) != 1L || is.na(threshold) || !is.finite(threshold)) {
        stop("a requested store has malformed p-value flag metadata; use candidate_source='full'", call. = FALSE)
      }
      if (pvalue_threshold > threshold) {
        stop("pvalue_threshold (", pvalue_threshold, ") is less selective than the store p-value flag threshold (",
             threshold, "); use candidate_source='full' for this threshold", call. = FALSE)
      }
      threshold
    }, numeric(1))
  }
  exact_order <- identical(pvalue_order, "require_exact")
  if (exact_order) {
    # Exact ordering needs CompreSSoR's stored rank domain (CompreSSoR#45).
    if (!fastmr_have_compressor_fn("read_pvalue_order")) {
      stop("exact p-value ordering is not available in the installed CompreSSoR (no read_pvalue_order()); resolve CompreSSoR#45 or use pvalue_order='reconstructed'", call. = FALSE)
    }
    for (i in seq_along(stores)) {
      domain <- fastmr_clump_default(stores[[i]]$manifest$domains, list())$pvalue_order
      if (!is.list(domain) || is.null(domain$threshold)) {
        stop("store '", labels[[i]], "' has no exact p-value ordering domain; rebuild it with pvalue_order=TRUE or use pvalue_order='reconstructed'", call. = FALSE)
      }
      if (pvalue_threshold > as.numeric(domain$threshold)) {
        stop("pvalue_threshold (", pvalue_threshold, ") exceeds the exact p-value ordering domain threshold (",
             as.numeric(domain$threshold), ") of store '", labels[[i]], "'", call. = FALSE)
      }
    }
  }
  data <- NULL
  rows <- NULL
  if (!is.null(flag_thresholds)) {
    # One read_candidates_batch(strategy = "pvalue_flag") pass: the flagged
    # rows' key, p (bit-identical to read_sumstats()) and, for exact ordering,
    # only the exact ranks of those rows.  The per-store fallback below decodes
    # the whole rank vector (~9M entries per store) to look up a few thousand.
    # The batch reader needs the flag's own threshold; the user threshold is
    # applied afterwards, as in the fallback, so membership stays the flag's.
    flag_rows <- NULL
    if (fastmr_have_flag_candidates_batch()) {
      # The flagged row ids are the request: the batch must return exactly
      # these rows.  A CompreSSoR build with "candidates_batch_rows_checked"
      # verifies that itself (every store's decoded rows must be identical to
      # the rows its own flag read selected, and their number must equal the
      # manifest's flagged-row count, or the batch stops), so the flag stream
      # is not decoded a second time here; the count, duplicate-row and
      # key/position checks below still run.  Older builds get the flagged row
      # ids read up front, and they are reused by the per-store fallback.
      if (!fastmr_have_checked_candidates_batch()) {
        flag_rows <- lapply(stores, fastmr_store_flag_rows, io_threads = io_threads)
      }
      got <- tryCatch(
        fastmr_read_candidates_batch(
          as.list(stats::setNames(paths, labels)), pvalue_threshold = unname(flag_thresholds),
          columns = c("key", "p_value", "chromosome", "base_pair_location", "z"),
          order = if (exact_order) "exact" else "none", threads = io_threads,
          strategy = "pvalue_flag"),
        error = function(e) e)
      # Trust the batch only if every store returned exactly its flagged
      # rows (CompreSSoR 0.7.0 could silently drop rows in mixed batches):
      # the manifest count (or, without the capability, the flagged row ids).
      problem <- if (inherits(got, "error")) {
        paste("read_candidates_batch() failed:", conditionMessage(got))
      } else {
        fastmr_flag_batch_problem(got, stores, labels, io_threads, flag_rows = flag_rows)
      }
      if (is.null(problem)) {
        data <- fastmr_candidates_from_batch(got, labels, pvalue_threshold, exact_order)
      } else {
        warning("batched p-value flag candidate read rejected (", problem,
                "); falling back to the per-store reader", call. = FALSE)
      }
    }
    if (is.null(data)) {
      rows <- if (!is.null(flag_rows)) flag_rows else {
        flag_reader <- getExportedValue("CompreSSoR", "read_pvalue_flag")
        lapply(stores, function(store) flag_reader(store, threads = io_threads))
      }
    }
  }
  columns <- c("chromosome", "base_pair_location", "effect_allele", "other_allele", "p_value", "z")
  # The aligned flag returns immutable zero-based row IDs.  The current
  # CompreSSoR batch reader accepts canonical keys, not row IDs, so use the
  # native single-store reader here; independent stores can still be decoded
  # concurrently.  This avoids a full variant-table scan just to translate
  # the flag rows back to keys.
  reader_threads <- if (length(paths) > 1L && io_threads > 1L) 1L else io_threads
  have_candidates <- is.null(rows) && fastmr_have_compressor_fn("read_candidates")
  reader <- function(i) {
    tryCatch({
      if (!is.null(rows)) {
        x <- CompreSSoR::read_sumstats(paths[[i]], variants = if (length(rows[[i]])) rows[[i]] else 0L,
                                       columns = columns, threads = reader_threads)
        if (!length(rows[[i]])) x <- x[0L, , drop = FALSE]
        if (nrow(x) != length(rows[[i]])) {
          stop("flagged-row read returned ", nrow(x), " of ", length(rows[[i]]), " rows")
        }
        attr(x, "fastmr_row") <- as.integer(rows[[i]])
        return(x)
      }
      # Full-store candidates: avoid passing seq.int(0, n - 1) as `variants`
      # (58 s / 6.4 GB per 10M-row store).  Same rows, same native order.
      # read_candidates() only *selects* rows (with a hair of slack, since its
      # erfc p can differ from read_sumstats() in the last ulp); the values
      # and the final p <= threshold test come from read_sumstats(), so the
      # result is identical to the full read.
      slack <- min(1, pvalue_threshold * (1 + 1e-9))
      sel <- NULL
      if (have_candidates) {
        sel <- tryCatch(
          CompreSSoR::read_candidates(paths[[i]],
                                      pvalue_threshold = slack,
                                      columns = "base_pair_location", order = "none",
                                      threads = reader_threads),
          error = function(e) NULL)
      }
      if (is.data.frame(sel) && "row" %in% names(sel)) {
        row <- as.integer(sel[["row"]])
      } else {
        pv <- CompreSSoR::read_sumstats(paths[[i]], columns = "p_value", threads = reader_threads)
        pv <- suppressWarnings(as.numeric(pv[["p_value"]]))
        row <- which(is.finite(pv) & pv <= slack) - 1L
      }
      x <- CompreSSoR::read_sumstats(paths[[i]], variants = if (length(row)) row else 0L,
                                     columns = columns, threads = reader_threads)
      if (!length(row)) x <- x[0L, , drop = FALSE]
      if (nrow(x) != length(row)) stop("candidate selection/read row mismatch")
      attr(x, "fastmr_row") <- row
      x
    }, error = function(e) stop("failed to read candidates from ", paths[[i]], ": ",
                               conditionMessage(e), call. = FALSE))
  }
  if (is.null(data) && is.null(flag_thresholds) && fastmr_have_one_pass_candidates() &&
      fastmr_have_checked_candidates_batch()) {
    # One read_candidates_batch() per exposure batch: candidate rows, key, p and
    # (for exact ordering) the exact rank in one block-selective pass per store,
    # with same-panel identity decoded once.  No p slack is needed because p is
    # bit-identical to the full read.
    got <- tryCatch(
      fastmr_read_candidates_batch(
        as.list(stats::setNames(paths, labels)), pvalue_threshold = pvalue_threshold,
        columns = c("key", "p_value", "chromosome", "base_pair_location", "z"),
        order = if (exact_order) "exact" else "none", threads = io_threads),
      error = function(e) NULL)
    if (!is.null(got)) {
      problem <- fastmr_full_batch_problem(got, labels)
      if (is.null(problem)) {
        data <- fastmr_candidates_from_batch(got, labels, pvalue_threshold, exact_order)
      } else {
        warning("batched candidate read rejected (", problem,
                "); falling back to the per-store reader", call. = FALSE)
      }
    }
  }
  if (is.null(data)) {
  pieces <- if (.Platform$OS.type != "windows" && io_threads > 1L && length(paths) > 1L) {
    parallel::mclapply(seq_along(paths), reader,
                       mc.cores = min(io_threads, length(paths)), mc.preschedule = TRUE)
  } else {
    lapply(seq_along(paths), reader)
  }
  names(pieces) <- labels
  data <- lapply(seq_along(pieces), function(i) {
    x <- pieces[[i]]
    if (!is.data.frame(x)) stop("compressed reader returned an invalid candidate table", call. = FALSE)
    if (!nrow(x)) return(data.frame(SNP = character(), id.exposure = character(),
                                    pval.exposure = numeric(), chr_name = character(),
                                    chrom_start = numeric(), .fastmr_abs_z = numeric(),
                                    stringsAsFactors = FALSE, check.names = FALSE))
    key <- CompreSSoR::compressor_variant_key(
      x[["chromosome"]], x[["base_pair_location"]],
      x[["other_allele"]], x[["effect_allele"]]
    )
    p <- suppressWarnings(as.numeric(x[["p_value"]]))
    keep <- is.finite(p) & p <= pvalue_threshold
    out <- data.frame(SNP = key[keep], id.exposure = labels[[i]], pval.exposure = p[keep],
                      chr_name = as.character(x[["chromosome"]][keep]),
                      chrom_start = as.numeric(x[["base_pair_location"]][keep]),
                      .fastmr_abs_z = fastmr_candidate_abs_z(x, keep),
                      stringsAsFactors = FALSE, check.names = FALSE)
    if (exact_order) {
      ranks <- CompreSSoR::read_pvalue_order(paths[[i]], as = "ranks", fallback = "error",
                                             threads = reader_threads)
      rank <- ranks[attr(x, "fastmr_row")[keep] + 1L]
      if (anyNA(rank) || any(rank <= 0L)) {
        stop("store '", labels[[i]], "' has candidates without an exact p-value rank; cannot guarantee exact clumping order",
             call. = FALSE)
      }
      out$pvalue_rank <- as.integer(rank)
    }
    out
  })
  }
  names(data) <- labels
  if (exact_order) {
    for (i in seq_along(data)) if (!"pvalue_rank" %in% names(data[[i]])) data[[i]]$pvalue_rank <- integer()
  }
  list(data = do.call(rbind, data), stores = stores,
       exact = exact_order,
       source = paste0(if (identical(candidate_source, "pvalue_flag")) "pvalue_flag" else "full_store",
                       if (exact_order) "_then_exact_rank_order" else "_then_reconstructed_p"))
}

#' Generate clumped instruments directly from Pcodec exposure stores
#'
#' Candidate variants are extracted from every exposure store and clumped
#' exactly with the strategy chosen by `partition`.  The default p-value flag
#' path reads only the rows marked by CompreSSoR's aligned flag domain (one
#' batched pass, including only those rows' exact ranks when
#' `pvalue_order = "require_exact"`), then keeps p <= `pvalue_threshold`.
#'
#' @param exposure_files Named character vector of Pcodec stores.
#' @param pvalue_threshold Candidate p-value threshold.
#' @param candidate_source `"pvalue_flag"` (default) or `"full"`.
#' @param pvalue_order `"reconstructed"` (default) or `"require_exact"`.
#'   The latter uses CompreSSoR's `read_pvalue_order()` rank domain when the
#'   installed CompreSSoR and each store provide it, and errors otherwise.
#' @param output Optional Parquet path for the clumped candidate table.
#' @param partition `"auto"` (default), `"graph"`, `"per_exposure"`,
#'   `"global"`, `"chromosome"`, or `"lead_row"`.  All return identical
#'   instruments.  `"auto"` ([fast_clump_data_auto()]) picks `"graph"` or
#'   `"per_exposure"` from the estimated number of candidate pairs within
#'   `clump_kb`, the number of exposures and `threads`, and records the choice
#'   in `diagnostics$auto`.  `"graph"` makes one PLINK2 all-pairs LD call over
#'   the candidate union plus a C++ greedy pass (lead-row fallback above
#'   `max_graph_pairs`, see [fast_clump_data_graph()]); `"per_exposure"` runs
#'   one PLINK2 `--clump` per exposure on a candidate-only reference, which is
#'   much faster when the pair graph is dense (permissive `clump_r2` and wide
#'   `clump_kb`), see [fast_clump_data_per_exposure()].  `"global"`,
#'   `"chromosome"` and `"lead_row"` are the earlier frontier strategies.
#' @param ... Arguments forwarded to the selected batched clumping function.
#' @return A list with `data`, named `instruments`, and `diagnostics`.
#' @export
fast_clump_compressed <- function(
    exposure_files, pvalue_threshold = 5e-8,
    candidate_source = c("pvalue_flag", "full"),
    pvalue_order = c("reconstructed", "require_exact"), output = NULL,
    partition = c("auto", "graph", "per_exposure", "global", "chromosome", "lead_row"), ...) {
  fastmr_require_compressor()
  paths <- fastmr_normalize_compressed_files(exposure_files, "exposure_files")
  pvalue_threshold <- fastmr_clump_number(pvalue_threshold, "pvalue_threshold", 0, 1)
  candidate_source <- match.arg(candidate_source)
  pvalue_order <- match.arg(pvalue_order)
  partition <- match.arg(partition)
  stores <- lapply(paths, CompreSSoR::open_compressor)
  invisible(lapply(stores, fastmr_validate_compressed_store))
  dots <- list(...)
  io_threads <- fastmr_positive_integer_scalar(
    fastmr_clump_default(dots$io_threads, 1L), "io_threads"
  )
  candidates <- fastmr_compressed_candidate_data(
    paths, names(paths), pvalue_threshold, candidate_source, pvalue_order,
    io_threads = io_threads
  )
  dots$io_threads <- NULL
  reference_manifest <- dots$reference_manifest
  clump_fun <- switch(partition,
                      auto = fast_clump_data_auto,
                      per_exposure = fast_clump_data_per_exposure,
                      global = fast_clump_data_batched,
                      chromosome = fast_clump_data_batched_chromosomal,
                      lead_row = fast_clump_data_lead_rows,
                      graph = fast_clump_data_graph)
  clumped <- do.call(clump_fun, c(list(dat = candidates$data), dots))
  clumped$data[[".fastmr_abs_z"]] <- NULL
  if (!is.null(output)) fast_write_parquet(clumped$data, output)
  clumped$diagnostics$compressed_input <- list(
    stores = unname(paths), pvalue_threshold = pvalue_threshold,
    candidate_source = candidates$source,
    pvalue_order = pvalue_order,
    partition = partition,
    pvalue_order_exact = candidates$exact,
    reference_manifest_md5 = if (!is.null(reference_manifest) && file.exists(reference_manifest))
      unname(tools::md5sum(reference_manifest)) else NULL
  )
  clumped
}
