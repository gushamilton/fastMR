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
            "--r2-unphased", "zs", "--ld-window-kb", format(clump_kb, trim = TRUE),
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
  dedup <- !duplicated(paste(as.character(dat$id.exposure), as.character(dat$SNP), sep = "\r"))
  dat <- dat[dedup, , drop = FALSE]
  p <- suppressWarnings(as.numeric(as.character(dat[[pcol]])))
  position <- fastmr_clump_position(dat)
  exposure_ids <- unique(as.character(dat$id.exposure))
  states <- lapply(exposure_ids, function(id) {
    ii <- which(as.character(dat$id.exposure) == id)
    ii <- ii[is.finite(p[ii]) & p[ii] <= clump_p1]
    ii <- ii[order(p[ii], as.character(dat$SNP[ii]), method = "radix")]
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
  dedup <- !duplicated(paste(as.character(dat$id.exposure), as.character(dat$SNP), sep = "\r"))
  dat <- dat[dedup, , drop = FALSE]
  p <- suppressWarnings(as.numeric(as.character(dat[[pcol]])))
  position <- fastmr_clump_position(dat)
  exposure_ids <- unique(as.character(dat$id.exposure))
  states <- lapply(exposure_ids, function(id) {
    ii <- which(as.character(dat$id.exposure) == id)
    ii <- ii[is.finite(p[ii]) & p[ii] <= clump_p1]
    ii <- ii[order(p[ii], as.character(dat$SNP[ii]), method = "radix")]
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

fastmr_compressed_candidate_data <- function(paths, labels, pvalue_threshold,
                                              candidate_source, pvalue_order,
                                              io_threads) {
  fastmr_require_compressor()
  stores <- lapply(paths, CompreSSoR::open_compressor)
  names(stores) <- labels
  if (identical(candidate_source, "pvalue_flag")) {
    if (!"read_pvalue_flag" %in% getNamespaceExports("CompreSSoR")) {
      stop("candidate_source='pvalue_flag' requires the current CompreSSoR read_pvalue_flag() API; see CompreSSoR#45", call. = FALSE)
    }
    flag_reader <- getExportedValue("CompreSSoR", "read_pvalue_flag")
    rows <- lapply(stores, function(store) {
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
      flag_reader(store, threads = io_threads)
    })
  } else {
    nrows <- vapply(stores, function(store) as.integer(fastmr_clump_default(store$manifest$n_rows, store$manifest$rows)), integer(1))
    rows <- lapply(nrows, function(n) seq.int(0L, n - 1L))
  }
  if (identical(pvalue_order, "require_exact")) {
    # Pcodec 0.5 reconstructs p-values from the stored Z stream.  The exact
    # order-preserving domain is tracked in CompreSSoR#45; do not claim exact
    # clumping order until that API has both manifest metadata and a reader.
    stop("exact p-value ordering is not available in the current CompreSSoR reader; resolve CompreSSoR#45 or use pvalue_order='reconstructed'", call. = FALSE)
  }
  columns <- c("chromosome", "base_pair_location", "effect_allele", "other_allele", "p_value")
  # The aligned flag returns immutable zero-based row IDs.  The current
  # CompreSSoR batch reader accepts canonical keys, not row IDs, so use the
  # native single-store reader here; independent stores can still be decoded
  # concurrently.  This avoids a full variant-table scan just to translate
  # the flag rows back to keys.
  reader_threads <- if (length(paths) > 1L && io_threads > 1L) 1L else io_threads
  reader <- function(i) {
    tryCatch(
      CompreSSoR::read_sumstats(paths[[i]], variants = rows[[i]], columns = columns,
                                threads = reader_threads),
      error = function(e) stop("failed to read candidates from ", paths[[i]], ": ",
                               conditionMessage(e), call. = FALSE)
    )
  }
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
                                    chrom_start = numeric(), stringsAsFactors = FALSE))
    key <- CompreSSoR::compressor_variant_key(
      x[["chromosome"]], x[["base_pair_location"]],
      x[["other_allele"]], x[["effect_allele"]]
    )
    p <- suppressWarnings(as.numeric(x[["p_value"]]))
    keep <- is.finite(p) & p <= pvalue_threshold
    data.frame(SNP = key[keep], id.exposure = labels[[i]], pval.exposure = p[keep],
               chr_name = as.character(x[["chromosome"]][keep]),
               chrom_start = as.numeric(x[["base_pair_location"]][keep]),
               stringsAsFactors = FALSE)
  })
  names(data) <- labels
  list(data = do.call(rbind, data), stores = stores,
       exact = FALSE,
       source = if (identical(candidate_source, "pvalue_flag")) "pvalue_flag_then_reconstructed_p" else "full_store_reconstructed_p")
}

#' Generate clumped instruments directly from Pcodec exposure stores
#'
#' Candidate variants are extracted from every exposure store, then one
#' exposure-grouped PLINK2 frontier is used to clump them.  The default
#' p-value flag path reads only the rows marked by CompreSSoR's aligned flag
#' domain; p-values are still reconstructed from stored Z values until the
#' exact ordering domain in CompreSSoR#45 is available.
#'
#' @param exposure_files Named character vector of Pcodec stores.
#' @param pvalue_threshold Candidate p-value threshold.
#' @param candidate_source `"pvalue_flag"` (default) or `"full"`.
#' @param pvalue_order `"reconstructed"` (default) or `"require_exact"`.
#' @param output Optional Parquet path for the clumped candidate table.
#' @param partition `"global"` (default), `"chromosome"`, or `"lead_row"`.
#'   Chromosome partitioning avoids an unnecessarily large cross-chromosome
#'   target union; lead-row mode shares only the LD row relevant to each
#'   current lead and is preferable when exposures overlap substantially.
#' @param ... Arguments forwarded to the selected batched clumping function.
#' @return A list with `data`, named `instruments`, and `diagnostics`.
#' @export
fast_clump_compressed <- function(
    exposure_files, pvalue_threshold = 5e-8,
    candidate_source = c("pvalue_flag", "full"),
    pvalue_order = c("reconstructed", "require_exact"), output = NULL,
    partition = c("global", "chromosome", "lead_row"), ...) {
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
                      global = fast_clump_data_batched,
                      chromosome = fast_clump_data_batched_chromosomal,
                      lead_row = fast_clump_data_lead_rows)
  clumped <- do.call(clump_fun, c(list(dat = candidates$data), dots))
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
