#' Fast local harmonisation of exposure and outcome summary statistics
#'
#' This is the local, dependency-free analogue of
#' [TwoSampleMR::harmonise_data()]. It implements all three native `action`
#' policies, all four native allele-information cases (2-2, 2-1, 1-2, and
#' 1-1), strand swaps, palindromic-frequency checks, and TwoSampleMR's local
#' indel recoding rules. Proxy lookup and remote extraction remain outside this
#' local summary-statistics layer.
#'
#' @param exposure_dat Exposure data with TwoSampleMR allele and beta columns.
#' @param outcome_dat Outcome data with TwoSampleMR allele and beta columns.
#' @param action Integer 1, 2, or 3 controlling palindromic handling: retain
#'   palindromes, infer strand from allele frequencies, or remove palindromes.
#' @param tolerance Frequency tolerance around 0.5 for palindromic SNPs.
#' @return A harmonised data frame with `mr_keep`, `remove`, `palindromic`, and
#'   `ambiguous` columns, plus a compact `log` attribute. Rows without an
#'   effect allele in either dataset follow native TwoSampleMR behavior and are
#'   omitted from the harmonised result.
#' @export
fast_harmonise_data <- function(exposure_dat, outcome_dat, action = 2,
                                tolerance = 0.08) {
  if (!is.data.frame(exposure_dat) || !is.data.frame(outcome_dat)) {
    stop("exposure_dat and outcome_dat must be data frames", call. = FALSE)
  }
  if (length(action) < 1L || any(!is.finite(action)) || any(!action %in% 1:3)) {
    stop("action must contain only 1, 2, or 3", call. = FALSE)
  }
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 || tolerance >= 0.5) {
    stop("tolerance must be a finite value between 0 and 0.5", call. = FALSE)
  }
  required <- function(dat, type) {
    needed <- c("SNP", "beta", "se", "effect_allele", "other_allele")
    names_needed <- paste0(needed[-1L], ".", type)
    missing <- setdiff(c("SNP", names_needed), names(dat))
    if (length(missing)) {
      stop("missing required ", type, " column(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
  }
  required(exposure_dat, "exposure")
  required(outcome_dat, "outcome")

  add_optional <- function(dat, column, value = NA_real_) {
    if (!column %in% names(dat)) dat[[column]] <- value
    dat
  }
  exposure_dat <- add_optional(exposure_dat, "eaf.exposure")
  outcome_dat <- add_optional(outcome_dat, "eaf.outcome")
  if (!"id.exposure" %in% names(exposure_dat)) exposure_dat$id.exposure <- "exposure"
  if (!"id.outcome" %in% names(outcome_dat)) outcome_dat$id.outcome <- "outcome"

  merged <- merge(outcome_dat, exposure_dat, by = "SNP", sort = FALSE)
  if (!nrow(merged)) {
    attr(merged, "log") <- data.frame(candidate_variants = nrow(exposure_dat),
                                       variants_absent_from_reference = nrow(exposure_dat),
                                       total_variants = 0L, total_variants_for_mr = 0L)
    return(merged)
  }

  # Match harmonise_cleanup_variables() in TwoSampleMR. In particular, NR is
  # a missing outcome frequency and an empty outcome other allele is missing.
  allele_cols <- c("effect_allele.exposure", "other_allele.exposure",
                   "effect_allele.outcome", "other_allele.outcome")
  for (column in allele_cols) merged[[column]] <- toupper(as.character(merged[[column]]))
  merged$beta.exposure <- suppressWarnings(as.numeric(merged$beta.exposure))
  merged$beta.outcome <- suppressWarnings(as.numeric(merged$beta.outcome))

  merged$eaf.exposure <- suppressWarnings(as.numeric(as.character(merged$eaf.exposure)))
  outcome_eaf <- as.character(merged$eaf.outcome)
  outcome_eaf[outcome_eaf %in% c("NR", "NR ")] <- NA_character_
  merged$eaf.outcome <- suppressWarnings(as.numeric(outcome_eaf))

  # TwoSampleMR partitions rows into four cases. Rows with no effect allele
  # are not handled by the native harmonise_* helpers and are therefore
  # omitted here as well.
  a1 <- merged$effect_allele.exposure
  a2 <- merged$other_allele.exposure
  b1 <- merged$effect_allele.outcome
  b2 <- merged$other_allele.outcome
  i22 <- !is.na(a1) & !is.na(a2) & !is.na(b1) & !is.na(b2)
  i21 <- !is.na(a1) & !is.na(a2) & !is.na(b1) & is.na(b2)
  i12 <- !is.na(a1) & is.na(a2) & !is.na(b1) & !is.na(b2)
  i11 <- !is.na(a1) & is.na(a2) & !is.na(b1) & is.na(b2)
  valid <- i22 | i21 | i12 | i11
  if (!any(valid)) {
    empty <- merged[FALSE, , drop = FALSE]
    attr(empty, "log") <- data.frame(candidate_variants = nrow(exposure_dat),
                                      variants_absent_from_reference = nrow(exposure_dat),
                                      total_variants = 0L, total_variants_for_mr = 0L)
    return(empty)
  }
  merged <- merged[valid, , drop = FALSE]
  i22 <- i22[valid]; i21 <- i21[valid]; i12 <- i12[valid]; i11 <- i11[valid]
  a1 <- merged$effect_allele.exposure
  a2 <- merged$other_allele.exposure
  b1 <- merged$effect_allele.outcome
  b2 <- merged$other_allele.outcome
  beta_a <- merged$beta.exposure
  beta_b <- merged$beta.outcome
  eaf_a <- merged$eaf.exposure
  eaf_b <- merged$eaf.outcome
  if (length(action) == 1L) action <- rep(action, length(unique(merged$id.outcome)))
  if (length(action) != length(unique(merged$id.outcome))) {
    stop("action must have length 1 or one value per unique id.outcome", call. = FALSE)
  }
  action_map <- setNames(action, unique(as.character(merged$id.outcome)))
  row_action <- unname(action_map[as.character(merged$id.outcome)])
  remove <- rep(FALSE, nrow(merged))
  palindromic <- rep(FALSE, nrow(merged))
  ambiguous <- rep(FALSE, nrow(merged))
  switched <- rep(FALSE, nrow(merged))
  flipped_basic <- rep(FALSE, nrow(merged))
  flipped_palindrome <- rep(FALSE, nrow(merged))
  flip_alleles <- function(x) chartr("ACGTacgt", "TGCAtgca", x)
  is_palindromic <- function(x, y) (x == "T" & y == "A") | (x == "A" & y == "T") |
    (x == "G" & y == "C") | (x == "C" & y == "G")
  minf <- 0.5 - tolerance
  maxf <- 0.5 + tolerance

  # 2-2: complete alleles in both studies, including native indel recoding.
  if (any(i22)) {
    ii <- which(i22)
    aa1 <- a1[ii]; aa2 <- a2[ii]; bb1 <- b1[ii]; bb2 <- b2[ii]
    nca1 <- nchar(aa1); nca2 <- nchar(aa2); ncb1 <- nchar(bb1); ncb2 <- nchar(bb2)
    indel <- nca1 > 1 | nca2 > 1 | aa1 %in% c("D", "I")
    r <- rep(TRUE, length(ii))
    z <- indel & nca1 > nca2 & bb1 == "I" & bb2 == "D"; bb1[z] <- aa1[z]; bb2[z] <- aa2[z]
    z <- indel & nca1 < nca2 & bb1 == "I" & bb2 == "D"; bb1[z] <- aa2[z]; bb2[z] <- aa1[z]
    z <- indel & nca1 > nca2 & bb1 == "D" & bb2 == "I"; bb1[z] <- aa2[z]; bb2[z] <- aa1[z]
    z <- indel & nca1 < nca2 & bb1 == "D" & bb2 == "I"; bb1[z] <- aa1[z]; bb2[z] <- aa2[z]
    z <- indel & ncb1 > ncb2 & aa1 == "I" & aa2 == "D"; aa1[z] <- bb1[z]; aa2[z] <- bb2[z]
    z <- indel & ncb1 < ncb2 & aa1 == "I" & aa2 == "D"; aa2[z] <- bb1[z]; aa1[z] <- bb2[z]
    z <- indel & ncb1 > ncb2 & aa1 == "D" & aa2 == "I"; aa2[z] <- bb1[z]; aa1[z] <- bb2[z]
    z <- indel & ncb1 < ncb2 & aa1 == "D" & aa2 == "I"; aa1[z] <- bb1[z]; aa2[z] <- bb2[z]
    r[indel & nca1 > 1 & nca1 == nca2 & bb1 %in% c("D", "I")] <- FALSE
    r[indel & ncb1 > 1 & ncb1 == ncb2 & aa1 %in% c("D", "I")] <- FALSE
    r[aa1 == aa2] <- FALSE; r[bb1 == bb2] <- FALSE
    status <- aa1 == bb1 & aa2 == bb2
    z <- aa1 == bb2 & aa2 == bb1
    switched[ii[z]] <- TRUE
    tmp <- bb1[z]; bb1[z] <- bb2[z]; bb2[z] <- tmp
    beta_b[ii[z]] <- -beta_b[ii[z]]; eaf_b[ii[z]] <- 1 - eaf_b[ii[z]]
    status <- aa1 == bb1 & aa2 == bb2
    pal <- is_palindromic(aa1, aa2)
    z <- !pal & !status
    bb1[z] <- flip_alleles(bb1[z]); bb2[z] <- flip_alleles(bb2[z]); flipped_basic[ii[z]] <- TRUE
    status <- aa1 == bb1 & aa2 == bb2
    z <- !pal & !status & aa1 == bb2 & aa2 == bb1
    switched[ii[z]] <- TRUE
    tmp <- bb1[z]; bb1[z] <- bb2[z]; bb2[z] <- tmp
    beta_b[ii[z]] <- -beta_b[ii[z]]; eaf_b[ii[z]] <- 1 - eaf_b[ii[z]]
    status <- aa1 == bb1 & aa2 == bb2
    rem <- !status; rem[!r] <- TRUE
    fa <- eaf_a[ii]; fb <- eaf_b[ii]; fa[is.na(fa)] <- 0.5; fb[is.na(fb)] <- 0.5
    amb <- ((fa > minf & fa < maxf) | (fb > minf & fb < maxf)) & pal
    z <- ((fa < 0.5 & fb > 0.5) | (fa > 0.5 & fb < 0.5)) & pal & !rem & row_action[ii] == 2L
    if (any(z)) {
      beta_b[ii[z]] <- -beta_b[ii[z]]; eaf_b[ii[z]] <- 1 - eaf_b[ii[z]]
      flipped_palindrome[ii[z]] <- TRUE
    }
    a1[ii] <- aa1; a2[ii] <- aa2; b1[ii] <- bb1; b2[ii] <- bb2
    remove[ii] <- rem; palindromic[ii] <- pal; ambiguous[ii] <- amb
  }

  # 2-1: both exposure alleles and only the outcome effect allele.
  if (any(i21)) {
    ii <- which(i21)
    aa1 <- a1[ii]; aa2 <- a2[ii]; bb1 <- b1[ii]; bb2 <- rep(NA_character_, length(ii))
    nca1 <- nchar(aa1); nca2 <- nchar(aa2); ncb1 <- nchar(bb1)
    indel <- nca1 > 1 | nca2 > 1 | aa1 %in% c("D", "I")
    r <- rep(TRUE, length(ii))
    z <- indel & nca1 > nca2 & bb1 == "I"; bb1[z] <- aa1[z]; bb2[z] <- aa2[z]
    z <- indel & nca1 < nca2 & bb1 == "I"; bb1[z] <- aa2[z]; bb2[z] <- aa1[z]
    z <- indel & nca1 > nca2 & bb1 == "D"; bb1[z] <- aa2[z]; bb2[z] <- aa1[z]
    z <- indel & nca1 < nca2 & bb1 == "D"; bb1[z] <- aa1[z]; bb2[z] <- aa2[z]
    r[indel & aa1 == "I" & aa2 == "D"] <- FALSE
    r[indel & aa1 == "D" & aa2 == "I"] <- FALSE
    r[indel & nca1 > 1 & nca1 == nca2 & bb1 %in% c("D", "I")] <- FALSE
    r[aa1 == aa2] <- FALSE
    pal <- is_palindromic(aa1, aa2); rem <- pal
    fa <- eaf_a[ii]; fb <- eaf_b[ii]; fa[is.na(fa)] <- 0.5; fb[is.na(fb)] <- 0.5
    similar1 <- (fa < minf & fb < minf) | (fa > maxf & fb > maxf)
    amb <- rep(FALSE, length(ii)); status <- aa1 == bb1
    amb[status & !similar1] <- TRUE; bb2[status] <- aa2[status]
    z <- aa2 == bb1
    switched[ii[z]] <- TRUE
    similar2 <- (fa < minf & fb > maxf) | (fa > maxf & fb < minf)
    amb[z & !similar2] <- TRUE
    bb2[z] <- bb1[z]; bb1[z] <- aa1[z]
    beta_b[ii[z]] <- -beta_b[ii[z]]; eaf_b[ii[z]] <- 1 - eaf_b[ii[z]]
    z <- aa1 != bb1 & aa2 != bb1
    amb[z] <- TRUE; bb1[z] <- flip_alleles(bb1[z]); flipped_basic[ii[z]] <- TRUE
    status <- aa1 == bb1; bb2[status] <- aa2[status]
    z <- aa2 == bb1
    bb2[z] <- bb1[z]; bb1[z] <- aa1[z]
    beta_b[ii[z]] <- -beta_b[ii[z]]; eaf_b[ii[z]] <- 1 - eaf_b[ii[z]]
    a1[ii] <- aa1; a2[ii] <- aa2; b1[ii] <- bb1; b2[ii] <- bb2
    remove[ii] <- rem | !r; palindromic[ii] <- pal; ambiguous[ii] <- amb | pal
  }

  # 1-2: only the exposure effect allele and both outcome alleles.
  if (any(i12)) {
    ii <- which(i12)
    aa1 <- a1[ii]; aa2 <- rep(NA_character_, length(ii)); bb1 <- b1[ii]; bb2 <- b2[ii]
    ncb1 <- nchar(bb1); ncb2 <- nchar(bb2); nca1 <- nchar(aa1)
    indel <- ncb1 > 1 | ncb2 > 1 | bb1 %in% c("D", "I")
    r <- rep(TRUE, length(ii))
    z <- indel & ncb1 > ncb2 & aa1 == "I"; aa1[z] <- bb1[z]; aa2[z] <- bb2[z]
    z <- indel & ncb1 < ncb2 & aa1 == "I"; aa2[z] <- bb1[z]; aa1[z] <- bb2[z]
    z <- indel & ncb1 > ncb2 & aa1 == "D"; aa2[z] <- bb1[z]; aa1[z] <- bb2[z]
    z <- indel & ncb1 < ncb2 & aa1 == "D"; aa1[z] <- bb1[z]; aa2[z] <- bb2[z]
    r[indel & bb1 == "I" & bb2 == "D"] <- FALSE
    r[indel & bb1 == "D" & bb2 == "I"] <- FALSE
    r[indel & ncb1 > 1 & ncb1 == ncb2 & aa1 %in% c("D", "I")] <- FALSE
    r[indel & bb1 == bb2] <- FALSE
    pal <- is_palindromic(bb1, bb2); rem <- pal
    fa <- eaf_a[ii]; fb <- eaf_b[ii]; fa[is.na(fa)] <- 0.5; fb[is.na(fb)] <- 0.5
    similar1 <- (fa < minf & fb < minf) | (fa > maxf & fb > maxf)
    amb <- rep(FALSE, length(ii)); status <- aa1 == bb1
    amb[status & !similar1] <- TRUE; aa2[status] <- bb2[status]
    z <- aa1 == bb2
    switched[ii[z]] <- TRUE
    similar2 <- (fa < minf & fb > maxf) | (fa > maxf & fb < minf)
    amb[z & !similar2] <- TRUE
    aa2[z] <- aa1[z]; aa1[z] <- bb1[z]
    beta_a[ii[z]] <- -beta_a[ii[z]]; eaf_a[ii[z]] <- 1 - eaf_a[ii[z]]
    z <- aa1 != bb1 & aa1 != bb2
    amb[z] <- TRUE; aa1[z] <- flip_alleles(aa1[z]); flipped_basic[ii[z]] <- TRUE
    status <- aa1 == bb1; aa2[status] <- bb2[status]
    z <- bb2 == aa1
    bb2[z] <- bb1[z]; bb1[z] <- aa1[z]
    beta_b[ii[z]] <- -beta_b[ii[z]]; eaf_b[ii[z]] <- 1 - eaf_b[ii[z]]
    a1[ii] <- aa1; a2[ii] <- aa2; b1[ii] <- bb1; b2[ii] <- bb2
    remove[ii] <- rem | !r; palindromic[ii] <- pal; ambiguous[ii] <- amb | pal
  }

  # 1-1: only the effect allele is available in both studies.
  if (any(i11)) {
    ii <- which(i11)
    status <- a1[ii] == b1[ii]
    fa <- eaf_a[ii]; fb <- eaf_b[ii]; fa[is.na(fa)] <- 0.5; fb[is.na(fb)] <- 0.5
    similar <- (fa < minf & fb < minf) | (fa > maxf & fb > maxf)
    remove[ii] <- !status
    ambiguous[ii] <- status & !similar
    a2[ii] <- NA_character_; b2[ii] <- NA_character_
  }

  merged$effect_allele.exposure <- a1
  merged$other_allele.exposure <- a2
  merged$effect_allele.outcome <- b1
  merged$other_allele.outcome <- b2
  merged$beta.exposure <- beta_a
  merged$beta.outcome <- beta_b
  merged$eaf.exposure <- eaf_a
  merged$eaf.outcome <- eaf_b
  merged$remove <- remove
  merged$palindromic <- palindromic
  merged$ambiguous <- ambiguous

  merged$mr_keep <- !remove
  merged$mr_keep[row_action == 2L] <- merged$mr_keep[row_action == 2L] & !ambiguous[row_action == 2L]
  merged$mr_keep[row_action == 3L] <- merged$mr_keep[row_action == 3L] & !palindromic[row_action == 3L] & !ambiguous[row_action == 3L]
  complete_mr <- stats::complete.cases(merged[, c("beta.exposure", "beta.outcome", "se.exposure", "se.outcome"), drop = FALSE])
  merged$mr_keep[!complete_mr] <- FALSE
  if (!"samplesize.outcome" %in% names(merged)) merged$samplesize.outcome <- NA_real_
  attr(merged, "log") <- data.frame(
    candidate_variants = nrow(exposure_dat),
    variants_absent_from_reference = nrow(exposure_dat) - nrow(merged),
    total_variants = nrow(merged),
    total_variants_for_mr = sum(merged$mr_keep, na.rm = TRUE),
    switched_alleles = sum(switched),
    flipped_alleles_basic = sum(flipped_basic),
    flipped_alleles_palindrome = sum(flipped_palindrome),
    incompatible_alleles = sum(remove, na.rm = TRUE),
    ambiguous_alleles = sum(ambiguous, na.rm = TRUE)
  )
  merged
}

#' Fast local LD clumping
#'
#' With `bfile` and `plink_bin`, delegates LD calculation to local PLINK.
#' With `ld_matrix`, performs deterministic greedy index-SNP clumping in R
#' without network access. The matrix row names (or `ld_snps`) must identify
#' the SNPs and be in the same order as its rows and columns.
#'
#' @param dat Data frame containing `SNP`, optionally `pval.exposure`, and
#'   optionally `id.exposure`, `chr_name`, and `chrom_start`.
#' @param clump_kb Maximum physical distance between index and secondary SNPs.
#' @param clump_r2 Minimum LD r-squared for removal.
#' @param clump_p1 Maximum p-value for index SNPs.
#' @param ld_matrix Optional local correlation matrix.
#' @param ld_snps Optional SNP names for `ld_matrix`.
#' @param bfile Optional PLINK binary-prefix for local clumping.
#' @param plink_bin Optional PLINK executable path.
#' @return `dat` filtered to retained index SNPs.
#' @export
fast_clump_data <- function(dat, clump_kb = 10000, clump_r2 = 0.001,
                            clump_p1 = 1, ld_matrix = NULL, ld_snps = NULL,
                            bfile = NULL, plink_bin = NULL) {
  if (!is.data.frame(dat) || !"SNP" %in% names(dat)) stop("dat must contain SNP", call. = FALSE)
  if (length(clump_kb) != 1L || !is.finite(clump_kb) || clump_kb < 0 ||
      length(clump_r2) != 1L || !is.finite(clump_r2) || clump_r2 < 0 || clump_r2 > 1 ||
      length(clump_p1) != 1L || !is.finite(clump_p1) || clump_p1 < 0 || clump_p1 > 1) {
    stop("clump_kb, clump_r2, and clump_p1 must be finite values in range", call. = FALSE)
  }
  pval_column <- if ("pval.exposure" %in% names(dat)) "pval.exposure" else if ("pval.outcome" %in% names(dat)) "pval.outcome" else NULL
  if (is.null(pval_column)) {
    dat$pval.exposure <- clump_p1
    pval_column <- "pval.exposure"
  }
  if (!"id.exposure" %in% names(dat)) dat$id.exposure <- "exposure"
  if (!is.null(bfile)) {
    if (is.null(plink_bin)) plink_bin <- Sys.which("plink")
    if (!nzchar(plink_bin)) stop("PLINK executable not found; provide plink_bin", call. = FALSE)
    pieces <- lapply(split(seq_len(nrow(dat)), dat$id.exposure, drop = TRUE), function(index) {
      if (length(index) <= 1L) return(index)
      stem <- tempfile("fastMR_clump_")
      input <- paste0(stem, ".txt")
      write.table(data.frame(SNP = dat$SNP[index], P = dat[[pval_column]][index]), input,
                  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
      on.exit(unlink(c(input, paste0(stem, ".clumped"), paste0(stem, ".log"), paste0(stem, ".nosex"))), add = TRUE)
      status <- system2(plink_bin, c("--bfile", bfile, "--clump", input,
        "--clump-p1", as.character(clump_p1), "--clump-r2", as.character(clump_r2),
        "--clump-kb", as.character(clump_kb), "--out", stem), stdout = FALSE, stderr = FALSE)
      if (!identical(status, 0L) || !file.exists(paste0(stem, ".clumped"))) {
        stop("PLINK clumping failed for exposure ", unique(dat$id.exposure[index]), call. = FALSE)
      }
      clumped <- tryCatch(read.table(paste0(stem, ".clumped"), header = TRUE, stringsAsFactors = FALSE),
                          error = function(e) data.frame(SNP = character()))
      index[dat$SNP[index] %in% clumped$SNP]
    })
    return(dat[sort(unlist(pieces, use.names = FALSE)), , drop = FALSE])
  }
  if (is.null(ld_matrix)) stop("provide ld_matrix for dependency-free local clumping or bfile for PLINK", call. = FALSE)
  ld_matrix <- as.matrix(ld_matrix)
  if (!is.numeric(ld_matrix) || nrow(ld_matrix) != ncol(ld_matrix)) stop("ld_matrix must be a numeric square matrix", call. = FALSE)
  if (is.null(ld_snps)) ld_snps <- rownames(ld_matrix)
  if (is.null(ld_snps) || length(ld_snps) != nrow(ld_matrix)) stop("ld_snps or row names are required for ld_matrix", call. = FALSE)
  ld_snps <- as.character(ld_snps)
  if (anyDuplicated(ld_snps)) stop("ld_snps must be unique", call. = FALSE)
  positions <- match(as.character(dat$SNP), ld_snps)
  chr <- if ("chr_name" %in% names(dat)) as.character(dat$chr_name) else rep(NA_character_, nrow(dat))
  bp <- if ("chrom_start" %in% names(dat)) suppressWarnings(as.numeric(dat$chrom_start)) else rep(NA_real_, nrow(dat))
  p <- suppressWarnings(as.numeric(dat[[pval_column]]))
  keep <- logical(nrow(dat))
  for (group in split(seq_len(nrow(dat)), dat$id.exposure, drop = TRUE)) {
    candidates <- group[!is.na(positions[group]) & is.finite(p[group]) & p[group] <= clump_p1]
    candidates <- candidates[order(p[candidates], as.character(dat$SNP[candidates]))]
    selected <- integer()
    for (i in candidates) {
      if (!length(selected)) {
        keep[i] <- TRUE; selected <- c(selected, i); next
      }
      same_chr <- is.na(chr[i]) | is.na(chr[selected]) | chr[i] == chr[selected]
      close <- is.na(bp[i]) | is.na(bp[selected]) | abs(bp[i] - bp[selected]) <= clump_kb * 1000
      r2 <- ld_matrix[positions[i], positions[selected]]^2
      if (!any(same_chr & close & is.finite(r2) & r2 >= clump_r2)) {
        keep[i] <- TRUE; selected <- c(selected, i)
      }
    }
  }
  dat[keep, , drop = FALSE]
}
