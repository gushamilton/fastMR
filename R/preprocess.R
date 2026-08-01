#' Fast local harmonisation of exposure and outcome summary statistics
#'
#' This is the local, dependency-free analogue of
#' [TwoSampleMR::harmonise_data()]. It handles common bi-allelic SNPs,
#' complements, strand swaps, palindromic-frequency checks, and the standard
#' `action` policies. Indel recoding and proxy lookup remain outside this
#' local summary-statistics layer.
#'
#' @param exposure_dat Exposure data with TwoSampleMR allele and beta columns.
#' @param outcome_dat Outcome data with TwoSampleMR allele and beta columns.
#' @param action Integer 1, 2, or 3 controlling palindromic handling.
#' @param tolerance Frequency tolerance around 0.5 for palindromic SNPs.
#' @return A harmonised data frame with `mr_keep`, `remove`, `palindromic`, and
#'   `ambiguous` columns, plus a compact `log` attribute.
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
  # Normalize the allele strings once; avoid a per-SNP R loop in the hot path.
  allele_cols <- c("effect_allele.exposure", "other_allele.exposure",
                   "effect_allele.outcome", "other_allele.outcome")
  for (column in allele_cols) merged[[column]] <- toupper(as.character(merged[[column]]))
  merged$beta.exposure <- suppressWarnings(as.numeric(merged$beta.exposure))
  merged$beta.outcome <- suppressWarnings(as.numeric(merged$beta.outcome))
  merged$eaf.exposure <- suppressWarnings(as.numeric(merged$eaf.exposure))
  merged$eaf.outcome <- suppressWarnings(as.numeric(merged$eaf.outcome))

  complement <- c(A = "T", T = "A", C = "G", G = "C")
  comp <- function(x) unname(complement[x])
  a1 <- merged$effect_allele.exposure
  a2 <- merged$other_allele.exposure
  b1 <- merged$effect_allele.outcome
  b2 <- merged$other_allele.outcome
  beta_out <- merged$beta.outcome
  eaf_out <- merged$eaf.outcome

  complete <- !is.na(a1) & !is.na(a2) & !is.na(b1) & !is.na(b2)
  palindromic <- complete & ((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
    (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))
  direct <- complete & a1 == b1 & a2 == b2
  reverse <- complete & a1 != a2 & a1 == b2 & a2 == b1
  to_swap <- reverse
  if (any(to_swap)) {
    tmp <- b1[to_swap]; b1[to_swap] <- b2[to_swap]; b2[to_swap] <- tmp
    beta_out[to_swap] <- -beta_out[to_swap]
    eaf_out[to_swap] <- 1 - eaf_out[to_swap]
  }

  # Resolve non-palindromic strand complements before checking compatibility.
  complementable <- complete & !palindromic & !direct & !to_swap &
    b1 %in% names(complement) & b2 %in% names(complement)
  if (any(complementable)) {
    b1[complementable] <- comp(b1[complementable])
    b2[complementable] <- comp(b2[complementable])
    reverse_comp <- complementable & a1 == b2 & a2 == b1
    if (any(reverse_comp)) {
      tmp <- b1[reverse_comp]; b1[reverse_comp] <- b2[reverse_comp]; b2[reverse_comp] <- tmp
      beta_out[reverse_comp] <- -beta_out[reverse_comp]
      eaf_out[reverse_comp] <- 1 - eaf_out[reverse_comp]
    }
  }
  aligned <- complete & a1 == b1 & a2 == b2
  remove <- !aligned
  minf <- 0.5 - tolerance
  maxf <- 0.5 + tolerance
  f_a <- merged$eaf.exposure
  f_b <- eaf_out
  temp_a <- f_a; temp_b <- f_b
  temp_a[is.na(temp_a)] <- 0.5
  temp_b[is.na(temp_b)] <- 0.5
  ambiguous <- palindromic & ((temp_a > minf & temp_a < maxf) |
                              (temp_b > minf & temp_b < maxf))
  if (length(action) == 1L) action <- rep(action, length(unique(merged$id.outcome)))
  if (length(action) != length(unique(merged$id.outcome))) {
    stop("action must have length 1 or one value per unique id.outcome", call. = FALSE)
  }
  action_map <- setNames(action, unique(merged$id.outcome))
  row_action <- unname(action_map[as.character(merged$id.outcome)])
  if (any(row_action == 2L)) {
    opposite <- palindromic & !remove &
      row_action == 2L & ((temp_a < 0.5 & temp_b > 0.5) | (temp_a > 0.5 & temp_b < 0.5))
    beta_out[opposite] <- -beta_out[opposite]
    eaf_out[opposite] <- 1 - eaf_out[opposite]
  }
  merged$effect_allele.outcome <- b1
  merged$other_allele.outcome <- b2
  merged$beta.outcome <- beta_out
  merged$eaf.outcome <- eaf_out
  merged$remove <- remove
  merged$palindromic <- palindromic
  merged$ambiguous <- ambiguous
  merged$mr_keep <- !remove
  merged$mr_keep[row_action == 2L] <- merged$mr_keep[row_action == 2L] & !ambiguous[row_action == 2L]
  merged$mr_keep[row_action == 3L] <- merged$mr_keep[row_action == 3L] & !palindromic[row_action == 3L] & !ambiguous[row_action == 3L]
  attr(merged, "log") <- data.frame(
    candidate_variants = nrow(exposure_dat),
    variants_absent_from_reference = nrow(exposure_dat) - nrow(merged),
    total_variants = nrow(merged),
    total_variants_for_mr = sum(merged$mr_keep, na.rm = TRUE),
    switched_alleles = sum(to_swap | (complementable & a1 == b2 & a2 == b1), na.rm = TRUE),
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
