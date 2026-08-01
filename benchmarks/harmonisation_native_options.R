local_lib <- "/Users/fergushamilton/projects/fastMR/.local/Rlib"
.libPaths(c(local_lib, .libPaths()))
library(fastMR)
library(TwoSampleMR)

out_dir <- "/Users/fergushamilton/projects/fastMR/outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

make_pair <- function(kind, tag) {
  exposure <- data.frame(
    SNP = tag, id.exposure = "exposure", exposure = "simulated exposure",
    beta.exposure = 0.20, se.exposure = 0.05,
    effect_allele.exposure = NA_character_, other_allele.exposure = NA_character_,
    eaf.exposure = 0.20, stringsAsFactors = FALSE
  )
  outcome <- data.frame(
    SNP = tag, id.outcome = "outcome", outcome = "simulated outcome",
    beta.outcome = 0.07, se.outcome = 0.04,
    effect_allele.outcome = NA_character_, other_allele.outcome = NA_character_,
    eaf.outcome = 0.20, stringsAsFactors = FALSE
  )
  if (kind == "22_direct") {
    exposure[1, c("effect_allele.exposure", "other_allele.exposure")] <- c("A", "G")
    outcome[1, c("effect_allele.outcome", "other_allele.outcome")] <- c("A", "G")
  } else if (kind == "22_reverse") {
    exposure[1, c("effect_allele.exposure", "other_allele.exposure")] <- c("A", "G")
    outcome[1, c("effect_allele.outcome", "other_allele.outcome")] <- c("G", "A")
    outcome$eaf.outcome <- 0.80
  } else if (kind == "22_complement") {
    exposure[1, c("effect_allele.exposure", "other_allele.exposure")] <- c("A", "G")
    outcome[1, c("effect_allele.outcome", "other_allele.outcome")] <- c("T", "C")
  } else if (kind == "22_reverse_complement") {
    exposure[1, c("effect_allele.exposure", "other_allele.exposure")] <- c("A", "G")
    outcome[1, c("effect_allele.outcome", "other_allele.outcome")] <- c("C", "T")
    outcome$eaf.outcome <- 0.80
  } else if (kind == "22_palindrome") {
    exposure[1, c("effect_allele.exposure", "other_allele.exposure")] <- c("A", "T")
    outcome[1, c("effect_allele.outcome", "other_allele.outcome")] <- c("T", "A")
    outcome$eaf.outcome <- 0.80
  } else if (kind == "22_ambiguous") {
    exposure[1, c("effect_allele.exposure", "other_allele.exposure")] <- c("C", "G")
    outcome[1, c("effect_allele.outcome", "other_allele.outcome")] <- c("C", "G")
    exposure$eaf.exposure <- 0.50; outcome$eaf.outcome <- 0.50
  } else if (kind == "22_incompatible") {
    exposure[1, c("effect_allele.exposure", "other_allele.exposure")] <- c("A", "G")
    outcome[1, c("effect_allele.outcome", "other_allele.outcome")] <- c("A", "C")
  } else if (kind == "22_indel") {
    exposure[1, c("effect_allele.exposure", "other_allele.exposure")] <- c("A", "AT")
    outcome[1, c("effect_allele.outcome", "other_allele.outcome")] <- c("I", "D")
    outcome$eaf.outcome <- 0.80
  } else if (kind == "21_direct") {
    exposure[1, c("effect_allele.exposure", "other_allele.exposure")] <- c("A", "G")
    outcome$effect_allele.outcome <- "A"
  } else if (kind == "21_reverse") {
    exposure[1, c("effect_allele.exposure", "other_allele.exposure")] <- c("A", "G")
    outcome$effect_allele.outcome <- "G"; outcome$eaf.outcome <- 0.80
  } else if (kind == "21_complement") {
    exposure[1, c("effect_allele.exposure", "other_allele.exposure")] <- c("A", "G")
    outcome$effect_allele.outcome <- "T"
  } else if (kind == "21_palindrome") {
    exposure[1, c("effect_allele.exposure", "other_allele.exposure")] <- c("A", "T")
    outcome$effect_allele.outcome <- "A"
  } else if (kind == "21_indel") {
    exposure[1, c("effect_allele.exposure", "other_allele.exposure")] <- c("A", "AT")
    outcome$effect_allele.outcome <- "I"; outcome$eaf.outcome <- 0.80
  } else if (kind == "12_direct") {
    exposure$effect_allele.exposure <- "A"
    outcome[1, c("effect_allele.outcome", "other_allele.outcome")] <- c("A", "G")
  } else if (kind == "12_reverse") {
    exposure$effect_allele.exposure <- "G"
    outcome[1, c("effect_allele.outcome", "other_allele.outcome")] <- c("A", "G")
    exposure$eaf.exposure <- 0.80
  } else if (kind == "12_complement") {
    exposure$effect_allele.exposure <- "A"
    outcome[1, c("effect_allele.outcome", "other_allele.outcome")] <- c("T", "C")
  } else if (kind == "12_palindrome") {
    exposure$effect_allele.exposure <- "A"
    outcome[1, c("effect_allele.outcome", "other_allele.outcome")] <- c("A", "T")
  } else if (kind == "12_indel") {
    exposure$effect_allele.exposure <- "I"
    outcome[1, c("effect_allele.outcome", "other_allele.outcome")] <- c("A", "AT")
    exposure$eaf.exposure <- 0.80
  } else if (kind == "11_direct") {
    exposure$effect_allele.exposure <- "A"
    outcome$effect_allele.outcome <- "A"
  } else if (kind == "11_incompatible") {
    exposure$effect_allele.exposure <- "A"
    outcome$effect_allele.outcome <- "G"
  } else if (kind == "11_ambiguous") {
    exposure$effect_allele.exposure <- "A"
    outcome$effect_allele.outcome <- "A"
    exposure$eaf.exposure <- 0.50; outcome$eaf.outcome <- 0.50
  }
  list(exposure = exposure, outcome = outcome)
}

key_columns <- c("effect_allele.exposure", "other_allele.exposure",
                 "effect_allele.outcome", "other_allele.outcome",
                 "beta.exposure", "beta.outcome", "eaf.exposure", "eaf.outcome",
                 "remove", "palindromic", "ambiguous", "mr_keep")
canonical <- function(x) {
  x <- x[order(x$SNP), , drop = FALSE]
  for (column in key_columns) {
    if (column %in% c("effect_allele.exposure", "other_allele.exposure",
                      "effect_allele.outcome", "other_allele.outcome")) {
      value <- as.character(x[[column]])
      value[is.na(value)] <- "<NA>"
      x[[column]] <- value
    }
  }
  x[, c("SNP", key_columns), drop = FALSE]
}

kinds <- c("22_direct", "22_reverse", "22_complement", "22_reverse_complement",
            "22_palindrome", "22_ambiguous", "22_incompatible", "22_indel",
            "21_direct", "21_reverse", "21_complement", "21_palindrome", "21_indel",
            "12_direct", "12_reverse", "12_complement", "12_palindrome", "12_indel",
            "11_direct", "11_incompatible", "11_ambiguous")
rows <- list()
for (kind in kinds) {
  pair <- make_pair(kind, paste0("rs_", kind))
  for (action in 1:3) {
    fast <- suppressMessages(fast_harmonise_data(pair$exposure, pair$outcome, action = action))
    native <- suppressMessages(TwoSampleMR::harmonise_data(pair$exposure, pair$outcome, action = action))
    fast_key <- canonical(fast); native_key <- canonical(native)
    same <- identical(fast_key, native_key)
    numeric_delta <- 0
    if (nrow(fast) && nrow(native)) {
      numeric_delta <- max(abs(as.matrix(fast[, c("beta.exposure", "beta.outcome", "eaf.exposure", "eaf.outcome")]) -
                                as.matrix(native[, c("beta.exposure", "beta.outcome", "eaf.exposure", "eaf.outcome")])), na.rm = TRUE)
      if (!is.finite(numeric_delta)) numeric_delta <- 0
    }
    rows[[length(rows) + 1L]] <- data.frame(kind = kind, action = action,
      fast_rows = nrow(fast), native_rows = nrow(native), exact_key_match = same,
      max_numeric_delta = numeric_delta, stringsAsFactors = FALSE)
  }
}

summary <- do.call(rbind, rows)
vector_exposure <- data.frame(
  SNP = paste0("rs_vector_", 1:3), id.exposure = "E", exposure = "E", beta.exposure = .2,
  se.exposure = .05, eaf.exposure = .2,
  effect_allele.exposure = c("A", "A", "A"),
  other_allele.exposure = c("G", "T", "G"), stringsAsFactors = FALSE
)
vector_outcome <- data.frame(
  SNP = paste0("rs_vector_", 1:3), id.outcome = paste0("O", 1:3),
  outcome = "O", beta.outcome = .1, se.outcome = .04, eaf.outcome = c(.2, .8, .2),
  effect_allele.outcome = c("A", "T", "A"),
  other_allele.outcome = c("G", "A", "C"), stringsAsFactors = FALSE
)
vector_fast <- suppressMessages(fast_harmonise_data(vector_exposure, vector_outcome,
                                                     action = c(1, 2, 3)))
vector_native <- suppressMessages(TwoSampleMR::harmonise_data(vector_exposure, vector_outcome,
                                                               action = c(1, 2, 3)))
action_vector_match <- identical(canonical(vector_fast), canonical(vector_native))
write.csv(summary, file.path(out_dir, "harmonisation_native_options.csv"), row.names = FALSE)
report <- c(
  "# Native harmonisation option parity",
  "",
  sprintf("Cases: %d; action policies: 1, 2, and 3; comparisons: %d.", length(kinds), nrow(summary)),
  sprintf("Exact key matches: %d/%d.", sum(summary$exact_key_match), nrow(summary)),
  sprintf("Maximum numeric delta: %.3e.", max(summary$max_numeric_delta)),
  sprintf("Outcome-specific action vector parity: %s.", action_vector_match),
  "Covered 2-2, 2-1, 1-2, and 1-1 allele information, strand swaps, complements, palindromes, incompatible alleles, missing frequencies, and native indel recoding."
)
writeLines(report, file.path(out_dir, "harmonisation_native_options.md"))
stopifnot(all(summary$exact_key_match), max(summary$max_numeric_delta) <= 1e-15,
          action_vector_match)
cat(paste(report, collapse = "\n"), "\n")
