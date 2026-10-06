# Rscript equiv-compare.R <reference.rds> <candidate.rds> [more candidates...]
args <- commandArgs(TRUE)
ref <- readRDS(args[1])
for (f in args[-1]) {
  cand <- readRDS(f)
  stopifnot(identical(names(ref$results), names(cand$results)))
  same_value <- mapply(function(a, b) identical(a$value, b$value), ref$results, cand$results)
  same_state <- mapply(function(a, b) identical(a$state, b$state), ref$results, cand$results)
  # Thread invariance inside the candidate: every /tN case equals its /t1 twin.
  nm <- names(cand$results)
  tn <- grepl("/t[248](/|$)", nm)
  twin <- sub("/t[248](/|$)", "/t1\\1", nm[tn])
  thread_same <- mapply(function(a, b) identical(cand$results[[a]], cand$results[[b]]), nm[tn], twin)
  cat(sprintf("%s vs %s: %d cases; value identical %d/%d; .Random.seed identical %d/%d; thread twins identical %d/%d\n",
              basename(f), basename(args[1]), length(nm), sum(same_value), length(nm), sum(same_state), length(nm),
              sum(thread_same), length(thread_same)))
  cat("  candidate mode path counts:", paste(names(cand$counts), cand$counts, sep = "=", collapse = " "), "\n")
  bad <- nm[!(same_value & same_state)]
  if (length(bad)) { cat("  MISMATCH:\n"); print(head(bad, 30)) }
  if (any(!thread_same)) { cat("  THREAD MISMATCH:\n"); print(head(nm[tn][!thread_same], 30)) }
}
