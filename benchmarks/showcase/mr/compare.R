#!/usr/bin/env Rscript
# Agreement between arms for one (scenario, set, size): reference = tsmr1 (else tsmr8). Appends to --out (agreement table) and
# --boot (per-row SE pairs for bootstrap-based methods, only up to 1000 pairs) CSVs. Base R only.
a <- commandArgs(TRUE); o <- list(); for (x in a) { kv <- regmatches(x, regexec("^--([^=]+)=(.*)$", x))[[1]]; if (length(kv) == 3) o[[kv[2]]] <- kv[3] }
fs <- file.path(o$rds, sprintf("%s_%s_%s_%s.rds", o$scen, o$size, o$set, c("tsmr1", "tsmr8", "tsmr1_seedB", "fast1", "fast8", "p2c1")))
names(fs) <- c("tsmr1", "tsmr8", "tsmr1_seedB", "fast1", "fast8", "p2c1"); fs <- fs[file.exists(fs)]
if (length(fs) < 2) quit(save = "no")
res <- lapply(fs, readRDS)
ref <- if ("tsmr1" %in% names(res)) "tsmr1" else if ("tsmr8" %in% names(res)) "tsmr8" else names(res)[1]
BOOT <- c("Weighted median", "Simple mode", "Weighted mode", "Weighted mode (NOME)", "Simple mode (NOME)")
rows <- list(); brow <- list()
w <- function(d, ...) write.table(d, ..., sep = ",", row.names = FALSE, qmethod = "double")
pairs <- list(); for (n in setdiff(names(res), ref)) pairs[[length(pairs) + 1]] <- c(ref, n)
if (all(c("fast1", "fast8") %in% names(res))) pairs[[length(pairs) + 1]] <- c("fast1", "fast8")   # thread invariance
for (pr in pairs) {
  A <- res[[pr[1]]]; B <- res[[pr[2]]]
  ka <- paste(A$key, A$method, sep = "#"); kb <- paste(B$key, B$method, sep = "#"); m <- match(ka, kb); ok <- !is.na(m)
  num <- setdiff(names(A), c("key", "method")); num <- num[vapply(A[num], is.numeric, NA)]
  meths <- unique(A$method)
  if (length(num) == 0) {   # key-set comparison (clump instruments)
    rows[[length(rows) + 1]] <- data.frame(scenario = o$scen, size = o$size, method_set = o$set, replicate = o$rep, ref_arm = pr[1], arm = pr[2],
      method = "instruments", metric = "set_identity", n_ref = nrow(A), n_arm = nrow(B), n_matched = sum(ok), max_abs_diff = NA, median_abs_diff = NA,
      n_mismatch = (nrow(A) - sum(ok)) + (nrow(B) - sum(ok)), stringsAsFactors = FALSE)
    next }
  for (me in meths) for (cn in num) {
    i <- which(A$method == me & ok); if (!length(i)) next
    d <- abs(A[[cn]][i] - B[[cn]][m[i]]); isflag <- grepl("^flag_", cn)
    rows[[length(rows) + 1]] <- data.frame(scenario = o$scen, size = o$size, method_set = o$set, replicate = o$rep, ref_arm = pr[1], arm = pr[2],
      method = me, metric = cn, n_ref = sum(A$method == me), n_arm = sum(B$method == me), n_matched = length(i),
      max_abs_diff = if (isflag) NA else suppressWarnings(max(d, na.rm = TRUE)), median_abs_diff = if (isflag) NA else median(d, na.rm = TRUE),
      n_mismatch = if (isflag) sum(d > 0, na.rm = TRUE) else NA, stringsAsFactors = FALSE)
    if (cn == "se" && me %in% BOOT && as.numeric(o$size) <= 1000 && o$scen %in% c("many", "single"))
      brow[[length(brow) + 1]] <- data.frame(scenario = o$scen, size = o$size, method_set = o$set, replicate = o$rep, ref_arm = pr[1], arm = pr[2],
        key = A$key[i], method = me, se_ref = A$se[i], se_arm = B$se[m[i]], stringsAsFactors = FALSE)
  }
}
if (length(rows)) { d <- do.call(rbind, rows); w(d, o$out, append = file.exists(o$out), col.names = !file.exists(o$out)) }
if (length(brow) && !is.null(o$boot)) { d <- do.call(rbind, brow); w(d, o$boot, append = file.exists(o$boot), col.names = !file.exists(o$boot)) }
