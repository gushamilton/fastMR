# Old (v7: fastMR fb62057 + CompreSSoR f91bb8d, lib6) vs new (v8: fastMR 543cae2 + CompreSSoR d158622, lib7) tables.
# Reads the assembled replicate dirs (final10*, final25*), the 50x50 fastMR dirs, the storage and multi-trait CSVs.
# Writes CSVs to showcase/summary_<NEW_TAG>/ and prints them. Rscript compare_v8.R
# Env NEW_TAG (default v8; v9 = the fix/clump-flag-preread build) and AB_DIR (default results_ab).
suppressMessages(library(data.table))
SH <- "/user/work/fh6520/showcase"; E <- file.path(SH, "e2e"); NT <- Sys.getenv("NEW_TAG", "v8"); OUT <- file.path(SH, paste0("summary_", NT)); dir.create(OUT, FALSE)
wall_of <- function(tf) { x <- sub(".*: ", "", grep("Elapsed \\(wall", readLines(tf), value = TRUE)); p <- as.numeric(strsplit(x, ":")[[1]])
  sum(p * 60^(rev(seq_along(p)) - 1)) }
cells <- function(root) rbindlist(lapply(list.dirs(root, recursive = TRUE), function(d) {
  b <- basename(d); if (!grepl("^[0-9]+x[0-9]+_[A-Z0-9]+$", b) || !file.exists(file.path(d, "time.txt"))) return(NULL)
  m <- fread(file.path(d, "meta.csv"), colClasses = "character")
  st <- if (file.exists(file.path(d, "stages.csv"))) fread(file.path(d, "stages.csv")) else data.table(stage = character(), wall_s = numeric())
  data.table(size = sub("_.*", "", b), arm = sub(".*_", "", b), rep = basename(dirname(d)), wall = wall_of(file.path(d, "time.txt")),
             rss_mb = as.numeric(sub(".*: ", "", grep("Maximum resident", readLines(file.path(d, "time.txt")), value = TRUE))) / 1024,
             host = m$host, cpu = m$cpu, fastmr = m$fastmr, compressor = m$compressor, stages = list(st)) }), fill = TRUE)
old <- rbind(cells(file.path(E, "final10")), cells(file.path(E, "final25")), cells(file.path(E, "results_v7_traits50_v7"))[size == "50x50"])
new <- rbind(cells(file.path(E, paste0("final10_", NT))), cells(file.path(E, paste0("final25_", NT))), cells(file.path(E, sprintf("results_%s_traits50_v8", NT)))[size == "50x50"])
new <- new[arm %in% c("C1", "C8")]
cat("== CPU models seen\n"); print(rbind(old[, .(set = "old", arm, cpu)], new[, .(set = "new", arm, cpu)])[, .N, by = .(set, cpu)])
smr <- function(d) d[, .(n = .N, median = median(wall), min = min(wall), max = max(wall), rss_mb = median(rss_mb)), by = .(size, arm)]
tot <- merge(smr(old), smr(new), by = c("size", "arm"), all = TRUE, suffixes = c("_old", "_new"))
tot[, change_pct := 100 * (median_new - median_old) / median_old]
setorder(tot, size, arm); fwrite(tot, file.path(OUT, "e2e_totals_old_vs_new.csv")); cat("\n== e2e totals\n"); print(tot, digits = 4)
# speedups by allocation: fastMR arm (new and old) vs TSMR 1 core (A), TSMR 8 CPUs (A8), GWAS-VCF + TSMR (B, MR serial)
ref <- smr(old)[arm %in% c("A", "A8", "B"), .(size, arm, ref_s = median)]
sp <- rbindlist(lapply(c("old", "new"), function(s) { f <- smr(get(s))[arm %in% c("C1", "C8"), .(size, fm_arm = arm, fm_s = median)]
  merge(f, ref, by = "size", allow.cartesian = TRUE)[, .(build = s, size, fm_arm, vs = arm, fm_s, ref_s, speedup = ref_s / fm_s)] }))
sp[, vs := fcase(vs == "A", "TSMR (MR serial, 1 core)", vs == "A8", "TSMR (MR over 8 CPUs)", vs == "B", "GWAS-VCF + TSMR (MR serial)")]
fwrite(sp, file.path(OUT, "e2e_speedups.csv")); cat("\n== speedups\n"); print(dcast(sp, size + fm_arm + vs ~ build, value.var = "speedup"), digits = 4)
stg <- function(d, s) d[arm %in% c("C1", "C8"), rbindlist(Map(function(x, r) x[, rep := r], stages, rep)), by = .(size, arm)][
  , .(median = median(wall_s)), by = .(size, arm, stage)][, build := s]
S <- dcast(rbind(stg(old, "old"), stg(new, "new")), size + arm + stage ~ build, value.var = "median")
S[, change_pct := 100 * (new - old) / old]; fwrite(S, file.path(OUT, "e2e_stages_old_vs_new.csv")); cat("\n== fastMR stages\n"); print(S, digits = 4)
# agreement
ag <- function(root, s, what) rbindlist(lapply(1:3, function(r) { f <- file.path(root, sprintf("rep%d/agreement_%s_%s.csv", r, s, what))
  if (file.exists(f)) fread(f)[, rep := r] }), fill = TRUE)
A <- rbindlist(lapply(list(c("old", "final10", "10x10"), c("old", "final25", "25x25"), c("new", paste0("final10_", NT), "10x10"), c("new", paste0("final25_", NT), "25x25"), c("old", "final10", "1x1"), c("new", paste0("final10_", NT), "1x1")),
  function(z) { root <- file.path(E, z[2]); i <- ag(root, z[3], "instr"); e <- ag(root, z[3], "est")
    if (!nrow(i)) return(NULL)
    i1 <- i[rep == 1 & arm %in% c("C1", "C8"), .(instr_ref = sum(n_ref), instr_arm = sum(n_arm), all_jaccard_1 = all(jaccard == 1)), by = arm]
    e1 <- e[arm %in% c("C1", "C8"), .(rows_matched = sum(pairs) / uniqueN(rep), ivw_med_db_se = median(med_db_se[method == "Inverse variance weighted"]),
                                     ivw_max_db_se = max(max_db_se[method == "Inverse variance weighted"]),
                                     wald_rows = sum(pairs[method == "Wald ratio"]) / uniqueN(rep)), by = arm]
    merge(i1, e1, by = "arm")[, `:=`(build = z[1], size = z[3])] }))
setcolorder(A, c("build", "size")); fwrite(A, file.path(OUT, "agreement_old_vs_new.csv")); cat("\n== agreement (C arms vs TSMR arm A)\n"); print(A, digits = 4)
dz <- setNames(c("10x10", "25x25"), paste0(c("final10_", "final25_"), NT))
W <- rbindlist(lapply(names(dz), function(z) ag(file.path(E, z), dz[[z]], "wald")), fill = TRUE)
if (nrow(W)) { fwrite(W, file.path(OUT, paste0("agreement_wald_", NT, ".csv"))); cat("\n== single-SNP Wald ratio check (", NT, ")\n"); print(W, digits = 4) }
X <- rbindlist(lapply(names(dz), function(z) ag(file.path(E, z), dz[[z]], "exact")[, size := dz[[z]]]), fill = TRUE)
if (nrow(X)) { X <- X[rep == 1]; fwrite(X, file.path(OUT, paste0("agreement_exact_", NT, ".csv"))); cat("\n== exactness fast_mr vs TSMR on the same harmonised data (", NT, ", rep1)\n"); print(X, digits = 3) }
# scaling
for (t in c("v7", NT)) { f <- file.path(E, sprintf("scaling_check_%s.csv", t)); if (file.exists(f)) { cat("\n== scaling check", t, "\n"); print(fread(f), digits = 4) } }
# storage (single FinnGen GWAS, cpr)
SR <- file.path(SH, "storage/results")
so <- function(tag) rbindlist(lapply(1:3, function(r) { f <- file.path(SR, sprintf("%s%d_ops.csv", tag, r)); if (file.exists(f)) fread(f, colClasses = list(character = "note")) }), fill = TRUE)[format == "cpr"]
st7 <- so("v7"); st8 <- so("v8")
sm <- function(d, s) d[, .(build = s, n = .N, seconds = median(seconds), min = min(seconds), max = max(seconds), wall = median(wall_s), rss_mb = median(maxrss_mb),
                           bytes = paste(unique(na.omit(bytes)), collapse = "/"), n_rows = paste(unique(n_rows), collapse = "/"), chk = paste(unique(chk_s), collapse = "/")), by = .(op, threads, size)]
ST <- merge(sm(st7, "old")[, !"build"], sm(st8, "new")[, !"build"], by = c("op", "threads", "size"), suffixes = c("_old", "_new"))
ST[, change_pct := 100 * (seconds_new - seconds_old) / seconds_old]; fwrite(ST, file.path(OUT, "storage_cpr_old_vs_new.csv")); cat("\n== storage (cpr)\n"); print(ST[, .(op, threads, size, seconds_old, seconds_new, min_new, max_new, change_pct, bytes_old, bytes_new, chk_same = chk_old == chk_new)], digits = 4)
sb <- file.path(SR, "store_bytes_v7_vs_v8.csv"); if (file.exists(sb)) { cat("\n== prep store bytes v7 vs v8\n"); print(fread(sb)) }
# multi-trait: cpr v7 vs v8, against the tsv.gz rows of the full run
MR <- file.path(SH, "storage/multi/results")
mr <- function(tag) rbindlist(lapply(1:3, function(r) { f <- file.path(MR, sprintf("%s%d_multi.csv", tag, r)); if (file.exists(f)) fread(f) }), fill = TRUE)
mm <- function(d) d[status == "ok", .(n = .N, wall = median(wall_s), n_rows = paste(unique(n_rows), collapse = "/")), by = .(op, format, method, threads, ntraits, size)]
M7 <- mm(mr("repv7")); M8 <- mm(mr("repv8")); MT <- mm(mr("rep")[format == "tsv_gz"])
M <- merge(M7, M8, by = c("op", "format", "method", "threads", "ntraits", "size"), suffixes = c("_old", "_new"))
M[, change_pct := 100 * (wall_new - wall_old) / wall_old]
M <- merge(M, MT[, .(op, threads, ntraits, size, tsv_wall = wall, tsv_rows = n_rows)], by = c("op", "threads", "ntraits", "size"), all.x = TRUE)
M[, `:=`(speedup_vs_tsv_new = tsv_wall / wall_new, rows_match_tsv = n_rows_new == tsv_rows)]
fwrite(M, file.path(OUT, "multi_cpr_old_vs_new.csv")); cat("\n== multi-trait (cpr vs tsv.gz)\n"); print(M[, !c("n_old", "n_new")], digits = 4)
fp <- file.path(MR, c("repv7_multi_footprint.csv", "repv8_multi_footprint.csv")); if (all(file.exists(fp))) { cat("\n== multi footprint\n"); print(merge(fread(fp[1]), fread(fp[2]), by = c("format", "ntraits"), suffixes = c("_v7", "_v8"))) }
# same-node A/B control (e2e_ab.sbatch): C8 on both builds, interleaved in one job per replicate
AB <- file.path(E, Sys.getenv("AB_DIR", "results_ab"))
if (dir.exists(AB)) {
  ab <- rbindlist(lapply(list.dirs(AB, recursive = TRUE), function(d) { b <- basename(d)
    if (!grepl("^[0-9]+x[0-9]+_C8_v[0-9]+_r[0-9]$", b) || !file.exists(file.path(d, "stages.csv"))) return(NULL)
    s <- fread(file.path(d, "stages.csv"))
    data.table(size = sub("_.*", "", b), build = sub(".*_(v[0-9]+)_.*", "\\1", b), rep = basename(dirname(d)), wall = wall_of(file.path(d, "time.txt")),
               clump = s[stage == "clump", wall_s], clump_cpu = s[stage == "clump", cpu_s], extract = s[stage == "extract", wall_s], mr = s[stage == "mr", wall_s]) }))
  ABs <- ab[, .(cells = .N, wall_median = median(wall), wall_min = min(wall), wall_max = max(wall), clump = median(clump), clump_cpu = median(clump_cpu),
                extract = median(extract), mr = median(mr)), by = .(size, build)]
  setorder(ABs, size, build); fwrite(ABs, file.path(OUT, "e2e_ab_same_node.csv")); cat("\n== same-node A/B (C8)\n"); print(ABs, digits = 4)
}
