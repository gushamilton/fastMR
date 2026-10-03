#!/usr/bin/env Rscript
# Aggregate MR-compute results: summary CSVs + plots (SVG and PNG).
# usage: Rscript aggregate.R <results_dir> [<out_dir>] [dry]    (dry => reads *_dry files)
args <- commandArgs(TRUE); rd <- if (length(args)) args[1] else "."; od <- if (length(args) > 1) args[2] else file.path(rd, "summary")
dry <- length(args) > 2 && args[3] == "dry"; dir.create(od, FALSE, TRUE)
suppressPackageStartupMessages(library(ggplot2))
pat <- function(p) list.files(rd, pattern = sprintf("^%s_?rep[0-9]+%s\\.csv$", p, if (dry) "_dry" else ""), full.names = TRUE)
rf <- list.files(rd, pattern = if (dry) "^rep[0-9]+_dry\\.csv$" else "^rep[0-9]+\\.csv$", full.names = TRUE)
rd_all <- function(fs) if (length(fs)) do.call(rbind, lapply(fs, read.csv, stringsAsFactors = FALSE)) else NULL
d <- rd_all(rf); write.csv(d, file.path(od, "all_replicates.csv"), row.names = FALSE)
sfx <- if (dry) "_dry" else ""
ag <- rd_all(list.files(rd, pattern = sprintf("^agreement_rep[0-9]+%s\\.csv$", sfx), full.names = TRUE))
bs <- rd_all(list.files(rd, pattern = sprintf("^bootse_rep[0-9]+%s\\.csv$", sfx), full.names = TRUE))
d$tool <- ifelse(d$arm == "p2c1", "plink2 --clump (per exposure)", paste0(d$package, " x", d$threads))
d$ok <- d$status == "ok"

# ---- 1. timing summary (median, min-max over replicates; censored cells flagged, never extrapolated)
key <- c("scenario", "method_set", "arm", "tool", "size")
cells <- unique(d[, key])
summ <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  x <- merge(d, cells[i, ], by = key); ok <- x[x$ok, ]
  st <- table(x$status[!x$ok]); cens <- if (length(st)) paste(names(st), st, sep = " x", collapse = "; ") else ""
  data.frame(cells[i, ], n_ok = nrow(ok), n_total = nrow(x), wall_median = if (nrow(ok)) median(ok$wall_s) else NA,
    wall_min = if (nrow(ok)) min(ok$wall_s) else NA, wall_max = if (nrow(ok)) max(ok$wall_s) else NA,
    cpu_median = if (nrow(ok)) median(ok$cpu_s) else NA, rss_median_mb = if (nrow(ok)) median(ok$peak_rss_mb) else NA,
    censored = cens, stringsAsFactors = FALSE) }))
summ <- summ[order(summ$scenario, summ$method_set, summ$arm, summ$size), ]
write.csv(summ, file.path(od, "summary_times.csv"), row.names = FALSE)

# ---- 2. speedup table: 1 vs 1 core (tsmr1/fast1) and 8 vs 8 cores (tsmr8/fast8); fast8 only exists in bootstrap cells
sp <- function(a, b, label) {
  A <- summ[summ$arm == a, c("scenario", "method_set", "size", "wall_median", "censored")]; B <- summ[summ$arm == b, c("scenario", "method_set", "size", "wall_median")]
  names(A)[4:5] <- c("tsmr_s", "tsmr_censored"); names(B)[4] <- "fast_s"; m <- merge(A, B); m$comparison <- label
  m$speedup <- m$tsmr_s / m$fast_s; m }
spd <- rbind(sp("tsmr1", "fast1", "1 vs 1 core"), sp("tsmr8", "fast8", "8 vs 8 cores"), sp("tsmr8", "fast1", "TSMR x8 vs fastMR x1 (cross)"))
spd <- spd[order(spd$comparison, spd$scenario, spd$method_set, spd$size), ]
write.csv(spd, file.path(od, "speedup_table.csv"), row.names = FALSE)
# where TSMR is censored or skipped but fastMR ran, record the bound as ">" (no extrapolation)
cen <- summ[summ$censored != "" & grepl("TSMR|TwoSampleMR", summ$tool), ]; write.csv(cen, file.path(od, "censored_tsmr_cells.csv"), row.names = FALSE)

# ---- 3. agreement table
if (!is.null(ag)) {
  ag$arm_vs_ref <- paste0(ag$arm, " vs ", ag$ref_arm)
  g <- aggregate(cbind(max_abs_diff, median_abs_diff, n_mismatch, n_matched, n_ref, n_arm) ~ scenario + size + method_set + arm_vs_ref + method + metric, ag,
                 function(x) if (all(is.na(x))) NA else if (length(x) > 0) max(x, na.rm = TRUE) else NA, na.action = na.pass)
  names(g)[7:8] <- c("max_abs_diff_worst_rep", "median_abs_diff_worst_rep")
  write.csv(g[order(g$scenario, g$method_set, g$size, g$arm_vs_ref, g$method, g$metric), ], file.path(od, "agreement_table.csv"), row.names = FALSE)
  # headline checks: deterministic methods vs targets
  det <- g[g$metric %in% c("b", "se") & g$method %in% c("Inverse variance weighted", "MR Egger") & grepl("^fast. vs tsmr", g$arm_vs_ref), ]
  wm <- g[g$metric == "b" & g$method %in% c("Weighted median", "Simple mode", "Weighted mode") & grepl("^fast. vs tsmr", g$arm_vs_ref), ]
  cat("\nIVW/Egger b,se: worst max|diff| =", format(max(det$max_abs_diff_worst_rep, na.rm = TRUE)), "(target ~1e-12)\n")
  cat("WM/modes b: worst max|diff| =", format(max(wm$max_abs_diff_worst_rep, na.rm = TRUE)), "(target 0)\n")
}
# ---- 3b. bootstrap SE: fastMR-vs-TSMR relative SE differences against the TSMR seed A vs seed B null
if (!is.null(bs)) {
  bs$rel <- (bs$se_arm - bs$se_ref) / bs$se_ref
  bs$comparison <- ifelse(bs$arm == "tsmr1_seedB", "NULL: TSMR seed B vs seed A", paste0(bs$arm, " vs ", bs$ref_arm))
  q <- function(x) c(n = length(x), mean_rel = mean(x), sd_rel = sd(x), q025 = unname(quantile(x, .025)), q975 = unname(quantile(x, .975)), max_abs_rel = max(abs(x)))
  bsum <- do.call(rbind, lapply(split(bs, list(bs$scenario, bs$method_set, bs$method, bs$comparison), drop = TRUE), function(x)
    data.frame(x[1, c("scenario", "method_set", "method", "comparison")], t(q(x$rel)), row.names = NULL)))
  write.csv(bsum, file.path(od, "bootstrap_se_diff_summary.csv"), row.names = FALSE)
  write.csv(bs, file.path(od, "bootstrap_se_diff_rows.csv"), row.names = FALSE)
  pb <- ggplot(bs[bs$arm %in% c("fast1", "tsmr1_seedB") & bs$method_set == "default_k10", ], aes(rel, colour = comparison, fill = comparison)) +
    geom_density(alpha = .25) + facet_wrap(~method, scales = "free_y") + theme_bw() +
    labs(x = "relative SE difference (arm - TSMR seed A) / TSMR seed A", y = "density", title = "Bootstrap SE: fastMR vs TSMR, against TSMR seed A vs B null") + theme(legend.position = "bottom")
  for (e in c("svg", "png")) ggsave(file.path(od, paste0("bootstrap_se_null.", e)), pb, width = 9, height = 4, dpi = 200)
}

# ---- 4. plots: time vs size (log-log); censored = open triangle at the cap (actual cap hit) or open square (predicted, skipped)
CAP <- 1500
cens_pts <- d[!d$ok & d$status %in% c("> cap", "> cap (predicted)"), ]
cens_pts <- unique(cens_pts[, c("scenario", "method_set", "tool", "size", "status")]); cens_pts$wall_s <- rep(CAP, nrow(cens_pts))
okd <- summ[!is.na(summ$wall_median), ]
okd$grp <- paste(okd$tool)
plot_one <- function(sc, title, xlab, fname, sets = NULL) {
  x <- okd[okd$scenario == sc, ]; cp <- cens_pts[cens_pts$scenario == sc, ]; if (!nrow(x)) return(invisible())
  p <- ggplot(x, aes(size, wall_median, colour = tool, shape = tool, group = tool)) + geom_line() + geom_point(size = 2.2) +
    geom_errorbar(aes(ymin = wall_min, ymax = wall_max), width = 0, alpha = .5) +
    geom_hline(yintercept = CAP, linetype = "dashed", colour = "grey40") +
    scale_x_log10() + scale_y_log10() + facet_wrap(~method_set, scales = "free_x") + theme_bw() +
    labs(x = xlab, y = "wall time (s), median and range of replicates", title = title, shape = NULL, colour = NULL) + theme(legend.position = "bottom")
  if (nrow(cp)) p <- p + geom_point(data = cp, aes(size, wall_s, colour = tool, shape = status), inherit.aes = FALSE, size = 3, stroke = 1, fill = NA) +
    scale_shape_manual(values = c(setNames(c(16, 17, 15, 18, 3, 4)[seq_along(unique(x$tool))], unique(x$tool)), "> cap" = 2, "> cap (predicted)" = 0))
  for (e in c("svg", "png")) ggsave(file.path(od, paste0(fname, ".", e)), p, width = 10, height = 4.5, dpi = 200)
}
plot_one("many", "Many pairs: time vs number of pairs (dashed = 25 min cap)", "number of exposure-outcome pairs", "time_vs_pairs_many")
plot_one("single", "Single pair: time vs SNPs", "SNPs in the pair", "time_vs_snps_single")
plot_one("harmonise", "Harmonisation (single pair)", "rows", "time_harmonise")
plot_one("steiger", "Steiger filtering", "rows", "time_steiger")
plot_one("clump", "LD clumping", "number of exposures", "time_clump")
for (sc in c("heterogeneity", "pleiotropy", "harmonise_many")) plot_one(sc, sc, "pairs", paste0("time_", sc))
cat("wrote", od, "\n")
