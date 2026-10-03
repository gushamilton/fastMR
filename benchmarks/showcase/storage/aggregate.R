#!/usr/bin/env Rscript
# Aggregate per-replicate CSVs into summary CSVs and plots (SVG + PNG). All timings are warm.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- Sys.getenv("STORAGE_ROOT", "/user/work/fh6520/showcase/storage")
TAG <- Sys.getenv("RESULT_TAG", "rep")
RES <- file.path(ROOT, "results"); PL <- file.path(ROOT, "plots"); dir.create(PL, FALSE, TRUE)
files <- list.files(RES, pattern = sprintf("^%s[0-9]+_ops\\.csv$", TAG), full.names = TRUE)
stopifnot(length(files) > 0)
ops <- rbindlist(lapply(files, fread, colClasses = list(character = c("note", "status", "rep"))), fill = TRUE)
if (!"load_s" %in% names(ops)) ops[, load_s := NA_real_]
ops[, ok := status == "ok"]
nrep <- uniqueN(ops$rep)
key <- c("op", "format", "method", "threads", "size")
med <- function(x) if (length(x)) as.numeric(median(x)) else NA_real_
summ <- ops[, .(n_ok = sum(ok), n_rep = .N, status = if (all(ok)) "ok" else paste(unique(status[!ok]), collapse = "; "),
  seconds_median = med(seconds[ok]), seconds_min = suppressWarnings(min(seconds[ok])), seconds_max = suppressWarnings(max(seconds[ok])),
  wall_median = med(wall_s[ok]), load_s_median = med(load_s[ok]), cpu_median = med((user_s + sys_s)[ok]), maxrss_mb_median = med(maxrss_mb[ok]),
  perq_median_s = med(perq_median_s[ok]), n_rows = med(n_rows[ok]), timing = "warm"), by = key]
summ[!is.finite(seconds_min), c("seconds_min", "seconds_max") := NA_real_]
fwrite(summ, file.path(RES, sprintf("%s_summary_all_ops.csv", TAG)))
best <- summ[n_ok > 0, .SD[which.min(seconds_median)], by = .(op, format, threads, size)]
capped <- summ[n_ok == 0, .(method = paste(method, collapse = ","), seconds_median = NA_real_, status = status[1]), by = .(op, format, threads, size)]
best <- rbind(best, capped[!best, on = .(op, format, threads, size)], fill = TRUE)
fwrite(best, file.path(RES, sprintf("%s_summary_best_method.csv", TAG)))

sizes <- fread(file.path(RES, "sizes.csv"))
sizes[, bytes_per_variant := total_bytes / 10e6]
pw <- fread(file.path(RES, "prep_writes.csv"))
wr <- summ[op == "write", .(format, write_seconds_median = seconds_median, write_wall_median = wall_median, write_rss_mb = maxrss_mb_median)]
tab <- merge(sizes, wr, by = "format", all.x = TRUE)
tab <- merge(tab, pw[format == "vcf", .(format, write_seconds_median_prep = seconds)], by = "format", all.x = TRUE)
tab[format == "vcf", write_seconds_median := write_seconds_median_prep][, write_seconds_median_prep := NULL]
fwrite(tab, file.path(RES, sprintf("%s_summary_size_write.csv", TAG)))
fp <- CJ(format = sizes$format, traits = c(1000, 3000, 10000))[sizes, on = "format"]
fp[, projected_TB := total_bytes * traits / 1e12]
fp[, projected_GB := total_bytes * traits / 1e9]
fwrite(fp[, .(format, traits, bytes_per_variant, total_bytes_r1 = total_bytes, projected_GB, projected_TB)],
       file.path(RES, sprintf("%s_summary_footprint.csv", TAG)))

lab <- c(tsv_gz = "TSV.gz", tsv_bgzip_tabix = "bgzip TSV + tabix", vcf = "GWAS-VCF + tabix",
         parquet = "Parquet (zstd)", qparquet = "Parquet quantised f32", cpr = "CompreSSoR .cpr")
cols <- c("TSV.gz" = "#8c8c8c", "bgzip TSV + tabix" = "#b59a6a", "GWAS-VCF + tabix" = "#c4713b",
          "Parquet (zstd)" = "#3b78b0", "Parquet quantised f32" = "#7aa8d1", "CompreSSoR .cpr" = "#2a9d6a")
lv <- unname(lab)
fl <- function(f) factor(lab[f], levels = lv)
th <- theme_minimal(base_size = 12) + theme(legend.position = "bottom", plot.title.position = "plot")
save2 <- function(p, name, w, h) {
  ggsave(file.path(PL, paste0(TAG, "_", name, ".svg")), p, width = w, height = h, device = svglite::svglite)
  ggsave(file.path(PL, paste0(TAG, "_", name, ".png")), p, width = w, height = h, dpi = 200)
}
# 1. Pareto: size vs best full-read time
fr <- best[op == "fullread" & !is.na(seconds_median)]
fr <- merge(fr, sizes[, .(format, total_bytes)], by = "format")
fr[, `:=`(fmt = fl(format), threads_f = factor(threads, labels = paste(sort(unique(threads)), "thread(s)")))]
p1 <- ggplot(fr, aes(total_bytes / 1e6, seconds_median, colour = fmt, shape = threads_f)) +
  geom_point(size = 3.5) + scale_colour_manual(values = cols, name = NULL) + scale_shape_discrete(name = NULL) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "On-disk size, MB (log)", y = "Full read, s (log, warm, best method)",
       title = "Size versus full-read time, R1 (10M variants)") + th
save2(p1, "pareto_size_vs_fullread", 7.5, 5)
# 2. Region and lookup
rl <- best[op %in% c("region", "lookup")]
rl[, task := fifelse(op == "region", "Region: 100 x +/-500 kb (total)", sprintf("Lookup: %s keys", format(size, big.mark = ",")))]
rl[, task := factor(task, levels = unique(task[order(op, size)]))]
rl[, `:=`(fmt = fl(format), threads_f = factor(threads, labels = paste(sort(unique(threads)), "thread(s)")))]
rl[, capped := is.na(seconds_median)]
rl[capped == TRUE, seconds_median := 600]
p2 <- ggplot(rl, aes(fmt, seconds_median, colour = fmt, shape = threads_f)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.5)) +
  geom_text(data = rl[capped == TRUE], aes(label = "> cap"), position = position_dodge(width = 0.5), vjust = -1, size = 3, show.legend = FALSE) +
  facet_wrap(~task, scales = "free_y", nrow = 1) + scale_y_log10() +
  scale_colour_manual(values = cols, guide = "none") + scale_shape_discrete(name = NULL) +
  labs(x = NULL, y = "Total time, s (log, warm, best method)",
       title = "Region and variant-lookup latency, R1") + th + theme(axis.text.x = element_text(angle = 40, hjust = 1))
save2(p2, "region_lookup_bars", 13, 5.5)
# 3. Projected footprint
fp[, `:=`(fmt = fl(format), traits_f = factor(traits, labels = paste(format(sort(unique(traits)), big.mark = ","), "traits")))]
p3 <- ggplot(fp, aes(traits_f, projected_TB, fill = fmt)) + geom_col(position = "dodge") +
  scale_fill_manual(values = cols, name = NULL) + scale_y_log10() +
  labs(x = NULL, y = "Projected disk footprint, TB (log)", title = "Projected footprint from R1 bytes per trait") + th
save2(p3, "footprint_projection", 8, 5)
cat("replicates:", nrep, "\n"); print(tab); print(best[op == "fullread"])
