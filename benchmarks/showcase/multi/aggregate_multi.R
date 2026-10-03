# Aggregate the multi-trait suite: median (min-max) of 3 replicates per cell; VCF rows from the LP-fixed re-run.
suppressMessages(library(data.table))
RES <- "/user/work/fh6520/showcase/storage/multi/results"
rd <- function(tag) rbindlist(lapply(1:3, function(r) { f <- file.path(RES, sprintf("%s%d_multi.csv", tag, r)); if (file.exists(f)) fread(f) }), fill = TRUE)
d <- rbind(rd("rep")[format != "vcf"], rd("repvcf"), fill = TRUE)
agg <- d[, .(n_ok = sum(status == "ok"), n_rep = .N, status = if (all(status == "ok")) "ok" else paste(unique(status), collapse = ";"),
             wall_median = median(wall_s[status == "ok"]), wall_min = min(wall_s[status == "ok"]), wall_max = max(wall_s[status == "ok"]),
             seconds_median = median(seconds[status == "ok"]), cpu_median = median((user_s + sys_s)[status == "ok"]),
             maxrss_mb_median = median(maxrss_mb[status == "ok"]), n_rows = paste(unique(n_rows), collapse = "/")),
         by = .(op, format, method, threads, ntraits, size)]
fwrite(agg, file.path(RES, "multi_summary_all.csv"))
# best method per tool (VCF: -R vs -T), and for cpr keep batch and loop separately
best <- agg[n_ok == n_rep][, tool := fifelse(format == "cpr", paste0("cpr_", method), format)]
best <- best[order(wall_median), .SD[1], by = .(op, tool, threads, ntraits, size)]
setorder(best, op, size, ntraits, threads, tool)
fwrite(best, file.path(RES, "multi_summary_best.csv"))
fp <- rbind(fread(file.path(RES, "rep_multi_footprint.csv"))[format != "vcf"], fread(file.path(RES, "repvcf_multi_footprint.csv")))
fp[, `:=`(projected_1000_GB = bytes_per_trait * 1000 / 1e9, projected_3000_GB = bytes_per_trait * 3000 / 1e9)]
fwrite(fp, file.path(RES, "multi_summary_footprint.csv"))
# row agreement across formats (instrument / extracted row counts should match)
print(agg[, .(n_rows = paste(unique(n_rows), collapse = " | ")), by = .(op, ntraits, size)])
print(best[ntraits == 20, .(op, size, tool, threads, wall_median, wall_min, wall_max, maxrss_mb_median)])
