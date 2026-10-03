# End-to-end summary: median (min-max) wall/RSS per arm and size, and median per-stage wall. VCF arm (B) from the LP-fixed re-run.
suppressMessages(library(data.table))
E <- "/user/work/fh6520/showcase/e2e"; R <- file.path(E, "results")
cells <- rbindlist(lapply(1:3, function(r) rbind(fread(file.path(R, sprintf("rep%d/cells.csv", r)))[arm != "B"],
                                                 fread(file.path(E, sprintf("results_B/rep%d/cells.csv", r))))))
tot <- cells[, .(n = .N, status = paste(unique(status), collapse = ";"), wall_median = median(wall_s), wall_min = min(wall_s),
                 wall_max = max(wall_s), cpu_median = median(user_s + sys_s), rss_mb_median = median(maxrss_kb) / 1024), by = .(size, arm, threads)]
setorder(tot, size, wall_median)
st <- rbindlist(lapply(1:3, function(r) rbindlist(lapply(c("1x1", "10x10"), function(s) rbindlist(lapply(c("A", "A8", "B", "C1", "C8"), function(a) {
  f <- file.path(R, sprintf("rep%d/%s_%s/stages.csv", r, s, a)); if (file.exists(f)) fread(f)[, `:=`(rep = r, size = s, arm = a)] }), fill = TRUE)), fill = TRUE)), fill = TRUE)
stg <- st[, .(wall_median = median(wall_s)), by = .(size, arm, stage)]
wide <- dcast(stg, size + arm ~ stage, value.var = "wall_median")
conv <- fread(file.path(E, "../prep/traits/conversion_times.csv"))
fwrite(tot, file.path(R, "e2e_summary_total.csv")); fwrite(wide, file.path(R, "e2e_summary_stages.csv"))
print(tot); print(wide)
cat("\nconversion per trait (old build) median s: tsv", median(conv$tsv_write_s), " vcf", median(conv$vcf_total_s), " cpr", median(conv$cpr_total_s), "\n")
