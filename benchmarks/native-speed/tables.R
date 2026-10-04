# Rscript tables.R <results dir>: markdown tables from per-method-*.csv and batch8-*.csv.
d <- commandArgs(TRUE)[1]
pm <- do.call(rbind, lapply(list.files(d, "^per-method-.*csv$", full.names = TRUE), read.csv))
a <- aggregate(us_per_pair_elapsed ~ label + k + method, pm, min)
w <- reshape(a, idvar = c("k", "method"), timevar = "label", direction = "wide")
names(w) <- sub("us_per_pair_elapsed.", "", names(w))
w$speedup <- w$main / w$new
w <- w[order(w$k, match(w$method, c("wald_ratio","ivw","egger","weighted_median","weighted_mode","five"))), ]
cat("| k | method | main us/pair | PR us/pair | speed-up |\n|---|---|---|---|---|\n")
for (i in seq_len(nrow(w))) cat(sprintf("| %d | %s | %.0f | %.0f | %.2fx |\n", w$k[i], w$method[i], w$main[i], w$new[i], w$speedup[i]))
b <- do.call(rbind, lapply(list.files(d, "^batch8-.*csv$", full.names = TRUE), read.csv))
bb <- aggregate(wall_s ~ label + workload + pairs, b, min)
bw <- reshape(bb, idvar = c("workload", "pairs"), timevar = "label", direction = "wide")
names(bw) <- sub("wall_s.", "", names(bw))
cat("\n| workload | pairs | main | hull off, sync | hull off, overlap | hull on, sync | PR (hull + overlap) | PR speed-up |\n|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(bw))) cat(sprintf("| %s | %d | %.2f | %.2f | %.2f | %.2f | %.2f | %.2fx |\n", bw$workload[i], bw$pairs[i], bw$main[i], bw[["nohull-sync"]][i], bw[["nohull-overlap"]][i], bw[["hull-sync"]][i], bw$new[i], bw$main[i] / bw$new[i]))
