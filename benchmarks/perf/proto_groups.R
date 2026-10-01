# Prototype: replace the O(groups x rows) which() scan in fast_mr() and
# fastmr_diagnostic_groups() with one split(); check identical() output.
.libPaths(c("/Users/fh6520/projects/fastMR/.local/Rlib", .libPaths())); suppressMessages(library(fastMR))
src <- commandArgs(TRUE)[1]; if (is.na(src)) src <- "fastmr_perf_bench.R"
ns <- asNamespace("fastMR")
fm2 <- fast_mr
b <- deparse(body(fm2))
i <- grep("group_index <- which", b)
b[i] <- "    group_index <- group_rows[[i]]"; b <- b[-(i + 1L)]
j <- grep("rows <- vector\\(\"list\", nrow\\(groups\\)\\)", b)
b <- append(b, "  group_rows <- unname(split(seq_len(n), factor(paste(id.exp, id.out, sep = \"\\r\"), levels = unique(paste(id.exp, id.out, sep = \"\\r\")))))", after = j)
body(fm2) <- parse(text = b)[[1]]; environment(fm2) <- ns
dg2 <- get("fastmr_diagnostic_groups", ns)
b <- deparse(body(dg2)); i <- grep("index <- which", b)
b[i] <- "    index <- group_rows[[i]]"; b <- b[-(i + 1L)]
j <- grep("groups <- vector\\(\"list\", nrow\\(pairs\\)\\)", b)
b <- append(b, "  group_rows <- unname(split(seq_len(n), factor(paste(id.exp, id.out, sep = \"\\r\"), levels = unique(paste(id.exp, id.out, sep = \"\\r\")))))", after = j)
body(dg2) <- parse(text = b)[[1]]; environment(dg2) <- ns
# reuse the generator from the benchmark
e <- new.env(); lines <- readLines(src); s <- grep("^gen_long <-", lines); t <- grep("^gen_csr <-", lines)
eval(parse(text = lines[s:(t - 1)]), e)
for (G in c(4000L, 20000L)) {
  d <- e$gen_long(G, 20L)
  t1 <- system.time(a <- fast_mr(d, methods = "ivw", nboot = 0))[["elapsed"]]
  t2 <- system.time(b2 <- fm2(d, methods = "ivw", nboot = 0))[["elapsed"]]
  cat(sprintf("fast_mr ivw G=%d: current %.2fs, split-fix %.2fs, identical=%s\n", G, t1, t2, identical(a, b2)))
}
G <- 400L; d <- e$gen_long(G, 20L); m <- c("ivw","egger","weighted_median","weighted_mode")
t1 <- system.time(a <- fast_mr(d, methods = m, nboot = 100, seed = 1))[["elapsed"]]
t2 <- system.time(b2 <- fm2(d, methods = m, nboot = 100, seed = 1))[["elapsed"]]
cat(sprintf("fast_mr mixed seeded G=%d: current %.2fs, split-fix %.2fs, identical=%s\n", G, t1, t2, identical(a, b2)))
G <- 20000L; d <- e$gen_long(G, 10L, steiger = TRUE)
t1 <- system.time(a <- fast_mr_steiger_filtering(d))[["elapsed"]]
orig <- get("fastmr_diagnostic_groups", ns); unlockBinding("fastmr_diagnostic_groups", ns)
assign("fastmr_diagnostic_groups", dg2, ns)
t2 <- system.time(b2 <- fast_mr_steiger_filtering(d))[["elapsed"]]
Rprof("steiger_fixed.prof", interval = 0.005); invisible(fast_mr_steiger_filtering(d)); Rprof(NULL)
assign("fastmr_diagnostic_groups", orig, ns)
cat(sprintf("steiger_filtering pairs=%d: current %.2fs, split-fix %.2fs, identical=%s\n", G, t1, t2, identical(a, b2)))
Rprof("fastmr_fixed.prof", interval = 0.005); invisible(fm2(e$gen_long(20000L, 20L), methods = "ivw", nboot = 0)); Rprof(NULL)
print(head(summaryRprof("fastmr_fixed.prof")$by.self, 10)); print(head(summaryRprof("steiger_fixed.prof")$by.self, 10))
