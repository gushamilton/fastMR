# Per-stage scaling model for the end-to-end study: time = a + b*units, units by stage
# (per exposure: select + clump; per outcome: extract; per pair: MR + harmonise + Steiger; plus fixed overhead),
# i.e. a + b*n + c*n^2 overall for an n x n study. Fit on 1x1 and 10x10 (medians of 3), then predict larger
# studies and compare with the measured 25x25 (and 50x50) where they exist. A model, not a measurement.
# Env: SC_FIT  = assembled replicate root with <size>_<arm>/{stages.csv,time.txt} for 1x1 and 10x10 (rep1..3)
#      SC_OBS  = space-separated list of further roots with measured sizes (e.g. final25 dirs, 50x50 dirs)
#      SC_TAG  = suffix for the output CSVs (scaling_fit_<tag>.csv, scaling_predictions_<tag>.csv, scaling_check_<tag>.csv)
suppressMessages(library(data.table))
E2E <- "/user/work/fh6520/showcase/e2e"
FIT <- Sys.getenv("SC_FIT", file.path(E2E, "final10")); OBS <- strsplit(Sys.getenv("SC_OBS", file.path(E2E, "final25")), " ")[[1]]
TAG <- Sys.getenv("SC_TAG", "v7")
wall_of <- function(tf) { x <- sub(".*: ", "", grep("Elapsed \\(wall", readLines(tf), value = TRUE)); p <- as.numeric(strsplit(x, ":")[[1]])
  sum(p * 60^(rev(seq_along(p)) - 1)) }
rd <- function(root) rbindlist(lapply(list.dirs(root, recursive = TRUE), function(d) {
  b <- basename(d); if (!grepl("^[0-9]+x[0-9]+_[A-Z0-9]+$", b) || !file.exists(file.path(d, "stages.csv")) || !file.exists(file.path(d, "time.txt"))) return(NULL)
  s <- fread(file.path(d, "stages.csv"))
  rbind(s, data.table(stage = "wall", wall_s = wall_of(file.path(d, "time.txt")), cpu_s = NA), fill = TRUE)[
    , `:=`(size = sub("_.*", "", b), arm = sub(".*_", "", b), rep = basename(dirname(d)))] }), fill = TRUE)
widen <- function(d) {
  m <- d[, .(t = median(wall_s)), by = .(arm, size, stage)]
  w <- dcast(m, arm + size ~ stage, value.var = "t", fill = 0); w[, n := as.integer(sub("x.*", "", size))]
  for (cn in c("select", "clump", "extract", "mr", "harmonise", "steiger", "write")) if (!cn %in% names(w)) w[, (cn) := 0]
  w[, `:=`(expo = select + clump, outc = extract, pair = mr + harmonise + steiger)]
  w[, fixed := pmax(0, wall - expo - outc - pair - write)]; w }
w <- widen(rd(FIT))[n %in% c(1L, 10L)]
w <- w[arm %in% w[, .N, by = arm][N == 2, arm]]   # arms measured at both fit sizes (A8 has no 1x1 cell)
fit <- w[, {
  o <- order(n); {
  lin <- function(y, u) { y <- y[o]; u <- u[o]; b <- (y[2] - y[1]) / (u[2] - u[1]); c(a = y[1] - b * u[1], b = b) }
  e <- lin(expo, n); oc <- lin(outc, n); pr <- lin(pair, n^2)
  list(fixed = mean(fixed + write), expo_a = e[["a"]], expo_b = e[["b"]], outc_a = oc[["a"]], outc_b = oc[["b"]],
       pair_a = pr[["a"]], pair_b = pr[["b"]]) } }, by = arm]
pred <- function(f, n) with(f, fixed + expo_a + expo_b * n + outc_a + outc_b * n + pair_a + pair_b * n^2)
ns <- c(1, 2, 5, 10, 25, 50, 100, 250, 1000)
P <- fit[, .(n = ns, predicted_s = pred(.SD, ns)), by = arm]
cat("== fitted per-stage coefficients (seconds; b per exposure, per outcome, per pair)\n"); print(fit, digits = 4)
cat("\n== measured medians used for the fit\n"); print(w[, .(arm, n, wall, expo, outc, pair, fixed)], digits = 4)
cat("\n== MODEL predictions (not measurements)\n"); print(dcast(P, n ~ arm, value.var = "predicted_s"), digits = 4)
ob <- rbindlist(lapply(OBS[nzchar(OBS)], function(r) { x <- rd(r); if (nrow(x)) x[stage == "wall", .(observed_s = median(wall_s), reps = .N), by = .(arm, size)] }))
if (nrow(ob)) {
  ob[, n := as.integer(sub("x.*", "", size))]
  chk <- merge(ob, P, by = c("arm", "n"))[, pct_error := 100 * (predicted_s - observed_s) / observed_s]
  cat("\n== prediction vs measured\n"); print(chk, digits = 4); fwrite(chk, file.path(E2E, sprintf("scaling_check_%s.csv", TAG)))
}
fwrite(fit, file.path(E2E, sprintf("scaling_fit_%s.csv", TAG))); fwrite(P, file.path(E2E, sprintf("scaling_predictions_%s.csv", TAG)))
cat("\nprediction written", format(Sys.time()), "\n")
