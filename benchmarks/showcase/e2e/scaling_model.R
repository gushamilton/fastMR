# Per-stage scaling model for the end-to-end study: time = a + b*units, units by stage.
# Fit on 1x1 and 10x10 (medians of 3), then predict larger studies. A model, not a measurement.
suppressMessages(library(data.table))
E2E <- "/user/work/fh6520/showcase/e2e"
rd <- function(dir, arm, size) rbindlist(lapply(1:3, function(r) {
  f <- file.path(dir, sprintf("rep%d/%s_%s/stages.csv", r, size, arm)); c <- file.path(dir, sprintf("rep%d/cells.csv", r))
  if (!file.exists(f)) return(NULL)
  s <- fread(f)[, rep := r]; aa <- arm; ss <- size; w <- fread(c)[get("arm") == aa & get("size") == ss, wall_s]
  rbind(s, data.table(stage = "wall", wall_s = w, cpu_s = NA, rep = r), fill = TRUE) }))
src <- list(A = file.path(E2E, "results"), A8 = file.path(E2E, "results"), B = file.path(E2E, "results_B"),
            C1 = file.path(E2E, "results_C5"), C8 = file.path(E2E, "results_C5"))
d <- rbindlist(lapply(names(src), function(a) rbindlist(lapply(c("1x1", "10x10"), function(s) {
  x <- rd(src[[a]], a, s); if (is.null(x) || !nrow(x)) return(NULL); x[, `:=`(arm = a, size = s)] }))))
m <- d[, .(t = median(wall_s)), by = .(arm, size, stage)]
w <- dcast(m, arm + size ~ stage, value.var = "t", fill = 0)
w[, n := as.integer(sub("x.*", "", size))]
# Stage groups: per-exposure work (select + clump), per-outcome work (extract), per-pair work (MR [+ harmonise, Steiger]),
# and fixed overhead = wall minus the timed stages.
for (cn in c("select", "harmonise", "steiger", "write", "mr_total_incl_extract")) if (!cn %in% names(w)) w[, (cn) := 0]
w[, `:=`(expo = select + clump, outc = extract, pair = mr + harmonise + steiger)]
w[, fixed := pmax(0, wall - expo - outc - pair - write)]
fit <- w[, {
  o <- order(n); n1 <- n[o][1]; n2 <- n[o][2]
  lin <- function(y, u) { y <- y[o]; u <- u[o]; b <- (y[2] - y[1]) / (u[2] - u[1]); c(a = y[1] - b * u[1], b = b) }
  e <- lin(expo, n); oc <- lin(outc, n); pr <- lin(pair, n^2)
  list(fixed = mean(fixed + write), expo_a = e[["a"]], expo_b = e[["b"]], outc_a = oc[["a"]], outc_b = oc[["b"]],
       pair_a = pr[["a"]], pair_b = pr[["b"]])
}, by = arm]
pred <- function(f, n) with(f, fixed + expo_a + expo_b * n + outc_a + outc_b * n + pair_a + pair_b * n^2)
ns <- c(1, 2, 5, 10, 25, 50, 100, 250, 1000)
P <- fit[, .(n = ns, predicted_s = pred(.SD, ns)), by = arm]
cat("== fitted per-stage coefficients (seconds; b per exposure, per outcome, per pair)\n"); print(fit, digits = 3)
cat("\n== measured medians\n"); print(w[, .(arm, n, wall, expo, outc, pair, fixed)], digits = 3)
cat("\n== MODEL predictions (not measurements)\n"); print(dcast(P, n ~ arm, value.var = "predicted_s"), digits = 3)
fwrite(fit, file.path(E2E, "scaling_fit.csv")); fwrite(P, file.path(E2E, "scaling_predictions.csv"))
cat("\nprediction written", format(Sys.time()), "\n")
