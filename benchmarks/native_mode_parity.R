args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else getwd()
.libPaths(c(file.path(root, ".local", "Rlib"),
            "/Users/fergushamilton/projects/twosamplemr-fast/.local/Rlib", .libPaths()))
suppressPackageStartupMessages({ library(fastMR); library(TwoSampleMR) })
set.seed(20260815)
method_specs <- data.frame(
  fast = c("simple_mode", "weighted_mode"),
  native = c("mr_simple_mode", "mr_weighted_mode"),
  stringsAsFactors = FALSE
)
rows <- vector("list", 24L)
k <- 0L
for (panel in seq_len(12L)) {
  n <- sample(6:20, 1L)
  dat <- data.frame(
    SNP = paste0("rs", seq_len(n)),
    beta.exposure = rnorm(n, 0, 0.12),
    beta.outcome = rnorm(n, 0, 0.12),
    se.exposure = runif(n, 0.01, 0.08),
    se.outcome = runif(n, 0.01, 0.08),
    id.exposure = "E", id.outcome = "O", exposure = "E", outcome = "O",
    mr_keep = TRUE, stringsAsFactors = FALSE
  )
  for (method_index in seq_len(nrow(method_specs))) {
    k <- k + 1L
    seed <- 20260900L + panel
    fast <- fast_mr(dat, methods = method_specs$fast[method_index],
                    nboot = 100L, seed = seed)
    set.seed(seed)
    native <- suppressMessages(suppressWarnings(
      TwoSampleMR::mr(dat, method_list = method_specs$native[method_index],
                      parameters = modifyList(TwoSampleMR::default_parameters(),
                                               list(nboot = 100L)))))
    rows[[k]] <- data.frame(
      panel = panel, nsnp = n, method = method_specs$fast[method_index],
      abs_beta_delta = abs(fast$b - native$b),
      abs_se_delta = abs(fast$se - native$se),
      abs_pval_delta = abs(fast$pval - native$pval)
    )
  }
}
parity <- do.call(rbind, rows)
set.seed(20260816)
ne <- 3L; no <- 3L; nsnp <- 18L
x <- matrix(rnorm(ne * nsnp), nrow = ne)
y <- matrix(rnorm(no * nsnp), nrow = no)
sx <- matrix(runif(ne * nsnp, 0.01, 0.08), nrow = ne)
sy <- matrix(runif(no * nsnp, 0.01, 0.08), nrow = no)
serial <- fast_mr_grid(x, y, sx, sy, methods = c("simple_mode", "weighted_mode"),
                        nboot = 20L, seed = 20260816, threads = 1L)
parallel <- fast_mr_grid(x, y, sx, sy, methods = c("simple_mode", "weighted_mode"),
                         nboot = 20L, seed = 20260816, threads = 5L)
grid_delta <- max(vapply(c("b", "se", "pval"), function(name) {
  max(abs(serial[[name]] - parallel[[name]]), na.rm = TRUE)
}, numeric(1)))
summary <- data.frame(
  panels = nrow(parity), max_abs_beta_delta = max(parity$abs_beta_delta),
  max_abs_se_delta = max(parity$abs_se_delta), max_abs_pval_delta = max(parity$abs_pval_delta),
  max_grid_serial_vs_5thread_delta = grid_delta,
  parity_gate = max(parity$abs_beta_delta, parity$abs_se_delta, parity$abs_pval_delta) <= 1e-12,
  thread_gate = grid_delta == 0
)
dir.create(file.path(root, "outputs"), showWarnings = FALSE, recursive = TRUE)
write.csv(parity, file.path(root, "outputs", "native_mode_parity.csv"), row.names = FALSE)
write.csv(summary, file.path(root, "outputs", "native_mode_parity_summary.csv"), row.names = FALSE)
cat("panels", summary$panels, "max beta", summary$max_abs_beta_delta,
    "max se", summary$max_abs_se_delta, "max p", summary$max_abs_pval_delta,
    "max grid thread delta", summary$max_grid_serial_vs_5thread_delta,
    "parity gate", summary$parity_gate, "thread gate", summary$thread_gate, "\n")
