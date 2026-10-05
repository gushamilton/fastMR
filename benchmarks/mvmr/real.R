# Real-data batched MVMR: UKB-PPP trans exposures x 2,940 outcomes conditional
# on platelet traits (inputs prepared from F-platelet; see prep_real.R).
# Usage: Rscript real.R <lib> <stores_dir> <inputs.rds> <outdir> <run> [threads]
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[[1]], "/user/work/fh6520/r_packages", .libPaths()))
suppressPackageStartupMessages(library(fastMR))
stores <- args[[2]]; inp <- readRDS(args[[3]]); outdir <- args[[4]]; run <- args[[5]]
threads <- if (length(args) >= 6) as.integer(args[[6]]) else 8L
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outcome_labels <- inp$excl$id.outcome
outcomes <- stats::setNames(file.path(stores, outcome_labels), outcome_labels)
exposure_labels <- names(inp$own)
exposures <- stats::setNames(file.path(stores, exposure_labels), exposure_labels)
cpu <- sub(".*: *", "", grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1])
cat(run, "| exposures", length(exposures), "| outcomes", length(outcomes), "| threads", threads,
    "|", cpu, "| nproc", parallel::detectCores(), "\n")
cov_all <- inp$covs
call <- switch(run,
  M1 = list(covariates = cov_all["PLT"], instruments = inp$m1, method = "mvmr",
            se_model = "multiplicative_floored"),
  M1_tsmr = list(covariates = cov_all["PLT"], instruments = inp$m1, method = "mvmr"),
  M1_shared = list(covariates = cov_all["PLT"], instruments = inp$m1, method = "mvmr",
                   se_model = "multiplicative_floored", weights = "shared",
                   shared_tolerance = 0.25),
  M1_greedy = list(covariates = cov_all["PLT"], instruments = inp$own, method = "mvmr",
                   covariate_instruments = list(PLT = inp$gsets$PLT), ld_pairs = inp$ld,
                   se_model = "multiplicative_floored"),
  R1 = list(covariates = cov_all["PLT"], instruments = inp$own, method = "residualised",
            covariate_instruments = list(PLT = inp$gsets$PLT)),
  R2 = list(covariates = cov_all[c("PLT", "MPV")], instruments = inp$own, method = "residualised",
            covariate_instruments = list(PLT = inp$gPM, MPV = character())),
  R3 = list(covariates = c(cov_all[c("PLT", "MPV")], list(GP6 = inp$gp6$frame)),
            instruments = inp$own, method = "residualised",
            covariate_instruments = list(PLT = inp$gPM, MPV = character(), GP6 = inp$gp6$instruments)),
  stop("unknown run"))
gc(reset = TRUE)
started <- proc.time()[["elapsed"]]
res <- do.call(fast_mvmr_compressed, c(list(
  exposure_files = exposures, outcome_files = outcomes, outcome_exclude = inp$excl,
  strict = FALSE, output_format = "matrix", threads = threads, io_threads = threads,
  covariate_estimates = TRUE), call))
wall <- proc.time()[["elapsed"]] - started
meta <- attr(res, "mvmr_input")
mem <- gc()
tm <- meta$timing
cat(sprintf("%s: io %.1f s | design %.1f s | estimator %.1f s | total %.1f s | wall %.1f s | max mem %.0f MB\n",
            run, tm$io_seconds, tm$design_seconds, tm$estimator_seconds, tm$total_seconds, wall,
            sum(mem[, ncol(mem)])))
print(summary(as.vector(res$nsnp)))
saveRDS(list(b = res$b, se = res$se, nsnp = res$nsnp, conditional_F = res$conditional_F,
             covariates = res$covariates, diagnostics = meta$diagnostics,
             instruments = if (run %in% c("M1_greedy")) meta$instruments else NULL,
             covariate_fit = meta$covariate_fit, timing = tm, counts = meta$counts[c("union_snps", "masked_cells")],
             shared_fits = meta$shared_fits, shared_deviation = meta$shared_deviation,
             threads = threads, cpu = cpu, wall = wall),
        file.path(outdir, paste0(run, ".rds")))
# Kernel-only re-timing on the in-memory panel is in kernel_real.R.
cat("done", run, "\n")
