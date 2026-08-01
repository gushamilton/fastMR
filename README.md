# fastMR

`fastMR` is a small R package for exact summary-statistics Mendelian
randomization. It accepts the tidy columns used by TwoSampleMR and uses a
registered C++17 backend for the shared exposure/outcome grid that dominates
large MR scans.

The first release includes IVW (under-dispersion-corrected, fixed-effects, and
multiplicative random-effects variants), MR-Egger, MR-Egger bootstrap, simple,
weighted, and penalised weighted median, simple mode, weighted mode, unweighted
regression, sign concordance, Wald ratio, and basic multivariable IVW. The
default methods are exact; there is no approximate mode estimator in the
default API.

## Quick start

```r
library(fastMR)

dat <- data.frame(
  SNP = paste0("rs", 1:4),
  beta.exposure = c(.10, .12, .08, .15),
  beta.outcome = c(.02, .03, .01, .04),
  se.exposure = rep(.02, 4),
  se.outcome = rep(.01, 4),
  id.exposure = "exposure_1",
  id.outcome = "outcome_1"
)

fast_mr(dat, methods = c("ivw", "egger"), nboot = 100, seed = 20260801)
```

For a shared grid, matrices have one exposure/outcome per row and one SNP per
column. Results are exposure-major, outcome-minor, so the first outcome row is
returned for every method before the next outcome row.

```r
fast_mr_grid(
  exposure_beta, outcome_beta, exposure_se, outcome_se,
  methods = c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode"),
  nboot = 100, seed = 20260801, threads = 10
)
```

The C++ kernel copies R's column-major matrices once into contiguous
row-major side layouts. Bootstrap normal draws are shared across all grid
pairs, exposure denominators and outcome numerators are generated once per
grid side, and worker count is bounded by the number of pairs. If the build
toolchain has no OpenMP, the same kernel uses a serial or `std::thread`
fallback; no nested thread pools are created.

For grids containing only IVW variants (`ivw`, `ivw_fe`, and `ivw_mre`),
fastMR takes a batched BLAS path: the exposure/outcome cross-products for all
pairs are computed together, and tidy results are flattened in one pass. This
is the high-throughput path for large exposure-by-outcome IVW scans; its
compact native return layout avoids allocating one nested result list per
pair. The mixed-method path remains available when other estimators are
requested.

The exact simple/weighted mode kernel also reuses its FFT workspaces and plan
across bootstrap draws. This reduces allocation and setup cost while retaining
the native density-grid semantics.

Penalised weighted median follows TwoSampleMR's chi-square down-weighting and
two-stream bootstrap semantics, with `penk = 20` by default.

For a local preprocessing path, harmonise alleles without a network dependency
and then clump with either a supplied LD matrix or a local PLINK reference:

```r
harmonised <- fast_harmonise_data(exposure_dat, outcome_dat, action = 2)
independent <- fast_clump_data(harmonised, bfile = "reference/eur")
```

The PLINK path delegates LD calculation to the installed binary; the matrix
path is useful when the same LD panel is reused across many exposure grids.

## Parquet

Parquet is intentionally optional. With the optional `arrow` package installed:

```r
fast_mr_parquet("summary_stats.parquet", methods = "ivw")
dat <- fast_read_parquet("summary_stats.parquet")
```

Install Arrow only when needed: `install.packages("arrow")`. The core package
has no mandatory data.table, dplyr, RcppArmadillo, or database dependency.

## Results and compatibility

Results use the familiar columns `id.exposure`, `id.outcome`, `method`,
`nsnp`, `b`, `se`, and `pval`, with method-specific diagnostics such as `Q`,
`Q_df`, `Q_pval`, `intercept`, and `sigma` retained when available. Method
codes are accepted alongside common TwoSampleMR function names; inspect
`fastmr_method_registry()` for the complete registry.

`fastMR` ports the validated exact native-grid implementation developed in the
companion prototype and follows TwoSampleMR's public method semantics. It does
not bundle TwoSampleMR or its remote OpenGWAS/data-extraction stack.

## Validation and benchmarks

The durable 82-row IL6/CRP fixture is in `inst/extdata/`. Reproducible package
benchmarks and optimization-loop history are in `benchmarks/` and `outputs/`.
The benchmark separates cold startup, warm single-pair compute, shared-grid
compute, and optional Arrow read time; see [HANDOFF.md](HANDOFF.md) for the
latest exact commands, correctness gates, and measured results.

```sh
/opt/homebrew/bin/R CMD INSTALL --library=.local/Rlib .
/opt/homebrew/bin/Rscript -e 'testthat::test_dir("tests/testthat")'
/opt/homebrew/bin/R CMD check --no-manual --as-cran .
```
