# fastMR

<p align="center">
  <img src="inst/fastMR-logo.png" alt="fastMR logo" width="180">
</p>

<p align="center"><strong>Exact, compiled summary-statistics Mendelian randomization for large R workflows.</strong></p>

[![R-CMD-check](https://github.com/gushamilton/fastMR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gushamilton/fastMR/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

`fastMR` is a small R package for doing TwoSampleMR-style summary-statistics
Mendelian randomization quickly when many exposures and outcomes must be
scanned together. It keeps the familiar tidy R interface, but moves the
repeated numerical work into a registered C++17 backend and exposes a shared
matrix-grid path for large scans.

## Why fastMR?

The usual R workflow is excellent for individual analyses, but a 50 × 50
exposure/outcome grid repeats data-frame construction, grouping, validation,
and R-level function calls thousands of times. `fastMR` prepares the matrices
once, reuses shared work across pairs, and returns the same tidy result shape.

It is intended as a fast local computation layer, not as a replacement for
OpenGWAS extraction, study metadata, or the broader TwoSampleMR ecosystem.

## Installation

Install the development release from GitHub with:

```r
install.packages("remotes")
remotes::install_github("gushamilton/fastMR")
```

For a local checkout:

```sh
R CMD INSTALL .
```

The core package requires R (≥ 4.1), `Rcpp`, and a C++17 compiler. Arrow is
optional and is only needed for Parquet input.

## Quick start

`fast_mr()` accepts the standard summary-statistics columns used by
TwoSampleMR:

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

For a grid, put exposures and outcomes in rows and SNPs in columns:

```r
fast_mr_grid(
  exposure_beta, outcome_beta, exposure_se, outcome_se,
  methods = c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode"),
  nboot = 100, seed = 20260801, threads = 10
)
```

The grid result is exposure-major and outcome-minor, with one tidy row per
pair and method.

## Methods

The method registry includes:

- IVW, fixed-effects IVW, and multiplicative random-effects IVW
- MR-Egger and seeded MR-Egger bootstrap
- simple, weighted, and penalised weighted median
- simple and weighted mode
- unweighted regression, sign concordance, and Wald ratio
- basic multivariable IVW

Use `fastmr_method_registry()` for method codes and descriptions. Results use
the familiar `id.exposure`, `id.outcome`, `method`, `nsnp`, `b`, `se`, and
`pval` columns, with method-specific diagnostics retained where applicable.

## Local preprocessing

The package also provides dependency-light local preprocessing:

```r
harmonised <- fast_harmonise_data(exposure_dat, outcome_dat, action = 2)
independent <- fast_clump_data(harmonised, bfile = "reference/eur")
```

Harmonisation supports the native TwoSampleMR actions 1, 2, and 3, strand and
complement handling, palindromic-frequency checks, indel recoding, incomplete
allele information, and outcome-specific action vectors. Clumping can use an
in-memory LD matrix or delegate LD calculation to a locally installed PLINK
binary. Parquet readers are available through optional `arrow` support:

```r
fast_mr_parquet("summary_stats.parquet", methods = "ivw")
```

## Benchmarks

Final Mac mini audit: native `TwoSampleMR` 0.7.9 comparator, five methods,
`nboot = 100`, and 10 fastMR threads for the grid run.

| Workload | fastMR | Native TwoSampleMR | Speedup |
|---|---:|---:|---:|
| Local BMI → CRP, 104 SNPs, median of 5 runs | 0.005 s | 0.127 s | **25.4×** |
| Independent simulation, 1 × 1, 400 SNPs | 0.008 s | 0.165 s | **20.6×** |
| Independent simulation, 50 × 50, 2,500 pairs, 400 SNPs | 2.030 s | 391.764 s | **193.0×** |

The 50 × 50 run processed about 1,232 pairs/s with fastMR versus 6.4
pairs/s natively. These are observed timings on one Mac mini, not a promise
for every machine or workload; the largest gains come from the shared-grid
implementation and reduced R-level allocation, with threading contributing
additional benefit for bootstrap-heavy methods.

The final tables and reproducible audit script are in
[`outputs/final_package_audit/`](outputs/final_package_audit/) and
[`benchmarks/final_package_audit.R`](benchmarks/final_package_audit.R).

### Adversarial 25 × 25 proof run

The package was also checked on a deliberately different, heavy-tailed and
heteroskedastic simulation: 25 exposures × 25 outcomes, 82 SNPs, all five
main methods, and `nboot = 1,000`. The run produced all 625 × 5 = 3,125
expected rows, repeated exactly across 1-thread and 10-thread fastMR runs,
and compared every pair with native TwoSampleMR 0.7.9. It also exercised
one- and two-SNP inputs, invalid/non-finite inputs, and deliberately flipped
alleles through the local harmoniser. The measured result is recorded in
[`outputs/adversarial_25x25_nboot1000/`](outputs/adversarial_25x25_nboot1000/)
and can be reproduced with
[`benchmarks/adversarial_25x25_nboot1000.R`](benchmarks/adversarial_25x25_nboot1000.R).

Measured on the Mac mini: fastMR took 3.630 s at 10 threads and 20.177 s at
one thread; native TwoSampleMR took 824.585 s. That is a 227.16× native
speedup, with a 5.56× fastMR thread multiplier and a maximum native point
estimate beta difference of `4.53e-12`. IVW and Egger beta/SE/p-value parity
was at floating-point noise; bootstrap SE/p-value distributions for the
median and mode methods are reported per method in `method_summary.csv`.

Correctness checks found exact harmonisation-key parity for local BMI → CRP,
floating-point parity for deterministic IVW/Egger estimates, and a maximum
beta difference of `6.9e-11` across the independent 50 × 50 simulation.
Bootstrap SE and p-value differences are reported separately because the two
implementations use independently sampled bootstrap streams.

## Design notes and limits

- Use `fast_mr_grid()` when scanning many exposure/outcome pairs; for one pair,
  R/data-frame overhead becomes a larger share of elapsed time.
- Pre-harmonise and clump once, then reuse the resulting matrices across scans.
- Increase `threads` for large grids, but do not expect linear scaling on tiny
  workloads.
- Named exposure/outcome grid matrices must have identical SNP column names in
  identical order; unnamed matrices are treated as already aligned.
- Kept tidy rows require finite beta values and positive standard errors;
  duplicate SNPs within each exposure/outcome group are collapsed to the first
  row, and non-empty SNP identifiers are required. Repeated SNPs across
  outcomes are restored after clumping.
- Parallel grid workers perform numerical kernels only; R probability values are
  filled after the workers join, keeping the threaded path safe for R.
- PLINK remains the preferred route when a large external LD reference panel
  is available. The matrix clumper is intended for local, reusable LD data.
- The package does not download GWAS data or reproduce TwoSampleMR's remote
  extraction stack.

## Development and verification

```sh
R CMD INSTALL --library=.local/Rlib .
Rscript -e 'testthat::test_local(".")'
_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual --as-cran .
```

The package has a 127-test suite, native harmonisation audits, simulation
tyre-kick tests, adversarial thread checks, and native TwoSampleMR benchmarks.
The `R CMD check` command above deliberately allows the optional Arrow
dependency to be absent; install Arrow first if Parquet checks are required.
See [`NEWS.md`](NEWS.md) and [`HANDOFF.md`](HANDOFF.md) for the optimization
history and detailed validation record.

## Logo

The mark is an original `MR` identity: a cobalt italic `MR` monogram under a
single silver orbital swoosh, finished with a coral acceleration tip. The
monogram identifies Mendelian randomization and the R ecosystem; the swoosh
represents a causal trajectory through the data; and the coral tip signals
fast movement through a large MR grid. It deliberately avoids copying the
official R Project mark.

## License

MIT © Fergus Hamilton.
