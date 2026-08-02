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

All timings below were measured on the same Mac mini against native
`TwoSampleMR` 0.7.9. Times are seconds; speedups are native/fastMR. The
results are observed measurements, not guarantees for every machine. Raw
CSV files and the scripts that generated them are linked in each table.

### End-to-end audits

| Workload | Shape / SNPs | Methods | `nboot` / threads | fastMR | Native TSMR | Speedup | Parity / checks |
|---|---|---:|---:|---:|---:|---:|---|
| Local BMI → CRP | 1 × 1 / 104 SNPs | 5 | 100 / 1 | 0.005 s | 0.127 s | **25.4×** | Harmonisation keys exact |
| Independent simulation | 1 × 1 / 400 SNPs | 5 | 100 / 1 | 0.008 s | 0.165 s | **20.6×** | Max beta delta `2.78e-16` |
| Independent heavy-tailed simulation | 50 × 50 / 400 SNPs, 2,500 pairs | 5 | 100 / 10 | 2.030 s | 391.764 s | **193.0×** | Max beta delta `6.9e-11` |
| Adversarial heavy-tailed simulation | 25 × 25 / 82 SNPs, 625 pairs / 3,125 rows | 5 | 1,000 / 10 | 3.630 s | 824.585 s | **227.2×** | 1-thread/10-thread exact; max beta delta `4.53e-12` |

The 50 × 50 run processed about 1,232 pairs/s with fastMR versus 6.4
pairs/s natively. The adversarial run took 20.177 s with one fastMR thread,
so the measured thread multiplier was 5.56×. The 25 × 25 run also exercised
one- and two-SNP inputs, invalid/non-finite inputs, deliberately flipped
alleles, and all five main methods. See the
[final package audit](outputs/final_package_audit/) and the
[25 × 25 audit](outputs/adversarial_25x25_nboot1000/), reproduced by
[`benchmarks/final_package_audit.R`](benchmarks/final_package_audit.R) and
[`benchmarks/adversarial_25x25_nboot1000.R`](benchmarks/adversarial_25x25_nboot1000.R).

### Per-method native comparison across five grid shapes

This is the complete five-method grid benchmark: balanced 50 × 50,
1 × 250, 250 × 1, 10 × 100, and 100 × 10; 400 SNPs, `nboot = 100`, and
five fastMR threads. Ranges are the minimum–maximum across those shapes.

| Method | fastMR time | Native TSMR time | Speedup | Max abs Δ beta | Max median abs Δ SE | Max abs Δ p |
|---|---:|---:|---:|---:|---:|---:|
| IVW | 0.001–0.001 s | 1.498–14.083 s | **1,498–14,083×** | `3.25e-18` | `4.34e-19` | `7.77e-16` |
| MR-Egger | 0.001–0.002 s | 1.581–14.494 s | **1,581–7,247×** | `8.67e-18` | `1.73e-18` | `9.44e-16` |
| Weighted median | 0.010–0.068 s | 1.826–17.296 s | **166.0–254.4×** | `4.34e-19` | `1.11e-4` | `0.148` |
| Simple mode | 0.156–1.286 s | 18.038–173.085 s | **114.9–134.6×** | `2.66e-16` | `2.15e-3` | `0.044` |
| Weighted mode | 0.158–1.290 s | 17.988–176.129 s | **85.3–136.5×** | `2.52e-16` | `2.31e-4` | `0.270` |

The full 25-row result is in
[`outputs/native_tsmr_method_benchmark.csv`](outputs/native_tsmr_method_benchmark.csv).
The additional registered methods were benchmarked across the same five
shapes as follows:

| Additional method | fastMR time | Native TSMR time | Speedup | Max abs Δ beta | Max median abs Δ SE | Max abs Δ p | Raw results |
|---|---:|---:|---:|---:|---:|---:|---|
| Simple median | 0.138–1.345 s | 1.753–17.365 s | **12.17–13.04×** | `3.30e-17` | `5.83e-4` | `0.142` | [CSV](outputs/native_tsmr_simple_median_benchmark.csv) |
| Penalised weighted median | 0.016–0.143 s | 2.068–20.544 s | **122.5–159.5×** | `6.51e-19` | `1.75e-4` | `0.154` | [CSV](outputs/native_tsmr_penalised_median_benchmark.csv) |
| MR-Egger bootstrap | 0.139–1.380 s | 1.843–18.293 s | **12.84–13.35×** | `6.37e-4` | `1.14e-4` | `0.110` | [CSV](outputs/native_tsmr_egger_bootstrap_benchmark.csv) |
| Unweighted regression | 0.132–1.307 s | 1.390–13.910 s | **10.45–10.66×** | `4.34e-19` | — | `0` | [CSV](outputs/native_tsmr_uwr_sign_benchmark.csv) |
| Sign concordance | 0.128–1.240 s | 1.457–13.708 s | **10.83–11.29×** | `0` | — | `0` | [CSV](outputs/native_tsmr_uwr_sign_benchmark.csv) |

Bootstrap SE and p-value differences are reported as distributions because
native TwoSampleMR and fastMR use independently sampled bootstrap streams;
the deterministic IVW/Egger estimates agree to floating-point precision.

### Scaling and implementation benchmarks

| Component | Workload | Result | Interpretation | Raw results |
|---|---|---:|---|---|
| Compact IVW kernel | 50 × 50 through 1,000 × 1,000; 82 SNPs | 0.00026–0.109 s; 9.17–10.59 million pairs/s | Linear large-grid throughput before tidy-frame allocation | [CSV](outputs/ivw_large_grid_scaling.csv) |
| Batched BLAS IVW path | 2,500 pairs; threads 1, 5, 10 | 0.00088 s; **1,492×** vs pre-BLAS scalar path | Same result and speed across tested thread counts | [CSV](outputs/ivw_algorithmic_benchmark.csv) |
| Five-method grid, no bootstrap | 2,500 pairs; 1 → 10 threads | 0.082 → 0.017 s; **4.82×** | Threading helps even without bootstrap work | [CSV](outputs/iteration_scaling.csv) |
| Five-method grid, `nboot = 100` | 2,500 pairs; 1 → 10 threads | 7.885 → 1.503 s; **5.25×** | Bootstrap-heavy work scales substantially with threads | [CSV](outputs/iteration_scaling.csv) |
| Mode workspace reuse | 2,500 pairs; simple / weighted mode | 1.886 → 1.012 s / 1.910 → 1.011 s | **1.86× / 1.89×** over the previous workspace path | [CSV](outputs/mode_optimization_benchmark.csv) |
| Local LD-matrix clumping | 1,500 pair rows | 0.037 s | Dependency-light local clumper | [CSV](outputs/iteration_scaling.csv) |
| Local harmonisation | 20,000 rows | 0.015 s vs 0.107 s native; **7.13×** | Same kept rows and beta values | [CSV](outputs/iteration_harmonise_benchmark.csv) |

### Correctness and adversarial checks

| Check | Coverage | Result |
|---|---|---|
| 25 × 25, `nboot = 1,000` proof run | 3,125 expected rows; 1 vs 10 threads; native comparison | All rows present; thread delta `0`; max native beta delta `4.53e-12` |
| Native harmonisation options | 63 combinations across actions 1, 2, and 3, including strand, palindromic, indel, missing-frequency, and incomplete-allele cases | **63/63 exact key matches**, max numeric delta `0` |
| Simulated harmonisation actions | 60 overlapping SNPs under actions 1, 2, and 3 | Zero keep/remove/ambiguous mismatches; max beta delta `0` |
| Unequal SNP sets | 47 exposure SNPs, 34 outcome SNPs | Expected 27-row overlap produced |
| Duplicate SNP regression | Five main methods; duplicate rows and repeated p-values | All `b`, `se`, `pval`, Q, and Egger diagnostic deltas `0` |
| Randomised threaded IVW | 100 adversarial panels | All scalar and threaded gates passed with delta `0` |
| Mixed threaded grids | 150 adversarial panels | All exactness and finite-result gates passed |
| Penalised median edge panels | 300 adversarial panels | All exactness and finite-result gates passed |

The proof-run details are in
[`outputs/adversarial_25x25_nboot1000/`](outputs/adversarial_25x25_nboot1000/);
the harmonisation and simulation files are in [`outputs/`](outputs/).

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

The package has a 134-test suite, native harmonisation audits, simulation
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
