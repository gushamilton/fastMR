# fastMR handoff

Status: complete first package implementation in `/Users/fergushamilton/projects/fastMR`.
The prototype and TwoSampleMR projects were read-only references and were not
modified.

## Delivered

* Conventional `fastMR` R package metadata, MIT license, generated Rd docs,
  registered Rcpp routines, tests, fixture, benchmark history, README, NEWS,
  and local git-ready layout.
* `fast_mr()` for TwoSampleMR-style tidy data, grouping by optional
  `id.exposure`/`id.outcome`, with `id.exposure`, `id.outcome`, `method`,
  `nsnp`, `b`, `se`, and `pval` plus method diagnostics.
* `fast_mr_grid()` with exposure-major/outcome-minor ordering. The exact C++17
  kernel copies R column-major matrices once into contiguous row-major storage,
  shares seeded bootstrap layouts and exposure/outcome grid-side layouts, and
  bounds worker count. It uses OpenMP when the toolchain exposes it and a
  portable bounded `std::thread`/serial fallback otherwise; Ferguss-Mini's
  Apple Clang build used the fallback.
* Exact default methods: IVW, fixed-effects IVW, multiplicative random-effects
  IVW, MR-Egger, weighted median, simple mode, weighted mode, and Wald ratio.
  The mode implementation matches native R `stats::density` semantics with the
  same weighted binning, extended FFT range, interpolation, and 512-point grid;
  no approximate mode estimate is routed through the default API. Seeded
  bootstraps use R RNG draw order and preserve the caller's RNG state.
* Basic multivariable IVW and optional `arrow::read_parquet()` wrappers. Arrow
  remains in Suggests and has a targeted error when unavailable.

The 82-row IL6/CRP fixture and provenance are in `inst/extdata/`. The C++
kernel is a direct R port of the validated shared-grid implementation in the
prototype; no upstream TwoSampleMR source or remote data-extraction stack is
bundled.

## Final verification

Commands used the required explicit paths:

```sh
/opt/homebrew/bin/R CMD INSTALL --library=/Users/fergushamilton/projects/fastMR/.local/Rlib .
/opt/homebrew/bin/Rscript -e 'lib <- "/Users/fergushamilton/projects/fastMR/.local/Rlib"; .libPaths(c(lib, .libPaths())); testthat::test_local(".")'
/opt/homebrew/bin/R CMD build .
/opt/homebrew/bin/R CMD check --no-manual --as-cran fastMR_0.1.0.tar.gz
```

Results: local install succeeded; testthat passed 36 expectations with one expected
skip because Arrow is installed; the tarball check completed with **0 errors,
0 warnings, and 3 notes**. The notes are the planned GitHub URL returning 404
until a remote is created, the full MIT text being reported as a non-DCF
license stub by this R build, and pandoc being absent for README/NEWS checks.
No GitHub remote was created or pushed.

## Optimization history

The complete per-loop hypothesis, exact command, wall time, throughput,
memory field, and correctness gate are in:

* `outputs/optimization_history.csv`
* `outputs/optimization_history.json`
* `outputs/optimization_history.md`

All 13 recorded rows passed their correctness gates. The package-side measured
loops on the 82-row fixture / 2,500-pair grid (`nboot=100`, seed `20260801`)
were:

| loop | wall seconds | pairs/s |
|---|---:|---:|
| single pair, five methods | 0.008 | 125 |
| IVW + Egger single pair | 0.001 | 1,000 |
| grid, nboot 0, threads 1 | 2.367 | 1,056 |
| grid, nboot 0, threads 4 | 2.349 | 1,064 |
| grid, nboot 0, threads 10 | 2.287 | 1,093 |
| grid, nboot 100, threads 1 | 16.725 | 149 |
| grid, nboot 100, threads 4 | 6.543 | 382 |
| grid, nboot 100, threads 10 | 5.089 | 491 |
| Arrow Parquet read | 0.001 | 82,000 rows/s |
| seeded thread correctness rerun | 5.103 | 490 |

Primary native-R comparison: outputs/native_tsmr_grid_benchmark.csv records 4.945 seconds for fastMR versus 338.919 seconds for the standard TwoSampleMR::mr() workflow on the same 82-row IL6, 2,500-pair, five-method, nboot=100 workload: a 68.538x speedup. The maximum absolute beta difference across all five methods and all 2,500 pairs was 3.47e-18. IVW, Egger, weighted median, simple mode, and weighted mode each matched native TwoSampleMR to machine precision in single-method seeded randomized parity tests (12 panels; maximum combined beta/SE/p-value difference 1.01e-12). The full-grid bootstrap SEs are intentionally independent Monte Carlo draws at nboot=100, so their median absolute differences were 9.28e-05 for weighted median, 8.89e-04 for simple mode, and 2.52e-04 for weighted mode; these are recorded in outputs/native_tsmr_grid_parity.csv rather than misreported as deterministic estimator errors.

The 5-thread adversarial gate passed 20 randomized grids plus near-zero, negative, duplicate-ratio, and single-SNP edge panels with a maximum serial-versus-5-thread difference of 0. The native-R mode point audit passed 40 randomized panels with maximum absolute beta difference 1.05e-13. The package test suite passed 36 expectations with one expected Arrow skip, and the built tarball check completed with 0 errors, 0 warnings, and 3 notes.

The thread correctness gate reported a maximum seeded difference of `0` for
`threads=1`, `4`, and `10`. Bootstrap values are deterministic within fastMR;
their seeded standard errors need not be bitwise identical to NumPy or
TwoSampleMR because those implementations use different random-number
streams.

## Per-method shape benchmarks

`outputs/native_tsmr_method_benchmark.csv` and `.md` compare each default MR
method separately with native `TwoSampleMR::mr()` across five grid shapes:
50x50, 1x250, 250x1, 10x100, and 100x10. Each run used the 82-row IL6/CRP
fixture, `nboot=100`, seed `20260802`, and five fastMR worker threads. The
measured speedup ranges were:

| method | speedup range |
|---|---:|
| IVW | 10.08-10.82x |
| MR-Egger | 8.67-9.85x |
| weighted median | 12.06-13.21x |
| simple mode | 45.12-48.99x |
| weighted mode | 44.82-48.28x |

All deterministic beta differences were at machine precision (maximum
`8e-18` for IVW/Egger/median and approximately `2.7e-16` for mode point
estimates). Bootstrap SE and p-value differences are Monte Carlo differences:
the two implementations generate independent bootstrap streams at `nboot=100`.
The grid kernel now also skips unused stochastic buffers for deterministic
method-only calls; on the balanced 50x50 grid this reduced IVW and Egger to
1.298 and 1.436 seconds respectively at five threads.

## Remaining limitations

The current package is intentionally local-summary-statistics only: it does
not implement OpenGWAS extraction, harmonisation, clumping, or the broader
TwoSampleMR method catalogue. The package URL is a planned destination only;
creating or pushing a remote was deliberately left to the owner.
