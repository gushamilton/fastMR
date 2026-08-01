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
  The validated native mode scan returns the baseline 512-grid coordinate; no
  approximate mode estimate is routed through the default API.
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

Results: local install succeeded; testthat passed 35 tests with one expected
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

Measured speedups: threads 10 versus 1 were 3.29x for the exact bootstrap
grid; the package grid was 18.67x faster than the existing Python reference
(95.032 s) at its comparable worker setting. The validated prototype native
artifact remains faster at 1.663 s / 1,503 pairs/s, so this first Rcpp port is
3.06x slower than that existing native artifact while keeping the same exact
kernel design and R-facing workflow. Against installed TwoSampleMR 0.7.9,
the warm five-method single-pair loop measured 0.008 s versus 0.139 s
(approximately 17.4x; startup is not included).

The thread correctness gate reported a maximum seeded difference of `0` for
`threads=1`, `4`, and `10`. Bootstrap values are deterministic within fastMR;
their seeded standard errors need not be bitwise identical to NumPy or
TwoSampleMR because those implementations use different random-number
streams.

## Remaining limitations

The current package is intentionally local-summary-statistics only: it does
not implement OpenGWAS extraction, harmonisation, clumping, or the broader
TwoSampleMR method catalogue. The package URL is a planned destination only;
creating or pushing a remote was deliberately left to the owner.
