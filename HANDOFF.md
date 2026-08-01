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
  portable bounded `std::thread`/serial fallback otherwise. The Mac mini uses
  the portable path because its R build does not expose standard OpenMP flags;
  other toolchains retain OpenMP automatically when R exposes them.
* IVW-only grids (`ivw`, `ivw_fe`, `ivw_mre`) use a batched BLAS kernel for
  numerator/denominator cross-products across all exposure/outcome pairs.
  Grid result tidying is also flattened into one data-frame construction rather
  than one data frame per pair.
* Exact default methods: IVW, fixed-effects IVW, multiplicative random-effects
  IVW, MR-Egger, weighted median, simple mode, weighted mode, and Wald ratio.
  The mode implementation matches native R `stats::density` semantics with the
  same weighted binning, extended FFT range, interpolation, and 512-point grid;
  no approximate mode estimate is routed through the default API. Seeded
  bootstraps use R RNG draw order and preserve the caller's RNG state.
* Basic multivariable IVW and optional `arrow::read_parquet()` wrappers. Arrow
  remains in Suggests and has a targeted error when unavailable.
* Simple median, `fast_harmonise_data()` for local allele alignment, and
  `fast_clump_data()` for either an in-memory LD matrix or a local PLINK
  reference panel. The preprocessing layer has no mandatory new dependency.

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

Results: local install succeeded; the testthat suite passed with one expected
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
| grid, nboot 0, threads 1 | 1.258 | 1,987 |
| grid, nboot 0, threads 4 | 1.152 | 2,170 |
| grid, nboot 0, threads 10 | 1.148 | 2,178 |
| grid, nboot 100, threads 1 | 16.587 | 151 |
| grid, nboot 100, threads 4 | 5.348 | 467 |
| grid, nboot 100, threads 10 | 3.864 | 647 |
| Arrow Parquet read | 0.001 | 82,000 rows/s |
| seeded thread correctness rerun | 5.103 | 490 |

Primary native-R comparison: outputs/native_tsmr_grid_benchmark.csv records 4.945 seconds for fastMR versus 338.919 seconds for the standard TwoSampleMR::mr() workflow on the same 82-row IL6, 2,500-pair, five-method, nboot=100 workload: a 68.538x speedup. The maximum absolute beta difference across all five methods and all 2,500 pairs was 3.47e-18. IVW, Egger, weighted median, simple mode, and weighted mode each matched native TwoSampleMR to machine precision in single-method seeded randomized parity tests (12 panels; maximum combined beta/SE/p-value difference 1.01e-12). The full-grid bootstrap SEs are intentionally independent Monte Carlo draws at nboot=100, so their median absolute differences were 9.28e-05 for weighted median, 8.89e-04 for simple mode, and 2.52e-04 for weighted mode; these are recorded in outputs/native_tsmr_grid_parity.csv rather than misreported as deterministic estimator errors.

The 5-thread adversarial gate passed 20 randomized grids plus near-zero, negative, duplicate-ratio, and single-SNP edge panels with a maximum serial-versus-5-thread difference of 0. The native-R mode point audit passed 40 randomized panels with maximum absolute beta difference 1.05e-13. The package test suite passed with one expected Arrow skip, and the built tarball check completed with 0 errors, 0 warnings, and 3 notes.

The thread correctness gate reported a maximum seeded difference of `0` for
`threads=1`, `4`, and `10`. Bootstrap values are deterministic within fastMR;
their seeded standard errors need not be bitwise identical to NumPy or
TwoSampleMR because those implementations use different random-number
streams.

The IVW-specific algorithmic benchmark is recorded in
`outputs/ivw_algorithmic_benchmark.csv` and `.md`. On the same 2,500-pair
fixture, the new BLAS plus one-pass tidy path took a median 0.235-0.237 seconds
over five warm repeats at 1, 5, and 10 requested threads, versus 1.313 seconds
for the measured pre-BLAS scalar fastMR path: 5.54-5.59x faster. The
implementation gain is therefore present before thread scaling.

The 100-panel adversarial IVW gate is recorded in
`outputs/adversarial_ivw_threads.csv`. It included zero, negative, near-zero,
and all-zero exposure rows. The maximum direct-versus-scalar beta/SE/p-value
delta was `2.08e-14`, with zero scalar-gate failures; serial-versus-five-thread
results were bitwise identical in all 100 panels.

## Implementation-versus-thread scaling

The direct `system.time()` scaling run is recorded in
`outputs/iteration_scaling.csv` and `.md`. On the 2,500-pair IL6 grid with all
five main methods, `nboot=0` took 1.470 seconds at one thread and 1.197 seconds
at ten threads. At `nboot=100`, one thread took 17.106 seconds and ten threads
took 3.887 seconds (4.40x thread scaling). Against the prior native
TwoSampleMR reference of 338.919 seconds, this is approximately 19.8x faster
serially and 87.2x faster at ten threads; the remaining gain is implementation
and shared-grid work, not just threading.

The local harmonisation path processed 20,000 synthetic strand/complement/
palindrome cases in 0.015 seconds versus 0.107 seconds for native
`TwoSampleMR::harmonise_data()` (7.13x), with zero beta, allele-frequency, or
`mr_keep` mismatches. Local LD-matrix clumping retained the expected index SNPs
from a 1,500-variant test panel in about 0.038 seconds. The PLINK path is a
direct local `system2()` delegation and was not timed on this Mini because no
PLINK executable is installed there.

The newly added Simple median is benchmarked separately in
`outputs/native_tsmr_simple_median_benchmark.csv` and `.md`: across the same
five grid shapes it ran 12.17-13.04x faster than native TwoSampleMR at five
threads, with a maximum beta difference of `3.30e-17`. Its bootstrap SE and
p-value differences are Monte Carlo differences from independent streams.

MR-Egger bootstrap is now also implemented as `egger_bootstrap`. Its seeded
single-pair parity audit over 12 randomized panels reached maximum beta, SE,
and p-value deltas of `2.22e-16`, `2.78e-17`, and `0`; intercept deltas were
below `1.1e-17`. The grid reuses raw exposure/outcome bootstrap draws across
pairs and remains bitwise deterministic across one and five threads. The
five-shape native benchmark is recorded in
`outputs/native_tsmr_egger_bootstrap_benchmark.csv` and `.md`; it measured
12.84-13.35x speedups at five threads. Its maximum grid beta delta was
`6.37e-4` because native TwoSampleMR and fastMR use independent per-pair
versus shared-grid bootstrap streams; the seeded single-pair estimator parity
above is the exact implementation gate.

Unweighted regression (`uwr`) and the exact binomial sign concordance test
(`sign`) are also implemented. Over 20 randomized single-pair panels their
maximum native beta/SE/p-value deltas were `8.33e-17`, `8.88e-16`, and
`2.22e-16` for UWR, and zero, zero, and `2.22e-16` for sign. Their five-shape
grid speedups and machine-precision point parity are recorded in
`outputs/native_tsmr_uwr_sign_benchmark.csv` and `.md`.

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
The grid kernel also skips unused stochastic buffers for deterministic
method-only calls. The dedicated IVW benchmark below measures the new BLAS
path separately; the older mixed-method measurements remain useful for
estimating the cost of requesting all five methods together.

`outputs/native_tsmr_ivw_benchmark.csv` and `.md` compare the new IVW-only
path with native TwoSampleMR across the same five grid shapes. At five
requested threads, speedups were 54.93-60.67x, with maximum beta deltas below
`3.3e-18`, maximum SE deltas below `3.1e-18`, and maximum p-value deltas below
`8e-16`.

## Remaining limitations

The current package is intentionally local-summary-statistics only: it does
not implement OpenGWAS extraction, proxy lookup, or the broader TwoSampleMR
method catalogue. The package URL is a planned destination only; creating or
pushing a remote was deliberately left to the owner. Harmonisation currently
targets common bi-allelic SNPs; indel recoding and proxy handling remain out of
scope for the local fast path.
