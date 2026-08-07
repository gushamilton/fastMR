# fastMR 0.1.9

- Updates the optional compressed-input integration for CompreSSoR 0.5's
  native Pcodec stores and strict prepared-input contract. The retired Python
  codec runtime is no longer required.
- Adds safe Zstandard-compressed Parquet result output through
  `fast_write_parquet()` and the `output` argument of `fast_mr()`,
  `fast_mr_grid()`, and `fast_mr_compressed()`.
- Adds opt-in `fast_clump_data_batched()` and `fast_clump_compressed()` APIs:
  exposure-specific greedy clumping is retained while PLINK2 LD queries are
  shared across all current exposure leads. Bounded work limits, reference
  manifest provenance, and explicit reconstructed-p-value labelling prevent
  silent approximation on large Pcodec runs.
- Adds tested internal masked and CSR IVW kernels for sparse exposure panels;
  adds exported `fast_mr_masked_ivw()` and `fast_mr_sparse_ivw()` wrappers with
  explicit masks/CSR contracts, duplicate-index checks, and output/native
  workspace bounds, plus a generic parity benchmark.

# fastMR 0.1.8

- Uses CompreSSoR's persistent wrapped-Pcodec reader and exact canonical keys
  for direct compressed-input MR.
- Adds a shared-instrument native IVW grid shortcut while retaining the
  established pair-specific seed stream for bootstrap-dependent methods.
- Separates explicit ten-read, optimized same-store deduplication, Tabix, and
  TSV.gz paths in the real FinnGen benchmark and applies the same IVW grid
  estimator to every fair-comparison path.

# fastMR 0.1.7

- Batches every compressed exposure and outcome read through one CompreSSoR
  process, removing repeated Python startup and reusing identical requests.
- Adds a corrected full-FinnGen 5 x 5 benchmark: median 0.140 seconds from
  Pcodec, 0.196 seconds from VCF.gz plus Tabix, and 18.175 seconds from ten
  TSV.gz scans, with all 25 REF/ALT keys and IVW results checked.

# fastMR 0.1.6

- Adds `fast_read_compressed()` and `fast_mr_compressed()` for indexed,
  canonical-key MR directly from self-contained CompreSSoR stores.
- Reads each exposure only for its instruments and each outcome once for the
  union, with optional file-level parallelism and explicit overlap counts.
- Fails clearly on unsupported backends, duplicate keys, corrupt parallel
  reads, and invalid statistics; non-strict runs report every omitted value.

# fastMR 0.1.5

* Added fast heterogeneity, MR-Egger pleiotropy, single-SNP, and leave-one-out
  utilities with TwoSampleMR-compatible tidy output.
* Added Steiger directionality testing and per-SNP Steiger filtering, including
  quantitative-trait, SD-scaled, and log-odds metadata paths. The Mac mini
  parity audit matches native TwoSampleMR exactly on the IL6 fixture.
* Added matrix-form multivariable MR with shared or exposure-specific
  instruments and an optional intercept, with native `mv_ivw` parity.
* Vectorized single-SNP Wald ratios and leave-one-out regression diagnostics;
  the Mac mini audit matches native TwoSampleMR row-for-row at floating-point
  precision, with leave-one-out about 2.9x faster and single-SNP about 2x faster
  on the 82-SNP fixture.

# fastMR 0.1.4

* Added a regression test proving that exact duplicate SNP rows leave all
  point estimates, standard errors, p-values, Q statistics, and Egger
  intercept diagnostics unchanged after deduplication.

# fastMR 0.1.3

* Duplicate SNP rows are now collapsed once per MR pair/clumping exposure,
  while repeated outcome rows are restored after clumping. Repeated p-values
  are treated as ordinary metadata.

# fastMR 0.1.2

* Added the original `MR` monogram logo and a reproducible adversarial
  25 × 25, `nboot = 1,000` validation run against native TwoSampleMR.
* Hardened the public API against duplicate methods/SNPs, malformed kept
  rows, ambiguous grouped IDs, mismatched named grid columns, and standalone
  mode-`phi` reporting. Grid Egger bootstrap now retains exact-zero exposure
  effects.
* Removed R probability-API calls from parallel workers; threaded grids now
  perform native numerical work in parallel and populate R p-values safely on
  the main thread.

# fastMR 0.1.1

* Refreshed the package identity with a distinctive R-inspired swoosh logo
  for GitHub, documentation, and package distribution.

# fastMR 0.1.0

* Added compact mixed-method grid results, preserving diagnostics while
  avoiding thousands of nested Rcpp lists. The full 50x50 IL6 five-method
  benchmark now runs in 1.564 seconds versus 339.139 seconds for native
  TwoSampleMR (216.8x faster).
* Reused weighted-median sort order across bootstrap draws and across the
  simple/weighted median pair. The current native five-shape benchmark gives
  166-254x speedups for weighted median, with a 50x50 grid at 0.068 seconds.
* Mixed-method grids now source all IVW, fixed-effects IVW, and multiplicative
  random-effects IVW rows from the BLAS batch, avoiding a second cross-product
  pass when IVW is requested alongside other methods.
* Added large-grid scaling coverage: the raw compact IVW kernel processes a
  1,000x1,000 grid (one million pairs, 82 SNPs) in about 0.102 seconds on the
  Mac mini, or about 9.8 million pairs per second before tidy-frame allocation.
* Added a compact all-IVW grid return path after the BLAS cross-products. On
  the Mac mini, the 50x50 one-method grid averages 0.00088 seconds over 100
  warm calls, 1,492x faster than the measured pre-BLAS scalar path; native
  TwoSampleMR shape speedups are 1,867-14,682x with machine-precision parity.
* Added native-compatible penalised weighted median with configurable `penk`
  and exact two-stream bootstrap parity; its five-shape benchmark reaches
  122-160x speedups over native TwoSampleMR.
* Reused exact mode-kernel FFT workspaces/plans and median-selection storage
  across bootstrap draws, reducing the 50x50 mode grid by about 7-8% without
  changing the native density semantics.
* Cached per-stage FFT twiddle factors in the exact mode transform. The
  50x50/nboot=100 mode grid is now 1.86-1.89x faster than the previous
  workspace implementation, with native parity unchanged.
* Added a batched BLAS IVW grid path for `ivw`, `ivw_fe`, and `ivw_mre`, plus
  one-pass flattening of grid results. On the Mac mini this makes a 50x50
  IVW grid about 5.6x faster than the prior scalar fastMR path and about
  55-61x faster than native TwoSampleMR across the tested grid shapes.
* Added Simple median, local dependency-free harmonisation, and local LD
  clumping through either a supplied LD matrix or a PLINK binary reference.
* Added seeded MR-Egger bootstrap with shared grid-side normal draws.
* Added unweighted regression and exact sign concordance methods.
* Initial GitHub-ready package with exact IVW, MR-Egger, weighted median,
  simple mode, weighted mode, Wald ratio, and basic multivariable IVW.
* Added a registered Rcpp C++17 shared-grid kernel with bounded parallelism and
  a serial fallback when OpenMP is unavailable.
* Added optional Arrow Parquet readers and a reproducible IL6/CRP benchmark.
* Added a reproducible simulation and harmonisation tyre-kick suite covering
  swapped, complemented, reverse-complemented, palindromic, incompatible,
  unequal-SNP, empty-overlap, and 7x11-grid cases; it matched native
  TwoSampleMR with zero harmonisation flag mismatches and a maximum
  representative grid beta delta of 1.33e-15.
* Extended local harmonisation to the full native `harmonise_data()` behavior:
  actions 1/2/3, 2-2/2-1/1-2/1-1 allele information, indel recoding, missing
  frequencies, outcome-specific action vectors, and native `mr_keep` handling.
  A 63-comparison audit against TwoSampleMR 0.7.9 matched every key field
  exactly.
* Completed a five-round Mac-mini optimization cycle: flat mixed-grid result
  storage and fused/linearized exact-mode scans were retained, while `-O3`
  and static thread chunks were rejected for portability or regression. The
  final 2,500-pair five-method workload improved from 1.544 to 1.503 seconds
  at ten threads.
