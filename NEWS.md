# fastMR (development)

- `fast_mr_compressed()` extraction: stores are opened and validated
  `io_threads` at a time instead of one by one, and the batched reader returns
  each row's numeric identity (global position and substitution code, as
  declared by the store manifest's `compressor_variant_identity_v1` encoding)
  instead of decoded chromosome/allele strings. Each row then takes the
  variant key of the requested key with the same identity code, so no key
  string is rebuilt per row. Results, counts, errors and warnings are
  unchanged; stores whose manifest does not declare that encoding use the
  previous string path.
- `fast_mr_compressed()`'s pairwise path (every non-IVW method set, or
  `estimator = "pairwise"`) and its `steiger = TRUE` pass now assemble the
  harmonised table with one vectorised pair index (one `match()` per store,
  no per-pair data frames or `do.call(rbind)`), then make the same single
  batched `fast_mr()` call. Results, counts, Steiger rows, errors and row
  order are identical to the per-pair loop; the estimator step now scales to
  millions of pairs. `strict = FALSE` omission warnings (pairwise and sparse
  IVW paths) list at most `getOption("fastMR.warning_pairs", 20)` entries
  and then the total count, instead of one string naming every pair; the
  per-pair detail stays in `attr(result, "compressed_input")$counts`, and
  `options(fastMR.warning_pairs = Inf)` restores the full listing.
- Faster mode bootstraps and threaded bootstrap batches, with every result
  (and the final `.Random.seed`) identical to before:
  - Mode densities for pairs with up to 1000 ratios first try a "hull" path that
    convolves and scans only the grid cells spanning the occupied bins, with a
    Gaussian-recurrence kernel (one exact `exp()` per 16 distances). Its guard
    carries a written error bound (see `src/fastmr.cpp`); any draw it cannot
    certify falls back to the existing direct/FFT paths.
  - Threaded bootstrap batches are double-buffered: the main thread draws the
    next batch's normals while the workers compute the current batch, and
    fills the previous batch's p-values meanwhile. R's RNG is still consumed
    only on the main thread, in the serial order.

- Single-instrument pairs (`nsnp = 1`): the IVW estimators (`"ivw"`,
  `"ivw_fe"`, `"ivw_mre"`) return the Wald ratio, `b = by / bx` and
  `se = se_y / |bx|`, exactly as TwoSampleMR's `mr()` reports it through
  `mr_wald_ratio()`, with `Q` and `sigma` `NA`. This holds for `fast_mr()`,
  `fast_mr_grid()` (tidy and compact), `fast_mr_sparse_ivw()`,
  `fast_mr_masked_ivw()` and every `fast_mr_compressed()` path. They used to
  return no estimate (NA/NaN), which left every single-instrument exposure (a
  quarter of the UKB-PPP cis exposures) without an IVW result. Other methods stay
  `NA` at `nsnp = 1`, as in TwoSampleMR. An exact fit with two or more
  instruments is unaffected (it keeps the fixed-effect se).
- Missing native results are now R's `NA` rather than `NaN`.

Correctness fixes from an adversarial review:

- `fast_mr()`: a `seed` for which `seed + (number of pairs) - 1` exceeds
  `.Machine$integer.max` is rejected up front (pair i is seeded with
  `seed + i - 1`); it used to fail part-way through the run. Every entry point
  now rejects seeds outside `[-.Machine$integer.max, .Machine$integer.max]`.
- `fast_clump_compressed()` candidate reads: the batched p-value flag read must
  return exactly each store's flagged row ids (not just as many rows), and
  every key's position must equal its `base_pair_location`; the full-store
  batch path checks unique row ids and keys too. Any mismatch falls back to the
  per-store reader with a warning. The per-store flagged-row reader now stops
  when `read_sumstats()` returns fewer rows than requested.
- Clumping (graph, per-exposure, auto, batched and lead-row partitions):
  candidates with equal p are now ordered by larger |z| before SNP ID. p
  underflows to 0 above |z| ~ 38 (CompreSSoR reconstructed p and many cis-pQTL
  files), and the lead used to be the lexicographically first SNP. |z| comes
  from `beta.exposure / se.exposure` when present, and from the stores' `z`
  for `fast_clump_compressed()`. Results are unchanged when p has no ties.
- Clumping: a repeated (exposure, SNP) pair is ordered by its smallest p
  (it was the first row's p); every row of a retained pair is still returned.
  This also applies to `fast_clump_data()`.
- `fast_clump_data_graph()`, `fast_clump_data_per_exposure()` and
  `fast_clump_data_auto()` count eligible candidates absent from the LD
  reference (the graph partition with one PLINK2 `--write-snplist` query),
  warn with the count, and stop when every candidate is absent (usually a
  SNP-ID scheme mismatch). New argument `absent = c("keep", "drop")`: `"keep"`
  (default) preserves the old results, `"drop"` removes them as TwoSampleMR
  does. `diagnostics$absent_from_reference` reports the count.
- `fast_clump_data_per_exposure()` warns when candidate positions or
  chromosome labels disagree with the reference (e.g. a different genome
  build) instead of delegating to the graph partition silently.
- `fast_harmonise_data()`: an exposure with identical alleles (e.g. A/A) is
  marked `remove` when the outcome has two alleles, and its outcome effect is
  never flipped (the outcome beta used to be negated). This matches
  TwoSampleMR, which also keeps (unflipped, subject to its ambiguity rules) a
  row whose outcome has only an effect allele.
- `fast_mr_compressed()`: an exposure with an empty instrument set is dropped
  with a warning when `strict = FALSE` (an error with `strict = TRUE`, and
  when every set is empty); it used to abort in both modes. The sparse IVW
  memory check now includes the outcome-by-union-instrument matrices and
  per-pair count matrices it builds, against
  `getOption("fastMR.sparse_ivw_max_memory_mb", 8192)` MiB (above it the
  pairwise path runs). The documentation now states that the shared-grid
  fast path is used even with `estimator = "pairwise"`.
- `fast_mr_grid()`: the OpenMP pair loop uses a 64-bit index, so grids with
  more than 2^31 - 1 pairs no longer overflow it.

- `fast_clump_compressed()` (`candidate_source = "pvalue_flag"`) checks the
  batched `read_candidates_batch(strategy = "pvalue_flag")` result against each
  store's flagged-row count (from the store manifest, or the flag stream when
  the manifest lacks it) and falls back to the per-store reader, with a
  warning, on any mismatch or error. CompreSSoR 0.7.0 could silently drop
  flagged rows in batches that mixed variant sets, which lost every instrument
  for some exposures. The full-store batch path, which has no count to check
  against, is used only with a CompreSSoR that reports the
  `"candidates_batch_rows_checked"` capability (>= 0.7.1).
- IVW on an exact fit (residual standard error 0, e.g. a self-pair with
  outcome = exposure) returns the fixed-effect standard error, as TwoSampleMR
  and fastMR <= 0.1.9's sparse kernel did, instead of se 0 and p NA (`"ivw"`)
  or se NA (`"ivw_fe"`). This applies to `fast_mr()`, the shared-grid and
  sparse IVW kernels (and so `fast_mr_compressed()`) and leave-one-out IVW.
  The sparse kernel had adopted the 0-se formula in d53afb8 (0.1.10).

# fastMR 0.2.0


- `fast_clump_compressed()` now defaults to `partition = "auto"`
  (`fast_clump_data_auto()`), which picks between the all-pairs graph and the
  new `fast_clump_data_per_exposure()` (one PLINK2 `--clump` per exposure on a
  candidate-only `--extract --make-pgen` subset, P = exact greedy rank,
  window/r2 arguments translated to the graph's inclusive comparisons, leads
  certified against the graph's `--r2-phased` LD) from the estimated
  candidate pair count. Instruments are identical to the graph
  and lead-row strategies; the choice is recorded in `diagnostics$auto`.
- `candidate_source = "pvalue_flag"` reads candidates with one
  `CompreSSoR::read_candidates_batch(strategy = "pvalue_flag")` pass, decoding
  only the flagged rows' exact ranks instead of each store's full rank vector.
  Membership is still the writer-time flag, filtered to `pvalue_threshold`.

- Graph clumping now asks PLINK2 for an uncompressed `.vcor` (about 30 bytes
  per edge, written to the work directory and deleted after parsing) and
  parses it in C++ straight into vertex ids; zstd is no longer needed for
  graph clumping. Unknown IDs, short lines and unrecognised headers error.
- Adds `fast_mr_steiger_r2()` as a composable vectorized primitive for
  continuous beta/SE/sample-size, standardized beta/EAF, and binary log-odds
  models. Scalar inputs recycle deterministically, invalid rows carry explicit
  validity/reason fields, and binary prevalence is never assumed.
- Updates `fast_mr_steiger_filtering()` to use the explicit R-squared models,
  compute continuous-trait R-squared without requiring p-values, and report
  why R-squared estimates or Steiger p-values are unavailable while preserving
  the established result columns.
- Adds optional pair-specific SNP filtering to `fast_mr_sparse_ivw()` through
  an outcome-by-concatenated-CSR-entry `pair_snp_keep` matrix. `NULL` preserves
  the existing path, `NA` is rejected, and returned `nsnp` is post-filter.
- Synchronizes package and citation metadata and excludes benchmark-result
  placeholders from source-package builds.

- Performance batch (streaming pipeline):
  - `fast_mr(threads = k)` now runs bootstrap methods (medians, penalised
    weighted median, modes, Egger bootstrap) on `k` threads. Standard normals
    are pre-drawn on the main thread in exactly the per-pair order (the
    caller's stream when `seed = NULL`, `set.seed(seed + i - 1)` per pair
    otherwise) into a reused native buffer of at most
    `options(fastMR.bootstrap_batch_draws)` draws (default 2^23, 64 MB), and
    p-values are computed serially. With `threads = 1`, or for a group whose
    draws alone exceed that budget, draws stream straight into the bootstrap
    layouts exactly as before, with no buffer (so such large groups run
    serially; raise the option to parallelise them at the cost of memory).
    Output and the post-call RNG state are byte-identical to the previous
    serial implementation for every thread count.
  - `fast_mr()` groups rows once and batches non-RNG methods in a single native
    call; output is byte-identical to the previous implementation (attribute
    order, compact `row.names`).
  - `fast_mr_steiger_filtering()` and diagnostic grouping are vectorised;
    supplied `effective_n`/`rsq_valid`/`rsq_reason` columns are still
    overwritten as before.
  - `fast_mr_heterogeneity()`, `fast_mr_pleiotropy_test()`,
    `fast_mr_singlesnp()`, `fast_mr_leaveoneout()` and
    `fast_mr_directionality_test()` no longer loop over pairs in R: each makes
    one batched `fast_mr()` call (Egger leave-one-out adds one native call
    over drop-one fits, with no expanded copy of the data) and builds its
    output column-wise. Output is identical to the per-pair code, including
    row order, `row.names` and attribute order. `threads` now takes effect:
    RNG-free `fast_mr()` groups run on `threads` workers, with identical
    results for every thread count. Egger leave-one-out of a pair with a
    single SNP used to error when a sample-size column was present; it now
    reports `NA` for that row's `samplesize`.
  - `fast_mr_sparse_ivw()` checks CSR rows for duplicate SNP indices in one
    vectorised pass.
  - `fast_mr_sparse_ivw()` gains `steiger_exposure_rsq`, `steiger_outcome_rsq`
    and `pair_snp_drop` so Steiger and drop masks are applied inside the
    kernel; Q uses a stable two-pass computation and `se` matches `fast_mr()`.
  - `fast_mr_compressed()` gains `estimator = c("auto", "pairwise")`. With the
    `"auto"` default, IVW-only runs with non-shared instrument sets use the
    sparse CSR kernel; counts, errors and row order match `"pairwise"` and
    estimates agree to within 1e-14 relative (usually bit-identical). The path
    taken is recorded as `estimator_path` in the `compressed_input` attribute.
  - `fast_mr_grid()` gains `return = c("tidy", "compact", "none")` and
    `chunk_pairs` for compact and streamed Parquet output
    (`fastmr_grid_chunk()` accessor); `fast_write_parquet()` gains
    `chunk_pairs`. Streamed IVW-only grids are computed in exposure blocks and
    agree with the tidy result to about 1e-15 relative.
  - New `partition = "graph"` clumping (`fast_clump_data_graph()`,
    `max_graph_pairs`): one exact PLINK2 all-pairs call per chromosome, with ID-only
    LD columns, streamed block parsing of the LD table, and empty-graph handling
    when no candidates are in the reference. Also adds
    `fast_clump_data_lead_rows()` and `fast_clump_data_batched_chromosomal()`.
  - Source builds now exclude `slurm/`.

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
