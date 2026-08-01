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
