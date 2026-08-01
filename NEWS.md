# fastMR 0.1.0

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
