# fastMR 0.1.0

* Added Simple median, local dependency-free harmonisation, and local LD
  clumping through either a supplied LD matrix or a PLINK binary reference.
* Added seeded MR-Egger bootstrap with shared grid-side normal draws.
* Initial GitHub-ready package with exact IVW, MR-Egger, weighted median,
  simple mode, weighted mode, Wald ratio, and basic multivariable IVW.
* Added a registered Rcpp C++17 shared-grid kernel with bounded parallelism and
  a serial fallback when OpenMP is unavailable.
* Added optional Arrow Parquet readers and a reproducible IL6/CRP benchmark.
