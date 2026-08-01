# IVW algorithmic optimisation benchmark

IL6 fixture: 82 SNPs; 50×50 grid; nboot=0; 100 warm repeats per thread count.
The scalar baseline is the measured pre-BLAS fastMR path on the same Mac mini and command (1.313 s).
The new path batches IVW numerator/denominator cross-products through BLAS and flattens tidy grid output once.

| threads | median fastMR s | min fastMR s | pre-BLAS scalar s | speedup |
|---:|---:|---:|---:|---:|
| 1 | 0.001 | 0.001 | 1.313 | 1492.05x |
| 5 | 0.001 | 0.001 | 1.313 | 1492.05x |
| 10 | 0.001 | 0.001 | 1.313 | 1492.05x |
