# IVW algorithmic optimisation benchmark

IL6 fixture: 82 SNPs; 50×50 grid; nboot=0; five warm repeats per thread count.
The scalar baseline is the measured pre-BLAS fastMR path on the same Mac mini and command (1.313 s).
The new path batches IVW numerator/denominator cross-products through BLAS and flattens tidy grid output once.

| threads | median fastMR s | min fastMR s | pre-BLAS scalar s | speedup |
|---:|---:|---:|---:|---:|
| 1 | 0.237 | 0.231 | 1.313 | 5.54x |
| 5 | 0.235 | 0.230 | 1.313 | 5.59x |
| 10 | 0.235 | 0.228 | 1.313 | 5.59x |
