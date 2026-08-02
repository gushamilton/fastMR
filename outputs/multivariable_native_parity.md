# fastMR multivariable MR versus native TwoSampleMR

TwoSampleMR 0.7.9; deterministic 800-SNP × 4-exposure simulation; timings are medians of 20 fast and 5 native calls.

| component | fastMR s | TwoSampleMR s | speedup | rows | key/nsnp match | max beta delta | max SE delta | max p delta |
|---|---:|---:|---:|---:|---|---:|---:|---:|
| mv_ivw | 0.002000 | 0.013000 | 6.50x | 4/4 | TRUE/TRUE | 1.665e-16 | 6.505e-19 | 0.000e+00 |
| mv_multiple_shared | 0.002000 | 0.002000 | 1.00x | 4/4 | TRUE/TRUE | 4.441e-16 | 7.589e-19 | 0.000e+00 |
| mv_multiple_specific | 0.001000 | 0.002000 | 2.00x | 4/4 | TRUE/TRUE | 1.665e-16 | 6.505e-19 | 0.000e+00 |
| mv_multiple_intercept | 0.001000 | 0.002000 | 2.00x | 4/4 | TRUE/TRUE | 2.498e-16 | 1.193e-18 | 0.000e+00 |
