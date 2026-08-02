# fastMR diagnostics versus native TwoSampleMR

TwoSampleMR 0.7.9; IL6 fixture: 82 SNPs; fastMR threads=5; timings are medians of 20 fast and 5 native calls.

| component | fastMR s | TwoSampleMR s | speedup | rows | row/key match | max primary delta | max secondary delta | max p delta |
|---|---:|---:|---:|---:|---|---:|---:|---:|
| heterogeneity | 0.002000 | 0.004000 | 2.00x | 2/2 | TRUE | 5.684e-14 | 0.000e+00 | 2.746e-35 |
| egger_pleiotropy | 0.001000 | 0.001000 | 1.00x | 1/1 | TRUE | 0.000e+00 | 2.168e-19 | 6.939e-17 |
| single_snp | 0.002500 | 0.004000 | 1.60x | 84/84 | TRUE/TRUE | 4.337e-18 | 1.735e-18 | 3.886e-16 |
| leave_one_out | 0.007000 | 0.019000 | 2.71x | 83/83 | TRUE/TRUE | 1.626e-18 | 6.939e-18 | 3.331e-16 |
| directionality | 0.000500 | 0.001000 | 2.00x | 1/1 | TRUE | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| steiger_filtering | 0.001000 | 0.001000 | 1.00x | 82/82 | TRUE/TRUE | 0.000e+00 | 0.000e+00 | 0.000e+00 |
