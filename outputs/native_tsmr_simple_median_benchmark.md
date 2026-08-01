# Native TwoSampleMR Simple median benchmark

IL6 fixture: 82 SNPs; nboot=100; fastMR threads=5; seed=20260805.

| scenario | pairs | fastMR s | TwoSampleMR s | speedup | max beta delta | median SE delta |
|---|---:|---:|---:|---:|---:|---:|
| balanced_50x50 | 2500 | 1.345 | 17.365 | 12.91x | 3.296e-17 | 5.295e-04 |
| one_exposure_250_outcomes | 250 | 0.144 | 1.753 | 12.17x | 3.296e-17 | 5.831e-04 |
| one_outcome_250_exposures | 250 | 0.138 | 1.762 | 12.77x | 3.296e-17 | 5.831e-04 |
| ten_exposures_100_outcomes | 1000 | 0.561 | 7.313 | 13.04x | 3.296e-17 | 5.374e-04 |
| hundred_exposures_ten_outcomes | 1000 | 0.586 | 7.260 | 12.39x | 3.296e-17 | 5.374e-04 |
