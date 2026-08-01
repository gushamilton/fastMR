# Native TwoSampleMR MR-Egger bootstrap benchmark

IL6 fixture: 82 SNPs; nboot=100; fastMR threads=5; seed=20260806.

| scenario | pairs | fastMR s | TwoSampleMR s | speedup | max beta delta | median SE delta |
|---|---:|---:|---:|---:|---:|---:|
| balanced_50x50 | 2500 | 1.380 | 18.293 | 13.26x | 6.374e-04 | 1.143e-04 |
| one_exposure_250_outcomes | 250 | 0.140 | 1.869 | 13.35x | 6.374e-04 | 1.112e-04 |
| one_outcome_250_exposures | 250 | 0.139 | 1.843 | 13.26x | 6.374e-04 | 1.112e-04 |
| ten_exposures_100_outcomes | 1000 | 0.574 | 7.441 | 12.96x | 6.374e-04 | 1.106e-04 |
| hundred_exposures_ten_outcomes | 1000 | 0.577 | 7.406 | 12.84x | 6.374e-04 | 1.106e-04 |
