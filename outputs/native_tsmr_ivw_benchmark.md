# Native TwoSampleMR IVW benchmark

IL6 fixture: 82 SNPs; nboot=100; fastMR threads=5; seed=20260802.

| scenario | pairs | fastMR s | TwoSampleMR s | speedup | max beta delta | max SE delta |
|---|---:|---:|---:|---:|---:|---:|
| balanced_50x50 | 2500 | 0.263 | 14.446 | 54.93x | 3.253e-18 | 3.036e-18 |
| one_exposure_250_outcomes | 250 | 0.024 | 1.447 | 60.29x | 1.952e-18 | 2.168e-18 |
| one_outcome_250_exposures | 250 | 0.025 | 1.427 | 57.08x | 1.735e-18 | 1.301e-18 |
| ten_exposures_100_outcomes | 1000 | 0.096 | 5.713 | 59.51x | 2.819e-18 | 3.036e-18 |
| hundred_exposures_ten_outcomes | 1000 | 0.095 | 5.764 | 60.67x | 2.602e-18 | 2.168e-18 |
