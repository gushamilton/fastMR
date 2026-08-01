# Native TwoSampleMR exact mode benchmark

IL6 fixture: 82 SNPs; nboot=100; fastMR threads=5; seed=20260818.

| scenario | method | pairs | fastMR s | TwoSampleMR s | speedup | max beta delta | max SE delta |
|---|---|---:|---:|---:|---:|---:|---:|
| balanced_50x50 | simple_mode | 2500 | 2.231 | 179.943 | 80.66x | 2.455e-16 | 7.050e-03 |
| balanced_50x50 | weighted_mode | 2500 | 2.423 | 178.288 | 73.58x | 2.520e-16 | 1.016e-03 |
| one_exposure_250_outcomes | simple_mode | 250 | 0.224 | 17.522 | 78.22x | 1.448e-16 | 5.910e-03 |
| one_exposure_250_outcomes | weighted_mode | 250 | 0.225 | 17.612 | 78.28x | 2.290e-16 | 9.642e-04 |
| one_outcome_250_exposures | simple_mode | 250 | 0.224 | 17.541 | 78.31x | 2.463e-16 | 5.372e-03 |
| one_outcome_250_exposures | weighted_mode | 250 | 0.224 | 17.535 | 78.28x | 2.515e-16 | 7.125e-04 |
| ten_exposures_100_outcomes | simple_mode | 1000 | 0.888 | 70.108 | 78.95x | 2.663e-16 | 6.395e-03 |
| ten_exposures_100_outcomes | weighted_mode | 1000 | 0.893 | 70.037 | 78.43x | 2.177e-16 | 1.035e-03 |
| hundred_exposures_ten_outcomes | simple_mode | 1000 | 0.885 | 70.397 | 79.54x | 2.472e-16 | 6.086e-03 |
| hundred_exposures_ten_outcomes | weighted_mode | 1000 | 0.904 | 69.998 | 77.43x | 2.520e-16 | 7.662e-04 |
