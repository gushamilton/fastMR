# Native TwoSampleMR method benchmarks

IL6 fixture: 82 SNPs; nboot=100; fastMR threads=5; seed=20260802.
The native comparator is the installed TwoSampleMR R workflow, timed one method at a time.
Bootstrap SE differences are Monte Carlo differences when comparing a shared fastMR grid layout with native per-pair streams.

| scenario | method | pairs | fastMR s | TwoSampleMR s | speedup | max beta delta | median SE delta |
|---|---|---:|---:|---:|---:|---:|---:|
| balanced_50x50 | ivw | 2500 | 0.001 | 14.083 | 14083.00x | 3.253e-18 | 4.337e-19 |
| balanced_50x50 | egger | 2500 | 0.002 | 14.494 | 7247.00x | 8.674e-18 | 4.337e-19 |
| balanced_50x50 | weighted_median | 2500 | 0.068 | 17.296 | 254.35x | 4.337e-19 | 1.108e-04 |
| balanced_50x50 | simple_mode | 2500 | 1.286 | 173.085 | 134.59x | 2.455e-16 | 2.147e-03 |
| balanced_50x50 | weighted_mode | 2500 | 1.290 | 176.129 | 136.53x | 2.520e-16 | 2.310e-04 |
| one_exposure_250_outcomes | ivw | 250 | 0.001 | 1.498 | 1498.00x | 1.952e-18 | 4.337e-19 |
| one_exposure_250_outcomes | egger | 250 | 0.001 | 1.806 | 1806.00x | 7.373e-18 | 1.735e-18 |
| one_exposure_250_outcomes | weighted_median | 250 | 0.011 | 1.826 | 166.00x | 4.337e-19 | 1.041e-04 |
| one_exposure_250_outcomes | simple_mode | 250 | 0.157 | 18.038 | 114.89x | 1.448e-16 | 1.573e-03 |
| one_exposure_250_outcomes | weighted_mode | 250 | 0.158 | 18.015 | 114.02x | 2.290e-16 | 2.224e-04 |
| one_outcome_250_exposures | ivw | 250 | 0.001 | 1.628 | 1628.00x | 1.735e-18 | 4.337e-19 |
| one_outcome_250_exposures | egger | 250 | 0.001 | 1.581 | 1581.00x | 7.373e-18 | 4.337e-19 |
| one_outcome_250_exposures | weighted_median | 250 | 0.010 | 1.911 | 191.10x | 4.337e-19 | 1.098e-04 |
| one_outcome_250_exposures | simple_mode | 250 | 0.156 | 18.386 | 117.86x | 2.463e-16 | 1.866e-03 |
| one_outcome_250_exposures | weighted_mode | 250 | 0.211 | 17.988 | 85.25x | 2.515e-16 | 1.744e-04 |
| ten_exposures_100_outcomes | ivw | 1000 | 0.001 | 6.001 | 6001.00x | 2.819e-18 | 4.337e-19 |
| ten_exposures_100_outcomes | egger | 1000 | 0.001 | 6.164 | 6164.00x | 6.939e-18 | 8.674e-19 |
| ten_exposures_100_outcomes | weighted_median | 1000 | 0.033 | 7.248 | 219.64x | 4.337e-19 | 1.091e-04 |
| ten_exposures_100_outcomes | simple_mode | 1000 | 0.581 | 72.206 | 124.28x | 2.663e-16 | 2.135e-03 |
| ten_exposures_100_outcomes | weighted_mode | 1000 | 0.580 | 72.436 | 124.89x | 2.177e-16 | 2.177e-04 |
| hundred_exposures_ten_outcomes | ivw | 1000 | 0.001 | 5.974 | 5974.00x | 2.602e-18 | 4.337e-19 |
| hundred_exposures_ten_outcomes | egger | 1000 | 0.001 | 6.186 | 6186.00x | 7.806e-18 | 4.337e-19 |
| hundred_exposures_ten_outcomes | weighted_median | 1000 | 0.030 | 7.260 | 242.00x | 4.337e-19 | 1.043e-04 |
| hundred_exposures_ten_outcomes | simple_mode | 1000 | 0.574 | 72.146 | 125.69x | 2.472e-16 | 1.852e-03 |
| hundred_exposures_ten_outcomes | weighted_mode | 1000 | 0.618 | 72.508 | 117.33x | 2.520e-16 | 1.867e-04 |
