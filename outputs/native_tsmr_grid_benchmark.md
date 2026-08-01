# Native TwoSampleMR grid benchmark

Fixture: 82 rows; grid: 2500 pairs; methods: 5; nboot: 100; seed: 20260801.

| implementation | wall seconds | pairs/s |
|---|---:|---:|
| fastMR exact R/C++ (threads=10) | 1.564000 | 1598.465 |
| TwoSampleMR native R workflow | 339.139000 | 7.372 |

Speedup: 216.841x.
Maximum absolute IVW beta difference: 3.253e-18.
