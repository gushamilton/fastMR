# Native TwoSampleMR grid benchmark

Fixture: 82 rows; grid: 2500 pairs; methods: 5; nboot: 100; seed: 20260801.

| implementation | wall seconds | pairs/s |
|---|---:|---:|
| fastMR exact R/C++ (threads=10) | 7.946000 | 314.624 |
| TwoSampleMR native R workflow | 339.105000 | 7.372 |

Speedup: 42.676x.
Maximum absolute IVW beta difference: 3.253e-18.
