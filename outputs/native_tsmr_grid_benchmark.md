# Native TwoSampleMR grid benchmark

Fixture: 82 rows; grid: 2500 pairs; methods: 5; nboot: 100; seed: 20260801.

| implementation | wall seconds | pairs/s |
|---|---:|---:|
| fastMR exact R/C++ (threads=10) | 4.945000 | 505.561 |
| TwoSampleMR native R workflow | 338.919000 | 7.376 |

Speedup: 68.538x.
Maximum absolute IVW beta difference: 3.469e-18.
