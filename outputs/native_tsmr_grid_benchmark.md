# Native TwoSampleMR grid benchmark

Fixture: 82 rows; grid: 2500 pairs; methods: 5; nboot: 100; seed: 20260801.

| implementation | wall seconds | pairs/s |
|---|---:|---:|
| fastMR exact R/C++ (threads=10) | 1.580000 | 1582.278 |
| TwoSampleMR native R workflow | 334.830000 | 7.466 |

Speedup: 211.918x.
Maximum absolute IVW beta difference: 3.469e-18.
