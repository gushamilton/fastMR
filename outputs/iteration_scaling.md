# fastMR implementation and scaling iteration

Prior native TwoSampleMR 50x50 reference: 339.139 seconds at nboot=100.

| component | nboot | threads | median seconds | minimum seconds | rows/pairs |
|---|---:|---:|---:|---:|---:|
| MR grid | 0 | 1 | 0.0840 | 0.0830 | 2500 |
| MR grid | 0 | 2 | 0.0430 | 0.0430 | 2500 |
| MR grid | 0 | 5 | 0.0230 | 0.0220 | 2500 |
| MR grid | 0 | 10 | 0.0170 | 0.0170 | 2500 |
| MR grid | 100 | 1 | 8.0630 | 8.0610 | 2500 |
| MR grid | 100 | 2 | 4.0860 | 4.0830 | 2500 |
| MR grid | 100 | 5 | 2.0810 | 2.0780 | 2500 |
| MR grid | 100 | 10 | 1.5440 | 1.5400 | 2500 |
| local LD-matrix clump | NA | NA | 0.0370 | 0.0360 | 1500 |
