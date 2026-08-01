# fastMR implementation and scaling iteration

Prior native TwoSampleMR 50x50 reference: 338.919 seconds at nboot=100.

| component | nboot | threads | median seconds | minimum seconds | rows/pairs |
|---|---:|---:|---:|---:|---:|
| MR grid | 0 | 1 | 0.0850 | 0.0830 | 2500 |
| MR grid | 0 | 2 | 0.0440 | 0.0440 | 2500 |
| MR grid | 0 | 5 | 0.0240 | 0.0230 | 2500 |
| MR grid | 0 | 10 | 0.0170 | 0.0170 | 2500 |
| MR grid | 100 | 1 | 8.0660 | 8.0210 | 2500 |
| MR grid | 100 | 2 | 4.0790 | 4.0770 | 2500 |
| MR grid | 100 | 5 | 2.0780 | 2.0770 | 2500 |
| MR grid | 100 | 10 | 1.5330 | 1.5310 | 2500 |
| local LD-matrix clump | NA | NA | 0.0400 | 0.0370 | 1500 |
