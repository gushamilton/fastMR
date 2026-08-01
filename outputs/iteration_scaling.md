# fastMR implementation and scaling iteration

Prior native TwoSampleMR 50x50 reference: 338.919 seconds at nboot=100.

| component | nboot | threads | median seconds | minimum seconds | rows/pairs |
|---|---:|---:|---:|---:|---:|
| MR grid | 0 | 1 | 0.5640 | 0.5600 | 2500 |
| MR grid | 0 | 2 | 0.4920 | 0.4870 | 2500 |
| MR grid | 0 | 5 | 0.4460 | 0.4450 | 2500 |
| MR grid | 0 | 10 | 0.4350 | 0.4330 | 2500 |
| MR grid | 100 | 1 | 15.4960 | 15.4810 | 2500 |
| MR grid | 100 | 2 | 8.1450 | 8.1400 | 2500 |
| MR grid | 100 | 5 | 4.0180 | 3.9640 | 2500 |
| MR grid | 100 | 10 | 2.9740 | 2.9640 | 2500 |
| local LD-matrix clump | NA | NA | 0.0360 | 0.0350 | 1500 |
