# Raw compact IVW grid scaling

IL6/CRP 82-SNP fixture; nboot=0; five requested native workers.
This measures the compact native grid before constructing the tidy R data frame.

| exposures | outcomes | pairs | repeats | seconds | pairs/s | layout |
|---:|---:|---:|---:|---:|---:|---|
| 50 | 50 | 2500 | 100 | 0.000260 | 9615385 | 1x2500 |
| 100 | 100 | 10000 | 100 | 0.000970 | 10309278 | 1x10000 |
| 250 | 250 | 62500 | 10 | 0.005900 | 10593220 | 1x62500 |
| 500 | 500 | 250000 | 5 | 0.024000 | 10416667 | 1x250000 |
| 1000 | 1000 | 1000000 | 3 | 0.109000 | 9174312 | 1x1000000 |
