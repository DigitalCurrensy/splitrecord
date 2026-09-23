# SPLITRECORD

SPLITRECORD is for a hydrologist who already holds two official records of one basin.

The residual is z(A) minus z(B). Each z-score subtracts the sample mean and divides by the sample standard deviation. The sample standard deviation uses divisor n-1. Lengths must already match. A length under 2 raises ValueError ("not enough"). A sample standard deviation of 0 raises ValueError ("no spread"). Unequal lengths raise ValueError ("align first").

Sen slope is the median pairwise slope of that residual. It needs at least two points.

mann_kendall returns S only. S is up steps minus down steps. S is a count, not a significance test.

The p function is a normal approximation. It returns None when n < 8. With no ties, var = n(n-1)(2n+5)/18. For each tie group of size t > 1, subtract t(t-1)(2t+5)/18. That variance is then multiplied by the Hamed-Rao factor, which is n/n*. If n < 3, the factor is 1. Otherwise the series is detrended by subtracting Sen's slope times the index, and the detrended series is ranked. Ties share the average rank. Ranks start at 1. For every lag i from 1 through n-1, rho is the Pearson correlation of those ranks with the same ranks shifted by i. If either side has no variation, rho is 0. Rho is clamped to [-0.999, 0.999]. Rho is kept only when abs(rho) is greater than 1.95996398454/sqrt(n). Otherwise that lag is treated as 0. The threshold is the two-sided 5% normal bound. The factor is 1 + (2 / (n*(n-1)*(n-2))) times the sum over every lag i of (n-i)*(n-i-1)*(n-i-2)*rho_i. The sum is over every lag, not a single lag and not a cap at 3. If the corrected variance is not positive, p is None. The continuity correction is z = (S - sign(S)) / sqrt(var), and sign is 0 when S is 0. p = erfc(|z| / sqrt(2)). This p is not a certificate. The worked file is not a basin study.

Do not average the records. Do not fetch granules. Do not sign a permit.

Copyright 2026 Digital Currensy Inc. License Apache-2.0. LICENSE is unmodified. Copyright notice is in NOTICE and the file headers.

```
PYTHONPATH=src python -m unittest tests.test_kernel
PYTHONPATH=src python -m splitrecord examples/left.csv examples/right.csv
```
