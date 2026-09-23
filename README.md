# SPLITRECORD

SPLITRECORD is for a hydrologist who already holds two official records of one basin.

The residual is z(A) minus z(B). Each z-score subtracts the sample mean and divides by the sample standard deviation. The sample standard deviation uses divisor n-1. Lengths must already match. A length under 2 raises ValueError ("not enough"). A sample standard deviation of 0 raises ValueError ("no spread"). Unequal lengths raise ValueError ("align first").

Sen slope is the median pairwise slope of that residual. It needs at least two points.

mann_kendall returns S only. S is up steps minus down steps. S is a count, not a significance test.

The p function is a normal approximation. It returns None when n < 8. With no ties, var = n(n-1)(2n+5)/18. For each tie group of size t > 1, subtract t(t-1)(2t+5)/18. Lag-1 is the Pearson correlation of the series with itself shifted by one, clamped to [-0.999, 0.999]. If there is no variation, lag-1 is 0. Only a positive lag-1 is used, and only as an AR(1) inflation, not the full Hamed-Rao sum. Then n_eff = n * (1 - r1) / (1 + r1). If n_eff < 8, p is withheld. If n_eff >= 8, var is multiplied by n / n_eff. A lag-1 that is not positive does not change var. Negative lag-1 does not shrink the variance. The continuity correction is z = (S - sign(S)) / sqrt(var), and p = erfc(|z| / sqrt(2)). A straight worked example prints p=dependent because the months move together. S and the Sen slope are still printed. Not a basin study. Not a certificate.

The worked months in examples/ are not a customer file and not InSAR or GRACE.

Do not average the records. Do not fetch granules. Do not sign a permit.

Copyright 2026 Digital Currensy Inc. License Apache-2.0. LICENSE is unmodified. Copyright notice is in NOTICE and the file headers.

```
PYTHONPATH=src python -m unittest tests.test_kernel
PYTHONPATH=src python -m splitrecord examples/left.csv examples/right.csv
```
