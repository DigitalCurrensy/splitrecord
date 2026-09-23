# SPLITRECORD

SPLITRECORD is for a hydrologist who already holds two official records of one basin.

The residual is z(A) minus z(B). Each z-score subtracts the sample mean and divides by the sample standard deviation. The sample standard deviation uses divisor n-1. Lengths must already match. A length under 2 raises ValueError ("not enough"). A sample standard deviation of 0 raises ValueError ("no spread"). Unequal lengths raise ValueError ("align first").

Sen slope is the median pairwise slope of that residual. It needs at least two points.

mann_kendall returns S only. S is up steps minus down steps. S is a count, not a significance test.

The p function is a normal approximation for n >= 8 under an independence assumption. It applies a continuity correction, z = (S - sign(S)) / sqrt(var), with var = n(n-1)(2n+5)/18, and p = erfc(|z| / sqrt(2)). It does not correct autocorrelation. It does not correct ties. It returns None when n < 8. It is not a certificate and not a basin study.

The worked months in examples/ are not a customer file and not InSAR or GRACE.

Do not average the records. Do not fetch granules. Do not sign a permit.

Copyright 2026 Digital Currensy Inc. License Apache-2.0. LICENSE is unmodified. Copyright notice is in NOTICE and the file headers.

```
PYTHONPATH=src python -m unittest tests.test_kernel
PYTHONPATH=src python -m splitrecord examples/left.csv examples/right.csv
```
