# SPLITRECORD

SPLITRECORD is for a hydrologist who already holds two official records of one basin.

The residual is z(A) minus z(B). Each z-score subtracts the sample mean and divides by the sample standard deviation. The sample standard deviation uses divisor n-1. Lengths must already match. A length under 2 raises ValueError ("not enough"). A sample standard deviation of 0 raises ValueError ("no spread"). Unequal lengths raise ValueError ("align first").

The size of the drift is the Theil-Sen slope, also called Sen's slope. It is the median of the pairwise slopes (v_j - v_i) / (j - i). An odd count takes the middle slope. An even count averages the two middle slopes. It needs at least two points. The unit is per row, not per date. There is no intercept and no confidence interval. The estimator does not assume a normal distribution. Its known breakdown point is 1 - 1/sqrt(2), about 29%. This library computes the median. It does not run a separate breakdown trial.

mann_kendall returns S only. S is up steps minus down steps. A tie adds nothing. S is a count, not a slope and not a p-value.

The p function is a normal approximation. It returns None when n < 8. The exact small-sample distribution is not computed. With no ties, var = n(n-1)(2n+5)/18. For each tie group of size t > 1, subtract t(t-1)(2t+5)/18. That variance is then multiplied by the Hamed-Rao factor, which is n/n*. If n < 3, the factor is 1. Otherwise the series is detrended by subtracting Sen's slope times the index, and the detrended series is ranked. Ties share the average rank. Ranks start at 1. For every lag i from 1 through n-1, rho uses one mean of the whole rank series. The numerator sums (rank_t - mean) times (rank_t+i - mean). The denominator sums (rank_t - mean) squared over every rank. If the ranks do not vary, rho is 0. Rho is clamped to [-0.999, 0.999]. Rho is kept only when abs(rho) is greater than 1.95996398454/sqrt(n). Otherwise that lag is treated as 0. The threshold is the two-sided 5% normal bound. The factor is 1 + (2 / (n*(n-1)*(n-2))) times the sum over every lag i of (n-i)*(n-i-1)*(n-i-2)*rho_i. The sum is over every lag, not a single lag and not a cap at 3. If the corrected variance is not positive, p is None. The continuity correction is z = (S - sign(S)) / sqrt(var), and sign is 0 when S is 0. p = erfc(|z| / sqrt(2)). This p is not a certificate. The worked file is not a basin study.

Do not average the records. Do not fetch granules. Do not sign a permit.


The printed line names every field.

```
rows=12 residual=z(A)-z(B) theil_sen_z_per_row=... mann_kendall_S=... variance=hamed-rao p=...
```

`variance=hamed-rao` means the tie-corrected variance was multiplied by the Hamed-Rao factor. That factor is for autocorrelation in one series. It does not remove a seasonal cycle. A January value is still compared with a July value. Under 8 rows, `p` is `short`. A non-positive corrected variance prints `p=dependent`.

A seasonal test is a separate command. `--seasons N` says row 0 is season 0 of year 0, row 1 is season 1 of that year, and so on. A monthly file uses `--seasons 12`. S is the sum of the within-season Mann-Kendall counts. The Theil-Sen slope is pooled inside seasons, so the unit is z per year, not z per row. The variance is the sum of the seasonal variances. Hamed-Rao is not applied on top of it.

```
rows=12 residual=z(A)-z(B) seasons=2 theil_sen_z_per_year=... seasonal_S=... variance=seasonal p=...
```

`--covariance` replaces that variance with the Hirsch-Slack estimator. The off-diagonal term is (K + 4 * sum of rank products - n * (n+1)^2) / 3. It needs a complete table, one value for every season of every year. A short last year raises ValueError ("uneven"). The estimator allows months inside one year to move together. It does not correct a strong correlation from one year to the next. Hirsch and Slack reported that the test is not reliable for very persistent series, or for a record of about five years.

```
pip install -e .
python -m unittest tests.test_kernel
python -m splitrecord examples/left.csv examples/right.csv
python -m splitrecord examples/left.csv examples/right.csv --seasons 12
```

Copyright 2026 Digital Currensy Inc. License Apache-2.0. LICENSE is unmodified. Copyright notice is in NOTICE and the file headers.
