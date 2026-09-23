# SPLITRECORD

SPLITRECORD is for a hydrologist who already holds two official records of one basin.

The residual is z(A) minus z(B). Each z-score subtracts the sample mean and divides by the sample standard deviation. The sample standard deviation uses divisor n-1. Lengths must already match. A length under 2 raises ValueError ("not enough"). A sample standard deviation of 0 raises ValueError ("no spread"). Unequal lengths raise ValueError ("align first").

The size of the drift is the Theil-Sen slope, also called Sen's slope. It is the median of the pairwise slopes (v_j - v_i) / (j - i). An odd count takes the middle slope. An even count averages the two middle slopes. It needs at least two points. The unit is per row, not per date. There is no intercept and no confidence interval. The estimator does not assume a normal distribution. Its known breakdown point is 1 - 1/sqrt(2), about 29%. This library computes the median. It does not run a separate breakdown trial.

The 95% interval is not the median. Sort the same pairwise slopes. Let k be how many there are, and let V be the tie-corrected Mann-Kendall variance with no Hamed-Rao factor. C = 1.95996398454 × √V. The lower rank is round((k − C) / 2). The upper rank is round((k + C) / 2 + 1). Ranks start at 1. A half rounds to even. The slopes at those two ranks are sen95_lo and sen95_hi, and the line says sen95=normal. Under 8 rows with no ties, that normal rank is not used. The line says sen95=exact. q is the largest integer with P(C < q) <= 0.025, where C counts upward pairs and every ordering is equally likely. The ranks are q and N+1-q. A tie under 8 rows is short. Ranks outside the list are wide. hamed95 uses that same rank rule with V multiplied by n/n*. It is not Sen's interval.

mann_kendall returns S only. S is up steps minus down steps. A tie adds nothing. S is a count, not a slope and not a p-value.

Kendall's tau-b is S divided by the square root of (n(n-1)/2 - sum t(t-1)/2) times n(n-1)/2. The sum is over tie groups. The row index has no ties. If every value is tied, the line prints tau=tied.

The p function is a normal approximation. It returns None when n < 8. The exact small-sample distribution is not computed. With no ties, var = n(n-1)(2n+5)/18. For each tie group of size t > 1, subtract t(t-1)(2t+5)/18. That variance is then multiplied by the Hamed-Rao factor, which is n/n*. If n < 3, the factor is 1. Otherwise the series is detrended by subtracting Sen's slope times the index, and the detrended series is ranked. Ties share the average rank. Ranks start at 1. For every lag i from 1 through n-1, rho uses one mean of the whole rank series. The numerator sums (rank_t - mean) times (rank_t+i - mean). The denominator sums (rank_t - mean) squared over every rank. If the ranks do not vary, rho is 0. Rho is clamped to [-0.999, 0.999]. Rho is kept only when abs(rho) is greater than 1.95996398454/sqrt(n). Otherwise that lag is treated as 0. The threshold is the two-sided 5% normal bound. The factor is 1 + (2 / (n*(n-1)*(n-2))) times the sum over every lag i of (n-i)*(n-i-1)*(n-i-2)*rho_i. The sum is over every lag, not a single lag and not a cap at 3. The weight (n−i)(n−i−1)(n−i−2) is the count, from the 1998 derivation, of how the signs inside S covary at lag i. The fraction 2/(n(n−1)(n−2)) puts that count on the scale of the independent variance, so the correction is a multiplier rather than a second variance formula. A positive autocorrelation makes the multiplier larger than 1. A negative one makes it smaller. If the multiplier is not positive, z and p are dependent. If the corrected variance is not positive, p is None. The continuity correction is z = (S - sign(S)) / sqrt(var), and sign is 0 when S is 0. p = erfc(|z| / sqrt(2)). This p is not a certificate. The worked file is not a basin study.

Do not average the records. Do not fetch granules. Do not sign a permit.


The printed line names every field.

```
rows=12 residual=z(A)-z(B) theil_sen_z_per_row=0.5547001962 sen95=normal sen95_lo=0.5547001962 sen95_hi=0.5547001962 hamed95_lo=0.5547001962 hamed95_hi=0.5547001962 mann_kendall_S=66 tau=1 var=212.6666667 n_over_nstar=1 z=4.457215629 variance=hamed-rao p=8.30311e-06
```

`var` is the tie-corrected variance times `n_over_nstar`. `z` is (S - sign(S)) / sqrt(var). `p` is erfc(|z| / sqrt(2)). On this file the factor is 1: after the Sen slope is removed, the ranks do not vary, so no lag is kept. A repeating series is different. For `5, 5, 5, 0` repeated three times, S is -9, the ordinary variance is 117, `n_over_nstar` is 1.339393939, and `var` is 156.7090909. Yue and Wang's lag-1-only factor is not this number. The sum still runs over every lag. On 1, 1.4, 1.1, 2.2, 1.8, 2.9, 2.4, 3.6, 3.0, 4.1, 3.7, 4.8 the ordinary variance is 212.6666667 and n/n* is 0.08583916084, so the corrected variance is 18.25512821. There are 66 pairwise slopes. Sen's C picks ranks 19 and 48, which are 0.25 and 0.38. The corrected C picks ranks 29 and 38, which are 0.3166666667 and 0.3454545455.

`variance=hamed-rao` means the tie-corrected variance was multiplied by the Hamed-Rao factor. That factor is for autocorrelation in one series. It does not remove a seasonal cycle. A January value is still compared with a July value. Under 8 rows, `p` is `short`. A non-positive corrected variance prints `p=dependent`.

A seasonal test is a separate command. `--seasons N` says row 0 is season 0 of year 0, row 1 is season 1 of that year, and so on. A monthly file uses `--seasons 12`. S is the sum of the within-season Mann-Kendall counts. The Theil-Sen slope is pooled inside seasons, so the unit is z per year, not z per row. The variance is the sum of the seasonal variances. Hamed-Rao is not applied on top of it.

```
rows=12 residual=z(A)-z(B) seasons=2 theil_sen_z_per_year=... seasonal_S=... variance=seasonal p=...
```

`--covariance` replaces that variance with the Hirsch-Slack estimator. The off-diagonal term is (K + 4 * sum of rank products - n * (n+1)^2) / 3. It needs a complete table, one value for every season of every year. A short last year raises ValueError ("uneven"). The estimator allows months inside one year to move together. It does not correct a strong correlation from one year to the next. Hirsch and Slack reported that the test is not reliable for very persistent series, or for a record of about five years.

`--prewhiten` is a third correction, and it is not stacked on the other two. It is the trend-free pre-whitening of Yue, Pilon, Phinney, and Cavadias (2002). The Theil-Sen slope is removed, lag-1 is estimated on that remainder with one mean and the full sum of squares, the remainder is whitened, and the slope is added back. The tested series is one row shorter. `variance=ordinary` means the tie-corrected Mann-Kendall variance with no Hamed-Rao factor. Combining `--prewhiten` with `--seasons` or `--covariance` raises `separate`.

von Storch pre-whitening is not implemented. That method removes lag-1 before removing the slope, and part of a real trend leaves with it.

## Pre-whitening example

`examples/pw_left.csv` is `0 1 2 3 4 5 6 7 8`. `examples/pw_right.csv` is `0 2 1 3 2 4 3 5 4`. Nine rows. The row index is the time, starting at 0. The rows are equally spaced. There is no date column.

Sen's slope is the median of all 36 pairwise slopes, `36 = 9×8/2`. Each slope is `(v_j − v_i) / (j − i)`. Thirty-six is even, so the median is the average of the 18th and 19th slopes after sorting. In this file those two slopes are the same number, `0.04892060565`. That is `removed_sen`.

The four steps on the residual `z(A)−z(B)`:

1. Remove that slope. Row `t` becomes `residual_t − removed_sen × t`. The remainder alternates `0.225955` and `-0.722729`.
2. Lag-1 of that remainder is `-0.888889`. One mean, full sum of squares. It is negative because each row is the opposite of the row before it.
3. Whiten from the second row: `remainder_t − r1 × remainder_{t−1}`. One row is used up.
4. Add the same slope back at the same index: whitened_t `+ removed_sen × t`.

```
python -m splitrecord examples/pw_left.csv examples/pw_right.csv --prewhiten
```

```
rows=9 residual=z(A)-z(B) series=trend-free-prewhiten whitened_rows=8 removed_sen=0.04892060565 r1=-0.888889 theil_sen_z_per_row=0.04892060565 sen95=normal sen95_lo=0.02783875459 sen95_hi=0.0840570241 mann_kendall_S=22 tau=0.7857142857 var=65.33333333 z=2.598076211 variance=ordinary p=0.00937477
```

`removed_sen` is the slope that was taken off. `theil_sen_z_per_row` is the Sen slope of the blended series. They match in this file. They do not match in every file. `variance=ordinary` is not Hamed-Rao.

```
pip install -e .
python -m unittest tests.test_kernel
python -m splitrecord examples/left.csv examples/right.csv
python -m splitrecord examples/left.csv examples/right.csv --seasons 12
python -m splitrecord examples/pw_left.csv examples/pw_right.csv --prewhiten
```

Copyright 2026 Digital Currensy Inc. License Apache-2.0. LICENSE is unmodified. Copyright notice is in NOTICE and the file headers.
