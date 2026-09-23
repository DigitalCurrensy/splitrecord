# Copyright 2026 Digital Currensy Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Small, owned statistics. Replace later; do not wrap a mystery package as the product."""

from __future__ import annotations

import math
from typing import Sequence


def zscores(values: Sequence[float]) -> list[float]:
    """Sample z-scores.

    Subtract the sample mean and divide by the sample standard deviation.
    The standard deviation is the square root of the sum of squared
    deviations divided by n-1. Requires len(values) >= 2. A standard
    deviation of 0 raises ValueError("no spread").
    """
    n = len(values)
    if n < 2:
        raise ValueError("not enough")
    mean = sum(values) / n
    variance = sum((x - mean) ** 2 for x in values) / (n - 1)
    if variance == 0.0:
        raise ValueError("no spread")
    sigma = math.sqrt(variance)
    return [(x - mean) / sigma for x in values]


def residual(a: Sequence[float], b: Sequence[float]) -> list[float]:
    """Return z(A) minus z(B), elementwise.

    Unequal lengths raise ValueError("align first"). Fewer than two
    points raise ValueError("not enough"). A constant series raises
    ValueError("no spread").
    """
    if len(a) != len(b):
        raise ValueError("align first")
    if len(a) < 2:
        raise ValueError("not enough")
    za = zscores(a)
    zb = zscores(b)
    return [x - y for x, y in zip(za, zb)]


def sen_slope(values: Sequence[float]) -> float:
    """Median pairwise slope. Requires at least two points."""
    n = len(values)
    slopes = []
    for i in range(n):
        for j in range(i + 1, n):
            denom = j - i
            slopes.append((values[j] - values[i]) / denom)
    slopes.sort()
    m = len(slopes)
    if m == 0:
        raise ValueError("not enough")
    if m % 2:
        return slopes[m // 2]
    return 0.5 * (slopes[m // 2 - 1] + slopes[m // 2])


def mann_kendall(values: Sequence[float]) -> int:
    """Return S, up steps minus down steps.

    S counts the sign of later minus earlier over every pair. Ties add
    nothing. S is a count, not a p-value and not a significance test.
    """
    s = 0
    n = len(values)
    for i in range(n):
        for j in range(i + 1, n):
            d = values[j] - values[i]
            s += 1 if d > 0 else -1 if d < 0 else 0
    return s


def mann_kendall_variance(n: int, ties: Sequence[int] = ()) -> float:
    """Variance of Mann-Kendall S.

    With no ties, var = n(n-1)(2n+5)/18. For each tie-group size t > 1,
    subtract t(t-1)(2t+5)/18. A call with only n is that no-tie formula.
    """
    var = n * (n - 1) * (2 * n + 5) / 18
    for t in ties:
        if t > 1:
            var -= t * (t - 1) * (2 * t + 5) / 18
    return var


def tie_counts(values: Sequence[float]) -> list[int]:
    """Return the sizes of values that occur more than once."""
    counts: dict[float, int] = {}
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    return [size for size in counts.values() if size > 1]


def lag1(values: Sequence[float]) -> float:
    """Pearson correlation of the series with itself shifted by one.

    Clamped to [-0.999, 0.999]. If there is no variation, return 0.0.
    """
    n = len(values)
    if n < 2:
        return 0.0
    left = values[:-1]
    right = values[1:]
    m = n - 1
    mean_left = sum(left) / m
    mean_right = sum(right) / m
    cov = 0.0
    var_left = 0.0
    var_right = 0.0
    for a, b in zip(left, right):
        da = a - mean_left
        db = b - mean_right
        cov += da * db
        var_left += da * da
        var_right += db * db
    if var_left == 0.0 or var_right == 0.0:
        return 0.0
    r = cov / math.sqrt(var_left * var_right)
    if r > 0.999:
        return 0.999
    if r < -0.999:
        return -0.999
    return r


def average_ranks(values: Sequence[float]) -> list[float]:
    """Ranks starting at 1. Tied values share the average rank."""
    n = len(values)
    order = sorted(range(n), key=lambda i: values[i])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and values[order[j + 1]] == values[order[i]]:
            j += 1
        average = ((i + 1) + (j + 1)) / 2.0
        for k in range(i, j + 1):
            ranks[order[k]] = average
        i = j + 1
    return ranks


def detrend(values: Sequence[float]) -> list[float]:
    """Subtract Sen's slope times the zero-based index.

    A series shorter than 2 raises ValueError("not enough").
    """
    if len(values) < 2:
        raise ValueError("not enough")
    slope = sen_slope(values)
    return [float(value) - slope * index for index, value in enumerate(values)]


def rank_autocorr(ranks: Sequence[float], lag: int) -> float:
    """Pearson correlation of ranks[:-lag] with ranks[lag:].

    If either side has no variation, return 0.0. Clamp to [-0.999, 0.999].
    """
    n = len(ranks)
    if lag <= 0 or lag >= n:
        return 0.0
    left = ranks[:-lag]
    right = ranks[lag:]
    m = len(left)
    if m == 0:
        return 0.0
    mean_left = sum(left) / m
    mean_right = sum(right) / m
    cov = 0.0
    var_left = 0.0
    var_right = 0.0
    for a, b in zip(left, right):
        da = a - mean_left
        db = b - mean_right
        cov += da * db
        var_left += da * da
        var_right += db * db
    if var_left == 0.0 or var_right == 0.0:
        return 0.0
    r = cov / math.sqrt(var_left * var_right)
    if r > 0.999:
        return 0.999
    if r < -0.999:
        return -0.999
    return r


def hamed_rao_factor(values: Sequence[float]) -> float:
    """Hamed and Rao (1998) variance factor n/n*.

    If n < 3, return 1.0. Otherwise detrend with Sen's slope, rank that
    series, and for every lag i from 1 through n-1 compute the rank
    autocorrelation. Keep rho only when abs(rho) > 1.95996398454 / sqrt(n);
    otherwise treat that lag as 0. The threshold is the two-sided 5% normal
    bound. The factor is

        1 + (2 / (n*(n-1)*(n-2))) * sum_i (n-i)*(n-i-1)*(n-i-2)*rho_i

    summed over every lag. The sum is not limited to the first lag and is not capped at 3.
    """
    n = len(values)
    if n < 3:
        return 1.0
    ranks = average_ranks(detrend(values))
    threshold = 1.95996398454 / math.sqrt(n)
    total = 0.0
    for i in range(1, n):
        rho = rank_autocorr(ranks, i)
        if abs(rho) <= threshold:
            rho = 0.0
        total += (n - i) * (n - i - 1) * (n - i - 2) * rho
    return 1.0 + (2.0 / (n * (n - 1) * (n - 2))) * total


def mann_kendall_p(values: Sequence[float]) -> float | None:
    """Two-sided normal approximation to Mann-Kendall S.

    Returns None when n < 8. var is the tie-corrected Mann-Kendall variance
    multiplied by the Hamed-Rao factor over every lag. If that corrected
    variance is not positive, return None. A clean straight line has no
    leftover rank correlation after Sen detrending, so the factor is 1 and
    p is the ordinary normal approximation. Then
    z = (S - sign(S)) / sqrt(var), with sign 0 when S is 0, and
    p = erfc(|z| / sqrt(2)). This p-value is not a certificate.
    """
    n = len(values)
    if n < 8:
        return None
    s = mann_kendall(values)
    var = mann_kendall_variance(n, tie_counts(values)) * hamed_rao_factor(values)
    if var <= 0.0:
        return None
    sign = (s > 0) - (s < 0)
    z = (s - sign) / math.sqrt(var)
    return math.erfc(abs(z) / math.sqrt(2))
