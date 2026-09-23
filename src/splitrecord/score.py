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


def mann_kendall_variance(n: int) -> float:
    """Variance of Mann-Kendall S with no tie correction: n(n-1)(2n+5)/18."""
    return n * (n - 1) * (2 * n + 5) / 18


def mann_kendall_p(values: Sequence[float]) -> float | None:
    """Two-sided normal approximation to Mann-Kendall S.

    Returns None when n < 8, because a short series is not a trend test.
    Otherwise z = (S - sign(S)) / sqrt(var) with the continuity correction,
    and p = erfc(|z| / sqrt(2)). var is mann_kendall_variance(n). Observations
    are assumed independent. Autocorrelation is not corrected, and ties are
    not corrected. This p-value is not a certificate.
    """
    n = len(values)
    if n < 8:
        return None
    s = mann_kendall(values)
    var = mann_kendall_variance(n)
    sign = (s > 0) - (s < 0)
    z = (s - sign) / math.sqrt(var)
    return math.erfc(abs(z) / math.sqrt(2))
