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
    _finite(values)
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


# Two-sided 5% standard-normal bound, Φ^{-1}(0.975).
Z_95 = 1.95996398454


def pairwise_slopes(
    values: Sequence[float], times: Sequence[float] | None = None
) -> list[float]:
    """Every (v_j - v_i) / (t_j - t_i), sorted ascending.

    Times default to the row index, so equally spaced rows use j - i.
    A repeated time raises ValueError("tied time"). A decreasing time
    raises ValueError("time order"). Requires at least two finite points.
    """
    n = len(values)
    if n < 2:
        raise ValueError("not enough")
    if times is None:
        times = list(range(n))
    if len(times) != n:
        raise ValueError("align first")
    _finite(values)
    _finite(list(times))
    slopes = []
    for i in range(n):
        for j in range(i + 1, n):
            gap = float(times[j]) - float(times[i])
            if gap == 0.0:
                raise ValueError("tied time")
            if gap < 0.0:
                raise ValueError("time order")
            slopes.append((values[j] - values[i]) / gap)
    slopes.sort()
    return slopes


def sen_slope(values: Sequence[float], times: Sequence[float] | None = None) -> float:
    """Median of the pairwise slopes. The Sen (1968) estimator.

    An odd count takes the middle slope. An even count averages the two
    middle slopes. There is no intercept. times defaults to the row index.
    """
    slopes = pairwise_slopes(values, times)
    m = len(slopes)
    if m % 2:
        return slopes[m // 2]
    return 0.5 * (slopes[m // 2 - 1] + slopes[m // 2])



def _inversion_counts(n: int) -> list[int]:
    """Number of permutations of n with each inversion count. Index is the count."""
    cap = n * (n - 1) // 2
    prev = [0] * (cap + 1)
    prev[0] = 1
    for m in range(1, n):
        # Adding the next item (m+1 items total) creates 0..m new inversions.
        width = (m + 1) * m // 2
        cur = [0] * (width + 1)
        for inversions, count in enumerate(prev):
            if count == 0:
                continue
            for extra in range(m + 1):
                cur[inversions + extra] += count
        prev = cur
    return prev


def exact_slope_ranks(n: int, alpha: float = 0.025) -> tuple[int, int] | None:
    """1-based ranks for the exact no-tie Sen interval.

    C is the number of upward pairs, N minus the inversion count.
    q is the largest integer with P(C < q) <= alpha. The ranks are q and
    N + 1 - q. None when no such q exists.
    """
    if n < 3:
        return None
    counts = _inversion_counts(n)
    total = float(sum(counts))
    pairs = n * (n - 1) // 2
    best: int | None = None
    running = 0
    # counts[k] is P-weight of C = pairs - k, walking C from 0 upward.
    for concordant in range(0, pairs + 1):
        inversions = pairs - concordant
        running += counts[inversions]
        # P(C < concordant + 1) == running / total after this concordant is included? 
        # After adding C=concordant, running/total = P(C <= concordant) = P(C < concordant+1)
        if running / total <= alpha:
            best = concordant + 1
        else:
            break
    if best is None:
        return None
    return best, pairs + 1 - best


def sen_exact_limits(
    values: Sequence[float], times: Sequence[float] | None = None
) -> tuple[str, str]:
    """Exact no-tie 95% Sen interval. Ties return ("short", "short").

    This is not the normal approximation. n under 3 returns short.
    times defaults to the row index. The ranks do not depend on the times.
    The slope values do.
    """
    n = len(values)
    if n < 3 or tie_counts(values):
        return "short", "short"
    ranks = exact_slope_ranks(n)
    if ranks is None:
        return "wide", "wide"
    slopes = pairwise_slopes(values, times)
    lo, hi = ranks
    if lo < 1 or hi > len(slopes) or lo > hi:
        return "wide", "wide"
    return f"{slopes[lo - 1]:.10g}", f"{slopes[hi - 1]:.10g}"


def sen_limits(values: Sequence[float], variance: float) -> tuple[str, str]:
    """95% limits from rounded ranks of the sorted slopes.

    C = Z_95 * sqrt(variance). With k slopes,
    rank_lo = round((k - C) / 2) and rank_up = round((k + C) / 2 + 1).
    Ranks are 1-based. A half rounds to even. The slopes at those ranks
    are the limits. This is the index rule used by the trend package.
    It does not interpolate. gilbert_limits does.

    The caller chooses variance. The tie-corrected Mann-Kendall variance
    is Sen's interval. That same variance times the Hamed-Rao factor is
    not Sen's interval. n < 8 returns ("short", "short"). A variance that
    is not positive returns ("dependent", "dependent"). Ranks outside 1..k
    return ("wide", "wide").
    """
    n = len(values)
    if n < 8:
        return "short", "short"
    if not math.isfinite(variance) or variance <= 0.0:
        return "dependent", "dependent"
    slopes = pairwise_slopes(values)
    k = len(slopes)
    c = Z_95 * math.sqrt(variance)
    rank_lo = round((k - c) / 2.0)
    rank_up = round((k + c) / 2.0 + 1.0)
    if rank_lo < 1 or rank_up > k or rank_lo > rank_up:
        return "wide", "wide"
    return f"{slopes[rank_lo - 1]:.10g}", f"{slopes[rank_up - 1]:.10g}"



def gilbert_limits(values: Sequence[float], variance: float) -> tuple[str, str]:
    """Interpolated 95% Sen interval.

    C = Z_95 * sqrt(variance). With k slopes,
    M1 = (k - C) / 2 and M2 = (k + C) / 2. Each limit is the straight
    line between the slopes at floor(M) and ceiling(M). Ranks are 1-based.
    An integer M is that one slope. This is the interpolation in the
    mannkendall package, which cites Gilbert (1987). It is not a page
    copied from that book, and it is not the rounded rank in sen_limits.

    n < 8 returns ("short", "short"). A variance that is not positive
    returns ("dependent", "dependent"). A rank outside 1..k returns
    ("wide", "wide").
    """
    n = len(values)
    if n < 8:
        return "short", "short"
    if not math.isfinite(variance) or variance <= 0.0:
        return "dependent", "dependent"
    slopes = pairwise_slopes(values)
    k = len(slopes)
    c = Z_95 * math.sqrt(variance)

    def at(rank: float) -> float | None:
        if rank < 1.0 or rank > k:
            return None
        lower = math.floor(rank)
        upper = math.ceil(rank)
        left = slopes[lower - 1]
        right = slopes[upper - 1]
        return left + (right - left) * (rank - lower)

    lo = at((k - c) / 2.0)
    hi = at((k + c) / 2.0)
    if lo is None or hi is None:
        return "wide", "wide"
    return f"{lo:.10g}", f"{hi:.10g}"


def mann_kendall(values: Sequence[float]) -> int:
    """Return S, up steps minus down steps.

    S counts the sign of later minus earlier over every pair. Ties add
    nothing. S is a count, not a p-value and not a significance test.
    """
    _finite(values)
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
    """Lag-1 autocorrelation. One mean, full sum of squares.

    numer = sum (x_t - mean) * (x_{t+1} - mean) over t = 0..n-2.
    denom = sum (x_t - mean)^2 over the whole series.
    This is not a Pearson correlation of the two windows. No variation
    returns 0. The result is clamped to [-0.999, 0.999].
    """
    n = len(values)
    if n < 2:
        return 0.0
    mean = sum(values) / n
    denom = 0.0
    for value in values:
        d = value - mean
        denom += d * d
    if denom == 0.0:
        return 0.0
    numer = 0.0
    for t in range(n - 1):
        numer += (values[t] - mean) * (values[t + 1] - mean)
    r = numer / denom
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


def _finite(values: Sequence[float]) -> None:
    for value in values:
        if not math.isfinite(value):
            raise ValueError("bad number")


def rank_autocorr(ranks: Sequence[float], lag: int) -> float:
    """Autocorrelation of ranks at one lag.

    One mean, taken over the whole rank series. The numerator is the sum of
    (rank_t - mean) * (rank_{t+lag} - mean) for t = 1..n-lag. The denominator
    is the sum of (rank_t - mean)^2 over the whole series. That is the
    autocorrelation in Hamed and Rao, not a Pearson correlation of the two
    windows. No variation returns 0. The result is clamped to [-0.999, 0.999].
    """
    n = len(ranks)
    if lag <= 0 or lag >= n:
        return 0.0
    mean = sum(ranks) / n
    denom = 0.0
    for rank in ranks:
        d = rank - mean
        denom += d * d
    if denom == 0.0:
        return 0.0
    numer = 0.0
    for t in range(n - lag):
        numer += (ranks[t] - mean) * (ranks[t + lag] - mean)
    r = numer / denom
    if r > 0.999:
        return 0.999
    if r < -0.999:
        return -0.999
    return r


def hamed_rao_factor(values: Sequence[float]) -> float:
    """Hamed and Rao (1998) variance factor n/n*.

    If n < 3, return 1.0. Otherwise detrend with Sen's slope, rank that
    series, and for every lag i from 1 through n-1 compute the rank
    autocorrelation. rho is the full-series autocorrelation defined in rank_autocorr. Keep rho only when abs(rho) > 1.95996398454 / sqrt(n);
    otherwise treat that lag as 0. The threshold is the two-sided 5% normal
    bound. The factor is

        1 + (2 / (n*(n-1)*(n-2))) * sum_i (n-i)*(n-i-1)*(n-i-2)*rho_i

    summed over every lag. The sum is not limited to the first lag and is not capped at 3.
    """
    n = len(values)
    if n < 3:
        return 1.0
    ranks = average_ranks(detrend(values))
    threshold = Z_95 / math.sqrt(n)
    total = 0.0
    for i in range(1, n):
        rho = rank_autocorr(ranks, i)
        if abs(rho) <= threshold:
            rho = 0.0
        total += (n - i) * (n - i - 1) * (n - i - 2) * rho
    return 1.0 + (2.0 / (n * (n - 1) * (n - 2))) * total


def kendall_tau(values: Sequence[float]) -> float | None:
    """Kendall's tau-b.

    The row index has no ties, so

        tau = S / sqrt( (n(n-1)/2 - sum t(t-1)/2) * n(n-1)/2 )

    The sum is over tie-group sizes t > 1. None means every value is tied,
    so the denominator is 0. Tau is not a p-value.
    """
    n = len(values)
    if n < 2:
        raise ValueError("not enough")
    _finite(values)
    n0 = n * (n - 1) / 2.0
    tie_pairs = sum(t * (t - 1) / 2.0 for t in tie_counts(values))
    comparable = n0 - tie_pairs
    if comparable <= 0.0:
        return None
    return mann_kendall(values) / math.sqrt(comparable * n0)


def mann_kendall_z(
    values: Sequence[float], *, hamed: bool = True
) -> tuple[int, float | None, float, float, float | None, float | None]:
    """Return S, tau, corrected variance, n/n*, z, and p.

    n/n* is the Hamed-Rao factor when hamed is true, otherwise 1.
    The corrected variance is the tie-corrected variance times that factor.
    z and p are None when n < 8 or the corrected variance is not positive.
    z = (S - sign(S)) / sqrt(var), and sign is 0 when S is 0.
    p = erfc(|z| / sqrt(2)). Yue and Wang's lag-1 formula is not this factor.
    """
    n = len(values)
    s = mann_kendall(values)
    tau = kendall_tau(values) if n >= 2 else None
    base = mann_kendall_variance(n, tie_counts(values))
    factor = hamed_rao_factor(values) if hamed else 1.0
    var = base * factor
    if n < 8 or var <= 0.0:
        return s, tau, var, factor, None, None
    sign = (s > 0) - (s < 0)
    z = (s - sign) / math.sqrt(var)
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return s, tau, var, factor, z, p


def mann_kendall_p(values: Sequence[float]) -> float | None:
    """Two-sided p from mann_kendall_z with the Hamed-Rao factor.

    Returns None when n < 8 or the corrected variance is not positive.
    This p-value is not a certificate.
    """
    return mann_kendall_z(values, hamed=True)[5]


def season_groups(values: Sequence[float], period: int) -> list[list[float]]:
    """Split a chronological series into seasons.

    Row 0 is season 0 of the first year. Row 1 is season 1 of that same year.
    A monthly file is January, February, and so on, then January again.
    period must be >= 2. A non-finite value raises ValueError("bad number").
    """
    if period < 2:
        raise ValueError("bad season")
    _finite(values)
    groups: list[list[float]] = [[] for _ in range(period)]
    for index, value in enumerate(values):
        groups[index % period].append(float(value))
    return groups


def seasonal_s(values: Sequence[float], period: int) -> int:
    """Seasonal Mann-Kendall S. Sum of the within-season S values.

    Each season is compared only with later rows of the same season.
    A different season is not a pair. Requires at least two years,
    meaning len(values) >= 2 * period. Otherwise ValueError("not enough").
    """
    if len(values) < 2 * period:
        raise ValueError("not enough")
    return sum(mann_kendall(group) for group in season_groups(values, period))


def seasonal_sen_slope(values: Sequence[float], period: int) -> float:
    """Theil-Sen slope pooled inside seasons. Unit is z per year-step.

    Pairwise slopes are formed only inside one season. The denominator is
    the year gap, not the row gap. The result is the median of that pool.
    An even count averages the two middle slopes. Requires two years.
    """
    if len(values) < 2 * period:
        raise ValueError("not enough")
    slopes: list[float] = []
    for group in season_groups(values, period):
        n = len(group)
        for i in range(n):
            for j in range(i + 1, n):
                slopes.append((group[j] - group[i]) / (j - i))
    if not slopes:
        raise ValueError("not enough")
    slopes.sort()
    m = len(slopes)
    if m % 2:
        return slopes[m // 2]
    return 0.5 * (slopes[m // 2 - 1] + slopes[m // 2])


def seasonal_variance(values: Sequence[float], period: int, covariance: bool = False) -> float:
    """Variance of seasonal Mann-Kendall S.

    covariance False is Hirsch, Slack, and Smith (1982): the sum of the
    tie-corrected variances of the seasons. Seasons are treated as independent.
    Hamed-Rao is not applied. Those are different corrections.

    covariance True is Hirsch and Slack (1984). Off-diagonal terms are the
    Dietz and Killeen estimator

        (K_gh + 4 * sum_i R_ig * R_ih - n * (n + 1)^2) / 3

    K_gh is the sum of sign products between the two seasons across years.
    R is the within-season rank, ties sharing the average rank. The diagonal
    stays the ordinary tie-corrected variance. This form needs a complete
    years-by-seasons table. If len(values) is not a multiple of period,
    raise ValueError("uneven"). The estimator assumes dependence inside a
    year, not a strong correlation from one year to the next.
    """
    if len(values) < 2 * period:
        raise ValueError("not enough")
    groups = season_groups(values, period)
    if not covariance:
        return sum(mann_kendall_variance(len(group), tie_counts(group)) for group in groups)
    if len(values) % period != 0:
        raise ValueError("uneven")
    years = len(values) // period
    ranks = [average_ranks(group) for group in groups]
    total = 0.0
    for g in range(period):
        for h in range(period):
            if g == h:
                total += mann_kendall_variance(years, tie_counts(groups[g]))
                continue
            k = 0
            for i in range(years - 1):
                for j in range(i + 1, years):
                    dg = groups[g][j] - groups[g][i]
                    dh = groups[h][j] - groups[h][i]
                    sg = 1 if dg > 0 else -1 if dg < 0 else 0
                    sh = 1 if dh > 0 else -1 if dh < 0 else 0
                    k += sg * sh
            sum_rr = sum(ranks[g][i] * ranks[h][i] for i in range(years))
            total += (k + 4.0 * sum_rr - years * (years + 1) ** 2) / 3.0
    return total


def seasonal_zp(
    values: Sequence[float], period: int, covariance: bool = False
) -> tuple[float | None, float | None, float]:
    """Return z, p, and the seasonal variance.

    z and p are None when len(values) < 8 or the variance is not positive.
    The Hamed-Rao factor is not multiplied in.
    """
    var = seasonal_variance(values, period, covariance=covariance)
    if len(values) < 8 or var <= 0.0:
        return None, None, var
    s = seasonal_s(values, period)
    sign = (s > 0) - (s < 0)
    z = (s - sign) / math.sqrt(var)
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return z, p, var


def seasonal_p(values: Sequence[float], period: int, covariance: bool = False) -> float | None:
    """Two-sided p from seasonal_zp. None when the normal approximation is refused."""
    return seasonal_zp(values, period, covariance=covariance)[1]


def mann_kendall_p_ordinary(values: Sequence[float]) -> float | None:
    """Two-sided p from mann_kendall_z with no Hamed-Rao factor."""
    return mann_kendall_z(values, hamed=False)[5]


def trend_free_prewhiten(values: Sequence[float]) -> tuple[list[float], float, float]:
    """Yue, Pilon, Phinney, and Cavadias (2002) trend-free pre-whitening.

    Returns the blended series, the lag-1 of the detrended remainder, and
    the Theil-Sen slope that was removed. Time is the row index, starting
    at 0. The blended series has length n-1. Requires len >= 3.

    Sen's slope is the median of every pairwise slope (v_j - v_i) / (j - i).
    An even count averages the two middle slopes. That is the slope removed
    here. The slope of the blended series is computed again by the caller
    and can differ.

    von Storch pre-whitening, which removes lag-1 before removing the slope,
    is not this function. That order deletes part of the trend. Hamed-Rao
    is not applied to the result. The caller uses mann_kendall_p_ordinary.
    """
    if len(values) < 3:
        raise ValueError("not enough")
    _finite(values)
    beta = sen_slope(values)
    detrended = [float(value) - beta * index for index, value in enumerate(values)]
    r1 = lag1(detrended)
    blended: list[float] = []
    for t in range(1, len(values)):
        whitened = detrended[t] - r1 * detrended[t - 1]
        blended.append(whitened + beta * t)
    return blended, r1, beta
