"""Small, owned statistics. Replace later; do not wrap a mystery package as the product."""
from typing import Sequence


def residual(a: Sequence[float], b: Sequence[float]) -> list[float]:
    if len(a) != len(b):
        raise ValueError("align first")
    return [x - y for x, y in zip(a, b)]


def sen_slope(values: Sequence[float]) -> float:
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
    s = 0
    n = len(values)
    for i in range(n):
        for j in range(i + 1, n):
            d = values[j] - values[i]
            s += 1 if d > 0 else -1 if d < 0 else 0
    return s
