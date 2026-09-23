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

"""Score two one-column CSV files. No header. One number per row."""

from __future__ import annotations

import math
import sys

from .score import (
    mann_kendall,
    mann_kendall_p,
    mann_kendall_p_ordinary,
    residual,
    seasonal_p,
    seasonal_s,
    seasonal_sen_slope,
    sen_slope,
    trend_free_prewhiten,
)


def read_column(path: str) -> list[float]:
    """Read a single numeric column. A bad row raises ValueError."""
    values: list[float] = []
    with open(path, encoding="utf-8") as handle:
        for lineno, line in enumerate(handle, start=1):
            text = line.strip()
            if text == "":
                raise ValueError(f"{path}:{lineno}: malformed row")
            try:
                value = float(text)
            except ValueError:
                raise ValueError(f"{path}:{lineno}: malformed row") from None
            if not math.isfinite(value):
                raise ValueError(f"{path}:{lineno}: malformed row")
            values.append(value)
    return values


def _parse(args: list[str]) -> tuple[str, str, int | None, bool, bool]:
    files: list[str] = []
    seasons: int | None = None
    covariance = False
    prewhiten = False
    i = 0
    while i < len(args):
        token = args[i]
        if token == "--seasons":
            if i + 1 >= len(args):
                raise ValueError("usage")
            try:
                seasons = int(args[i + 1])
            except ValueError:
                raise ValueError("bad season") from None
            if seasons < 2:
                raise ValueError("bad season")
            i += 2
            continue
        if token == "--covariance":
            covariance = True
            i += 1
            continue
        if token == "--prewhiten":
            prewhiten = True
            i += 1
            continue
        files.append(token)
        i += 1
    if len(files) != 2 or (covariance and seasons is None):
        raise ValueError("usage")
    if prewhiten and (seasons is not None or covariance):
        raise ValueError("separate")
    return files[0], files[1], seasons, covariance, prewhiten


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    try:
        left_path, right_path, seasons, covariance, prewhiten = _parse(args)
    except ValueError as exc:
        if str(exc) == "usage":
            print(
                "usage: python -m splitrecord LEFT.csv RIGHT.csv [--seasons N] [--covariance] [--prewhiten]",
                file=sys.stderr,
            )
            return 2
        print(exc, file=sys.stderr)
        return 1
    try:
        series = residual(read_column(left_path), read_column(right_path))
        whitened: list[float] | None = None
        r1 = 0.0
        if prewhiten:
            whitened, r1 = trend_free_prewhiten(series)
            slope = sen_slope(whitened)
            s = mann_kendall(whitened)
            p = mann_kendall_p_ordinary(whitened)
            tested = whitened
        elif seasons is None:
            slope = sen_slope(series)
            s = mann_kendall(series)
            p = mann_kendall_p(series)
            tested = series
        else:
            slope = seasonal_sen_slope(series, seasons)
            s = seasonal_s(series, seasons)
            p = seasonal_p(series, seasons, covariance=covariance)
            tested = series
    except (OSError, UnicodeError, ValueError) as exc:
        print(exc, file=sys.stderr)
        return 1
    if len(tested) < 8:
        p_text = "short"
    elif p is None:
        p_text = "dependent"
    else:
        p_text = f"{p:.6g}"
    if prewhiten:
        assert whitened is not None
        print(
            f"rows={len(series)} residual=z(A)-z(B) series=trend-free-prewhiten "
            f"whitened_rows={len(whitened)} r1={r1:.6g} "
            f"theil_sen_z_per_row={slope:.10g} mann_kendall_S={s} "
            f"variance=ordinary p={p_text}"
        )
    elif seasons is None:
        print(
            f"rows={len(series)} residual=z(A)-z(B) "
            f"theil_sen_z_per_row={slope:.10g} mann_kendall_S={s} "
            f"variance=hamed-rao p={p_text}"
        )
    else:
        variance = "hirsch-slack" if covariance else "seasonal"
        print(
            f"rows={len(series)} residual=z(A)-z(B) seasons={seasons} "
            f"theil_sen_z_per_year={slope:.10g} seasonal_S={s} "
            f"variance={variance} p={p_text}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
