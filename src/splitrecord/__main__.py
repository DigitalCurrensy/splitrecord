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

from .score import mann_kendall, mann_kendall_p, residual, sen_slope


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


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 2:
        print("usage: python -m splitrecord LEFT.csv RIGHT.csv", file=sys.stderr)
        return 2
    try:
        left = read_column(args[0])
        right = read_column(args[1])
        series = residual(left, right)
        slope = sen_slope(series)
        s = mann_kendall(series)
        p = mann_kendall_p(series)
    except (OSError, UnicodeError, ValueError) as exc:
        print(exc, file=sys.stderr)
        return 1
    p_text = "short" if p is None else f"{p:.6g}"
    print(f"n={len(series)} sen={slope:.10g} S={s} p={p_text}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
