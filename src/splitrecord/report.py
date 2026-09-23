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

"""One labeled summary of the residual."""

from __future__ import annotations

from .score import mann_kendall, mann_kendall_p, residual, sen_slope

WORD_CAP = 80


def fnv1a_32(text: str) -> str:
    h = 2166136261
    for ch in text.encode("utf-8"):
        h ^= ch
        h = (h * 16777619) & 0xFFFFFFFF
    return f"{h:08x}"


def compile_report(left: list[float], right: list[float]) -> dict:
    r = residual(left, right)
    slope = sen_slope(r)
    s = mann_kendall(r)
    p = mann_kendall_p(r)
    if len(r) < 8:
        p_text = "short"
    elif p is None:
        p_text = "dependent"
    else:
        p_text = f"{p:.6g}"
    body = (
        f"Residual is z(A) minus z(B). rows {len(r)}. "
        f"Theil-Sen {slope:.4f} z per row. "
        f"Mann-Kendall S {s}. variance Hamed-Rao. p {p_text}."
    )
    words = len(body.split())
    return {
        "rows": len(r),
        "sen": slope,
        "mk": s,
        "variance": "hamed-rao",
        "p": p_text,
        "checksum": fnv1a_32(body),
        "body": body,
        "words": words,
        "fetched": False,
    }
