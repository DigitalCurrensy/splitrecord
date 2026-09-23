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

"""Cite and checksum. No granules. Counsel stays unsigned."""

from __future__ import annotations

from .score import mann_kendall, residual, sen_slope

OFFER = "Unsigned. Not an invoice."
WORD_CAP = 80
HUC8 = "18030012"
BASIN = "Tulare Lake"


def fnv1a_32(text: str) -> str:
    h = 2166136261
    for ch in text.encode("utf-8"):
        h ^= ch
        h = (h * 16777619) & 0xFFFFFFFF
    return f"{h:08x}"


def compile_report(left: list[float], right: list[float]) -> dict:
    r = residual(left, right)
    body = (
        f"SPLITRECORD. {BASIN} HUC8 {HUC8}. Residual named. "
        "Sen and Mann-Kendall on the fight, not a merged map. "
        "No granules. Not NASA-endorsed. Counsel unsigned."
    )
    words = len(body.split())
    return {
        "basin": BASIN,
        "huc8": HUC8,
        "sen": sen_slope(r),
        "mk": mann_kendall(r),
        "checksum": fnv1a_32(body),
        "body": body,
        "words": words,
        "offer": OFFER,
        "fetched": False,
    }
