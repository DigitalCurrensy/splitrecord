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

"""One JSON record. Not WaterML. Not a USGS response. Not a stamp."""

from __future__ import annotations

import json

ABSENT = ["stamp", "measured_months", "laz_points"]
# ok and pass are gate words. scored is a finished statistical line.
# path is a walk that reached the goal. None of them is a keep.
PASS = {"ok", "pass", "scored", "path"}


def file_word(words: list[str]) -> str:
    for word in words:
        if word not in PASS:
            return word
    return words[0] if words else "missing"


def record(desk: str, word: str, formula: str, rows: list[dict[str, str]]) -> dict[str, object]:
    return {
        "absent": ABSENT,
        "desk": desk,
        "exit": 0 if word in PASS else 1,
        "formula": formula,
        "keep": False,
        "rows": rows,
        "word": word,
    }


def emit(payload: dict[str, object]) -> int:
    print(json.dumps(payload, separators=(",", ":"), sort_keys=True))
    return int(payload["exit"])


def finish(desk: str, formula: str, lines: list[str], as_json: bool, words: list[str] | None = None) -> int:
    parsed = words if words is not None else [line.split()[0] for line in lines]
    word = file_word(parsed)
    rows = [{"line": line, "word": item} for line, item in zip(lines, parsed)]
    payload = record(desk, word, formula, rows)
    if as_json:
        return emit(payload)
    for line in lines:
        print(line)
    return int(payload["exit"])
