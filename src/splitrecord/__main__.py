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

"""Score two numeric columns. CSV, USGS RDB, WaterML, TimeseriesML, or a USGS JSON file."""

from __future__ import annotations

import json
import math
import sys
import xml.etree.ElementTree as ET

from .record import finish
from .score import (
    mann_kendall_variance,
    mann_kendall_z,
    residual,
    seasonal_s,
    seasonal_sen_slope,
    seasonal_zp,
    sen_exact_limits,
    gilbert_limits,
    sen_limits,
    tie_counts,
    sen_slope,
    trend_free_prewhiten,
)


def read_column(path: str) -> list[float]:
    """Read one numeric column from a file already on disk.

    CSV is one number per row. `.rdb` is the USGS tab file, column `dv_va`.
    `.xml` is WaterML 1.1 (`value` plus `dateTime`), one WaterML 2.0
    measurement time-value pair series (`MeasurementTVP`), or one
    TimeseriesML 1.0 measurement series (OGC 15-042r3, namespace
    `http://www.opengis.net/tsml/1.0`). A categorical series, a
    domain-range series, a nil, a second series, and TimeseriesML 1.2
    or 1.3 are refused.
    `.json` is either the legacy WaterServices tree
    (`value.timeSeries`) or an OGC FeatureCollection whose properties use
    `value`. This reader does not call USGS. A bad row raises ValueError.
    """
    lower = path.lower()
    if lower.endswith(".rdb"):
        return _read_rdb(path)
    if lower.endswith(".xml"):
        return _read_waterml(path)
    if lower.endswith(".json"):
        return _read_usgs_json(path)
    return _read_plain(path)


def _finite(path: str, raw: object) -> float:
    if isinstance(raw, bool) or raw is None:
        raise ValueError(f"{path}: malformed row")
    if isinstance(raw, (int, float)):
        value = float(raw)
    else:
        text = str(raw).strip()
        if text in {"", "Ice", "Ssn", "Eqp", "Rat", "Dis", "Mnt"}:
            raise ValueError(f"{path}: malformed row")
        try:
            value = float(text)
        except ValueError:
            raise ValueError(f"{path}: malformed row") from None
    if not math.isfinite(value):
        raise ValueError(f"{path}: malformed row")
    return value


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _read_plain(path: str) -> list[float]:
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


def _read_rdb(path: str) -> list[float]:
    names: list[str] | None = None
    values: list[float] = []
    with open(path, encoding="utf-8") as handle:
        for lineno, line in enumerate(handle, start=1):
            text = line.rstrip("\n")
            if text == "" or text.startswith("#"):
                continue
            parts = text.split("\t")
            if names is None:
                names = parts
                continue
            if all(part.endswith(("s", "n", "d")) and len(part) > 1 and part[:-1].isdigit() for part in parts):
                continue
            if "dv_va" not in names:
                raise ValueError(f"{path}:{lineno}: malformed row")
            raw = parts[names.index("dv_va")]
            if raw in {"", "Ice", "Ssn"}:
                raise ValueError(f"{path}:{lineno}: malformed row")
            try:
                value = float(raw)
            except ValueError:
                raise ValueError(f"{path}:{lineno}: malformed row") from None
            if not math.isfinite(value):
                raise ValueError(f"{path}:{lineno}: malformed row")
            values.append(value)
    return values


def _nil(child: ET.Element) -> bool:
    for key, raw in child.attrib.items():
        if key.rsplit("}", 1)[-1] == "nil" and str(raw).strip().lower() == "true":
            return True
    return False


def _tag_ns(tag: str) -> str:
    if tag.startswith("{"):
        return tag[1:].split("}", 1)[0]
    return ""


# TimeseriesML 1.2 is OGC 15-042r5. TimeseriesML 1.3 is OGC 15-042r6.
# Neither namespace is 1.0. The TVP element is still named MeasurementTVP.
_LATER_TSML = (
    "http://www.opengis.net/timeseriesml/1.2",
    "http://www.opengis.net/timeseriesml/1.3",
)
_TSML_10 = "http://www.opengis.net/tsml/1.0"
_SERIES = {"MeasurementTimeseries", "TimeseriesTVP", "Timeseries"}


def _series_count(root: ET.Element) -> int:
    """Count result series, not a series nested inside another."""
    parent: dict[ET.Element, ET.Element] = {}
    for node in root.iter():
        for child in list(node):
            parent[child] = node
    count = 0
    for node in root.iter():
        if _local(node.tag) not in _SERIES:
            continue
        owner = parent.get(node)
        nested = False
        while owner is not None:
            if _local(owner.tag) in _SERIES:
                nested = True
                break
            owner = parent.get(owner)
        if not nested:
            count += 1
    return count


def _type_tokens(root: ET.Element) -> list[str]:
    """Observation-type tokens. Only an element named `type` counts."""
    tokens: list[str] = []
    for node in root.iter():
        if _local(node.tag) != "type":
            continue
        if node.text and node.text.strip():
            tokens.append(node.text.strip())
        for key, raw in node.attrib.items():
            if _local(key) in {"href", "title"} and str(raw).strip():
                tokens.append(str(raw).strip())
    return tokens


def _read_waterml(path: str) -> list[float]:
    """Read one numeric column from WaterML or TimeseriesML already on disk.

    WaterML 2.0 Part 1 (OGC 10-126r4) and TimeseriesML 1.0 (OGC 15-042r3)
    both encode a measurement time-value pair as `MeasurementTVP`. Annex C
    of 15-042r3 maps that element onto the WaterML 2.0 element. The direct
    child named `value` is the column. A nested `value` (the uncertainty
    quantity) is not the column. `uom` is not converted. `time` is not the
    slope axis. Interpolation is not applied.

    The 1.0 result element is `TimeseriesTVP` in the schema. The collection
    example in that schema directory uses `Timeseries` for the same points.
    Both names are one series. Two of them are refused.

    Refused: `CategoricalTVP`, a domain-range document (`TimeseriesDomainRange`,
    `domainSet`, `rangeSet`), an observation type whose token says Categorical
    or DomainRange, `xsi:nil`, and a second series. TimeseriesML 1.2
    (`http://www.opengis.net/timeseriesml/1.2`) and 1.3
    (`http://www.opengis.net/timeseriesml/1.3`) are refused. 1.2 changes
    domain-range metadata and time periods. 1.3 adds `numberTimeSteps`.
    This function does not read those documents.

    WaterML 1.1 is the older CUAHSI document, not an OGC type. A `value`
    element that carries `dateTime` is that column. It is not read when a
    WaterML 2.0 or TimeseriesML result is present.
    """
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as exc:
        raise ValueError(f"{path}: malformed row") from exc
    namespaces = {_tag_ns(el.tag) for el in root.iter()}
    if any(ns in _LATER_TSML for ns in namespaces):
        raise ValueError(f"{path}: not timeseriesml 1.0")
    names = [_local(el.tag) for el in root.iter()]
    tokens = _type_tokens(root)
    if (
        "TimeseriesDomainRange" in names
        or "domainSet" in names
        or "rangeSet" in names
        or any("DomainRange" in token for token in tokens)
    ):
        raise ValueError(f"{path}: domain-range")
    if (
        "CategoricalTimeseries" in names
        or "CategoricalTVP" in names
        or any("Categorical" in token for token in tokens)
    ):
        raise ValueError(f"{path}: categorical timeseries")
    if _series_count(root) > 1:
        raise ValueError(f"{path}: more than one series")
    points = [el for el in root.iter() if _local(el.tag) == "MeasurementTVP"]
    if points:
        values: list[float] = []
        for point in points:
            kids = [child for child in list(point) if _local(child.tag) == "value"]
            if len(kids) != 1:
                raise ValueError(f"{path}: malformed row")
            child = kids[0]
            if _nil(child):
                raise ValueError(f"{path}: nil value")
            values.append(_finite(path, child.text))
        return values
    if _TSML_10 in namespaces:
        raise ValueError(f"{path}: malformed row")
    values = []
    for el in root.iter():
        if _local(el.tag) != "value":
            continue
        if not any(_local(key) == "dateTime" for key in el.attrib):
            continue
        values.append(_finite(path, el.text))
    if not values:
        raise ValueError(f"{path}: malformed row")
    return values


def _read_usgs_json(path: str) -> list[float]:
    try:
        with open(path, encoding="utf-8") as handle:
            payload = json.load(handle)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{path}: malformed row") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"{path}: malformed row")
    if payload.get("type") == "FeatureCollection":
        features = payload.get("features")
        if not isinstance(features, list) or not features:
            raise ValueError(f"{path}: malformed row")
        values: list[float] = []
        for feature in features:
            if not isinstance(feature, dict):
                raise ValueError(f"{path}: malformed row")
            props = feature.get("properties")
            if not isinstance(props, dict) or "value" not in props:
                raise ValueError(f"{path}: malformed row")
            values.append(_finite(path, props["value"]))
        return values
    series = payload.get("value")
    if isinstance(series, dict) and isinstance(series.get("timeSeries"), list):
        values = []
        for item in series["timeSeries"]:
            if not isinstance(item, dict):
                continue
            for block in item.get("values") or []:
                if not isinstance(block, dict):
                    continue
                for row in block.get("value") or []:
                    if isinstance(row, dict) and "dateTime" in row:
                        values.append(_finite(path, row.get("value")))
        if values:
            return values
    raise ValueError(f"{path}: malformed row")


def _parse(args: list[str]) -> tuple[str, str, int | None, bool, bool, bool]:
    files: list[str] = []
    seasons: int | None = None
    covariance = False
    prewhiten = False
    as_json = False
    i = 0
    while i < len(args):
        token = args[i]
        if token == "--json":
            as_json = True
            i += 1
            continue
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
    return files[0], files[1], seasons, covariance, prewhiten, as_json


def _sen_fields(values: list[float], ordinary: float, corrected: float | None = None) -> str:
    if len(values) < 8 and not tie_counts(values):
        lo, hi = sen_exact_limits(values)
        method = "exact"
        hlo, hhi = ("short", "short") if corrected is not None else (None, None)
    else:
        lo, hi = sen_limits(values, ordinary)
        method = "normal" if lo not in {"short", "dependent", "wide"} else lo
        if corrected is None or len(values) < 8:
            hlo, hhi = (None, None) if corrected is None else ("short", "short")
        else:
            hlo, hhi = sen_limits(values, corrected)
    body = f"sen95={method} sen95_lo={lo} sen95_hi={hi}"
    if hlo is not None:
        body = f"{body} hamed95_lo={hlo} hamed95_hi={hhi}"
    if method == "normal":
        glo, ghi = gilbert_limits(values, ordinary)
        body = f"{body} gilbert95_lo={glo} gilbert95_hi={ghi}"
    return body


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    try:
        left_path, right_path, seasons, covariance, prewhiten, as_json = _parse(args)
    except ValueError as exc:
        if str(exc) == "usage":
            print(
                "usage: python -m splitrecord LEFT RIGHT [--seasons N] [--covariance] [--prewhiten] [--json]",
                file=sys.stderr,
            )
            return 2
        print(exc, file=sys.stderr)
        return 1
    try:
        series = residual(read_column(left_path), read_column(right_path))
        whitened: list[float] | None = None
        r1 = 0.0
        removed = 0.0
        if prewhiten:
            whitened, r1, removed = trend_free_prewhiten(series)
            slope = sen_slope(whitened)
            s, tau, var, _factor, z, p = mann_kendall_z(whitened, hamed=False)
            tested = whitened
        elif seasons is None:
            slope = sen_slope(series)
            s, tau, var, factor, z, p = mann_kendall_z(series, hamed=True)
            ordinary = mann_kendall_variance(len(series), tie_counts(series))
            tested = series
        else:
            slope = seasonal_sen_slope(series, seasons)
            s = seasonal_s(series, seasons)
            z, p, var = seasonal_zp(series, seasons, covariance=covariance)
            tested = series
    except (OSError, UnicodeError, ValueError) as exc:
        print(exc, file=sys.stderr)
        return 1
    if len(tested) < 8:
        z_text = "short"
        p_text = "short"
    elif z is None or p is None:
        z_text = "dependent"
        p_text = "dependent"
    else:
        z_text = f"{z:.10g}"
        p_text = f"{p:.6g}"
    var_text = f"{var:.10g}"
    if prewhiten:
        assert whitened is not None
        tau_text = "tied" if tau is None else f"{tau:.10g}"
        line = (
            f"rows={len(series)} residual=z(A)-z(B) series=trend-free-prewhiten "
            f"whitened_rows={len(whitened)} removed_sen={removed:.10g} r1={r1:.6g} "
            f"theil_sen_z_per_row={slope:.10g} {_sen_fields(whitened, var)} "
            f"mann_kendall_S={s} tau={tau_text} "
            f"var={var_text} z={z_text} variance=ordinary p={p_text}"
        )
    elif seasons is None:
        tau_text = "tied" if tau is None else f"{tau:.10g}"
        line = (
            f"rows={len(series)} residual=z(A)-z(B) "
            f"theil_sen_z_per_row={slope:.10g} {_sen_fields(series, ordinary, var)} "
            f"mann_kendall_S={s} tau={tau_text} "
            f"var={var_text} n_over_nstar={factor:.10g} z={z_text} "
            f"variance=hamed-rao p={p_text}"
        )
    else:
        variance = "hirsch-slack" if covariance else "seasonal"
        line = (
            f"rows={len(series)} residual=z(A)-z(B) seasons={seasons} "
            f"theil_sen_z_per_year={slope:.10g} seasonal_S={s} "
            f"var={var_text} z={z_text} variance={variance} p={p_text}"
        )
    word = "scored" if z_text not in {"short", "dependent"} else z_text
    return finish("splitrecord", "z(A)-z(B), then Theil-Sen and Mann-Kendall", [line], as_json, [word])


if __name__ == "__main__":
    sys.exit(main())
