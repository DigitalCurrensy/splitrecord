"""Cite + FNV checksum + paid report. Field-locked briefs. No granules."""

from __future__ import annotations

from .score import mann_kendall, residual, sen_slope

OFFER = "paid report $4k–$12k. Not an invoice."
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
