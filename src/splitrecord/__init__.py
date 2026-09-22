"""SPLITRECORD — when two official records disagree, that disagreement is the product."""

from .report import compile_report, fnv1a_32
from .score import mann_kendall, residual, sen_slope

__all__ = ["compile_report", "fnv1a_32", "mann_kendall", "residual", "sen_slope"]
