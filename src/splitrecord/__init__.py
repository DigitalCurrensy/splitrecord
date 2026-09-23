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

"""Theil-Sen slope and Mann-Kendall test of z(A) minus z(B)."""

from .report import compile_report, fnv1a_32
from .score import (
    mann_kendall,
    mann_kendall_p,
    mann_kendall_p_ordinary,
    mann_kendall_variance,
    residual,
    seasonal_p,
    seasonal_s,
    seasonal_sen_slope,
    seasonal_variance,
    sen_slope,
    trend_free_prewhiten,
    zscores,
)

__all__ = [
    "compile_report",
    "fnv1a_32",
    "mann_kendall",
    "mann_kendall_p",
    "mann_kendall_variance",
    "residual",
    "seasonal_p",
    "seasonal_s",
    "seasonal_sen_slope",
    "seasonal_variance",
    "mann_kendall_p_ordinary",
    "sen_slope",
    "trend_free_prewhiten",
    "zscores",
]
