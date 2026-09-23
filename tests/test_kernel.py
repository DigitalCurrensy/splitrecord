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

from __future__ import annotations

import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from splitrecord.report import compile_report, fnv1a_32  # noqa: E402
from splitrecord.score import (  # noqa: E402
    mann_kendall,
    mann_kendall_p,
    mann_kendall_variance,
    residual,
    sen_slope,
    zscores,
)


class ScoreTests(unittest.TestCase):
    def test_length_mismatch_align_first(self) -> None:
        with self.assertRaises(ValueError) as ctx:
            residual([1.0, 2.0], [1.0])
        self.assertEqual(str(ctx.exception), "align first")

    def test_constant_series_no_spread(self) -> None:
        with self.assertRaises(ValueError) as ctx:
            zscores([5.0, 5.0, 5.0])
        self.assertEqual(str(ctx.exception), "no spread")
        with self.assertRaises(ValueError) as ctx:
            residual([1.0, 1.0, 1.0], [0.0, 1.0, 2.0])
        self.assertEqual(str(ctx.exception), "no spread")

    def test_hand_computed_z_residual(self) -> None:
        # A = [1, 2, 4]
        # mean = 7/3
        # sample variance = ((1-7/3)^2 + (2-7/3)^2 + (4-7/3)^2) / 2 = 7/3
        # sample sd = sqrt(7/3)
        # z(A) = (A - 7/3) / sqrt(7/3)
        # B = [0, 3, 6]
        # mean = 3, sample variance = 9, sample sd = 3
        # z(B) = [-1, 0, 1]
        # expected = z(A) - z(B)
        expected = (
            0.12712843905603033,
            -0.2182178902359925,
            0.09108945117996181,
        )
        got = residual([1.0, 2.0, 4.0], [0.0, 3.0, 6.0])
        self.assertEqual(len(got), len(expected))
        for g, e in zip(got, expected):
            self.assertLess(abs(g - e), 1e-9)
        self.assertEqual(zscores([0.0, 3.0, 6.0]), [-1.0, 0.0, 1.0])

    def test_sen_slope_hand_median(self) -> None:
        # Pairwise slopes of [1, 2, 4] are 1, 1.5, and 2. The median is 1.5.
        slope = sen_slope([1.0, 2.0, 4.0])
        self.assertGreater(slope, 0)
        self.assertEqual(slope, 1.5)

    def test_mann_kendall_s_sign(self) -> None:
        self.assertEqual(mann_kendall([1.0, 2.0, 4.0]), 3)
        self.assertGreater(mann_kendall([1.0, 2.0, 4.0]), 0)
        self.assertLess(mann_kendall([4.0, 2.0, 1.0]), 0)

    def test_p_none_when_short(self) -> None:
        self.assertIsNone(mann_kendall_p([1.0, 2.0, 3.0]))
        self.assertIsNone(mann_kendall_p([float(i) for i in range(7)]))
        self.assertIsNotNone(mann_kendall_p([float(i) for i in range(8)]))

    def test_p_long_series(self) -> None:
        increasing = [float(i) for i in range(10)]
        p_up = mann_kendall_p(increasing)
        self.assertIsNotNone(p_up)
        assert p_up is not None
        self.assertGreaterEqual(p_up, 0.0)
        self.assertLessEqual(p_up, 1.0)
        self.assertLess(p_up, 0.05)
        # Hand check of the normal approximation, not a call back into the helper.
        s = 10 * 9 // 2
        var = 10 * 9 * (2 * 10 + 5) / 18
        z = (s - 1) / math.sqrt(var)
        self.assertLess(abs(p_up - math.erfc(abs(z) / math.sqrt(2))), 1e-12)
        self.assertEqual(mann_kendall_variance(10), var)
        # Alternating series: no monotonic pattern.
        wobble = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
        p_flat = mann_kendall_p(wobble)
        self.assertIsNotNone(p_flat)
        assert p_flat is not None
        self.assertGreaterEqual(p_flat, 0.0)
        self.assertLessEqual(p_flat, 1.0)
        self.assertGreater(p_flat, 0.05)

    def test_empty_list_raises(self) -> None:
        for call in (
            lambda: residual([], []),
            lambda: zscores([]),
            lambda: sen_slope([]),
        ):
            with self.assertRaises(ValueError) as ctx:
                call()
            self.assertEqual(str(ctx.exception), "not enough")


class ReportTests(unittest.TestCase):
    def test_report_stays_under_80_and_checksums(self) -> None:
        rep = compile_report([1.0, 2.0, 3.0], [1.0, 1.5, 2.0])
        self.assertEqual(rep["huc8"], "18030012")
        self.assertLessEqual(rep["words"], 80)
        self.assertEqual(rep["checksum"], fnv1a_32(rep["body"]))
        self.assertFalse(rep["fetched"])
        self.assertNotIn("certificate", rep["body"].lower())
        self.assertIn("z(A) minus z(B)", rep["body"])


class CommandTests(unittest.TestCase):
    def test_examples_print_one_line(self) -> None:
        env = dict(os.environ, PYTHONPATH=str(SRC))
        proc = subprocess.run(
            [
                sys.executable,
                "-m",
                "splitrecord",
                str(ROOT / "examples" / "left.csv"),
                str(ROOT / "examples" / "right.csv"),
            ],
            cwd=ROOT,
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertNotIn("Traceback", proc.stderr)
        lines = proc.stdout.splitlines()
        self.assertEqual(len(lines), 1)
        self.assertRegex(lines[0], r"^n=12 sen=\S+ S=-?\d+ p=(short|\d)")

    def test_malformed_row_exits(self) -> None:
        env = dict(os.environ, PYTHONPATH=str(SRC))
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / "bad.csv"
            good = Path(tmp) / "good.csv"
            bad.write_text("1\nnope\n3\n", encoding="utf-8")
            good.write_text("1\n2\n3\n", encoding="utf-8")
            proc = subprocess.run(
                [sys.executable, "-m", "splitrecord", str(bad), str(good)],
                cwd=ROOT,
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )
        self.assertNotEqual(proc.returncode, 0)
        self.assertNotIn("Traceback", proc.stderr)
        self.assertNotIn("Traceback", proc.stdout)
        err_lines = [line for line in proc.stderr.splitlines() if line]
        self.assertEqual(len(err_lines), 1)
        self.assertIn("malformed row", err_lines[0])


if __name__ == "__main__":
    unittest.main()
