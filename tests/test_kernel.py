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
    hamed_rao_factor,
    rank_autocorr,
    mann_kendall,
    mann_kendall_p,
    mann_kendall_variance,
    lag1,
    mann_kendall_p_ordinary,
    residual,
    seasonal_p,
    seasonal_s,
    seasonal_sen_slope,
    seasonal_variance,
    sen_slope,
    tie_counts,
    trend_free_prewhiten,
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
        # Four points, six slopes: 1, 1, 1, 5/3, 2, 3. Even count, so the
        # median is the mean of the two middle slopes, 1 and 5/3.
        even = sen_slope([1.0, 2.0, 3.0, 6.0])
        self.assertAlmostEqual(even, (1.0 + 5.0 / 3.0) / 2.0)


    def test_mann_kendall_s_sign(self) -> None:
        self.assertEqual(mann_kendall([1.0, 2.0, 4.0]), 3)
        self.assertGreater(mann_kendall([1.0, 2.0, 4.0]), 0)
        self.assertLess(mann_kendall([4.0, 2.0, 1.0]), 0)

    def test_p_none_when_short(self) -> None:
        self.assertIsNone(mann_kendall_p([1.0, 2.0, 3.0]))
        self.assertIsNone(mann_kendall_p([float(i) for i in range(7)]))

    def test_p_long_series(self) -> None:
        self.assertEqual(mann_kendall_variance(10), 10 * 9 * 25 / 18)
        repeated = [1.0, 1.0, 1.0, 2.0, 3.0, 4.0]
        bare = mann_kendall_variance(len(repeated))
        term = 3 * (3 - 1) * (2 * 3 + 5) / 18
        corrected = mann_kendall_variance(len(repeated), tie_counts(repeated))
        self.assertEqual(tie_counts(repeated), [3])
        self.assertEqual(corrected, bare - term)
        self.assertLess(corrected, bare)
        straight = [float(i) for i in range(10)]
        self.assertAlmostEqual(hamed_rao_factor(straight), 1.0, delta=1e-9)
        p_line = mann_kendall_p(straight)
        self.assertIsNotNone(p_line)
        assert p_line is not None
        self.assertLess(p_line, 0.05)
        s = mann_kendall(straight)
        var = mann_kendall_variance(10)
        sign = (s > 0) - (s < 0)
        z = (s - sign) / math.sqrt(var)
        self.assertAlmostEqual(p_line, math.erfc(abs(z) / math.sqrt(2.0)))
        # Three identical shapes. After Sen detrending the ranks still repeat,
        # so at least one lag clears the normal bound and the factor is not 1.
        repeats = [5.0, 5.0, 5.0, 0.0] * 3
        self.assertEqual(len(repeats), 12)
        self.assertGreater(hamed_rao_factor(repeats), 1.0)
        # ranks 1,2,3,4. mean 2.5. denom 5. lag-1 numerator 1.25. rho 0.25.
        self.assertAlmostEqual(rank_autocorr([1.0, 2.0, 3.0, 4.0], 1), 0.25)
        with self.assertRaises(ValueError) as ctx:
            sen_slope([1.0, float("nan")])
        self.assertEqual(str(ctx.exception), "bad number")

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
        self.assertLessEqual(rep["words"], 80)
        self.assertEqual(rep["checksum"], fnv1a_32(rep["body"]))
        self.assertFalse(rep["fetched"])
        self.assertEqual(rep["variance"], "hamed-rao")
        self.assertNotIn("certificate", rep["body"].lower())
        self.assertNotIn("Tulare", rep["body"])
        self.assertIn("z(A) minus z(B)", rep["body"])
        self.assertIn("variance Hamed-Rao", rep["body"])


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
        self.assertRegex(
            lines[0],
            r"^rows=12 residual=z\(A\)-z\(B\) theil_sen_z_per_row=\S+ "
            r"mann_kendall_S=-?\d+ variance=hamed-rao p=\S+$",
        )
        s_token = next(part for part in lines[0].split() if part.startswith("mann_kendall_S="))
        s_value = int(s_token.split("=", 1)[1])
        self.assertNotEqual(s_value, 0)

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



class SeasonalTests(unittest.TestCase):
    def test_same_season_only(self) -> None:
        # Two seasons, three years, both rising: 1,1, 2,2, 3,3.
        series = [1.0, 1.0, 2.0, 2.0, 3.0, 3.0]
        self.assertEqual(seasonal_s(series, 2), 6)
        one = 3 * 2 * 11 / 18
        self.assertAlmostEqual(seasonal_variance(series, 2), 2 * one)
        self.assertAlmostEqual(seasonal_variance(series, 2, covariance=True), 4 * one)
        self.assertEqual(seasonal_sen_slope(series, 2), 1.0)
        self.assertIsNone(seasonal_p(series, 2))

    def test_opposite_seasons_cancel_in_the_covariance(self) -> None:
        series = [1.0, 3.0, 2.0, 2.0, 3.0, 1.0]
        self.assertEqual(seasonal_s(series, 2), 0)
        self.assertAlmostEqual(seasonal_variance(series, 2, covariance=True), 0.0)

    def test_uneven_table_refuses_covariance(self) -> None:
        with self.assertRaises(ValueError) as ctx:
            seasonal_variance([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0], 2, covariance=True)
        self.assertEqual(str(ctx.exception), "uneven")

    def test_hamed_rao_is_not_applied(self) -> None:
        series = [float(i) for i in range(12)]
        plain = mann_kendall_variance(6) * 2
        self.assertEqual(seasonal_variance(series, 2), plain)
        factor = hamed_rao_factor(series)
        if factor != 1.0:
            self.assertNotEqual(seasonal_variance(series, 2), plain * factor)

    def test_bad_season(self) -> None:
        with self.assertRaises(ValueError) as ctx:
            seasonal_s([1.0, 2.0, 3.0, 4.0], 1)
        self.assertEqual(str(ctx.exception), "bad season")




class SeasonalCommandTests(unittest.TestCase):
    def test_seasons_line_names_the_variance(self) -> None:
        env = dict(os.environ, PYTHONPATH=str(SRC))
        proc = subprocess.run(
            [
                sys.executable, "-m", "splitrecord",
                str(ROOT / "examples" / "left.csv"),
                str(ROOT / "examples" / "right.csv"),
                "--seasons", "2",
            ],
            cwd=ROOT, env=env, capture_output=True, text=True, check=False,
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        line = proc.stdout.strip()
        self.assertRegex(
            line,
            r"^rows=12 residual=z\(A\)-z\(B\) seasons=2 "
            r"theil_sen_z_per_year=\S+ seasonal_S=-?\d+ variance=seasonal p=\S+$",
        )
        cov = subprocess.run(
            [
                sys.executable, "-m", "splitrecord",
                str(ROOT / "examples" / "left.csv"),
                str(ROOT / "examples" / "right.csv"),
                "--seasons", "2", "--covariance",
            ],
            cwd=ROOT, env=env, capture_output=True, text=True, check=False,
        )
        self.assertEqual(cov.returncode, 0, cov.stderr)
        self.assertIn("variance=hirsch-slack", cov.stdout)



class PrewhitenTests(unittest.TestCase):
    def test_lag1_uses_one_mean(self) -> None:
        # mean 2.5, denom 5, numer 1.25. A two-window Pearson correlation is 1.
        self.assertAlmostEqual(lag1([1.0, 2.0, 3.0, 4.0]), 0.25)

    def test_line_keeps_its_slope(self) -> None:
        series = [float(i) for i in range(6)]
        whitened, r1, _beta = trend_free_prewhiten(series)
        self.assertEqual(r1, 0.0)
        self.assertEqual(len(whitened), 5)
        self.assertEqual(sen_slope(whitened), 1.0)

    def test_blend_matches_the_formula(self) -> None:
        series = [0.0, 1.0, 0.5, 1.5, 1.0, 2.0, 1.5, 2.5, 2.0]
        beta = sen_slope(series)
        detrended = [value - beta * index for index, value in enumerate(series)]
        r1 = lag1(detrended)
        whitened, got_r, _beta = trend_free_prewhiten(series)
        self.assertAlmostEqual(got_r, r1)
        self.assertEqual(len(whitened), len(series) - 1)
        for t in range(1, len(series)):
            expected = (detrended[t] - r1 * detrended[t - 1]) + beta * t
            self.assertAlmostEqual(whitened[t - 1], expected)
        ordinary = mann_kendall_variance(len(whitened), tie_counts(whitened))
        s = mann_kendall(whitened)
        sign = (s > 0) - (s < 0)
        p = math.erfc(abs((s - sign) / math.sqrt(ordinary)) / math.sqrt(2))
        self.assertAlmostEqual(mann_kendall_p_ordinary(whitened), p)

    def test_prewhiten_does_not_stack(self) -> None:
        env = dict(os.environ, PYTHONPATH=str(SRC))
        proc = subprocess.run(
            [sys.executable, "-m", "splitrecord", "a.csv", "b.csv", "--prewhiten", "--seasons", "2"],
            cwd=ROOT, env=env, capture_output=True, text=True, check=False,
        )
        self.assertEqual(proc.returncode, 1)
        self.assertIn("separate", proc.stderr)
        line = subprocess.run(
            [
                sys.executable, "-m", "splitrecord",
                str(ROOT / "examples" / "left.csv"),
                str(ROOT / "examples" / "right.csv"),
                "--prewhiten",
            ],
            cwd=ROOT, env=env, capture_output=True, text=True, check=False,
        )
        self.assertEqual(line.returncode, 0, line.stderr)
        self.assertRegex(
            line.stdout.strip(),
            r"^rows=12 residual=z\(A\)-z\(B\) series=trend-free-prewhiten "
            r"whitened_rows=11 removed_sen=\S+ r1=\S+ theil_sen_z_per_row=\S+ "
            r"mann_kendall_S=-?\d+ variance=ordinary p=\S+$",
        )
        self.assertNotIn("hamed-rao", line.stdout)



class PrewhitenExampleTests(unittest.TestCase):
    def test_pw_files_match_the_readme_line(self) -> None:
        env = dict(os.environ, PYTHONPATH=str(SRC))
        proc = subprocess.run(
            [
                sys.executable, "-m", "splitrecord",
                str(ROOT / "examples" / "pw_left.csv"),
                str(ROOT / "examples" / "pw_right.csv"),
                "--prewhiten",
            ],
            cwd=ROOT, env=env, capture_output=True, text=True, check=False,
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        line = proc.stdout.strip()
        self.assertEqual(
            line,
            "rows=9 residual=z(A)-z(B) series=trend-free-prewhiten "
            "whitened_rows=8 removed_sen=0.04892060565 r1=-0.888889 "
            "theil_sen_z_per_row=0.04892060565 mann_kendall_S=22 "
            "variance=ordinary p=0.00937477",
        )
        series = residual(
            [float(v) for v in (ROOT / "examples" / "pw_left.csv").read_text().split()],
            [float(v) for v in (ROOT / "examples" / "pw_right.csv").read_text().split()],
        )
        slopes = [
            (series[j] - series[i]) / (j - i)
            for i in range(len(series))
            for j in range(i + 1, len(series))
        ]
        self.assertEqual(len(slopes), 36)
        slopes.sort()
        median = 0.5 * (slopes[17] + slopes[18])
        self.assertAlmostEqual(median, sen_slope(series))
        self.assertIn(f"removed_sen={median:.10g}", line)


if __name__ == "__main__":
    unittest.main()
