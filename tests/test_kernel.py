from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from splitrecord.report import compile_report, fnv1a_32  # noqa: E402
from splitrecord.score import mann_kendall, residual, sen_slope  # noqa: E402


class ScoreTests(unittest.TestCase):
    def test_residual_and_trend(self) -> None:
        r = residual([1.0, 2.0, 4.0], [0.0, 0.0, 0.0])
        self.assertEqual(r, [1.0, 2.0, 4.0])
        self.assertGreater(sen_slope(r), 0)
        self.assertGreater(mann_kendall(r), 0)

    def test_align_first(self) -> None:
        with self.assertRaises(ValueError):
            residual([1.0], [1.0, 2.0])


class ReportTests(unittest.TestCase):
    def test_report_stays_under_80_and_checksums(self) -> None:
        rep = compile_report([1.0, 2.0, 3.0], [1.0, 1.5, 2.0])
        self.assertEqual(rep["huc8"], "18030012")
        self.assertLessEqual(rep["words"], 80)
        self.assertEqual(rep["checksum"], fnv1a_32(rep["body"]))
        self.assertFalse(rep["fetched"])


if __name__ == "__main__":
    unittest.main()
