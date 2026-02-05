from __future__ import annotations

import unittest
from pathlib import Path

from imz2anndata.cli import build_parser


class CLITests(unittest.TestCase):
    def test_cli_parser_reads_common_options(self):
        parser = build_parser()
        args = parser.parse_args(
            [
                "data/input.imzML",
                "data/output.h5ad",
                "--dataset-id",
                "case_a",
                "--mz-bin-width",
                "0.02",
                "--min-feature-occurrence",
                "3",
                "--min-intensity",
                "10",
                "--disable-peak-picking",
            ]
        )

        self.assertEqual(args.input_imzml, Path("data/input.imzML"))
        self.assertEqual(args.output_h5ad, Path("data/output.h5ad"))
        self.assertEqual(args.dataset_id, "case_a")
        self.assertEqual(args.mz_bin_width, 0.02)
        self.assertEqual(args.min_feature_occurrence, 3)
        self.assertEqual(args.min_intensity, 10)
        self.assertTrue(args.disable_peak_picking)


if __name__ == "__main__":
    unittest.main()
