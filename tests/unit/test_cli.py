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
                "--peak-picking-snr",
                "2.5",
                "--signal-to-noise-window",
                "150",
                "--min-intensity",
                "10",
                "--assess-spatial-filtering-only",
                "--enable-spatial-filtering",
                "--min-feature-pixels",
                "5",
                "--min-feature-total-intensity",
                "42",
                "--min-morans-i",
                "0.2",
                "--skip-morans-i",
                "--disable-peak-picking",
            ]
        )

        self.assertEqual(args.input_imzml, Path("data/input.imzML"))
        self.assertEqual(args.output_h5ad, Path("data/output.h5ad"))
        self.assertEqual(args.dataset_id, "case_a")
        self.assertEqual(args.mz_bin_width, 0.02)
        self.assertEqual(args.min_feature_occurrence, 3)
        self.assertEqual(args.peak_picking_snr, 2.5)
        self.assertEqual(args.signal_to_noise_window, 150.0)
        self.assertEqual(args.min_intensity, 10)
        self.assertTrue(args.assess_spatial_filtering_only)
        self.assertTrue(args.enable_spatial_filtering)
        self.assertEqual(args.min_feature_pixels, 5)
        self.assertEqual(args.min_feature_total_intensity, 42)
        self.assertEqual(args.min_morans_i, 0.2)
        self.assertTrue(args.skip_morans_i)
        self.assertTrue(args.disable_peak_picking)


if __name__ == "__main__":
    unittest.main()
