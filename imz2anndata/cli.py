from __future__ import annotations

import argparse
from pathlib import Path

from imz2anndata.config import AlignmentConfig, PeakPickingConfig, PipelineConfig
from imz2anndata.pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert imzML to AnnData")
    parser.add_argument("input_imzml", type=Path, help="Path to input .imzML")
    parser.add_argument("output_h5ad", type=Path, help="Path to output .h5ad")
    parser.add_argument("--dataset-id", default="imaging_dataset", help="Dataset identifier")
    parser.add_argument("--mz-bin-width", type=float, default=0.01, help="m/z bin width for alignment")
    parser.add_argument(
        "--min-feature-occurrence",
        type=int,
        default=3,
        help="Minimum number of spectra a feature must appear in",
    )
    parser.add_argument(
        "--min-intensity",
        type=float,
        default=0.0,
        help="Drop picked peaks below this intensity",
    )
    parser.add_argument(
        "--disable-peak-picking",
        action="store_true",
        help="Skip OpenMS peak picking and keep original spectra",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    config = PipelineConfig(
        input_imzml=args.input_imzml,
        output_h5ad=args.output_h5ad,
        dataset_id=args.dataset_id,
        peak_picking=PeakPickingConfig(min_intensity=args.min_intensity),
        alignment=AlignmentConfig(
            mz_bin_width=args.mz_bin_width,
            min_feature_occurrence=args.min_feature_occurrence,
        ),
        enable_peak_picking=not args.disable_peak_picking,
    )

    run_pipeline(config)


if __name__ == "__main__":
    main()
