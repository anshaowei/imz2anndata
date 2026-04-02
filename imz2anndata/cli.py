from __future__ import annotations

import argparse
from pathlib import Path

from imz2anndata.config import AlignmentConfig, PeakPickingConfig, PipelineConfig, SpatialFilterConfig
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
        "--peak-picking-snr",
        type=float,
        default=0.0,
        help="Signal-to-noise threshold for OpenMS PeakPickerHiRes",
    )
    parser.add_argument(
        "--signal-to-noise-window",
        type=float,
        default=200.0,
        help="Sliding window length used for signal-to-noise estimation during peak picking",
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
    parser.add_argument(
        "--enable-spatial-filtering",
        action="store_true",
        help="Filter features by prevalence, total intensity, and Moran's I",
    )
    parser.add_argument(
        "--assess-spatial-filtering-only",
        action="store_true",
        help="Compute spatial-filter statistics without subsetting the feature matrix",
    )
    parser.add_argument(
        "--min-feature-pixels",
        type=int,
        default=3,
        help="Minimum number of pixels in which a feature must be detected",
    )
    parser.add_argument(
        "--min-feature-total-intensity",
        type=float,
        default=0.0,
        help="Minimum summed intensity required for a feature to be retained",
    )
    parser.add_argument(
        "--min-morans-i",
        type=float,
        default=0.1,
        help="Minimum Moran's I score required for a feature to be retained",
    )
    parser.add_argument(
        "--skip-morans-i",
        action="store_true",
        help="Skip Moran's I computation and retain features using prevalence/intensity filters only",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    config = PipelineConfig(
        input_imzml=args.input_imzml,
        output_h5ad=args.output_h5ad,
        dataset_id=args.dataset_id,
        peak_picking=PeakPickingConfig(
            snr=args.peak_picking_snr,
            signal_to_noise_window=args.signal_to_noise_window,
            min_intensity=args.min_intensity,
        ),
        alignment=AlignmentConfig(
            mz_bin_width=args.mz_bin_width,
            min_feature_occurrence=args.min_feature_occurrence,
        ),
        spatial_filter=SpatialFilterConfig(
            enabled=args.enable_spatial_filtering,
            assess_only=args.assess_spatial_filtering_only,
            min_pixels=args.min_feature_pixels,
            min_total_intensity=args.min_feature_total_intensity,
            min_morans_i=args.min_morans_i,
            compute_morans_i=not args.skip_morans_i,
        ),
        enable_peak_picking=not args.disable_peak_picking,
    )

    run_pipeline(config)


if __name__ == "__main__":
    main()
