from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(slots=True)
class PeakPickingConfig:
    snr: float = 0.0
    signal_to_noise_window: float = 200.0
    min_intensity: float = 0.0


@dataclass(slots=True)
class AlignmentConfig:
    mz_bin_width: float = 0.01
    min_feature_occurrence: int = 3


@dataclass(slots=True)
class SpatialFilterConfig:
    enabled: bool = False
    assess_only: bool = False
    min_pixels: int = 3
    min_total_intensity: float = 0.0
    min_morans_i: float = 0.1
    connectivity: int = 8
    compute_morans_i: bool = True


@dataclass(slots=True)
class PipelineConfig:
    input_imzml: Path
    output_h5ad: Path
    peak_picking: PeakPickingConfig = field(default_factory=PeakPickingConfig)
    alignment: AlignmentConfig = field(default_factory=AlignmentConfig)
    spatial_filter: SpatialFilterConfig = field(default_factory=SpatialFilterConfig)
    dataset_id: str = "imaging_dataset"
    enable_peak_picking: bool = True
