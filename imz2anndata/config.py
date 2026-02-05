from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(slots=True)
class PeakPickingConfig:
    snr: float = 3.0
    signal_to_noise_window: float = 1.0
    min_intensity: float = 0.0


@dataclass(slots=True)
class AlignmentConfig:
    mz_bin_width: float = 0.01
    min_feature_occurrence: int = 3


@dataclass(slots=True)
class PipelineConfig:
    input_imzml: Path
    output_h5ad: Path
    peak_picking: PeakPickingConfig = field(default_factory=PeakPickingConfig)
    alignment: AlignmentConfig = field(default_factory=AlignmentConfig)
    dataset_id: str = "imaging_dataset"
    enable_peak_picking: bool = True
