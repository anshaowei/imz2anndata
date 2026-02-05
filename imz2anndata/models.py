from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse


@dataclass(slots=True)
class SpectrumSignal:
    mz: np.ndarray
    intensity: np.ndarray


@dataclass(slots=True)
class PixelCoordinate:
    x: int
    y: int
    z: int | None = None


@dataclass(slots=True)
class SpectrumRecord:
    spectrum_id: int
    coordinate: PixelCoordinate
    signal: SpectrumSignal


@dataclass(slots=True)
class FeatureTable:
    matrix: sparse.csr_matrix
    feature_mz: np.ndarray
