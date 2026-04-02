from __future__ import annotations

import warnings

import numpy as np
from pyopenms import MSSpectrum, PeakPickerHiRes

from imz2anndata.config import PeakPickingConfig
from imz2anndata.models import SpectrumSignal


class OpenMSPeakExtractor:
    def __init__(self, config: PeakPickingConfig) -> None:
        self._config = config
        self._picker = PeakPickerHiRes()
        params = self._picker.getParameters()
        params.setValue(b"signal_to_noise", float(config.snr))
        params.setValue(b"SignalToNoise:win_len", float(config.signal_to_noise_window))
        self._picker.setParameters(params)

    def _validate_signal(self, signal: SpectrumSignal) -> str | None:
        if signal.mz.ndim != 1 or signal.intensity.ndim != 1:
            return "expected 1D m/z and intensity arrays"
        if signal.mz.size != signal.intensity.size:
            return "m/z and intensity arrays must have the same length"
        if signal.mz.size < 3:
            return "PeakPickerHiRes requires at least three profile points"
        if not np.isfinite(signal.mz).all() or not np.isfinite(signal.intensity).all():
            return "m/z and intensity arrays must be finite"
        if np.any(signal.intensity < 0):
            return "intensity values must be non-negative"
        if np.any(np.diff(signal.mz) <= 0):
            return "m/z values must be strictly increasing"
        return None

    def extract(self, signal: SpectrumSignal) -> SpectrumSignal:
        if signal.mz.size == 0:
            return signal

        validation_error = self._validate_signal(signal)
        if validation_error is not None:
            warnings.warn(
                f"Skipping OpenMS peak picking for incompatible spectrum: {validation_error}",
                RuntimeWarning,
                stacklevel=2,
            )
            return signal

        input_spectrum = MSSpectrum()
        input_spectrum.set_peaks((signal.mz, signal.intensity.astype(np.float64)))

        picked = MSSpectrum()
        try:
            self._picker.pick(input_spectrum, picked)
            mz_values, intensity_values = picked.get_peaks()
            mz = np.asarray(mz_values, dtype=np.float64)
            intensity = np.asarray(intensity_values, dtype=np.float32)
        except Exception as exc:
            warnings.warn(
                f"OpenMS PeakPickerHiRes failed; returning the original spectrum. {exc}",
                RuntimeWarning,
                stacklevel=2,
            )
            mz = signal.mz
            intensity = signal.intensity

        if self._config.min_intensity > 0:
            keep = intensity >= self._config.min_intensity
            mz = mz[keep]
            intensity = intensity[keep]

        return SpectrumSignal(mz=mz, intensity=intensity)
