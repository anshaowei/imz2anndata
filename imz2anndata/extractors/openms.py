from __future__ import annotations

import numpy as np
from pyopenms import MSSpectrum, PeakPickerHiRes

from imz2anndata.config import PeakPickingConfig
from imz2anndata.models import SpectrumSignal


class OpenMSPeakExtractor:
    def __init__(self, config: PeakPickingConfig) -> None:
        self._config = config
        self._picker = PeakPickerHiRes()

    def extract(self, signal: SpectrumSignal) -> SpectrumSignal:
        if signal.mz.size == 0:
            return signal

        input_spectrum = MSSpectrum()
        input_spectrum.set_peaks((signal.mz, signal.intensity.astype(np.float64)))

        picked = MSSpectrum()
        try:
            self._picker.pick(input_spectrum, picked)
            mz_values, intensity_values = picked.get_peaks()
            mz = np.asarray(mz_values, dtype=np.float64)
            intensity = np.asarray(intensity_values, dtype=np.float32)
        except Exception:
            # Fallback path keeps pipeline usable when OpenMS fails on edge-case spectra.
            mz = signal.mz
            intensity = signal.intensity

        if self._config.min_intensity > 0:
            keep = intensity >= self._config.min_intensity
            mz = mz[keep]
            intensity = intensity[keep]

        return SpectrumSignal(mz=mz, intensity=intensity)
