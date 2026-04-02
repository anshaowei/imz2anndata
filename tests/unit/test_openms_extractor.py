from __future__ import annotations

import unittest
import warnings

import numpy as np

from imz2anndata.config import PeakPickingConfig
from imz2anndata.extractors.openms import OpenMSPeakExtractor
from imz2anndata.models import SpectrumSignal


class OpenMSPeakExtractorTests(unittest.TestCase):
    def test_config_defaults_match_openms_defaults(self):
        cfg = PeakPickingConfig()

        self.assertEqual(cfg.snr, 0.0)
        self.assertEqual(cfg.signal_to_noise_window, 200.0)

    def test_peak_picker_receives_configured_parameters(self):
        cfg = PeakPickingConfig(snr=7.5, signal_to_noise_window=42.0, min_intensity=10.0)
        extractor = OpenMSPeakExtractor(cfg)

        params = extractor._picker.getParameters()

        self.assertEqual(float(params.getValue(b"signal_to_noise")), 7.5)
        self.assertEqual(float(params.getValue(b"SignalToNoise:win_len")), 42.0)

    def test_unsorted_spectrum_is_reported_and_returned_unchanged(self):
        extractor = OpenMSPeakExtractor(PeakPickingConfig())
        signal = SpectrumSignal(
            mz=np.array([100.2, 100.0, 100.1], dtype=np.float64),
            intensity=np.array([10.0, 0.0, 5.0], dtype=np.float32),
        )

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = extractor.extract(signal)

        self.assertEqual(len(caught), 1)
        self.assertIn("strictly increasing", str(caught[0].message))
        np.testing.assert_array_equal(result.mz, signal.mz)
        np.testing.assert_array_equal(result.intensity, signal.intensity)


if __name__ == "__main__":
    unittest.main()
