from __future__ import annotations

import unittest
from pathlib import Path

import anndata as ad
import numpy as np

from imz2anndata.config import PipelineConfig
from imz2anndata.models import PixelCoordinate, SpectrumRecord, SpectrumSignal
from imz2anndata.pipeline import run_pipeline


class DummyExtractor:
    def __init__(self, *_args, **_kwargs):
        pass

    def extract(self, signal):
        return signal


class PipelineTests(unittest.TestCase):
    def test_run_pipeline_writes_h5ad(self):
        from unittest import mock

        records = [
            SpectrumRecord(
                spectrum_id=0,
                coordinate=PixelCoordinate(x=1, y=1),
                signal=SpectrumSignal(np.array([100.0]), np.array([10.0], dtype=np.float32)),
            ),
            SpectrumRecord(
                spectrum_id=1,
                coordinate=PixelCoordinate(x=1, y=2),
                signal=SpectrumSignal(
                    np.array([100.0, 200.0]),
                    np.array([4.0, 6.0], dtype=np.float32),
                ),
            ),
        ]

        with mock.patch("imz2anndata.pipeline.iter_imzml_spectra", return_value=iter(records)):
            with mock.patch("imz2anndata.pipeline.OpenMSPeakExtractor", DummyExtractor):
                output_h5ad = Path("/tmp/imz2anndata_pipeline_test.h5ad")
                cfg = PipelineConfig(
                    input_imzml=Path("dummy.imzML"),
                    output_h5ad=output_h5ad,
                    dataset_id="test",
                )
                adata = run_pipeline(cfg)

        self.assertTrue(output_h5ad.exists())
        loaded = ad.read_h5ad(output_h5ad)
        self.assertEqual(loaded.shape, adata.shape)
        self.assertEqual(loaded.uns["dataset_id"], "test")
        self.assertIn("tic_raw", loaded.obs.columns)
        self.assertIn("tic_processed", loaded.obs.columns)


if __name__ == "__main__":
    unittest.main()
