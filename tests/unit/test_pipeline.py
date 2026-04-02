from __future__ import annotations

import unittest
import warnings
from pathlib import Path

import anndata as ad
import numpy as np

from imz2anndata.config import AlignmentConfig, PipelineConfig, SpatialFilterConfig
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

        with mock.patch("imz2anndata.pipeline.load_imzml_records", return_value=(None, records)):
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

    def test_run_pipeline_skips_peak_picking_for_centroid_mode(self):
        from unittest import mock

        records = [
            SpectrumRecord(
                spectrum_id=0,
                coordinate=PixelCoordinate(x=1, y=1),
                signal=SpectrumSignal(np.array([100.0]), np.array([10.0], dtype=np.float32)),
            )
        ]

        with mock.patch("imz2anndata.pipeline.load_imzml_records", return_value=("centroid", records)):
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    output_h5ad = Path("/tmp/imz2anndata_pipeline_centroid_test.h5ad")
                    cfg = PipelineConfig(
                        input_imzml=Path("dummy.imzML"),
                        output_h5ad=output_h5ad,
                        dataset_id="test_centroid",
                    )
                    adata = run_pipeline(cfg)

        self.assertFalse(adata.uns["peak_picking_enabled"])
        self.assertTrue(adata.uns["peak_picking_requested"])
        self.assertEqual(adata.uns["peak_picking_strategy"], "auto_skipped_for_centroid")
        self.assertEqual(adata.uns["spectrum_mode"], "centroid")
        self.assertEqual(len(caught), 1)
        self.assertIn("skipping peak picking", str(caught[0].message).lower())

    def test_run_pipeline_applies_peak_picking_for_profile_mode(self):
        from unittest import mock

        records = [
            SpectrumRecord(
                spectrum_id=0,
                coordinate=PixelCoordinate(x=1, y=1),
                signal=SpectrumSignal(np.array([100.0]), np.array([10.0], dtype=np.float32)),
            )
        ]

        with mock.patch("imz2anndata.pipeline.load_imzml_records", return_value=("profile", records)):
            with mock.patch("imz2anndata.pipeline.OpenMSPeakExtractor", DummyExtractor):
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    output_h5ad = Path("/tmp/imz2anndata_pipeline_profile_test.h5ad")
                    cfg = PipelineConfig(
                        input_imzml=Path("dummy.imzML"),
                        output_h5ad=output_h5ad,
                        dataset_id="test_profile",
                    )
                    adata = run_pipeline(cfg)

        self.assertTrue(adata.uns["peak_picking_enabled"])
        self.assertTrue(adata.uns["peak_picking_requested"])
        self.assertEqual(adata.uns["peak_picking_strategy"], "auto_applied_for_profile")
        self.assertEqual(adata.uns["spectrum_mode"], "profile")
        self.assertEqual(len(caught), 1)
        self.assertIn("applying openms peak picking", str(caught[0].message).lower())

    def test_run_pipeline_skips_peak_picking_for_unknown_mode(self):
        from unittest import mock

        records = [
            SpectrumRecord(
                spectrum_id=0,
                coordinate=PixelCoordinate(x=1, y=1),
                signal=SpectrumSignal(np.array([100.0]), np.array([10.0], dtype=np.float32)),
            )
        ]

        with mock.patch("imz2anndata.pipeline.load_imzml_records", return_value=(None, records)):
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    output_h5ad = Path("/tmp/imz2anndata_pipeline_unknown_test.h5ad")
                    cfg = PipelineConfig(
                        input_imzml=Path("dummy.imzML"),
                        output_h5ad=output_h5ad,
                        dataset_id="test_unknown",
                    )
                    adata = run_pipeline(cfg)

        self.assertFalse(adata.uns["peak_picking_enabled"])
        self.assertTrue(adata.uns["peak_picking_requested"])
        self.assertEqual(adata.uns["peak_picking_strategy"], "auto_skipped_for_unknown_mode")
        self.assertEqual(adata.uns["spectrum_mode"], "unknown")
        self.assertEqual(len(caught), 1)
        self.assertIn("skipping peak picking", str(caught[0].message).lower())

    def test_run_pipeline_applies_spatial_feature_filtering(self):
        from unittest import mock

        records = [
            SpectrumRecord(
                spectrum_id=0,
                coordinate=PixelCoordinate(x=0, y=0),
                signal=SpectrumSignal(np.array([100.0, 300.0]), np.array([9.0, 1.0], dtype=np.float32)),
            ),
            SpectrumRecord(
                spectrum_id=1,
                coordinate=PixelCoordinate(x=1, y=0),
                signal=SpectrumSignal(np.array([100.0]), np.array([8.0], dtype=np.float32)),
            ),
            SpectrumRecord(
                spectrum_id=2,
                coordinate=PixelCoordinate(x=2, y=0),
                signal=SpectrumSignal(np.array([200.0]), np.array([9.0], dtype=np.float32)),
            ),
            SpectrumRecord(
                spectrum_id=3,
                coordinate=PixelCoordinate(x=3, y=0),
                signal=SpectrumSignal(np.array([200.0, 300.0]), np.array([8.0, 1.0], dtype=np.float32)),
            ),
        ]

        with mock.patch("imz2anndata.pipeline.load_imzml_records", return_value=(None, records)):
            with mock.patch("imz2anndata.pipeline.OpenMSPeakExtractor", DummyExtractor):
                output_h5ad = Path("/tmp/imz2anndata_pipeline_filtered_test.h5ad")
                cfg = PipelineConfig(
                    input_imzml=Path("dummy.imzML"),
                    output_h5ad=output_h5ad,
                    dataset_id="test_filtered",
                    alignment=AlignmentConfig(mz_bin_width=0.01, min_feature_occurrence=1),
                    spatial_filter=SpatialFilterConfig(
                        enabled=True,
                        min_pixels=2,
                        min_total_intensity=5.0,
                        min_morans_i=0.15,
                        connectivity=4,
                    ),
                )
                adata = run_pipeline(cfg)

        self.assertEqual(adata.n_vars, 2)
        self.assertTrue(np.allclose(adata.var["mz"].to_numpy(), np.array([100.0, 200.0])))
        self.assertIn("morans_i", adata.var.columns)
        self.assertTrue(adata.uns["spatial_filter"]["applied"])


if __name__ == "__main__":
    unittest.main()
