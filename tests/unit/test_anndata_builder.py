from __future__ import annotations

import unittest

import numpy as np
import scipy.sparse as sp

from imz2anndata.anndata_builder import build_anndata
from imz2anndata.models import FeatureTable, PixelCoordinate, SpectrumRecord, SpectrumSignal


class AnnDataBuilderTests(unittest.TestCase):
    def test_build_anndata_populates_fields(self):
        records = [
            SpectrumRecord(
                spectrum_id=0,
                coordinate=PixelCoordinate(x=10, y=20, z=1),
                signal=SpectrumSignal(np.array([100.0]), np.array([1.0], dtype=np.float32)),
            ),
            SpectrumRecord(
                spectrum_id=1,
                coordinate=PixelCoordinate(x=11, y=21, z=1),
                signal=SpectrumSignal(np.array([200.0]), np.array([2.0], dtype=np.float32)),
            ),
        ]
        table = FeatureTable(
            matrix=sp.csr_matrix(np.array([[1.0, 0.0], [0.0, 2.0]], dtype=np.float32)),
            feature_mz=np.array([100.0, 200.0], dtype=np.float64),
        )

        adata = build_anndata(records, table, dataset_id="demo")

        self.assertEqual(adata.shape, (2, 2))
        self.assertEqual(int(adata.obs.loc["pixel_0", "x"]), 10)
        self.assertEqual(int(adata.obs.loc["pixel_1", "y"]), 21)
        self.assertAlmostEqual(float(adata.var.loc["mz_100.000000", "mz"]), 100.0)
        self.assertEqual(adata.uns["dataset_id"], "demo")
        self.assertEqual(adata.uns["feature_space"], "binned_mz")
        self.assertIn("spatial", adata.obsm)
        self.assertTrue(np.array_equal(adata.obsm["spatial"], np.array([[10, 20], [11, 21]])))


if __name__ == "__main__":
    unittest.main()
