from __future__ import annotations

import unittest

import numpy as np
import scipy.sparse as sp

from imz2anndata.anndata_builder import build_anndata
from imz2anndata.models import FeatureTable, PixelCoordinate, SpectrumRecord, SpectrumSignal


class BoundaryTests(unittest.TestCase):
    def test_non_regular_coordinates_are_preserved(self):
        records = [
            SpectrumRecord(
                spectrum_id=0,
                coordinate=PixelCoordinate(x=1, y=1),
                signal=SpectrumSignal(np.array([100.0]), np.array([1.0], dtype=np.float32)),
            ),
            SpectrumRecord(
                spectrum_id=1,
                coordinate=PixelCoordinate(x=3, y=1),
                signal=SpectrumSignal(np.array([100.0]), np.array([2.0], dtype=np.float32)),
            ),
            SpectrumRecord(
                spectrum_id=2,
                coordinate=PixelCoordinate(x=10, y=4),
                signal=SpectrumSignal(np.array([200.0]), np.array([3.0], dtype=np.float32)),
            ),
        ]
        table = FeatureTable(
            matrix=sp.csr_matrix(np.array([[1.0, 0.0], [2.0, 0.0], [0.0, 3.0]], dtype=np.float32)),
            feature_mz=np.array([100.0, 200.0], dtype=np.float64),
        )

        adata = build_anndata(records, table, dataset_id="irregular")
        spatial = np.asarray(adata.obsm["spatial"])

        self.assertTrue(np.array_equal(spatial[:, 0], np.array([1, 3, 10])))
        self.assertTrue(np.array_equal(spatial[:, 1], np.array([1, 1, 4])))


if __name__ == "__main__":
    unittest.main()
