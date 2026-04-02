from __future__ import annotations

import unittest

import numpy as np

from imz2anndata.aligner import align_records
from imz2anndata.config import AlignmentConfig
from imz2anndata.models import PixelCoordinate, SpectrumRecord, SpectrumSignal


def make_record(spectrum_id: int, x: int, y: int, mz: list[float], intensity: list[float]) -> SpectrumRecord:
    return SpectrumRecord(
        spectrum_id=spectrum_id,
        coordinate=PixelCoordinate(x=x, y=y),
        signal=SpectrumSignal(
            mz=np.asarray(mz, dtype=np.float64),
            intensity=np.asarray(intensity, dtype=np.float32),
        ),
    )


class AlignerTests(unittest.TestCase):
    def test_align_records_builds_sparse_matrix(self):
        records = [
            make_record(0, 1, 1, mz=[100.001, 100.004, 200.0], intensity=[10.0, 5.0, 3.0]),
            make_record(1, 1, 2, mz=[100.004, 300.0], intensity=[7.0, 8.0]),
        ]

        table = align_records(records, AlignmentConfig(mz_bin_width=0.01, min_feature_occurrence=1))

        self.assertEqual(table.matrix.shape, (2, 3))
        self.assertTrue(np.allclose(table.feature_mz, np.array([100.0, 200.0, 300.0])))

        dense = table.matrix.toarray()
        self.assertAlmostEqual(float(dense[0, 0]), 15.0)
        self.assertAlmostEqual(float(dense[0, 1]), 3.0)
        self.assertAlmostEqual(float(dense[1, 2]), 8.0)

    def test_align_records_respects_min_feature_occurrence(self):
        records = [
            make_record(0, 1, 1, mz=[100.0, 200.0], intensity=[1.0, 2.0]),
            make_record(1, 1, 2, mz=[100.0, 300.0], intensity=[3.0, 4.0]),
        ]

        table = align_records(records, AlignmentConfig(mz_bin_width=0.01, min_feature_occurrence=2))

        self.assertEqual(table.matrix.shape, (2, 1))
        self.assertTrue(np.allclose(table.feature_mz, np.array([100.0])))
        self.assertTrue(np.allclose(table.matrix.toarray(), np.array([[1.0], [3.0]], dtype=np.float32)))


if __name__ == "__main__":
    unittest.main()
