from __future__ import annotations

import unittest

import anndata as ad
import numpy as np
import scipy.sparse as sp

from imz2anndata.config import SpatialFilterConfig
from imz2anndata.spatial_filtering import filter_spatial_features


class SpatialFilteringTests(unittest.TestCase):
    def test_filter_spatial_features_scores_and_filters_features(self):
        matrix = sp.csr_matrix(
            np.array(
                [
                    [5.0, 0.0, 1.0],
                    [4.0, 0.0, 0.0],
                    [0.0, 4.0, 0.0],
                    [0.0, 5.0, 1.0],
                ],
                dtype=np.float32,
            )
        )
        adata = ad.AnnData(X=matrix)
        adata.obs["x"] = [0, 1, 2, 3]
        adata.obs["y"] = [0, 0, 0, 0]
        adata.obsm["spatial"] = np.array([[0, 0], [1, 0], [2, 0], [3, 0]], dtype=np.int32)
        adata.var["mz"] = [100.0, 200.0, 300.0]

        filtered = filter_spatial_features(
            adata,
            SpatialFilterConfig(
                enabled=True,
                min_pixels=2,
                min_total_intensity=5.0,
                min_morans_i=0.15,
                connectivity=4,
            ),
        )

        self.assertEqual(filtered.n_vars, 2)
        self.assertTrue(np.allclose(filtered.var["mz"].to_numpy(), np.array([100.0, 200.0])))
        self.assertIn("morans_i", filtered.var.columns)
        self.assertTrue(np.all(filtered.var["pass_morans_filter"].to_numpy()))
        self.assertEqual(filtered.uns["spatial_filter"]["n_features_before"], 3)
        self.assertEqual(filtered.uns["spatial_filter"]["n_features_after"], 2)
        self.assertTrue(filtered.uns["spatial_filter"]["applied"])

    def test_filter_spatial_features_can_assess_without_subsetting(self):
        matrix = sp.csr_matrix(np.array([[1.0, 0.0], [0.0, 2.0]], dtype=np.float32))
        adata = ad.AnnData(X=matrix)
        adata.obsm["spatial"] = np.array([[0, 0], [1, 0]], dtype=np.int32)

        assessed = filter_spatial_features(
            adata,
            SpatialFilterConfig(
                enabled=False,
                assess_only=True,
                min_pixels=1,
                min_total_intensity=0.0,
                min_morans_i=-1.0,
            ),
        )

        self.assertEqual(assessed.shape, (2, 2))
        self.assertIn("selected", assessed.var.columns)
        self.assertEqual(assessed.uns["spatial_filter"]["n_features_before"], 2)
        self.assertEqual(assessed.uns["spatial_filter"]["n_features_after"], 2)

    def test_filter_spatial_features_can_skip_morans_i(self):
        matrix = sp.csr_matrix(
            np.array(
                [
                    [5.0, 0.0, 1.0],
                    [4.0, 0.0, 0.0],
                    [0.0, 4.0, 0.0],
                    [0.0, 5.0, 1.0],
                ],
                dtype=np.float32,
            )
        )
        adata = ad.AnnData(X=matrix)
        adata.obsm["spatial"] = np.array([[0, 0], [1, 0], [2, 0], [3, 0]], dtype=np.int32)

        assessed = filter_spatial_features(
            adata,
            SpatialFilterConfig(
                enabled=False,
                assess_only=True,
                min_pixels=2,
                min_total_intensity=5.0,
                min_morans_i=0.99,
                compute_morans_i=False,
            ),
        )

        self.assertTrue(np.isnan(assessed.var["morans_i"].to_numpy()).all())
        self.assertTrue(assessed.var["pass_morans_filter"].to_numpy().all())
        self.assertTrue(
            np.array_equal(
                assessed.var["selected"].to_numpy(),
                np.array([True, True, False]),
            )
        )
        self.assertFalse(assessed.uns["spatial_filter"]["morans_i_computed"])

    def test_filter_spatial_features_skips_assessment_when_disabled(self):
        matrix = sp.csr_matrix(np.array([[1.0, 0.0], [0.0, 2.0]], dtype=np.float32))
        adata = ad.AnnData(X=matrix)
        adata.obsm["spatial"] = np.array([[0, 0], [1, 0]], dtype=np.int32)

        assessed = filter_spatial_features(
            adata,
            SpatialFilterConfig(
                enabled=False,
                assess_only=False,
            ),
        )

        self.assertNotIn("morans_i", assessed.var.columns)
        self.assertTrue(np.array_equal(assessed.var["selected"].to_numpy(), np.array([True, True])))
        self.assertFalse(assessed.uns["spatial_filter"]["assessed"])
        self.assertFalse(assessed.uns["spatial_filter"]["applied"])


if __name__ == "__main__":
    unittest.main()
