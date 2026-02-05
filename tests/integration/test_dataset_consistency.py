from __future__ import annotations

import hashlib
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np
from pyimzml.ImzMLParser import ImzMLParser
from scipy import sparse

from imz2anndata.config import AlignmentConfig, PipelineConfig
from imz2anndata.pipeline import run_pipeline


def get_project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def find_imzml_files() -> list[Path]:
    data_dir = get_project_root() / "data"
    files = sorted(data_dir.glob("*.imzML"))
    files.extend(sorted((data_dir / "real").glob("*.imzML")))
    return files


def array_hash(values: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(values).tobytes()).hexdigest()


def coordinate_stats(coords: np.ndarray) -> tuple[tuple[int, int], tuple[int, int], int, int, tuple[int, ...], tuple[int, ...]]:
    xs = coords[:, 0]
    ys = coords[:, 1]
    x_unique = np.unique(xs)
    y_unique = np.unique(ys)
    x_steps = tuple(np.unique(np.diff(x_unique)).astype(int).tolist()) if x_unique.size > 1 else ()
    y_steps = tuple(np.unique(np.diff(y_unique)).astype(int).tolist()) if y_unique.size > 1 else ()
    return (
        (int(xs.min()), int(xs.max())),
        (int(ys.min()), int(ys.max())),
        int(x_unique.size),
        int(y_unique.size),
        x_steps,
        y_steps,
    )


def ion_image_reference(parser: ImzMLParser, target_mz: float, tol: float) -> np.ndarray:
    half = tol / 2.0
    ref = np.zeros(len(parser.coordinates), dtype=np.float64)
    for i in range(len(parser.coordinates)):
        mz, intensity = parser.getspectrum(i)
        left = np.searchsorted(mz, target_mz - half, side="left")
        right = np.searchsorted(mz, target_mz + half, side="right")
        ref[i] = float(np.sum(intensity[left:right]))
    return ref


def binned_peak_count(parser: ImzMLParser, mz_bin_width: float) -> np.ndarray:
    counts = np.zeros(len(parser.coordinates), dtype=np.int32)
    for i in range(len(parser.coordinates)):
        mz, _ = parser.getspectrum(i)
        bins = np.rint(np.asarray(mz, dtype=np.float64) / mz_bin_width).astype(np.int64)
        counts[i] = int(np.unique(bins).size)
    return counts


class DatasetConsistencyTests(unittest.TestCase):
    def test_core_consistency_and_reproducibility(self):
        imzml_files = find_imzml_files()
        if not imzml_files:
            self.skipTest("No .imzML file found under data/ or data/real/")

        for input_imzml in imzml_files:
            with self.subTest(dataset=input_imzml.name):
                output1 = Path(f"/tmp/{input_imzml.stem}_consistency_1.h5ad")
                output2 = Path(f"/tmp/{input_imzml.stem}_consistency_2.h5ad")

                config = PipelineConfig(
                    input_imzml=input_imzml,
                    output_h5ad=output1,
                    dataset_id=input_imzml.stem,
                    alignment=AlignmentConfig(mz_bin_width=0.01, min_feature_occurrence=1),
                    enable_peak_picking=False,
                )
                adata1 = run_pipeline(config)
                adata2 = run_pipeline(replace(config, output_h5ad=output2))

                parser = ImzMLParser(str(input_imzml))
                try:
                    coords = np.asarray([[int(c[0]), int(c[1])] for c in parser.coordinates], dtype=np.int32)

                    # Pixel count and expected matrix container
                    self.assertEqual(len(parser.coordinates), adata1.n_obs)
                    self.assertTrue(sparse.issparse(adata1.X))

                    # Spatial coordinates must be present and consistent
                    self.assertIn("spatial", adata1.obsm)
                    spatial = np.asarray(adata1.obsm["spatial"], dtype=np.int32)
                    self.assertTrue(np.array_equal(spatial, coords))
                    self.assertEqual(coordinate_stats(spatial), coordinate_stats(coords))

                    # TIC reconstruction should match raw TIC when peak picking is disabled.
                    raw_tic = np.array([float(np.sum(parser.getspectrum(i)[1])) for i in range(len(parser.coordinates))])
                    reconstructed_tic = np.asarray(adata1.X.sum(axis=1)).ravel().astype(np.float64)
                    rel_err = np.abs(reconstructed_tic - raw_tic) / np.maximum(raw_tic, 1e-12)
                    self.assertLess(float(np.quantile(rel_err, 0.99)), 1e-4)

                    # With min_feature_occurrence=1, nnz per pixel should equal unique m/z bins.
                    expected_binned_count = binned_peak_count(parser, mz_bin_width=config.alignment.mz_bin_width)
                    converted_peak_count = np.diff(adata1.X.indptr)
                    self.assertTrue(np.array_equal(converted_peak_count, expected_binned_count))

                    # AnnData structural fields
                    for required_obs in ["x", "y", "tic_raw", "tic_processed", "raw_peak_count", "processed_peak_count"]:
                        self.assertIn(required_obs, adata1.obs.columns)
                    self.assertIn("mz", adata1.var.columns)

                    # Reproducibility: shape + TIC hash + coordinate hash must match.
                    tic_hash_1 = array_hash(np.asarray(adata1.X.sum(axis=1)).ravel())
                    tic_hash_2 = array_hash(np.asarray(adata2.X.sum(axis=1)).ravel())
                    coord_hash_1 = array_hash(np.asarray(adata1.obsm["spatial"]))
                    coord_hash_2 = array_hash(np.asarray(adata2.obsm["spatial"]))
                    self.assertEqual(adata1.shape, adata2.shape)
                    self.assertEqual(tic_hash_1, tic_hash_2)
                    self.assertEqual(coord_hash_1, coord_hash_2)
                finally:
                    handle = getattr(parser, "m", None)
                    if handle is not None and hasattr(handle, "close"):
                        handle.close()

    def test_downstream_ion_images_against_reference(self):
        imzml_files = find_imzml_files()
        if not imzml_files:
            self.skipTest("No .imzML file found under data/ or data/real/")

        input_imzml = imzml_files[0]
        output_h5ad = Path(f"/tmp/{input_imzml.stem}_downstream.h5ad")
        adata = run_pipeline(
            PipelineConfig(
                input_imzml=input_imzml,
                output_h5ad=output_h5ad,
                dataset_id=f"{input_imzml.stem}_downstream",
                alignment=AlignmentConfig(mz_bin_width=0.01, min_feature_occurrence=1),
                enable_peak_picking=False,
            )
        )

        parser = ImzMLParser(str(input_imzml))
        try:
            feature_sums = np.asarray(adata.X.sum(axis=0)).ravel()
            top_cols = np.argsort(feature_sums)[-3:]

            for col in top_cols:
                target_mz = float(adata.var.iloc[int(col)]["mz"])
                converted = np.asarray(adata.X[:, int(col)].toarray()).ravel().astype(np.float64)
                reference = ion_image_reference(parser, target_mz=target_mz, tol=0.01)

                if np.std(reference) == 0 and np.std(converted) == 0:
                    corr = 1.0
                else:
                    corr = float(np.corrcoef(converted, reference)[0, 1])
                self.assertGreater(corr, 0.95)
        finally:
            handle = getattr(parser, "m", None)
            if handle is not None and hasattr(handle, "close"):
                handle.close()


if __name__ == "__main__":
    unittest.main()
