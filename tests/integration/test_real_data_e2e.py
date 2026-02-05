from __future__ import annotations

import unittest
from pathlib import Path

from imz2anndata.config import PipelineConfig
from imz2anndata.pipeline import run_pipeline


def get_project_root() -> Path:
    """Return the project root path."""
    return Path(__file__).resolve().parents[2]


def first_imzml_file() -> Path | None:
    project_root = get_project_root()
    data_path = project_root / "data"
    candidates = sorted(data_path.glob("*.imzML")) + sorted((data_path / "real").glob("*.imzML"))
    return candidates[0] if candidates else None


class RealDataE2ETests(unittest.TestCase):
    def test_pipeline_with_real_data(self):
        input_imzml = first_imzml_file()
        if input_imzml is None:
            self.skipTest("No .imzML file found under data/ or data/real/")

        output_h5ad = Path("/tmp/imz2anndata_real_e2e.h5ad")
        adata = run_pipeline(
            PipelineConfig(
                input_imzml=input_imzml,
                output_h5ad=output_h5ad,
                dataset_id="real_validation",
            )
        )

        self.assertGreater(adata.n_obs, 0)
        self.assertGreater(adata.n_vars, 0)
        self.assertTrue(output_h5ad.exists())


if __name__ == "__main__":
    unittest.main()
