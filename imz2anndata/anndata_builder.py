from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from imz2anndata.models import FeatureTable, SpectrumRecord


def build_anndata(records: list[SpectrumRecord], table: FeatureTable, dataset_id: str) -> ad.AnnData:
    obs = pd.DataFrame(
        {
            "spectrum_id": [r.spectrum_id for r in records],
            "x": [r.coordinate.x for r in records],
            "y": [r.coordinate.y for r in records],
            "z": [float(r.coordinate.z) if r.coordinate.z is not None else np.nan for r in records],
        }
    )
    obs.index = [f"pixel_{i}" for i in range(len(records))]

    var = pd.DataFrame({"mz": table.feature_mz})
    var.index = [f"mz_{mz:.6f}" for mz in table.feature_mz]

    adata = ad.AnnData(X=table.matrix, obs=obs, var=var)
    adata.obsm["spatial"] = obs[["x", "y"]].to_numpy()
    adata.uns["dataset_id"] = dataset_id
    adata.uns["feature_space"] = "binned_mz"
    return adata
