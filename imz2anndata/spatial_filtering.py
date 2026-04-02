from __future__ import annotations

from dataclasses import asdict

import anndata as ad
import numpy as np
from scipy import sparse

from imz2anndata.config import SpatialFilterConfig


def _build_adjacency(coords: np.ndarray, connectivity: int) -> sparse.csr_matrix:
    if connectivity not in {4, 8}:
        raise ValueError("connectivity must be 4 or 8")

    n_obs = coords.shape[0]
    if n_obs == 0:
        return sparse.csr_matrix((0, 0), dtype=np.float32)

    rows: list[int] = []
    cols: list[int] = []
    coord_to_index = {(int(x), int(y)): idx for idx, (x, y) in enumerate(coords.tolist())}
    offsets = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    if connectivity == 8:
        offsets.extend([(-1, -1), (-1, 1), (1, -1), (1, 1)])

    for idx, (x, y) in enumerate(coords.tolist()):
        for dx, dy in offsets:
            neighbor = coord_to_index.get((int(x + dx), int(y + dy)))
            if neighbor is not None:
                rows.append(idx)
                cols.append(neighbor)

    data = np.ones(len(rows), dtype=np.float32)
    adjacency = sparse.coo_matrix((data, (rows, cols)), shape=(n_obs, n_obs), dtype=np.float32)
    adjacency = adjacency.maximum(adjacency.T).tocsr()
    adjacency.setdiag(0)
    adjacency.eliminate_zeros()
    return adjacency


def _morans_i_for_matrix(matrix: sparse.csr_matrix, adjacency: sparse.csr_matrix) -> np.ndarray:
    n_obs, n_vars = matrix.shape
    if n_vars == 0:
        return np.array([], dtype=np.float64)

    s0 = float(adjacency.sum())
    if n_obs == 0 or s0 == 0.0:
        return np.zeros(n_vars, dtype=np.float64)

    scores = np.zeros(n_vars, dtype=np.float64)
    for col in range(n_vars):
        values = matrix[:, col].toarray().ravel().astype(np.float64, copy=False)
        centered = values - values.mean()
        denominator = float(centered @ centered)
        if denominator <= 0.0:
            scores[col] = 0.0
            continue
        numerator = float(centered @ adjacency.dot(centered))
        scores[col] = (n_obs / s0) * (numerator / denominator)
    return scores


def filter_spatial_features(adata: ad.AnnData, config: SpatialFilterConfig) -> ad.AnnData:
    should_assess = bool(config.enabled or config.assess_only)

    if not sparse.issparse(adata.X):
        matrix = sparse.csr_matrix(np.asarray(adata.X, dtype=np.float32))
    else:
        matrix = adata.X.tocsr().astype(np.float32, copy=False)

    if not should_assess:
        adata.var["selected"] = np.ones(adata.n_vars, dtype=bool)
        adata.uns["spatial_filter"] = {
            **asdict(config),
            "morans_i_computed": False,
            "n_features_before": int(adata.n_vars),
            "n_features_after": int(adata.n_vars),
            "applied": False,
            "assessed": False,
        }
        return adata

    matrix_csc = matrix.tocsc()
    detection_count = np.diff(matrix_csc.indptr).astype(np.int32, copy=False)
    total_intensity = np.asarray(matrix_csc.sum(axis=0)).ravel().astype(np.float64, copy=False)

    pass_prevalence = detection_count >= int(config.min_pixels)
    pass_intensity = total_intensity >= float(config.min_total_intensity)
    candidate_mask = pass_prevalence & pass_intensity

    morans_i = np.full(adata.n_vars, np.nan, dtype=np.float64)
    pass_morans = np.ones(adata.n_vars, dtype=bool)
    if config.compute_morans_i:
        if "spatial" not in adata.obsm:
            raise KeyError("AnnData must contain adata.obsm['spatial'] for spatial filtering")

        adjacency = _build_adjacency(np.asarray(adata.obsm["spatial"], dtype=np.int32), config.connectivity)
        pass_morans = np.zeros(adata.n_vars, dtype=bool)
        if np.any(candidate_mask):
            candidate_scores = _morans_i_for_matrix(matrix[:, candidate_mask].tocsr(), adjacency)
            morans_i[candidate_mask] = candidate_scores
        pass_morans[candidate_mask] = morans_i[candidate_mask] >= float(config.min_morans_i)

    keep_mask = pass_prevalence & pass_intensity & pass_morans

    adata.var["detection_count"] = detection_count
    adata.var["detection_rate"] = detection_count / max(adata.n_obs, 1)
    adata.var["total_intensity"] = total_intensity
    adata.var["morans_i"] = morans_i
    adata.var["pass_prevalence_filter"] = pass_prevalence
    adata.var["pass_intensity_filter"] = pass_intensity
    adata.var["pass_morans_filter"] = pass_morans
    adata.var["selected"] = keep_mask
    adata.uns["spatial_filter"] = {
        **asdict(config),
        "morans_i_computed": bool(config.compute_morans_i),
        "n_features_before": int(adata.n_vars),
        "n_features_after": int(np.count_nonzero(keep_mask)),
        "applied": False,
        "assessed": True,
    }

    if not config.enabled:
        return adata

    filtered = adata[:, keep_mask].copy()
    filtered.uns["spatial_filter"] = {
        **adata.uns["spatial_filter"],
        "applied": True,
    }
    return filtered
