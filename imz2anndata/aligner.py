from __future__ import annotations

from collections import Counter

import numpy as np
from scipy import sparse

from imz2anndata.config import AlignmentConfig
from imz2anndata.models import FeatureTable, SpectrumRecord


def _to_bin(mz: np.ndarray, mz_bin_width: float) -> np.ndarray:
    return np.rint(mz / mz_bin_width).astype(np.int64)


def align_records(records: list[SpectrumRecord], config: AlignmentConfig) -> FeatureTable:
    if not records:
        return FeatureTable(
            matrix=sparse.csr_matrix((0, 0), dtype=np.float32),
            feature_mz=np.array([], dtype=np.float64),
        )

    feature_counter: Counter[int] = Counter()
    for record in records:
        bins = _to_bin(record.signal.mz, config.mz_bin_width)
        feature_counter.update(np.unique(bins).tolist())

    kept_bins = [b for b, cnt in feature_counter.items() if cnt >= config.min_feature_occurrence]
    kept_bins.sort()
    if not kept_bins:
        return FeatureTable(
            matrix=sparse.csr_matrix((len(records), 0), dtype=np.float32),
            feature_mz=np.array([], dtype=np.float64),
        )

    kept_bins_arr = np.asarray(kept_bins, dtype=np.int64)
    row_chunks: list[np.ndarray] = []
    col_chunks: list[np.ndarray] = []
    data_chunks: list[np.ndarray] = []

    for row, record in enumerate(records):
        bins = _to_bin(record.signal.mz, config.mz_bin_width)
        if bins.size == 0:
            continue

        order = np.argsort(bins)
        bins_sorted = bins[order]
        intensity_sorted = record.signal.intensity[order]

        uniq_bins, uniq_start = np.unique(bins_sorted, return_index=True)
        aggregated_intensity = np.add.reduceat(intensity_sorted, uniq_start).astype(np.float32, copy=False)

        cols = np.searchsorted(kept_bins_arr, uniq_bins)
        in_bounds = cols < kept_bins_arr.size
        valid = np.zeros_like(in_bounds, dtype=bool)
        valid[in_bounds] = kept_bins_arr[cols[in_bounds]] == uniq_bins[in_bounds]
        if not np.any(valid):
            continue

        valid_cols = cols[valid].astype(np.int32, copy=False)
        valid_data = aggregated_intensity[valid]
        valid_rows = np.full(valid_cols.size, row, dtype=np.int32)

        row_chunks.append(valid_rows)
        col_chunks.append(valid_cols)
        data_chunks.append(valid_data)

    if data_chunks:
        rows = np.concatenate(row_chunks)
        cols = np.concatenate(col_chunks)
        data = np.concatenate(data_chunks)
    else:
        rows = np.array([], dtype=np.int32)
        cols = np.array([], dtype=np.int32)
        data = np.array([], dtype=np.float32)

    matrix = sparse.coo_matrix((data, (rows, cols)), shape=(len(records), kept_bins_arr.size), dtype=np.float32).tocsr()

    feature_mz = np.asarray(kept_bins, dtype=np.float64) * config.mz_bin_width
    return FeatureTable(matrix=matrix, feature_mz=feature_mz)
