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
    bin_to_col = {b: i for i, b in enumerate(kept_bins)}

    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []

    for row, record in enumerate(records):
        bins = _to_bin(record.signal.mz, config.mz_bin_width)
        for b, intensity in zip(bins, record.signal.intensity, strict=False):
            col = bin_to_col.get(int(b))
            if col is None:
                continue
            rows.append(row)
            cols.append(col)
            data.append(float(intensity))

    matrix = sparse.coo_matrix(
        (data, (rows, cols)),
        shape=(len(records), len(kept_bins)),
        dtype=np.float32,
    ).tocsr()
    matrix.sum_duplicates()

    feature_mz = np.asarray(kept_bins, dtype=np.float64) * config.mz_bin_width
    return FeatureTable(matrix=matrix, feature_mz=feature_mz)
