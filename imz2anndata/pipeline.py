from __future__ import annotations

import warnings

from imz2anndata.aligner import align_records
from imz2anndata.anndata_builder import build_anndata
from imz2anndata.config import PipelineConfig
from imz2anndata.extractors.openms import OpenMSPeakExtractor
from imz2anndata.io import load_imzml_records
from imz2anndata.models import SpectrumRecord
from imz2anndata.spatial_filtering import filter_spatial_features


class IdentityExtractor:
    def extract(self, signal):
        return signal


def run_pipeline(config: PipelineConfig):
    spectrum_mode, input_records = load_imzml_records(config.input_imzml)
    peak_picking_requested = config.enable_peak_picking
    peak_picking_effective = peak_picking_requested
    peak_picking_strategy = "requested" if peak_picking_requested else "disabled_by_user"

    if spectrum_mode == "centroid" and peak_picking_requested:
        warnings.warn(
            "Detected centroid-mode imzML input; skipping peak picking and using the centroid spectra directly.",
            RuntimeWarning,
            stacklevel=2,
        )
        peak_picking_effective = False
        peak_picking_strategy = "auto_skipped_for_centroid"
    elif spectrum_mode == "profile" and peak_picking_requested:
        warnings.warn(
            "Detected profile-mode imzML input; applying OpenMS peak picking before alignment.",
            RuntimeWarning,
            stacklevel=2,
        )
        peak_picking_strategy = "auto_applied_for_profile"
    elif spectrum_mode is None and peak_picking_requested:
        warnings.warn(
            "Could not detect imzML spectrum mode from metadata; skipping peak picking by default and using the original spectra.",
            RuntimeWarning,
            stacklevel=2,
        )
        peak_picking_effective = False
        peak_picking_strategy = "auto_skipped_for_unknown_mode"

    extractor = OpenMSPeakExtractor(config.peak_picking) if peak_picking_effective else IdentityExtractor()

    records: list[SpectrumRecord] = []
    raw_tic: list[float] = []
    raw_peak_count: list[int] = []
    extracted_tic: list[float] = []
    extracted_peak_count: list[int] = []

    for record in input_records:
        raw_tic.append(float(record.signal.intensity.sum()))
        raw_peak_count.append(int(record.signal.mz.size))
        record.signal = extractor.extract(record.signal)
        extracted_tic.append(float(record.signal.intensity.sum()))
        extracted_peak_count.append(int(record.signal.mz.size))
        records.append(record)

    feature_table = align_records(records, config.alignment)
    adata = build_anndata(records, feature_table, dataset_id=config.dataset_id)
    adata.obs["tic_raw"] = raw_tic
    adata.obs["tic_processed"] = extracted_tic
    adata.obs["raw_peak_count"] = raw_peak_count
    adata.obs["processed_peak_count"] = extracted_peak_count
    adata.uns["spectrum_mode"] = spectrum_mode or "unknown"
    adata.uns["peak_picking_requested"] = peak_picking_requested
    adata.uns["peak_picking_enabled"] = peak_picking_effective
    adata.uns["peak_picking_strategy"] = peak_picking_strategy
    adata = filter_spatial_features(adata, config.spatial_filter)
    adata.write_h5ad(config.output_h5ad)
    return adata
