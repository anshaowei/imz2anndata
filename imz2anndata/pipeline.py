from __future__ import annotations

from imz2anndata.aligner import align_records
from imz2anndata.anndata_builder import build_anndata
from imz2anndata.config import PipelineConfig
from imz2anndata.extractors.openms import OpenMSPeakExtractor
from imz2anndata.io import iter_imzml_spectra
from imz2anndata.models import SpectrumRecord


class IdentityExtractor:
    def extract(self, signal):
        return signal


def run_pipeline(config: PipelineConfig):
    extractor = OpenMSPeakExtractor(config.peak_picking) if config.enable_peak_picking else IdentityExtractor()

    records: list[SpectrumRecord] = []
    raw_tic: list[float] = []
    raw_peak_count: list[int] = []
    extracted_tic: list[float] = []
    extracted_peak_count: list[int] = []

    for record in iter_imzml_spectra(config.input_imzml):
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
    adata.uns["peak_picking_enabled"] = config.enable_peak_picking
    adata.write_h5ad(config.output_h5ad)
    return adata
