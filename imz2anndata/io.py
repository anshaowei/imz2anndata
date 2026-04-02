from __future__ import annotations

from pathlib import Path
from typing import Iterator

import numpy as np
from pyimzml.ImzMLParser import ImzMLParser

from imz2anndata.models import PixelCoordinate, SpectrumRecord, SpectrumSignal


def _detect_spectrum_mode_from_text(imzml_path: Path) -> str | None:
    text = imzml_path.read_text(encoding="utf-8", errors="ignore")
    if 'accession="MS:1000127"' in text or 'name="centroid spectrum"' in text:
        return "centroid"
    if 'accession="MS:1000128"' in text or 'name="profile spectrum"' in text:
        return "profile"
    return None


def detect_imzml_spectrum_mode(imzml_path: Path) -> str | None:
    if not imzml_path.exists():
        return None
    try:
        parser = ImzMLParser(str(imzml_path), ibd_file=None)
    except Exception:
        return _detect_spectrum_mode_from_text(imzml_path)

    try:
        if parser.spectrum_mode in {"centroid", "profile"}:
            return parser.spectrum_mode
    finally:
        ibd_handle = getattr(parser, "m", None)
        if ibd_handle is not None and hasattr(ibd_handle, "close"):
            ibd_handle.close()

    return _detect_spectrum_mode_from_text(imzml_path)


def _records_from_parser(parser: ImzMLParser) -> list[SpectrumRecord]:
    records: list[SpectrumRecord] = []
    for spectrum_id, (x, y, *rest) in enumerate(parser.coordinates):
        z = rest[0] if rest else None
        mz_values, intensity_values = parser.getspectrum(spectrum_id)
        records.append(
            SpectrumRecord(
                spectrum_id=spectrum_id,
                coordinate=PixelCoordinate(x=int(x), y=int(y), z=int(z) if z is not None else None),
                signal=SpectrumSignal(
                    mz=np.asarray(mz_values, dtype=np.float64),
                    intensity=np.asarray(intensity_values, dtype=np.float32),
                ),
            )
        )
    return records


def load_imzml_records(imzml_path: Path) -> tuple[str | None, list[SpectrumRecord]]:
    parser = ImzMLParser(str(imzml_path))
    try:
        spectrum_mode = parser.spectrum_mode if parser.spectrum_mode in {"centroid", "profile"} else None
        if spectrum_mode is None:
            spectrum_mode = _detect_spectrum_mode_from_text(imzml_path)
        return spectrum_mode, _records_from_parser(parser)
    finally:
        ibd_handle = getattr(parser, "m", None)
        if ibd_handle is not None and hasattr(ibd_handle, "close"):
            ibd_handle.close()


def iter_imzml_spectra(imzml_path: Path) -> Iterator[SpectrumRecord]:
    parser = ImzMLParser(str(imzml_path))
    try:
        yield from _records_from_parser(parser)
    finally:
        ibd_handle = getattr(parser, "m", None)
        if ibd_handle is not None and hasattr(ibd_handle, "close"):
            ibd_handle.close()
