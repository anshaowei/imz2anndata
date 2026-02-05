from __future__ import annotations

from pathlib import Path
from typing import Iterator

import numpy as np
from pyimzml.ImzMLParser import ImzMLParser

from imz2anndata.models import PixelCoordinate, SpectrumRecord, SpectrumSignal


def iter_imzml_spectra(imzml_path: Path) -> Iterator[SpectrumRecord]:
    parser = ImzMLParser(str(imzml_path))
    try:
        for spectrum_id, (x, y, *rest) in enumerate(parser.coordinates):
            z = rest[0] if rest else None
            mz_values, intensity_values = parser.getspectrum(spectrum_id)

            yield SpectrumRecord(
                spectrum_id=spectrum_id,
                coordinate=PixelCoordinate(x=int(x), y=int(y), z=int(z) if z is not None else None),
                signal=SpectrumSignal(
                    mz=np.asarray(mz_values, dtype=np.float64),
                    intensity=np.asarray(intensity_values, dtype=np.float32),
                ),
            )
    finally:
        ibd_handle = getattr(parser, "m", None)
        if ibd_handle is not None and hasattr(ibd_handle, "close"):
            ibd_handle.close()
