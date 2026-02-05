from __future__ import annotations

from typing import Protocol

from imz2anndata.models import SpectrumSignal


class SignalExtractor(Protocol):
    def extract(self, signal: SpectrumSignal) -> SpectrumSignal:
        ...
