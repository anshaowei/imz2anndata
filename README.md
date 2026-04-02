# imz2anndata

`imz2anndata` converts `imzML` mass spectrometry imaging datasets into `AnnData` (`.h5ad`) for downstream spatial omics analysis.

## What It Does

The conversion pipeline:

1. Reads spectra from `imzML` pixel by pixel.
2. Detects whether the input spectra are `centroid` or `profile`.
3. Applies OpenMS peak picking only for `profile` inputs.
4. Aligns spectra into a shared binned `m/z` feature space.
5. Builds an `AnnData` object with sparse intensities, pixel metadata, and spatial coordinates.

The output schema is summarized in [`docs/anndata_schema.md`](docs/anndata_schema.md).

## Installation

Recommended: Python 3.11 or newer.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install .
```

For local development:

```bash
pip install -r requirements-dev.txt
```

## CLI Usage

```bash
imz2anndata <input.imzML> <output.h5ad> \
  --dataset-id demo \
  --mz-bin-width 0.01 \
  --min-feature-occurrence 3 \
  --peak-picking-snr 0 \
  --signal-to-noise-window 200 \
  --min-intensity 0
```

You can also run:

```bash
python -m imz2anndata.cli <input.imzML> <output.h5ad>
```

Optional spatial filtering flags:

- `--enable-spatial-filtering` filters features by prevalence, total intensity, and Moran's I.
- `--assess-spatial-filtering-only` computes filter statistics without subsetting the matrix.

Peak-picking behavior is automatic by default:

- `centroid` input: peak picking is skipped.
- `profile` input: OpenMS `PeakPickerHiRes` is applied.
- unknown mode: peak picking is skipped and a warning is emitted.

## Python API

```python
from pathlib import Path

from imz2anndata import PipelineConfig, run_pipeline

adata = run_pipeline(
    PipelineConfig(
        input_imzml=Path("sample.imzML"),
        output_h5ad=Path("sample.h5ad"),
        dataset_id="sample",
    )
)
```

## Tests

Run the full test suite:

```bash
python -m unittest discover -s tests -p "test_*.py" -v
```

Run only integration tests:

```bash
python -m unittest discover -s tests/integration -p "test_*.py" -v
```

Integration tests look for local `.imzML` files under `data/raw/` and skip automatically when no test data is present.

## Local Data

Local validation data is not committed by default. If you want to exercise the integration tests with real inputs, place files under:

- `data/raw/sample.imzML`
- `data/raw/sample.ibd`

See [`data/README.md`](data/README.md) for the expected layout.
