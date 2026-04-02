# imz2anndata

`imz2anndata` converts `imzML` mass spectrometry imaging datasets into `AnnData` (`.h5ad`) for downstream spatial omics analysis.

Mass spectrometry imaging (MSI) provides spatially resolved molecular measurements, but standard `imzML` files are not directly optimized for matrix-based computational workflows. `imz2anndata` turns `imzML` data into an analysis-ready `AnnData` representation by parsing spectra with spatial coordinates, applying adaptive spectral preprocessing, aligning signals into a shared `m/z` feature space, and building a sparse pixel-by-feature matrix while preserving spatial structure and preprocessing provenance.

The workflow also supports spatial feature assessment and filtering based on abundance and spatial autocorrelation, making it easier to move from raw MSI exchange files to scalable downstream analysis in the Python scientific ecosystem.

## What It Does

The conversion pipeline:

1. Reads spectra from `imzML` pixel by pixel.
2. Detects whether the input spectra are `centroid` or `profile`.
3. Applies OpenMS peak picking only for `profile` inputs.
4. Aligns spectra into a shared binned `m/z` feature space.
5. Builds an `AnnData` object with sparse intensities, pixel metadata, and spatial coordinates.

The output schema is summarized in [`docs/anndata_schema.md`](docs/anndata_schema.md).

## Why AnnData

`AnnData` provides a compact, matrix-oriented representation that fits naturally into modern analysis workflows. In the converted object:

- `X` stores the sparse pixel-by-feature intensity matrix.
- `obs` stores pixel-level metadata such as coordinates and TIC summaries.
- `var` stores aligned `m/z` features.
- `obsm["spatial"]` preserves the tissue coordinate system for spatial downstream analysis.

This makes the converted MSI data easier to inspect, subset, benchmark, serialize, and integrate with existing Python tools.

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
pip install -r requirements.txt
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

## Examples

An end-to-end bladder MSI demo notebook is available at [`examples/bladder_msi_demo.ipynb`](examples/bladder_msi_demo.ipynb).

The notebook expects the local demo dataset files at:

- `data/Bladder-MSI.imzML`
- `data/Bladder-MSI.ibd`

The mouse bladder MSI dataset is accessible via the ProteomeXchange Consortium under dataset identifier `PXD001283`.

## Output Characteristics

The converted `AnnData` object is designed to preserve information needed for downstream spatial analysis while remaining compact and practical to store:

- sparse matrix storage for pixel-by-feature intensities
- retained spatial coordinates for each pixel
- recorded preprocessing decisions such as detected spectrum mode and peak-picking behavior
- optional feature-level spatial statistics for assessment or filtering

## Tests

Run the full test suite:

```bash
python -m unittest discover -s tests -p "test_*.py" -v
```

Run only integration tests:

```bash
python -m unittest discover -s tests/integration -p "test_*.py" -v
```

Integration tests look for local `.imzML` files under `data/` and `data/raw/`, and skip automatically when no test data is present.

## Local Data

The repository does not include the raw bladder MSI dataset. To run the demo notebook or the real-data integration tests, download the dataset separately and place the files at:

- `data/Bladder-MSI.imzML`
- `data/Bladder-MSI.ibd`

The mouse bladder MSI dataset is accessible via the ProteomeXchange Consortium under dataset identifier `PXD001283`.

The integration tests also support an optional alternative layout:

- `data/raw/sample.imzML`
- `data/raw/sample.ibd`

See [`data/README.md`](data/README.md) for the expected layout.
