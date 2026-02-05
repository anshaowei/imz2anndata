# imz2anndata

`imz2anndata` converts **imzML (spatial metabolomics / MSI)** data into **AnnData (`.h5ad`)**.

## Pipeline

1. Read spectra from `imzML` (pixel by pixel)
2. Extract signals with OpenMS peak picking
3. Align and merge spectra in a shared `m/z` feature space (binned)
4. Build `AnnData`
   - `X`: sparse intensity matrix (`pixel x feature`)
   - `obs`: pixel coordinates and spectrum metadata
   - `var`: feature `m/z`

## Installation

Recommended: Python 3.11+.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -r requirements.txt
```

Or install as a package (recommended for release usage):

```bash
pip install .
```

For development and testing:

```bash
pip install -r requirements-dev.txt
```

For notebook execution (if not already available in your environment):

```bash
pip install notebook
```

## CLI usage

```bash
python main.py <input.imzML> <output.h5ad> \
  --dataset-id demo \
  --mz-bin-width 0.01 \
  --min-feature-occurrence 3 \
  --min-intensity 0
```

Example:

```bash
python main.py data/real/sample.imzML results/sample.h5ad --dataset-id pilot
```

## Tests

### Unit and integration tests

- `tests/unit/test_aligner.py`: `m/z` alignment and sparse matrix construction
- `tests/unit/test_anndata_builder.py`: AnnData field integrity
- `tests/unit/test_pipeline.py`: pipeline output writing and mocked end-to-end behavior
- `tests/unit/test_cli.py`: CLI argument parsing
- `tests/unit/test_boundaries.py`: irregular-coordinate boundary behavior
- `tests/integration/test_real_data_e2e.py`: real-data end-to-end run
- `tests/integration/test_dataset_consistency.py`: dataset-based consistency, reproducibility, and ion-image reference checks

Run all tests:

```bash
python -m unittest discover -s tests -p "test_*.py"
```

Run only integration tests:

```bash
python -m unittest discover -s tests/integration -p "test_*.py"
```

### Notebook-based interactive tests

- `notebooks/01_smoke_test.ipynb`: quick synthetic test for matrix/AnnData behavior
- `notebooks/02_real_imzml_test.ipynb`: real-data validation with consistency metrics, TIC image, and ion images

Run notebooks:

```bash
jupyter lab
```

## Local test data layout

Place local validation files under `data/` (not committed by default):

- `data/real/sample.imzML`
- `data/real/sample.ibd`

See `data/README.md` for details.

For quick runs, `.imzML` files directly under `data/` are also auto-discovered.

## Feature filtering default

By default, `min_feature_occurrence=3`, meaning a feature bin must appear in at least
3 spectra/pixels to be kept in the output matrix. This reduces one-off single-pixel noise
features in typical MSI runs.
