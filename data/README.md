# Data Layout

Place local validation data here when running integration tests with real inputs.

Recommended structure:

- `data/raw/` for source `imzML` and paired `ibd` files
- `data/results/` for generated `h5ad` outputs or temporary summaries

Examples:

- `data/raw/sample.imzML`
- `data/raw/sample.ibd`
- `data/results/sample.h5ad`

The integration tests in `tests/integration/` automatically scan `data/raw/` for `.imzML` files and skip when none are present.
