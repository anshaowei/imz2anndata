# Changelog

All notable changes to this project are documented in this file.

## [0.1.0] - 2026-04-02

Initial public release.

### Added
- Core `imzML` to `AnnData` conversion pipeline with CLI and Python API support.
- Automatic spectrum-mode detection with conditional OpenMS peak picking for profile-mode inputs.
- `m/z` alignment, sparse feature-matrix construction, and AnnData export to `.h5ad`.
- Optional spatial feature filtering and assessment utilities.
- Unit and integration test suite.
- AnnData schema documentation in [`docs/anndata_schema.md`](docs/anndata_schema.md).
- Public example notebook in [`examples/bladder_msi_demo.ipynb`](examples/bladder_msi_demo.ipynb).
- Dataset usage and redistribution guidance in [`data/README.md`](data/README.md).

### Notes
- The repository includes example code and tests, but raw MSI datasets are not redistributed with the public release.
