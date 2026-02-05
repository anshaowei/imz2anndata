# Changelog

All notable changes to this project are documented in this file.

## [0.1.0] - 2026-02-05

### Added
- Packaging metadata via `pyproject.toml` and console entry point `imz2anndata`.
- Open-source license file (`MIT`).
- Initial `CHANGELOG.md`.
- CI workflow for unit/integration test execution with `unittest`.

### Changed
- Default `min_feature_occurrence` from `1` to `3` to reduce single-pixel noise features in default runs.
- Integration tests now resolve project-root `data/` paths robustly and close parser handles explicitly.
- Real-data notebook path handling stabilized to avoid local subfolder `data` confusion.
