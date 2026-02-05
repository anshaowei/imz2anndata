# Test Data Layout

Put local validation files here (not committed to git by default).

Recommended structure:

- `data/<dataset>.imzML` - direct placement is supported and auto-discovered by tests/notebooks.
- `data/<dataset>.ibd` - paired binary file required by imzML parser.
- `data/real/sample.imzML` - real MSI file for integration test.
- `data/real/sample.ibd` - paired binary file required by imzML parser.
- `data/demo/` - optional tiny demo files for quick checks.

The tests `tests/integration/test_real_data_e2e.py` and `tests/integration/test_dataset_consistency.py`
automatically discover `.imzML` files in both `data/` and `data/real/`.
