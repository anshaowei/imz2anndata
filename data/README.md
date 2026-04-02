# Data Files

This repository may include example analysis code that expects a local bladder MSI dataset under the `data/` directory, but the project is intended to be published without redistributing the raw dataset files.

## Demo Dataset

The demo notebook uses the following local paths:

- `data/Bladder-MSI.imzML`
- `data/Bladder-MSI.ibd`

The mouse bladder MSI dataset is accessible via the ProteomeXchange Consortium under dataset identifier `PXD001283`.

## Public Repository Policy

- Analysis example code and notebooks may be included in the public repository.
- Raw MSI data files are not required to be committed to the repository.
- To run the demo notebook locally, download the dataset separately and place the paired `imzML` and `ibd` files in this `data/` directory with the filenames shown above.

## Suggested Layout

Recommended local layout for users who want to reproduce the demo:

- `data/Bladder-MSI.imzML`
- `data/Bladder-MSI.ibd`
- `data/results/` for generated `h5ad` outputs

The integration tests in `tests/integration/` discover `.imzML` files from both `data/` and the optional `data/raw/` layout. They skip automatically when compatible local test data is not present.
