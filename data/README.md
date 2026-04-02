# Data Files

Raw MSI data files are not distributed with this repository.

## Demo Dataset

The demo notebook expects the following local paths:

- `data/Bladder-MSI.imzML`
- `data/Bladder-MSI.ibd`

The mouse bladder MSI dataset is accessible via the ProteomeXchange Consortium under dataset identifier `PXD001283`.

## Download And Placement

To run the demo notebook locally, download the dataset separately and place the paired `imzML` and `ibd` files in this directory using the filenames shown above.

## Suggested Layout

Recommended local layout:

- `data/Bladder-MSI.imzML`
- `data/Bladder-MSI.ibd`
- `data/results/` for generated `h5ad` outputs

The integration tests in `tests/integration/` also discover `.imzML` files from the optional `data/raw/` layout. They skip automatically when compatible local test data is not present.
