# AnnData Schema

`imz2anndata` writes a single `AnnData` object per input dataset.

## Matrix

- `X`: sparse intensity matrix with shape `n_pixels x n_features`

## Observations

Each row in `obs` represents one pixel spectrum.

- `spectrum_id`: sequential spectrum identifier from the parsed input
- `x`: pixel x coordinate
- `y`: pixel y coordinate
- `z`: pixel z coordinate when available, otherwise `NaN`
- `tic_raw`: total ion current before preprocessing
- `tic_processed`: total ion current after peak picking or pass-through preprocessing
- `raw_peak_count`: number of peaks in the original spectrum
- `processed_peak_count`: number of peaks after preprocessing

## Variables

Each row in `var` represents one aligned feature bin.

- `mz`: representative `m/z` value for the feature bin

When spatial filtering is enabled or assessed, additional feature-level statistics may also be present in `var`, such as selection flags or Moran's I values.

## Spatial Coordinates

- `obsm["spatial"]`: integer pixel coordinates with columns `x` and `y`

## Unstructured Metadata

- `uns["dataset_id"]`: dataset identifier supplied by the user
- `uns["feature_space"]`: currently `"binned_mz"`
- `uns["spectrum_mode"]`: detected input mode, one of `centroid`, `profile`, or `unknown`
- `uns["peak_picking_requested"]`: whether peak picking was requested by configuration
- `uns["peak_picking_enabled"]`: whether peak picking was actually applied
- `uns["peak_picking_strategy"]`: how the final peak-picking decision was made

When spatial filtering is used, `uns["spatial_filter"]` records the applied settings and summary statistics.
