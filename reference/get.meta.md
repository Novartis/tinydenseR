# Extract Sample-Level Metadata with Automatic Format Detection

Convenience wrapper that automatically detects input format (Seurat
v4/v5, SingleCellExperiment, or HDF5AnnData) and calls the appropriate
metadata extraction function. Simplifies workflow by eliminating need to
know object version.

## Usage

``` r
get.meta(.obj, .sample.var, .verbose = TRUE)
```

## Arguments

- .obj:

  Object containing cell-level metadata. Can be:

  - Seurat object (v4 or v5) - calls `get.meta.Seurat` or
    `get.meta.Seurat5`

  - SingleCellExperiment object - calls `get.meta.SCE`

  - HDF5AnnData object - calls `get.meta.HDF5AnnData`

- .sample.var:

  Character: column name in object metadata identifying samples.

- .verbose:

  Logical: print progress messages? Default TRUE.

## Value

Data frame with one row per sample, containing only sample-level
metadata columns. Rownames are sample IDs.

## Details

This function inspects the class of `.obj` to determine format:

- `Seurat`: Checks first assay for Assay5 class to distinguish v4 from
  v5

- `SingleCellExperiment`: Uses standard `colData` accessor

- `HDF5AnnData`: Uses `obs` accessor

Automatically filters to sample-level metadata by testing which columns
have consistent values within each sample. Cell-level columns (varying
within samples) are excluded with a warning.

## See also

[`get.meta.Seurat`](https://opensource.nibr.com/tinydenseR/reference/get.meta.Seurat.md),
[`get.meta.Seurat5`](https://opensource.nibr.com/tinydenseR/reference/get.meta.Seurat5.md),
[`get.meta.SCE`](https://opensource.nibr.com/tinydenseR/reference/get.meta.SCE.md),
[`get.meta.HDF5AnnData`](https://opensource.nibr.com/tinydenseR/reference/get.meta.HDF5AnnData.md),
[`get.cells`](https://opensource.nibr.com/tinydenseR/reference/get.cells.md)
for extracting count data

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Seurat object (auto-detects v4 vs v5)
meta <- get.meta(.obj = seurat.obj,
                 .sample.var = "sample_id")

# Example 2: SingleCellExperiment object
meta <- get.meta(.obj = sce.obj,
                 .sample.var = "Sample")

# Example 3: HDF5AnnData object
meta <- get.meta(.obj = adata,
                 .sample.var = "sample_id")

# Use with get.cells
cells <- get.cells(.exprs = seurat.obj,
                   .meta = meta,
                   .sample.var = "sample_id")
} # }
```
