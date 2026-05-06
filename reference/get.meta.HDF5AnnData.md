# Extract Sample-Level Metadata from HDF5AnnData Object

Extracts sample-wise metadata from an HDF5-backed AnnData object by
identifying which `obs` columns have consistent values within each
sample. Returns a data frame suitable for `setup.tdr.obj(.meta = ...)`.

## Usage

``` r
get.meta.HDF5AnnData(.h5ad.obj, .sample.var, .verbose = TRUE)
```

## Arguments

- .h5ad.obj:

  An HDF5AnnData object (created via
  `anndataR::read_h5ad(..., backend = "HDF5AnnData")`).

- .sample.var:

  Character: column name in `.h5ad.obj$obs` identifying samples.

- .verbose:

  Logical: print progress messages? Default TRUE.

## Value

Data frame with one row per sample, containing only sample-level
metadata columns. Rownames are sample IDs from `.sample.var`.

## Details

AnnData objects store cell-level annotations in `obs`. This function:

1.  Reads `obs` as a data frame

2.  Groups by `.sample.var`

3.  Identifies columns with unique values within each sample
    (sample-level)

4.  Excludes varying columns (cell-level) with warning

5.  Returns one row per sample

## See also

[`get.meta`](https://opensource.nibr.com/tinydenseR/reference/get.meta.md)
for automatic format detection,
[`get.meta.Seurat`](https://opensource.nibr.com/tinydenseR/reference/get.meta.Seurat.md)
for Seurat objects,
[`get.meta.SCE`](https://opensource.nibr.com/tinydenseR/reference/get.meta.SCE.md)
for SingleCellExperiment objects

## Examples

``` r
if (FALSE) { # \dontrun{
adata <- anndataR::read_h5ad("data.h5ad", backend = "HDF5AnnData")
meta <- get.meta.HDF5AnnData(.h5ad.obj = adata,
                             .sample.var = "sample_id")
} # }
```
