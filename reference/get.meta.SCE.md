# Extract Sample-Level Metadata from SingleCellExperiment Object

Extracts sample-wise metadata from a SingleCellExperiment object by
identifying which colData columns have consistent values within each
sample. Identical logic to Seurat metadata extraction but adapted for
SCE's `colData` structure. Returns a data frame suitable for
`setup.tdr.obj(.meta = ...)`.

## Usage

``` r
get.meta.SCE(.sce.obj, .sample.var, .verbose = TRUE)
```

## Arguments

- .sce.obj:

  SingleCellExperiment object with cell-level metadata in `colData`.

- .sample.var:

  Character: column name in `colData(.sce.obj)` identifying samples.

- .verbose:

  Logical: print progress messages? Default TRUE.

## Value

Data frame with one row per sample, containing only sample-level
metadata columns. Rownames are sample IDs from `.sample.var`.

## Details

SingleCellExperiment stores cell-level metadata in `colData`. This
function:

1.  Converts `colData` to data frame

2.  Groups by `.sample.var`

3.  Identifies columns with unique values within each sample
    (sample-level)

4.  Excludes varying columns (cell-level) with warning

5.  Returns one row per sample

Common sample-level variables include experimental conditions, batches,
donors. Common cell-level variables include QC metrics, clusters, cell
types.

## See also

[`get.meta`](https://opensource.nibr.com/tinydenseR/reference/get.meta.md)
for automatic format detection,
[`get.meta.Seurat`](https://opensource.nibr.com/tinydenseR/reference/get.meta.Seurat.md)
for Seurat objects,
[`get.cells.SCE`](https://opensource.nibr.com/tinydenseR/reference/get.cells.SCE.md)
for extracting count data

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract sample metadata from SCE object
meta <- get.meta.SCE(.sce.obj = sce.obj,
                     .sample.var = "sample_id")

# Use with get.cells
cells <- get.cells.SCE(.sce.obj = sce.obj,
                       .meta = meta,
                       .sample.var = "sample_id")
} # }
```
