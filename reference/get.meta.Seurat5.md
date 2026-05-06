# Extract Sample-Level Metadata from Seurat v5 Object

Extracts sample-wise metadata from a Seurat v5 object by identifying
which metadata columns have consistent values within each sample
(sample-level variables). Cell-level variables are excluded with a
warning. Returns a data frame suitable for `setup.tdr.obj(.meta = ...)`.

## Usage

``` r
get.meta.Seurat5(.seurat.obj, .sample.var, .verbose = TRUE)
```

## Arguments

- .seurat.obj:

  Seurat v5 object with cell-level metadata in `@meta.data`.

- .sample.var:

  Character: column name in `.seurat.obj@meta.data` identifying samples.

- .verbose:

  Logical: print progress messages? Default TRUE.

## Value

Data frame with one row per sample, containing only sample-level
metadata columns. Rownames are sample IDs from `.sample.var`.

## Details

This function automatically distinguishes sample-level from cell-level
metadata:

- **Sample-level**: Variables with identical values for all cells within
  each sample (e.g., Condition, Batch, Treatment). These are kept.

- **Cell-level**: Variables that vary between cells in the same sample
  (e.g., nCount_RNA, percent.mt, cluster). These are excluded with a
  warning.

The `.sample.var` column is always included. Only distinct rows are
returned (one per sample).

Useful when you've stored experimental design information in the Seurat
object and want to extract it for tinydenseR analysis.

## See also

[`get.meta`](https://opensource.nibr.com/tinydenseR/reference/get.meta.md)
for automatic format detection,
[`get.meta.Seurat`](https://opensource.nibr.com/tinydenseR/reference/get.meta.Seurat.md)
for Seurat v4,
[`get.cells.Seurat5`](https://opensource.nibr.com/tinydenseR/reference/get.cells.Seurat5.md)
for extracting count data

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming Seurat object has sample metadata
meta <- get.meta.Seurat5(.seurat.obj = seurat.obj,
                         .sample.var = "sample_id")

# Use with get.cells
cells <- get.cells.Seurat5(.seurat.obj = seurat.obj,
                           .meta = meta,
                           .sample.var = "sample_id")
} # }
```
