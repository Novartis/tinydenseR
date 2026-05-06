# Extract Sample-Level Metadata from Seurat v4 Object

Extracts sample-wise metadata from a Seurat v4 object by identifying
which metadata columns have consistent values within each sample.
Identical functionality to `get.meta.Seurat5` but for v4 objects.
Returns a data frame suitable for `setup.tdr.obj(.meta = ...)`.

## Usage

``` r
get.meta.Seurat(.seurat.obj, .sample.var, .verbose = TRUE)
```

## Arguments

- .seurat.obj:

  Seurat v4 object with cell-level metadata in `@meta.data`.

- .sample.var:

  Character: column name in `.seurat.obj@meta.data` identifying samples.

- .verbose:

  Logical: print progress messages? Default TRUE.

## Value

Data frame with one row per sample, containing only sample-level
metadata columns. Rownames are sample IDs from `.sample.var`.

## Details

Automatically distinguishes sample-level (consistent within samples)
from cell-level (varying within samples) metadata. Only sample-level
columns are retained. Cell-level columns are excluded with a warning if
`.verbose = TRUE`.

The algorithm groups by `.sample.var` and tests each column for
uniqueness within groups. Columns with multiple unique values within any
sample are classified as cell-level.

## See also

[`get.meta`](https://opensource.nibr.com/tinydenseR/reference/get.meta.md)
for automatic format detection,
[`get.meta.Seurat5`](https://opensource.nibr.com/tinydenseR/reference/get.meta.Seurat5.md)
for Seurat v5,
[`get.cells.Seurat`](https://opensource.nibr.com/tinydenseR/reference/get.cells.Seurat.md)
for extracting count data

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract sample metadata from Seurat v4 object
meta <- get.meta.Seurat(.seurat.obj = seurat.obj,
                        .sample.var = "sample_id")

# Use with get.cells
cells <- get.cells.Seurat(.seurat.obj = seurat.obj,
                          .meta = meta,
                          .sample.var = "sample_id")
} # }
```
