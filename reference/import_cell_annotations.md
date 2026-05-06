# Import multiple cell-level annotation columns as landmark-level celltyping solutions

Scans a cell-level metadata data.frame for categorical columns and
imports each as a named landmark-level celltyping solution via
[`celltyping`](https://opensource.nibr.com/tinydenseR/reference/celltyping.md)
Mode B. After importing, users can switch between solutions with
[`set_active_celltyping`](https://opensource.nibr.com/tinydenseR/reference/set_active_celltyping.md).

## Usage

``` r
import_cell_annotations(x, ...)

# S3 method for class 'TDRObj'
import_cell_annotations(
  x,
  .cell.meta,
  .sample.var,
  .celltype.vec = NULL,
  .verbose = TRUE,
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, or SingleCellExperiment object.

- ...:

  Additional arguments passed to methods.

- .cell.meta:

  A `data.frame` with rownames set to cell IDs (original IDs before the
  tinydenseR sample prefix). Must contain at least the `.sample.var`
  column.

- .sample.var:

  Character(1). Column name in `.cell.meta` identifying sample
  membership.

- .celltype.vec:

  Character(1) or `NULL`. If non-NULL, the name of a column in
  `.cell.meta` that should be set as the active `$ids` after import.

- .verbose:

  Logical (default `TRUE`). If `TRUE`, prints a message listing the
  imported columns.

## Value

The updated object with imported celltyping solutions stored in
`@landmark.annot$celltyping$<column_name>`.

## Details

This function is designed to be called after
[`get.graph`](https://opensource.nibr.com/tinydenseR/reference/get.graph.md)
and **before**
[`get.map`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
in the pipeline. Because
[`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
has not yet run,
[`.refresh_celltyping()`](https://opensource.nibr.com/tinydenseR/reference/dot-refresh_celltyping.md)
inside each
[`celltyping()`](https://opensource.nibr.com/tinydenseR/reference/celltyping.md)
call is a no-op, making the loop efficient.

## Column selection criteria

A column in `.cell.meta` is imported if it is:

- `character` or `factor` type (not numeric, logical, integer, etc.)

- Not the `.sample.var` column

- Has more than 1 unique non-NA value

- Has fewer unique non-NA values than total cells (excludes ID-like
  columns)

- Column name is not `"ids"` (reserved by the multi-solution system)

## See also

[`celltyping`](https://opensource.nibr.com/tinydenseR/reference/celltyping.md),
[`set_active_celltyping`](https://opensource.nibr.com/tinydenseR/reference/set_active_celltyping.md),
[`list_celltyping_solutions`](https://opensource.nibr.com/tinydenseR/reference/list_celltyping_solutions.md)

## Examples

``` r
if (FALSE) { # \dontrun{
tdr.obj <- setup.tdr.obj(.cells, .meta) |>
  get.landmarks() |>
  get.graph()

# Import all categorical cell-level columns
tdr.obj <- import_cell_annotations(tdr.obj,
  .cell.meta = cell_metadata,
  .sample.var = "sample_id",
  .celltype.vec = "cell_type_l1"
)

# Check what was imported
list_celltyping_solutions(tdr.obj)

# Switch to a different annotation
tdr.obj <- set_active_celltyping(tdr.obj, "cell_type_l2")
} # }
```
