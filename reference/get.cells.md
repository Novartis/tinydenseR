# Create .cells Object with Automatic Format Detection

Convenience wrapper that automatically detects input format and calls
the appropriate conversion function. Supports Seurat (v4/v5),
SingleCellExperiment, and list of count matrices. Simplifies workflow by
eliminating need to know which specific function to call.

## Usage

``` r
get.cells(
  .exprs,
  .meta,
  .sample.var = NULL,
  .assay = "RNA",
  .layer.pattern = "counts",
  .slot = "counts",
  .min.cells.per.sample = 10,
  .compress = FALSE,
  .verbose = TRUE
)
```

## Arguments

- .exprs:

  Expression data in one of three formats:

  - Named list of count matrices (calls `get.cells.list.mat`)

  - Seurat object v4 or v5 (calls `get.cells.Seurat` or
    `get.cells.Seurat5`)

  - SingleCellExperiment object (calls `get.cells.SCE`)

- .meta:

  Data frame with sample-level metadata. Required for list input.
  Rownames must match sample IDs. For Seurat/SCE, rownames must match
  values in `.sample.var` column.

- .sample.var:

  Character: metadata column identifying samples. Required for
  Seurat/SCE, ignored for list input.

- .assay:

  Character: assay name. Default "RNA" for Seurat, "counts" for SCE.
  Ignored for list input.

- .layer.pattern:

  Character: layer name pattern for Seurat v5. Default "counts". Ignored
  for other formats.

- .slot:

  Character: data slot for Seurat v4. Default "counts". Ignored for
  other formats.

- .min.cells.per.sample:

  Integer: minimum cells per sample. Default 10. Ignored for list input.

- .compress:

  Logical: compress RDS files? Default FALSE.

- .verbose:

  Logical: print progress messages? Default TRUE.

## Value

Named list of file paths to temporary RDS files, one per sample.

## Details

This function inspects the class of `.exprs` to determine format:

- `Seurat`: Checks for Assay5 class to distinguish v4 from v5

- `SingleCellExperiment`: Maps "RNA" assay to "counts" if needed

- `list`: Validates names match `.meta` rownames

For Seurat v5, automatically handles layered data structure. For SCE,
uses standard `assay()` accessor.

## See also

[`get.cells.list.mat`](https://opensource.nibr.com/tinydenseR/reference/get.cells.list.mat.md),
[`get.cells.Seurat`](https://opensource.nibr.com/tinydenseR/reference/get.cells.Seurat.md),
[`get.cells.Seurat5`](https://opensource.nibr.com/tinydenseR/reference/get.cells.Seurat5.md),
[`get.cells.SCE`](https://opensource.nibr.com/tinydenseR/reference/get.cells.SCE.md),
[`setup.tdr.obj`](https://opensource.nibr.com/tinydenseR/reference/setup.tdr.obj.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: List of count matrices
count.matrices <- list(A_R1 = counts.A_R1, B_R1 = counts.B_R1)
cells <- get.cells(.exprs = count.matrices,
                   .meta = sim_trajectory.meta[c("A_R1", "B_R1"), ])

# Example 2: Seurat object (auto-detects v4 vs v5)
cells <- get.cells(.exprs = seurat.obj,
                   .meta = sim_trajectory.meta,
                   .sample.var = "sample_id")

# Example 3: SingleCellExperiment object
cells <- get.cells(.exprs = sce.obj,
                   .meta = sim_trajectory.meta,
                   .sample.var = "Sample",
                   .assay = "counts")
} # }
```
