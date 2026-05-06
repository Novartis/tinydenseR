# Summarize gene expression patterns for genes of interest

Generates summaries of gene expression for specified genes of interest
(GOI), including detection status (positive/negative), cell counts, and
percentages across samples. Useful for identifying marker gene
expression patterns across clusters or cell types.

## Usage

``` r
goi.summary(x, ...)

# S3 method for class 'TDRObj'
goi.summary(
  x,
  .source = NULL,
  .goi,
  .id.idx = NULL,
  .id = NULL,
  .id.from = "clustering",
  .verbose = TRUE,
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object
  initialized with `setup.tdr.obj` and processed with `get.graph` and
  `get.map`. Must contain RNA assay data.

- ...:

  Additional arguments passed to methods.

- .source:

  The raw data object for non-file backends. `NULL` (default) for the
  files backend; otherwise a Seurat, SingleCellExperiment, or anndataR
  AnnData object. Used by `.get_sample_matrix()` to retrieve per-sample
  expression matrices.

- .goi:

  Character vector of gene names to summarize. All genes must exist in
  `.tdr.obj$raw.landmarks` column names.

- .id.idx:

  Integer index of a specific landmark to analyze. If NULL (default),
  uses `.id` parameter or all cells. Rarely used directly.

- .id:

  Character vector of cluster or cell type IDs to filter analysis. If
  NULL (default), includes all clusters/cell types. Use this to focus on
  specific populations (e.g., `c("cluster.01", "cluster.02")`).

- .id.from:

  Character specifying whether `.id` refers to "clustering" or
  "celltyping" identifiers. Default is "clustering".

- .verbose:

  Logical, whether to print progress messages. Default TRUE.

## Value

A named list with one element per gene in `.goi`. Each gene's element
contains:

- `$clustering`:

  List with `$ids`, `$cell.count`, and `$cell.perc` for cluster-level
  summaries with pos./neg. prefixes

- `$celltyping`:

  Same structure as clustering but for cell types (NULL if cell typing
  not performed)

- `$all`:

  Same structure but only pos./neg. categories without cluster/cell type
  breakdown

Cell counts are matrices with samples as rows and categories as columns.
Cell percentages show within-sample proportions.

## Details

For each gene of interest, the function:

1.  Determines which cells express the gene (pos.) vs don't express it
    (neg.)

2.  Creates three parallel summary structures:

    - `clustering` - summaries by cluster IDs with pos./neg. prefix

    - `celltyping` - summaries by cell type labels with pos./neg. prefix
      (if available)

    - `all` - summaries by pos./neg. status only

3.  Computes cell counts and percentages for each category per sample

The function operates on mapped cells (after `get.map`), allowing you to
optionally filter to specific clusters or cell types using `.id`
parameter.

## Examples

``` r
if (FALSE) { # \dontrun{
# Complete workflow with gene expression summary
lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks() |>
  get.graph() |>
  get.map()

# Summarize expression of marker genes
markers <- c("CD4", "CD8A", "CD19")
goi_results <- goi.summary(lm.cells, markers)

# Check CD4 expression by cluster
goi_results$CD4$clustering$cell.perc

# Focus analysis on specific clusters only
goi_results <- goi.summary(lm.cells, markers, 
                           .id = c("cluster.01", "cluster.02"))
} # }
```
