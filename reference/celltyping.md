# Manually assign cell type labels to clusters

This function allows manual annotation of clusters by mapping one or
more clusters to biologically meaningful cell type labels. This is
particularly useful in cytometry experiments, where there typically is a
larger number of cells and lineage marker-based hierarchical analysis is
common.

## Usage

``` r
celltyping(x, ...)

# S3 method for class 'TDRObj'
celltyping(x, .celltyping.map, .verbose = TRUE, .name = NULL, ...)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object.

- ...:

  Additional arguments passed to methods.

- .celltyping.map:

  Cell type assignments, supplied in one of two mutually exclusive
  formats:

  Mode A — cluster map (named `list`)

  :   List names are cell-type labels (e.g., `"CD4.T.cells"`), and each
      element is a character vector of cluster IDs (e.g.,
      `c("cluster.01", "cluster.02")`) that belong to that cell type.
      Every cluster in `.tdr.obj$landmark.annot$clustering$ids` must
      appear in exactly one element.

  Mode B — per-cell labels (named `character` vector)

  :   A named character vector where
      [`names()`](https://rdrr.io/r/base/names.html) are the original
      cell IDs (before the `paste0(sample, "_", ...)` prefix added by
      tinydenseR) and values are cell-type labels. The vector must cover
      **all** cells across all samples; landmark-only labels are
      extracted automatically. Duplicate names are not allowed.

- .verbose:

  Logical; if `TRUE`, print progress messages when refreshing downstream
  slots (default `TRUE`).

- .name:

  Character name for storing the solution (default `NULL`, which
  auto-generates `"manual.<timestamp>"` for Mode A or
  `"labels.<timestamp>"` for Mode B). Cannot be `"ids"`, which is
  reserved for the active solution.

## Value

The `.tdr.obj` with the following updated fields:

- `$landmark.annot$celltyping$ids`:

  Factor vector of cell type labels for each landmark (active solution)

- `$landmark.annot$celltyping[[.name]]`:

  Copy of `$ids`, stored as a named solution for multi-solution
  workflows

If
[`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
has already been run, the following slots are also refreshed (no
re-mapping required):

- `$cellmap$celltyping$ids`:

  Per-sample list of cell-level celltype assignments, re-derived from
  the existing fuzzy graph

- `$density$composition$celltyping$cell.count`:

  Samples x cell types count matrix

- `$density$composition$celltyping$cell.perc`:

  Samples x cell types percentage matrix

If
[`get.lm()`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md)
has been run, any existing `$results$lm[[model]]$trad$celltyping` fits
are invalidated (set to `NULL`) with a warning.

## Details

`celltyping()` can be called at any point in the pipeline — before or
after
[`get.map`](https://opensource.nibr.com/tinydenseR/reference/get.map.md).
When called **after**
[`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md),
it automatically refreshes all celltyping-dependent downstream slots
(cell-level IDs, composition matrices, and the summary heatmap) without
re-reading expression data or re-running the UMAP transform. Any
existing
[`get.lm()`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md)
results that depend on celltyping are invalidated with a warning so the
user can re-run
[`get.lm()`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md).

Manual cell type annotation is useful when:

- Domain expertise is needed for fine-grained annotations

- Custom groupings are required for specific analyses

- Working with cytometry data where marker combinations define cell
  types

Cell type labels assigned here will be used in downstream analyses:

- [`get.map`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
  propagates labels to all cells via nearest landmarks

- [`get.lm`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md)
  enables cell-type-level differential density testing in traditional
  analysis

- [`get.pbDE`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md)
  enables cell-type-specific pseudobulk differential expression analysis

## See also

[`get.map`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
for automatic cell typing with symphony reference objects,
[`lm.cluster`](https://opensource.nibr.com/tinydenseR/reference/lm.cluster.md)
for the clustering that produces cluster IDs used here,
[`import_cell_annotations`](https://opensource.nibr.com/tinydenseR/reference/import_cell_annotations.md)
for automatic multi-column ingestion of cell-level annotations from
metadata

## Examples

``` r
if (FALSE) { # \dontrun{
# After clustering with get.graph()
lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks() |>
  get.graph()

# Map clusters to cell types (before get.map — traditional order)
celltype_map <- list(
  "CD4.T.cells" = c("cluster.01", "cluster.02"),
  "CD8.T.cells" = c("cluster.03"),
  "B.cells" = c("cluster.04", "cluster.05"),
  "unknown" = c("cluster.06")
)
lm.cells <- celltyping(lm.cells, celltype_map)

# Late celltyping — after get.map (new supported workflow)
lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks() |>
  get.graph() |>
  get.map() |>
  celltyping(celltype_map)

# Mode B — per-cell labels (named character vector)
# cell_labels is a named character vector: names = original cell IDs,
# values = cell-type labels, covering ALL cells across all samples.
lm.cells <- celltyping(lm.cells, cell_labels)
} # }
```
