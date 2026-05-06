# Create a hierarchical subset TDRObj

Given a parent `TDRObj` that has been processed through
[`get.map`](https://opensource.nibr.com/tinydenseR/reference/get.map.md),
creates a new *child* `TDRObj` containing only the cells matching the
specified selection criteria. The child references the parent's
expression data **without copying**: for file-backed data, the same RDS
paths are reused with an index filter; for in-memory backends (Seurat,
SCE, matrix, cytoset), the index vectors are intersected.

## Usage

``` r
get.subset(x, ...)

# S3 method for class 'TDRObj'
get.subset(
  x,
  .source = NULL,
  .id = NULL,
  .id.idx = NULL,
  .id.from = "clustering",
  .label.confidence = NULL,
  .prop.landmarks = 0.1,
  .min.cells.per.sample = 10,
  .verbose = TRUE,
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, or SingleCellExperiment object that has been processed through
  [`get.map`](https://opensource.nibr.com/tinydenseR/reference/get.map.md).

- ...:

  Additional arguments (currently unused).

- .source:

  The raw data object for non-file backends (`NULL` for the files
  backend). When calling on a bare `TDRObj`, pass the original Seurat or
  SCE object here for seurat/sce backends. When calling via a dispatch
  wrapper (`get.subset.Seurat`), this is filled automatically.

- .id:

  Optional character vector of cluster or celltype labels to keep.
  Exactly one of `.id` or `.id.idx` must be provided.

- .id.idx:

  Optional integer vector of landmark indices. Cells are selected via
  fuzzy confidence-thresholded voting using the stored UMAP fuzzy
  graphs.

- .id.from:

  Character: `"clustering"` (default) or `"celltyping"`. Source of
  labels when `.id` is used.

- .label.confidence:

  Numeric scalar in `[0,1]`. Minimum confidence for fuzzy label
  assignment. Defaults to the parent's stored threshold
  (`@config\$label.confidence`), or 0.5 if not set.

- .prop.landmarks:

  Numeric between 0 and 1 specifying proportion of cells to use as
  landmarks when the child pipeline runs
  [`get.landmarks()`](https://opensource.nibr.com/tinydenseR/reference/get.landmarks.md).
  Default 0.1 (10 percent).

- .min.cells.per.sample:

  Integer. Samples with fewer qualifying cells than this threshold are
  excluded from the child (default 10).

- .verbose:

  Logical: print progress? Default `TRUE`.

## Value

A new
[`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
with:

- `@cells`:

  Modified per-sample cell references (no expression data copied).

- `@metadata`:

  Filtered to retained samples, with updated `n.cells`.

- `@config`:

  Inherited parameters (backend, assay.type, markers, harmony.var) plus
  provenance in `@config\$.subset.provenance`.

- All other slots:

  Empty – ready for
  [`get.landmarks()`](https://opensource.nibr.com/tinydenseR/reference/get.landmarks.md).

## Details

The child object is returned with empty analysis slots (assay,
landmarks, graphs, density, cellmap, results). The user must re-run the
full pipeline (`get.landmarks` \\\to\\ `get.graph` \\\to\\ `get.map`
\\\to\\ `get.lm` \\\to\\ ...) on the child to obtain subset-specific
results.

Cell selection uses the same machinery as
[`get.pbDE`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md):

- `.id`: select cells whose per-cell label (from
  `@cellmap\$clustering\$ids` or `@cellmap\$celltyping\$ids`) matches
  the specified identifiers.

- `.id.idx`: select cells via fuzzy confidence-thresholded voting from
  the specified landmark indices (using the UMAP fuzzy simplicial set
  stored in `@cellmap\$fuzzy.graphs`).

Nested subsetting is supported: calling `get.subset()` on a child object
creates a grandchild that still references the original expression data
without additional copies.

## Session-scoped

Subset objects are session-scoped. Serialization via
[`saveRDS`](https://rdrr.io/r/base/readRDS.html) is not supported for
child TDRObjs that reference in-memory backends (seurat, sce, matrix,
cyto). For the files backend, the child can be saved as long as the
parent's RDS files remain on disk.

## See also

[`get.map`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
(required predecessor),
[`get.pbDE`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md)
(uses the same `.id`/`.id.idx` selection semantics)

## Examples

``` r
if (FALSE) { # \dontrun{
# After running the full pipeline on a parent population
parent <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks() |> get.graph() |>
  celltyping(.celltyping.map = list(CD4T = c("1","2"), CD8T = "3")) |>
  get.map()

# Create subset of CD8 T cells
child <- get.subset(parent, .id = "CD8T", .id.from = "celltyping")

# Run full pipeline on the subset (new landmarks, graph, etc.)
child <- child |>
  get.landmarks() |> get.graph() |>
  celltyping(.celltyping.map = list(Tem = "1", Tn = "2")) |>
  get.map() |>
  get.lm(.design = design)

# Nested subsetting
grandchild <- get.subset(child, .id = "Tem", .id.from = "celltyping")
grandchild <- grandchild |> get.landmarks() |> get.graph() |> get.map()
} # }
```
