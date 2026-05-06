# Mapping cells to landmarks

Projects all cells onto the landmark graph to compute fuzzy graph edge
weights between cells and landmarks. In addition, transfers cluster/cell
type labels from landmarks to all cells.

## Usage

``` r
get.map(x, ...)

# S3 method for class 'TDRObj'
get.map(
  x,
  .source = NULL,
  .ref.obj = NULL,
  .celltype.col.name = "cell_type",
  .verbose = TRUE,
  .seed = 123,
  .label.confidence = 0.5,
  .cache.on.disk = TRUE,
  .cache.path = NULL,
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object with
  `$graph` component populated by `get.graph`.

- ...:

  Additional arguments passed to methods.

- .source:

  The raw data object for non-file backends. `NULL` (default) for the
  files backend; otherwise a Seurat, SingleCellExperiment, or anndataR
  AnnData object. Used by `.get_sample_matrix()` to retrieve per-sample
  expression matrices.

- .ref.obj:

  Optional Symphony reference object for cell type annotation. Must have
  `Z_corr` field (harmony-corrected embeddings) and metadata with cell
  type labels. Only compatible with RNA assays. Replaces any existing
  `$graph$celltyping`.

- .celltype.col.name:

  Column name in `.ref.obj$meta_data` containing cell type labels
  (default "cell_type"). Only relevant when `.ref.obj` is provided.

- .verbose:

  Logical for progress messages (default TRUE).

- .seed:

  Integer seed for reproducibility (default 123).

- .label.confidence:

  Numeric scalar in `[0,1]` controlling the minimum posterior confidence
  required to assign a cell to a landmark‑derived cluster/celltype
  label.

- .cache.on.disk:

  Logical (default TRUE). When `TRUE`, four large per-sample slots
  (`clustering$ids`, `celltyping$ids`, `nearest.lm`, `fuzzy.graphs`) are
  serialized to disk as uncompressed RDS files and stored in `@cellmap`
  as attributed path strings. Downstream accessors (e.g.\\ in
  `get.pbDE`, `goi.summary`) read them back lazily on a per-sample
  basis. Set to `FALSE` to keep everything in memory.

- .cache.path:

  Character scalar or `NULL` (default). When non-`NULL`, overrides the
  default [`tempdir()`](https://rdrr.io/r/base/tempfile.html)-based
  cache root so that cached `@cellmap` files are written to this
  user-specified directory instead. This is useful on HPC systems where
  session suspension or temporary-directory cleanup would otherwise
  invalidate the `@cellmap` slot. The directory is created if it does
  not exist. Ignored when `.cache.on.disk = FALSE`.

  Cache files are stored under the system temporary directory
  ([`tempdir()`](https://rdrr.io/r/base/tempfile.html)) and are
  automatically removed when the R session ends via a registered
  finalizer. This means the cache is **ephemeral** and never persists
  across R sessions. There are no implications for reproducibility since
  the cache only stores intermediate results that are recomputed
  deterministically.

  Override the cache root by setting `.cache.path` to a user-controlled
  directory (e.g.\\ a project-level path on shared storage). This is
  useful on HPC systems where session suspension or cleanup can delete
  the system [`tempdir()`](https://rdrr.io/r/base/tempfile.html),
  rendering the `@cellmap` slot unusable.

## Value

Updated `.tdr.obj` with `@density` containing:

- `raw`: Matrix of raw fuzzy graph density sums (landmarks × samples).
  Each entry is the sum of cell-landmark fuzzy edge weights before
  size-factor normalization.

- `norm`: Matrix of size-factor-normalized fuzzy densities (landmarks ×
  samples). Each sample column is divided by \\n_j / \bar{n}\\ where
  \\n_j\\ is the cell count for sample \\j\\.

- `log.norm`: Matrix of log2-transformed normalized densities (landmarks
  × samples): `log2(norm + 0.5)`. Used by
  [`get.lm()`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md)
  for linear modeling and
  [`get.embedding()`](https://opensource.nibr.com/tinydenseR/reference/get.embedding.md)
  for unsupervised embeddings.

- `size.factors`: Named numeric vector (length N) of per-sample size
  factors used to normalize `raw` into `norm`.

- `clustering$ids`: List of named character vectors (one per sample)
  with cluster assignments for all cells.

- `clustering$cell.count`: Matrix (samples × clusters) of cell counts
  per cluster per sample. Used for "traditional" compositional
  statistics.

- `clustering$cell.perc`: Matrix (samples × clusters) of percentage of
  cells per cluster per sample. Used for "traditional" compositional
  statistics.

- `celltyping$ids`: List of named character vectors (one per sample)
  with cell type assignments (only if celltyping available or `.ref.obj`
  provided).

- `celltyping$cell.count`: Matrix (samples × cell types) of cell counts
  per cell type per sample. Used for "traditional" compositional
  statistics.

- `celltyping$cell.perc`: Matrix (samples × cell types) of percentage of
  cells per cell type per sample. Used for "traditional" compositional
  statistics.

- `nearest.lm`: List of matrices (one per sample) with nearest landmark
  indices for all cells from UMAP transform.

When `.cache.on.disk = TRUE`, the four cell-level slots above are stored
as attributed path strings (with `schema_v` and `bytes` attributes) in
`@cellmap` rather than in-memory objects. The cache root is stored in
`@config$.cache.root`. Use
[`tdr_cache_cleanup()`](https://opensource.nibr.com/tinydenseR/reference/tdr_cache_cleanup.md)
to remove cached files. If `.ref.obj` provided, also updates
`@landmark.annot$celltyping$ids` (factor of cell type assignments for
landmarks) and stores a named copy under
`@landmark.annot$celltyping[[.celltype.col.name]]`.

## Details

**Workflow Overview:**

For each sample, the function:

1.  Loads expression data and normalizes (size factors for RNA, marker
    subset for cytometry)

2.  Projects to PCA/Harmony space (matching landmark processing)

3.  Uses landmark UMAP model to compute fuzzy graph (cell-landmark edge
    weights) and find nearest landmarks

4.  Assigns clusters/cell types by confidence-thresholded voting

5.  Aggregates fuzzy graph edge weights into landmark densities per
    sample

**Label transfer confidence model:**

- Without `.ref.obj`: label confidence is the normalized fuzzy-mass
  ratio, \\\mathrm{conf}(c,\ell)=\sum\_{m\in\ell}w\_{c,m}/\sum_m
  w\_{c,m}\\, where \\w\_{c,m}\\ are UMAP-derived cell-landmark
  connection strengths.

- With `.ref.obj`: label confidence is kNN voting frequency in reference
  space, \\\mathrm{conf}(c,\ell)=N\_{c,\ell}/k\\ with \\k=10\\ nearest
  neighbors.

- In both modes, a label is accepted only if \\\mathrm{conf}(c,\ell) \ge
  {.label.confidence}\\; otherwise the cell is labeled
  `"..low.confidence.."`.

**Reference-Based Cell Typing:**

When `.ref.obj` is provided:

- Expression is mapped to reference via Symphony

- Cell types assigned by kNN voting (k = 10) in reference embedding

- Landmark cell types updated and used for visualization/statistics

- Overwrites the active celltyping solution (`$ids`). Previously stored
  named solutions (from
  [`celltyping()`](https://opensource.nibr.com/tinydenseR/reference/celltyping.md)
  or
  [`import_cell_annotations()`](https://opensource.nibr.com/tinydenseR/reference/import_cell_annotations.md))
  are preserved and can be restored via
  [`set_active_celltyping()`](https://opensource.nibr.com/tinydenseR/reference/set_active_celltyping.md)

**Fuzzy Graph Densities:**

The `norm` density matrix quantifies how strongly each landmark is
connected to cells in each sample, after size-factor normalization. High
values indicate the landmark's neighborhood is enriched in that sample.
The `raw` matrix stores the pre-normalization sums, enabling users to
explore alternative normalizations. This forms the basis for
differential density testing in `get.lm`.

## See also

[`get.graph`](https://opensource.nibr.com/tinydenseR/reference/get.graph.md),
[`get.lm`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md),
[`celltyping`](https://opensource.nibr.com/tinydenseR/reference/celltyping.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Complete workflow with mapping
lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks(.nHVG = 500) |>
  get.graph() |>
  get.map()

# Use Symphony reference for cell typing (RNA data only)
ref <- readRDS("pbmc_reference.rds")
lm.cells <- get.map(lm.cells, 
                    .ref.obj = ref, 
                    .celltype.col.name = "cell_type")
} # }
```
