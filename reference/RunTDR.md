# Run the full tinydenseR pipeline

S3 generic that executes the complete tinydenseR analysis pipeline:
landmark selection, graph construction, optional cell typing, and
mapping.

## Usage

``` r
RunTDR(x, ...)

# Default S3 method
RunTDR(
  x,
  .celltype.vec = NULL,
  .celltype.vec.overwrite = FALSE,
  .verbose = TRUE,
  .seed = 123,
  ...
)

# S3 method for class 'Seurat'
RunTDR(
  x,
  .sample.var,
  .assay = "RNA",
  .layer = "counts",
  .assay.type = "RNA",
  .harmony.var = NULL,
  .markers = NULL,
  .celltype.vec = NULL,
  .min.cells.per.sample = 10,
  .verbose = TRUE,
  .seed = 123,
  .prop.landmarks = 0.1,
  .n.threads = if (is.hpc()) {
     max(RhpcBLASctl::blas_get_num_procs(),
    RhpcBLASctl::omp_get_num_procs(), RhpcBLASctl::omp_get_max_threads(), na.rm = TRUE)

    } else {
     parallel::detectCores(logical = TRUE)
 },
  ...
)

# S3 method for class 'SingleCellExperiment'
RunTDR(
  x,
  .sample.var,
  .assay = "counts",
  .assay.type = "RNA",
  .bpcells.dir = NULL,
  .harmony.var = NULL,
  .markers = NULL,
  .celltype.vec = NULL,
  .min.cells.per.sample = 10,
  .verbose = TRUE,
  .seed = 123,
  .prop.landmarks = 0.1,
  .n.threads = if (is.hpc()) {
     max(RhpcBLASctl::blas_get_num_procs(),
    RhpcBLASctl::omp_get_num_procs(), RhpcBLASctl::omp_get_max_threads(), na.rm = TRUE)

    } else {
     parallel::detectCores(logical = TRUE)
 },
  ...
)

# S3 method for class 'HDF5AnnData'
RunTDR(
  x,
  .sample.var,
  .assay.type = "RNA",
  .h5ad.group = NULL,
  .bpcells.dir = NULL,
  .harmony.var = NULL,
  .markers = NULL,
  .celltype.vec = NULL,
  .min.cells.per.sample = 10,
  .verbose = TRUE,
  .seed = 123,
  .prop.landmarks = 0.1,
  .n.threads = if (is.hpc()) {
     max(RhpcBLASctl::blas_get_num_procs(),
    RhpcBLASctl::omp_get_num_procs(), RhpcBLASctl::omp_get_max_threads(), na.rm = TRUE)

    } else {
     parallel::detectCores(logical = TRUE)
 },
  ...
)

# S3 method for class 'character'
RunTDR(
  x,
  .sample.var,
  .assay.type = "RNA",
  .h5ad.group = NULL,
  .bpcells.dir = NULL,
  .harmony.var = NULL,
  .markers = NULL,
  .celltype.vec = NULL,
  .min.cells.per.sample = 10,
  .verbose = TRUE,
  .seed = 123,
  .prop.landmarks = 0.1,
  .n.threads = if (is.hpc()) {
     max(RhpcBLASctl::blas_get_num_procs(),
    RhpcBLASctl::omp_get_num_procs(), RhpcBLASctl::omp_get_max_threads(), na.rm = TRUE)

    } else {
     parallel::detectCores(logical = TRUE)
 },
  ...
)

# S3 method for class 'cytoset'
RunTDR(
  x,
  .sample.var,
  .assay.type = "cyto",
  .harmony.var = NULL,
  .markers = NULL,
  .celltype.vec = NULL,
  .min.cells.per.sample = 10,
  .verbose = TRUE,
  .seed = 123,
  .prop.landmarks = 0.1,
  .n.threads = if (is.hpc()) {
     max(RhpcBLASctl::blas_get_num_procs(),
    RhpcBLASctl::omp_get_num_procs(), RhpcBLASctl::omp_get_max_threads(), na.rm = TRUE)

    } else {
     parallel::detectCores(logical = TRUE)
 },
  ...
)

# S3 method for class 'flowSet'
RunTDR(
  x,
  .sample.var,
  .assay.type = "cyto",
  .harmony.var = NULL,
  .markers = NULL,
  .celltype.vec = NULL,
  .min.cells.per.sample = 10,
  .verbose = TRUE,
  .seed = 123,
  .prop.landmarks = 0.1,
  .n.threads = if (is.hpc()) {
     max(RhpcBLASctl::blas_get_num_procs(),
    RhpcBLASctl::omp_get_num_procs(), RhpcBLASctl::omp_get_max_threads(), na.rm = TRUE)

    } else {
     parallel::detectCores(logical = TRUE)
 },
  ...
)

# S3 method for class 'dgCMatrix'
RunTDR(x, .cell.meta, ...)

# S3 method for class 'DelayedMatrix'
RunTDR(x, .cell.meta, .bpcells.dir = NULL, ...)

# S3 method for class 'IterableMatrix'
RunTDR(x, .cell.meta, ...)
```

## Arguments

- x:

  A `DelayedMatrix` (features × cells for RNA).

- ...:

  Additional arguments passed to `.run_tdr_matrix` (e.g. `.sample.var`,
  `.assay.type`, `.harmony.var`, `.markers`, `.celltype.vec`,
  `.min.cells.per.sample`, `.verbose`, `.seed`, `.prop.landmarks`,
  `.n.threads`).

- .celltype.vec:

  Character(1). Column name in `.cell.meta` containing per-cell type
  labels, or `NULL`.

- .celltype.vec.overwrite:

  Logical. If `TRUE`, replace an existing `.celltype.vec` stored in the
  TDRObj config.

- .verbose:

  Logical. Print progress messages.

- .seed:

  Integer. Random seed.

- .sample.var:

  Character(1). Column name in `.cell.meta` identifying sample
  membership.

- .assay:

  Character(1). Name of the assay in `assayNames(x)`.

- .layer:

  Character(1). Layer within the assay (e.g. `"counts"`).

- .assay.type:

  Character. `"RNA"` or `"cyto"`.

- .harmony.var:

  Character vector of batch variable column names in sample-level
  metadata, or `NULL`.

- .markers:

  Character vector of marker names (required for cyto).

- .min.cells.per.sample:

  Integer. Minimum cells for a sample to be included.

- .prop.landmarks:

  Numeric in (0, 1\]. Proportion of cells as landmarks.

- .n.threads:

  Integer. Number of threads.

- .bpcells.dir:

  Character(1) or `NULL`. Directory path for the BPCells on-disk matrix.
  If `NULL` (default), uses a temporary directory
  ([`tempdir()`](https://rdrr.io/r/base/tempfile.html)) that is cleaned
  up on session end. If a path is given and already contains a valid
  BPCells matrix, the conversion step is skipped (cache hit).

- .h5ad.group:

  Character(1) or `NULL`. HDF5 group path containing the count matrix
  (e.g. `"/layers/counts"` or `"/X"`). If `NULL` (default), auto-detects
  by probing `/layers/counts` first, then falling back to `/X`.

- .cell.meta:

  A `data.frame` of per-cell metadata. Must have one row per cell with
  rownames matching cell IDs in `x` (`colnames` for RNA, `rownames` for
  cyto).

## Value

Depends on the method; see individual method documentation.

The updated
[`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md).

The Seurat object `x` with the TDRObj stored in
`Misc(x, slot = "tdr.obj")`.

The `SingleCellExperiment` `x` with the TDRObj stored in
`S4Vectors::metadata(x)$tdr.obj`.

A
[`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md).

A
[`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md).

The updated
[`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md).

The updated
[`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md).

A
[`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md).

A
[`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md).

## Details

`RunTDR` is a convenience wrapper around the explicit step-by-step
workflow
(`get.meta → get.cells → setup.tdr.obj → get.landmarks → get.graph → get.map → get.embedding`).
For equivalent inputs, seed, and arguments the two paths produce
**numerically identical** results.

## Methods (by class)

- `RunTDR(default)`: Run the pipeline on an existing TDRObj

  Executes
  `get.landmarks → get.graph → celltyping → get.map → get.embedding` on
  a pre-built
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md).

- `RunTDR(Seurat)`: Run the pipeline directly on a Seurat object

  Builds a
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
  from a Seurat object and executes the full pipeline. The finished
  TDRObj is stored in `SeuratObject::Misc(x, slot = "tdr.obj")`.

  All categorical cell-level columns in `x@meta.data` are automatically
  imported as named celltyping solutions via
  [`import_cell_annotations`](https://opensource.nibr.com/tinydenseR/reference/import_cell_annotations.md).
  The `.celltype.vec` column (if specified) is set as the active
  annotation.

- `RunTDR(SingleCellExperiment)`: Run the pipeline directly on a
  SingleCellExperiment

  Builds a
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
  from a `SingleCellExperiment` and executes the full pipeline. The
  finished TDRObj is stored in `S4Vectors::metadata(x)$tdr.obj`.

  If the assay is a `DelayedMatrix` (e.g., HDF5-backed), it is converted
  to a BPCells on-disk `IterableMatrix` and routed through
  `.run_tdr_matrix()` for efficient lazy access. If the assay is an
  in-memory matrix (e.g., `dgCMatrix`), the existing SCE backend path is
  used.

- `RunTDR(HDF5AnnData)`: Run the pipeline on an HDF5AnnData object

  Builds a
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
  from an HDF5-backed AnnData object. The expression matrix is converted
  to a BPCells on-disk directory for efficient lazy access, then the
  proven `IterableMatrix` pipeline is used for all downstream steps.

  Metadata (obs) is read via `anndataR`; the expression matrix is read
  and converted via `BPCells`.

- `RunTDR(character)`: Run the pipeline directly from an h5ad file path

  Reads metadata and gene/cell names from the h5ad file using `rhdf5`,
  opens the expression matrix with `BPCells`, and delegates to the
  internal `IterableMatrix` pipeline.

  This method supports **all** H5AD format versions, including pre-0.8.0
  files generated by older versions of Python anndata that are not
  supported by the `anndataR` package.

- `RunTDR(cytoset)`: Run the pipeline on a flowWorkspace cytoset

  Builds a
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
  from a `cytoset` object (one FCS sample per `cytoframe`) and executes
  the full pipeline.

- `RunTDR(flowSet)`: Run the pipeline on a flowCore flowSet

  Builds a
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
  from a `flowSet` object (one FCS sample per `flowFrame`) and executes
  the full pipeline. Requires only flowCore (not flowWorkspace).

- `RunTDR(dgCMatrix)`: Run the pipeline on a sparse matrix (dgCMatrix)

  Builds a
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
  from a `dgCMatrix` and per-cell metadata, then executes the full
  pipeline.

- `RunTDR(DelayedMatrix)`: Run the pipeline on a DelayedMatrix

  Converts the `DelayedMatrix` to a BPCells on-disk `IterableMatrix` for
  efficient lazy access, then delegates to the proven `IterableMatrix`
  pipeline via `.run_tdr_matrix()`.

- `RunTDR(IterableMatrix)`: Run the pipeline on a BPCells IterableMatrix
