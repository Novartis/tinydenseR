# Multi-Backend Workflows with tinydenseR

## Introduction

`tinydenseR` supports multiple single-cell data formats through a
unified S3 dispatch system. Regardless of whether your data lives in a
Seurat object, a SingleCellExperiment, an HDF5-backed AnnData file, or a
bare sparse matrix, the same analysis functions work identically.

The supported backends are:

| Backend | Input class | Result container |
|:---|:---|:---|
| **TDRObj** (direct) | `TDRObj` | `TDRObj` |
| **Seurat** | `Seurat` | `Seurat` (TDRObj in `Misc` slot) |
| **SingleCellExperiment** | `SingleCellExperiment` | `SingleCellExperiment` (TDRObj in `metadata`) |
| **H5AD file path** | `character` (path to `.h5ad`) | `TDRObj` (bare) |
| **HDF5AnnData** | `HDF5AnnData` (via `anndataR`, modern files only) | `TDRObj` (bare) |
| **Sparse matrix** | `dgCMatrix` | `TDRObj` (bare) |
| **On-disk matrix** | `IterableMatrix` (BPCells) | `TDRObj` (bare) |
| **Delayed matrix** | `DelayedMatrix` | `TDRObj` (bare) |
| **Cytometry** | `flowSet` (flowCore) or `cytoset` (flowWorkspace) | `TDRObj` (bare) |

This means you can adopt `tinydenseR` without converting your data into
a specific format first.

## Entry Points: RunTDR()

[`RunTDR()`](https://opensource.nibr.com/tinydenseR/reference/RunTDR.md)
is the main entry point. It accepts different input types and runs the
full pipeline (landmark selection, graph construction, optional cell
typing, mapping, and embedding). The result is returned in a form that
matches the input container when possible.

### Seurat

``` r

library(tinydenseR)
library(SeuratObject)

seurat_result <- RunTDR(seurat_obj,
                        .sample.var = "sample_id",
                        .assay = "RNA",
                        .layer = "counts",
                        .assay.type = "RNA")
```

The returned object is the same Seurat object with the `TDRObj` stored
in `Misc(seurat_result, slot = "tdr.obj")`.

### SingleCellExperiment

``` r

library(SingleCellExperiment)

sce_result <- RunTDR(sce_obj,
                     .sample.var = "sample_id",
                     .assay = "counts",
                     .assay.type = "RNA")
```

The returned object is the same SCE with the `TDRObj` stored in
`S4Vectors::metadata(sce_result)$tdr.obj`.

### H5AD File Path (recommended)

The simplest way to work with `.h5ad` files is to pass the file path
directly. This supports **all** H5AD format versions, including
pre-0.8.0 files, using `rhdf5` for metadata and `BPCells` for the
expression matrix.

``` r

tdr <- RunTDR("data.h5ad",
              .sample.var = "sample_id",
              .assay.type = "RNA")
```

### HDF5AnnData (modern files only)

Alternatively, if you already have an `HDF5AnnData` object from
`anndataR`, you can pass it directly. Note that `anndataR` only supports
H5AD files created with anndata \>= 0.8.0.

``` r

library(anndataR)

adata <- read_h5ad("data.h5ad", as = "HDF5AnnData")
h5ad_result <- RunTDR(adata,
                      .sample.var = "sample_id",
                      .assay.type = "RNA")
```

Both H5AD entry points convert the expression matrix to a BPCells
on-disk format for efficient lazy access. The result is a bare `TDRObj`.

### Sparse matrix (dgCMatrix)

``` r

tdr <- RunTDR(dgc_matrix,
              .cell.meta = cell_metadata,
              .sample.var = "sample_id",
              .assay.type = "RNA")
```

### On-disk matrix (BPCells IterableMatrix)

``` r

tdr <- RunTDR(BPCells_matrix,
              .cell.meta = cell_metadata,
              .sample.var = "sample_id",
              .assay.type = "RNA")
```

For both sparse and on-disk matrix inputs, `.cell.meta` must be a
`data.frame` with one row per cell and rownames matching cell IDs. The
result is a bare `TDRObj`.

## Downstream Analysis: Identical API

Once
[`RunTDR()`](https://opensource.nibr.com/tinydenseR/reference/RunTDR.md)
has been called, the same pipeline code works regardless of the input
class. S3 dispatch routes each function call through the appropriate
wrapper automatically.

``` r

# These work identically whether result is Seurat or SCE:
result <- result |>
  get.lm(.design = design)

# Plotting works the same way
plotPCA(result)
plotBeeswarm(result, .coefs = "ConditionB")
```

This uniformity extends to all analysis and plotting functions in the
package:
[`get.pbDE()`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md),
[`get.plsD()`](https://opensource.nibr.com/tinydenseR/reference/get.plsD.md),
[`plotUMAP()`](https://opensource.nibr.com/tinydenseR/reference/plotUMAP.md),
[`plotHeatmap()`](https://opensource.nibr.com/tinydenseR/reference/plotHeatmap.md),
and so on.

## How It Works: GetTDR / SetTDR

The dispatch mechanism follows a simple three-step pattern. When you
call an analysis function on a container object (e.g., a Seurat object),
the S3 method:

1.  **Extracts** the `TDRObj` with
    [`GetTDR()`](https://opensource.nibr.com/tinydenseR/reference/GetTDR.md)
2.  **Runs** the computation on the `TDRObj`
3.  **Stores** the updated `TDRObj` back with
    [`SetTDR()`](https://opensource.nibr.com/tinydenseR/reference/SetTDR.md)

&nbsp;

    User calls:  get.lm(seurat_obj, .design = design)
                    |
                    v
    S3 dispatch --> get.lm.Seurat(x, ...)
                    |
                    +-- tdr <- GetTDR(x)            # Extract TDRObj from Misc slot
                    +-- tdr <- get.lm.TDRObj(tdr, ...)  # Run computation
                    +-- SetTDR(x, tdr)              # Store updated TDRObj back
                    |
                    v
    Returns:     updated seurat_obj (with results inside)

Where the `TDRObj` is stored depends on the container:

| Container | [`GetTDR()`](https://opensource.nibr.com/tinydenseR/reference/GetTDR.md) reads from | [`SetTDR()`](https://opensource.nibr.com/tinydenseR/reference/SetTDR.md) writes to |
|:---|:---|:---|
| Seurat | `Misc(x, slot = "tdr.obj")` | `Misc(x, slot = "tdr.obj")` |
| SingleCellExperiment | `S4Vectors::metadata(x)$tdr.obj` | `S4Vectors::metadata(x)$tdr.obj` |
| TDRObj | returns `x` directly | returns `tdr` directly |

Note: H5AD and matrix-class inputs always produce a bare `TDRObj`, so
`GetTDR`/`SetTDR` are not needed for those backends.

## Extracting the TDRObj

You can always extract the `TDRObj` directly to inspect results or work
with them programmatically:

``` r

# Extract from any container
tdr <- GetTDR(seurat_obj)

# Access linear model results
tdr$results$lm$default$fit$coefficients

# Access landmark PCA coordinates
head(tdr$landmark.embed$pca$coord)

# Access clustering assignments
tdr$landmark.annot$clustering$ids
```

The `$` accessor on a `TDRObj` maps to the underlying S4 slots, so
`tdr$results` is equivalent to `tdr@results`.

## Matrix-Class Backends

The `dgCMatrix`, `DelayedMatrix`, `IterableMatrix`, `flowSet`, and
`cytoset` backends produce a bare `TDRObj` directly — they do not wrap
the result back into a container object. After calling
[`RunTDR()`](https://opensource.nibr.com/tinydenseR/reference/RunTDR.md),
you work with the `TDRObj` directly:

``` r

# From a sparse matrix
tdr <- RunTDR(dgc_matrix,
              .cell.meta = cell_metadata,
              .sample.var = "sample_id",
              .assay.type = "RNA")

# All downstream functions work directly on the TDRObj
tdr <- tdr |>
  get.lm(.design = design)

plotPCA(tdr)
plotBeeswarm(tdr, .coefs = "ConditionB")
```

For `DelayedMatrix` inputs (e.g., HDF5-backed `SingleCellExperiment`
assays), the data is automatically converted to a BPCells on-disk format
for efficient lazy access before running the pipeline. You can control
the conversion directory with the `.bpcells.dir` argument.

The `flowSet` and `cytoset` backends (from `flowCore` and
`flowWorkspace`, respectively) are designed for flow, mass, and spectral
cytometry data and require `.assay.type = "cyto"` along with a
`.markers` argument specifying the channels to use. Both formats use the
same internal backend (`"cyto"`), so downstream behaviour is identical.
