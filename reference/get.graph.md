# Graph embedding of landmarks

Builds k-NN graph and trains UMAP model for landmark cells. This creates
the graph structure and transformation model required for generating
fuzzy density matrices in `get.map`. Clustering is performed as a
secondary step for visualization and interpretation, but is not the
primary objective.

## Usage

``` r
get.graph(x, ...)

# S3 method for class 'TDRObj'
get.graph(
  x,
  .k = 20,
  .scale = if (!is.null(x@integration$harmony.obj) || x@config$assay.type == "RNA") FALSE
    else TRUE,
  .verbose = TRUE,
  .seed = 123,
  .cl.method = "snn",
  .cl.resolution.parameter = 0.8,
  .small.size = floor(x = nrow(x = x@assay$expr)/200),
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object.

- ...:

  Additional arguments passed to methods.

- .k:

  Integer number of nearest neighbors for graph construction (default
  20). Higher values create more connected graphs and smoother UMAP
  embeddings. Used for k-NN graph construction and passed to
  [`uwot::umap`](https://jlmelville.github.io/uwot/reference/umap.html)
  as `n_neighbors`.

- .scale:

  Logical indicating whether to scale features before UMAP. Defaults to
  FALSE when Harmony batch correction is active (any assay type) or for
  RNA assays (PCA already scaled). Defaults to TRUE for cytometry
  without Harmony (raw marker input to UMAP).

- .verbose:

  Logical for progress messages (default TRUE).

- .seed:

  Integer seed for reproducibility of UMAP and clustering (default 123).

- .cl.method:

  Character specifying clustering method: "snn" (shared nearest
  neighbors via Jaccard similarity) or "fgraph" (fuzzy graph from UMAP).
  Default "snn".

- .cl.resolution.parameter:

  Numeric controlling cluster granularity (default 0.8). Higher values
  produce more fine-grained clusters.

- .small.size:

  Integer threshold for straggler clusters (default 0.5% of landmarks).
  Clusters smaller than this are absorbed into nearest large cluster.

## Value

Updated `.tdr.obj` with `$graph` component containing:

- `uwot`: UMAP model with multiple components:

  - `$embedding`: 2D UMAP coordinates for visualization

  - `$nn$euclidean$idx`: k-NN indices matrix (landmarks × k neighbors)

  - `$nn$euclidean$dist`: k-NN distance matrix

  - `$fgraph`: Fuzzy graph (landmark-landmark probabilistic edges) for
    optional clustering

  - `$model`: Trained UMAP model for projecting query cells in `get.map`

- `adj.matrix`: Sparse k-NN adjacency matrix (non-symmetric).

- `snn`: Shared nearest neighbor graph via Jaccard similarity.

- `LE`: Laplacian Eigenmap components computed from symmetrized k-NN
  graph, including `$W.sym` (symmetrized adjacency), `$L` (Laplacian),
  `$Disqrt` (degree matrix), `$vectors`, `$values`, `$converged`,
  `$nontriv` (non-trivial eigenvalue indices), `$elbow` (selected
  dimensionality), `$keep` (retained eigenvectors), and `$embed` (final
  normalized embedding). Used for clustering initialization only.

- `clustering`: List containing `$ids` (factor of cluster assignments),
  `$median.exprs` (matrix of mean expression per cluster), and
  `$pheatmap` (heatmap object).

## Details

The function orchestrates graph-based analysis in the following order:

**1. k-NN Graph Construction:** Builds k-nearest neighbor graph from
PCA/Harmony embeddings and converts it to multiple representations for
different downstream applications:

- `adj.matrix`: Non-symmetric sparse adjacency matrix (used to compute
  SNN graph)

- `snn`: Shared nearest neighbor graph via Jaccard similarity (used for
  clustering)

- `W.sym`: Symmetrized binary adjacency (used for Laplacian Eigenmap
  spectral analysis)

**2. UMAP Model Training:** Trains UMAP transformation on landmark cells
with `ret_model=TRUE` to enable projection of query cells later. The
UMAP model includes the fuzzy graph (landmark-landmark probabilistic
edges) which is used for optional clustering and, after projection in
`get.map`, generates cell-landmark edge weights essential for computing
fuzzy density matrices. The UMAP model is the primary output of
`get.graph`.

**3. Laplacian Eigenmap (LE):** Computes spectral embedding by solving
the generalized eigenvalue problem of the graph Laplacian.
Dimensionality is automatically selected via elbow detection on
eigenvalue spectrum. The LE embedding is currently used only for
warm-start initialization in clustering.

**4. Clustering (Secondary):** Applies Leiden algorithm to identify
communities for visualization and interpretation. Uses SNN graph
(Jaccard similarity) or UMAP fuzzy graph depending on `.cl.method`.
Small "straggler" clusters are absorbed into neighbors. While useful for
exploration, clustering is not required for downstream statistical
analysis, unless traditional analysis is requested later in `get.lm`.

## See also

[`setup.tdr.obj`](https://opensource.nibr.com/tinydenseR/reference/setup.tdr.obj.md),
[`get.landmarks`](https://opensource.nibr.com/tinydenseR/reference/get.landmarks.md),
[`lm.cluster`](https://opensource.nibr.com/tinydenseR/reference/lm.cluster.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Typical workflow after landmark selection
lm.cells <- setup.tdr.obj(.cells = .cells, 
                          .meta = .meta,
                          .assay.type = "RNA") |>
  get.landmarks(.nHVG = 500, .nPC = 3) |>
  get.graph(.k = 10)

# Higher resolution for finer clusters
lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks() |>
  get.graph(.cl.resolution.parameter = 200)

# Use fuzzy graph for clustering instead of SNN
lm.cells <- get.graph(lm.cells, .cl.method = "fgraph")
} # }
```
