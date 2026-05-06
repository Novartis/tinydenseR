# Leiden clustering of landmarks

Applies Leiden community detection to partition landmarks into clusters
based on either the SNN graph or UMAP fuzzy graph. Clusters are used for
visualization and as a coarser grouping for statistical testing. This is
a wrapper around `leiden.cluster`.

## Usage

``` r
lm.cluster(x, ...)

# S3 method for class 'TDRObj'
lm.cluster(
  x,
  .cl.method = "snn",
  .cl.resolution.parameter = 0.8,
  .seed = 123,
  .verbose = TRUE,
  .small.size = 3,
  .column.name = NULL,
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

- .cl.method:

  Character specifying clustering method: "snn" (shared nearest
  neighbors) or "fgraph" (fuzzy graph from UMAP). Default "snn".

- .cl.resolution.parameter:

  Numeric resolution for Leiden clustering (default 0.8). Internally
  scaled by 1e-3 for CPM objective. Higher values yield finer/more
  clusters.

- .seed:

  Integer seed for reproducibility (default 123).

- .verbose:

  Logical for progress messages (default TRUE).

- .small.size:

  Integer threshold for straggler absorption (default 3). Clusters
  smaller than this are merged into most connected neighbors.

- .column.name:

  Character name for storing the clustering solution (default `NULL`,
  which auto-generates
  `paste0("leiden.res.", .cl.resolution.parameter)`). Cannot be `"ids"`,
  which is reserved for the active solution.

## Value

Updated `.tdr.obj` with `@landmark.annot$clustering` updated:

- `$ids`: Factor of cluster assignments for each landmark (active
  solution)

- `$[[.column.name]]`: Same factor, stored as a named solution for
  multi-solution workflows

## Details

The fuzzy graph method ("fgraph") uses probabilistic UMAP edges (pruned
with tolerance 1/20 to remove weak connections) and may produce more
visually coherent clusters in UMAP space. The SNN method uses Jaccard
similarity of k-NN overlaps and often gives more stable results for
downstream statistics.

## See also

[`leiden.cluster`](https://opensource.nibr.com/tinydenseR/reference/leiden.cluster.md),
[`get.graph`](https://opensource.nibr.com/tinydenseR/reference/get.graph.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# After graph construction
lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks() |>
  get.graph()

# Extract clustering (called internally by get.graph)
lm.cells <- lm.cluster(lm.cells)

# Higher resolution clustering
lm.cells <- lm.cluster(lm.cells, .cl.resolution.parameter = 200)

# Use fuzzy graph instead of SNN
lm.cells <- lm.cluster(lm.cells, .cl.method = "fgraph")
} # }
```
