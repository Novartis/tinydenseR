# Recluster landmarks

Re-runs Leiden community detection with new parameters and refreshes all
clustering-dependent downstream slots (cell-level IDs, composition
matrices, traditional fits). This is the recommended way to explore
different resolutions after the initial
[`get.graph()`](https://opensource.nibr.com/tinydenseR/reference/get.graph.md)
call.

## Usage

``` r
# S3 method for class 'Seurat'
recluster(x, ...)

# S3 method for class 'SingleCellExperiment'
recluster(x, ...)

recluster(x, ...)

# S3 method for class 'TDRObj'
recluster(
  x,
  .cl.resolution.parameter = 0.8,
  .cl.method = "snn",
  .small.size = 3,
  .column.name = NULL,
  .seed = 123,
  .verbose = TRUE,
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object with
  [`get.graph()`](https://opensource.nibr.com/tinydenseR/reference/get.graph.md)
  already run.

- ...:

  Additional arguments passed to methods.

- .cl.resolution.parameter:

  Numeric resolution for Leiden CPM (default 0.8). Internally scaled by
  1e-3. Higher = finer clusters.

- .cl.method:

  Character: `"snn"` or `"fgraph"` (default `"snn"`).

- .small.size:

  Integer threshold for straggler absorption (default 3).

- .column.name:

  Character name for storing the solution (default
  `paste0("leiden.res.", .cl.resolution.parameter)`). Cannot be `"ids"`.

- .seed:

  Integer for reproducibility (default 123).

- .verbose:

  Logical (default TRUE).

## Value

Updated object with new active clustering and refreshed downstream
slots.
