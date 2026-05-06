# Marker Gene/Protein Identification (Deprecated)

**\[deprecated\]**

`get.markerDE()` is deprecated. Use
[`get.pbDE`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md)`(.mode = "marker")`
instead. This function dispatches to `get.pbDE` with `.mode = "marker"`.

## Usage

``` r
get.markerDE(x, ...)

# S3 method for class 'TDRObj'
get.markerDE(
  x,
  .source = NULL,
  .geneset.ls = NULL,
  .id1.idx = NULL,
  .id2.idx = NULL,
  .id1 = NULL,
  .id2 = "..all.other.landmarks..",
  .id.from = "clustering",
  .model.name = "default",
  .comparison.name = NULL,
  .force.recalc = FALSE,
  .label.confidence = 0.5,
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object
  processed through
  [`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md).

- ...:

  Additional arguments passed to `get.pbDE`.

- .source:

  The raw data object (or NULL for file backend).

- .geneset.ls:

  Optional named list of character vectors for GSVA. RNA only.

- .id1.idx:

  Optional integer vector of landmark indices for group 1.

- .id2.idx:

  Optional integer vector of landmark indices for group 2.

- .id1:

  Character vector of cluster/celltype IDs for group 1.

- .id2:

  Character vector of IDs for group 2. Default
  `"..all.other.landmarks.."`.

- .id.from:

  Character: "clustering" or "celltyping".

- .model.name:

  Character: model name (default "default").

- .comparison.name:

  Character: comparison name, passed as `.result.name`.

- .force.recalc:

  Logical: recalculate even if results exist? Default FALSE.

- .label.confidence:

  Numeric (0-1): confidence threshold. Default 0.5.

## Value

The modified `x` with results stored in `$markerDE`.

## See also

[`get.pbDE`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md)
