# Deprecated: Use get.pbDE(.mode = "marker") instead

`get.marker()` is deprecated. Use
[`get.pbDE`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md)`(.mode = "marker")`
instead.

## Usage

``` r
get.marker(
  .tdr.obj,
  .geneset.ls = NULL,
  .id1.idx = NULL,
  .id2.idx = NULL,
  .id1 = NULL,
  .id2 = "..all.other.landmarks..",
  .id.from = "clustering",
  .label.confidence = 0.5
)
```

## Arguments

- .tdr.obj:

  A tinydenseR object processed through
  [`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md).

- .geneset.ls:

  Optional named list of character vectors for GSVA.

- .id1.idx:

  Optional landmark indices for group 1.

- .id2.idx:

  Optional landmark indices for group 2.

- .id1:

  Cluster/celltype IDs for group 1.

- .id2:

  Reference group IDs. Default `"..all.other.landmarks.."`.

- .id.from:

  "clustering" or "celltyping".

- .label.confidence:

  Numeric (0-1): minimum confidence for cell assignment (default: 0.5).

## Value

A .tdr.obj with results stored in .tdr.obj\$markerDE.

## See also

[`get.pbDE`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md)
