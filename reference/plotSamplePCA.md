# Plot Sample PCA (Deprecated)

**\[deprecated\]**

This function is deprecated. Please use
[`plotSampleEmbedding`](https://opensource.nibr.com/tinydenseR/reference/plotSampleEmbedding.md)
with `.embedding = "pca"` instead.

## Usage

``` r
plotSamplePCA(x, ...)

# S3 method for class 'TDRObj'
plotSamplePCA(
  x,
  .labels.from = colnames(x = .tdr.obj@metadata)[1],
  .cat.feature.color = Color.Palette[1, 1:5],
  .point.size = 1,
  .panel.size = 2,
  .midpoint = if (is.numeric(x = .tdr.obj@metadata[[.labels.from]])) stats::median(x =
    .tdr.obj@metadata[[.labels.from]]) else NA,
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object
  processed through
  [`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
  and
  [`get.embedding()`](https://opensource.nibr.com/tinydenseR/reference/get.embedding.md).

- ...:

  Additional arguments passed to methods.

- .labels.from:

  Character specifying metadata column for coloring points.

- .cat.feature.color:

  Character vector of colors for categorical labels.

- .point.size:

  Numeric point size (default 1).

- .panel.size:

  Numeric panel width/height in inches (default 2).

- .midpoint:

  Numeric midpoint for diverging color scale (continuous labels only).

## Value

A `ggplot` object showing PC1 vs PC2 of sample-level density profiles.

## See also

[`plotSampleEmbedding`](https://opensource.nibr.com/tinydenseR/reference/plotSampleEmbedding.md)
for the recommended replacement
