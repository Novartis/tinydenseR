# Plot Differential Expression Analysis Results (Deprecated)

**\[deprecated\]**

`plotDEA()` has been renamed to
[`plotPbDE()`](https://opensource.nibr.com/tinydenseR/reference/plotPbDE.md)
for clarity. This function is provided for backward compatibility and
will be removed in a future version.

## Usage

``` r
plotDEA(x, ...)

# S3 method for class 'TDRObj'
plotDEA(
  x,
  .dea.obj,
  .coefs = colnames(x = .dea.obj$coefficients),
  .order.by = "clustering",
  .markers = NULL,
  .q = 0.1,
  .row.space.scaler = 0.2,
  .col.space.scaler = 0.065,
  .label.substr.rm = "",
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object.

- .dea.obj:

  Differential expression results from
  [`get.dea()`](https://opensource.nibr.com/tinydenseR/reference/get.dea.md)
  or
  [`get.pbDE()`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md).

- .coefs:

  Character vector of coefficient names to plot.

- .order.by:

  Character: "clustering" or "celltyping".

- .markers:

  Character vector of feature names.

- .q:

  Numeric adjusted p-value threshold.

- .row.space.scaler:

  Numeric row height scaling.

- .col.space.scaler:

  Numeric column width scaling.

- .label.substr.rm:

  Character substring to remove from labels.

## Value

A `ggplot` heatmap.

## See also

[`plotPbDE`](https://opensource.nibr.com/tinydenseR/reference/plotPbDE.md)
