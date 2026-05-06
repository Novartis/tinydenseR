# Plot Traditional Statistics

Visualizes results from traditional cluster/cell type-level differential
abundance testing. Shows effect sizes as heatmap with significance
markers, providing a complementary view to landmark-based analysis for
easier interpretation at the population level.

## Usage

``` r
plotTradStats(x, ...)

# S3 method for class 'TDRObj'
plotTradStats(
  x,
  .model.name = "default",
  .split.by = "clustering",
  .coefs = NULL,
  .q = 0.1,
  .row.space.scaler = 0.2,
  .col.space.scaler = 0.07,
  .label.substr.rm = "",
  .order.ids = FALSE,
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
  [`get.lm()`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md).
  Statistical results should be stored in `.tdr.obj$map$lm`.

- ...:

  Additional arguments passed to methods.

- .model.name:

  Character string naming which model fit to use from `.tdr.obj$map$lm`
  (default "default"). Must match a name used in
  `get.lm(.model.name = ...)`.

- .split.by:

  Character: "clustering" (default) or "celltyping" - which population
  grouping to use.

- .coefs:

  Character vector of coefficient names to plot. Defaults to all
  coefficients from traditional model.

- .q:

  Numeric q-value threshold for significance stars (default 0.1).

- .row.space.scaler:

  Numeric scaling for row height (default 0.2 inches per population).

- .col.space.scaler:

  Numeric scaling for column width (default 0.07).

- .label.substr.rm:

  Character substring to remove from labels (default "").

- .order.ids:

  Logical whether to order IDs based on dendrogram order (default
  FALSE).

## Value

A `patchwork` object combining two plots: (1) a dot plot showing mean
cell percentages per population, and (2) a heatmap showing log fold
changes colored by magnitude with significance indicated by asterisks.

## Details

Traditional analysis tests for DA at the cluster/celltype level by
aggregating cell counts per sample, then using standard models
(typically edgeR). This provides:

- Easier biological interpretation (known populations vs abstract
  landmarks)

- Comparison to landmark-based results for validation

- Detection of broad shifts affecting entire populations

Note: Traditional analysis is less sensitive to subtle within-cluster
variation that landmark-based methods can detect.

## See also

[`get.lm`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md),
[`plotTradPerc`](https://opensource.nibr.com/tinydenseR/reference/plotTradPerc.md),
[`plotBeeswarm`](https://opensource.nibr.com/tinydenseR/reference/plotBeeswarm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# After get.lm with traditional analysis
lm.obj <- get.lm(lm.obj, .design = design)

# Heatmap of cluster-level changes
plotTradStats(lm.obj, .split.by = "clustering")

# Cell type-level changes (using a different model)
plotTradStats(lm.obj, .model.name = "full", .split.by = "celltyping")
} # }
```
