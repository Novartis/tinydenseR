# Plot UMAP Plot UMAP

Visualizes landmarks in UMAP space with flexible coloring by features,
clusters, or statistical results. Supports both continuous (e.g., gene
expression, fold changes) and categorical (e.g., clusters, cell types)
overlays. Optional interactive hover shows landmark feature signatures.
Identical interface to
[`plotPCA()`](https://opensource.nibr.com/tinydenseR/reference/plotPCA.md)
for easy comparison.

## Usage

``` r
plotUMAP(x, ...)

# S3 method for class 'TDRObj'
plotUMAP(
  x,
  .feature = .tdr.obj@landmark.annot$clustering$ids,
  .cat.feature.color = Color.Palette[1, 1:5],
  .panel.size = 2,
  .midpoint = NULL,
  .plot.title = "",
  .color.label = "",
  .legend.position = "right",
  .point.size = 0.1,
  .seed = 123,
  .hover.stats = "none",
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object with
  UMAP computed via
  [`get.graph()`](https://opensource.nibr.com/tinydenseR/reference/get.graph.md).

- ...:

  Additional arguments passed to methods.

- .feature:

  Numeric vector or factor of length `nrow(.tdr.obj$landmarks)` to color
  points by. Defaults to cluster IDs from
  `$landmark.annot$clustering$ids`. Can be:

  - Numeric: gene expression, fold changes, statistics (uses diverging
    color scale)

  - Factor: clusters, cell types, conditions (uses discrete color
    palette)

- .cat.feature.color:

  Character vector of colors for categorical features. Default uses
  `Color.Palette[1,1:5]` (5-color palette, interpolated to match factor
  levels).

- .panel.size:

  Numeric panel width/height in inches (default 2).

- .midpoint:

  Numeric midpoint for diverging color scale (continuous features only).
  Defaults to median of `.feature`.

- .plot.title:

  Character plot title (default "").

- .color.label:

  Character legend title (default "").

- .legend.position:

  Character legend position: "right" (default), "left", "top", "bottom",
  or "none".

- .point.size:

  Numeric point size (default 0.1).

- .seed:

  Integer random seed for plot point ordering (default 123).

- .hover.stats:

  Character specifying hover information: "none" (default) or "marker"
  (shows landmark feature signatures from
  [`get.features()`](https://opensource.nibr.com/tinydenseR/reference/get.features.md)).
  Requires ggiraph.

## Value

Plot object (class depends on interactivity):

- Static ggplot:

  `.hover.stats = "none"` returns a `ggplot` object

- Interactive girafe:

  `.hover.stats != "none"` returns a `girafe` object (if ggiraph
  installed), otherwise falls back to static `ggplot` with warning

## Note

Requires
[`get.graph()`](https://opensource.nibr.com/tinydenseR/reference/get.graph.md)
to have been run. Interactive hover features require the ggiraph
package. Install with `install.packages("ggiraph")`.

## See also

[`plotPCA`](https://opensource.nibr.com/tinydenseR/reference/plotPCA.md)
for PCA visualization,
[`get.graph`](https://opensource.nibr.com/tinydenseR/reference/get.graph.md),
[`get.features`](https://opensource.nibr.com/tinydenseR/reference/get.features.md)
for hover feature computation

## Examples

``` r
if (FALSE) { # \dontrun{
# After graph construction
lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks() |>
  get.graph()

# Basic UMAP
plotUMAP(lm.cells, .point.size = 1)

# Color by gene expression
plotUMAP(lm.cells, .feature = lm.cells$landmarks[,"CD4"], .color.label = "CD4")

# Color by fold change from get.lm()
condition.stats <- get.lm(lm.cells, .design = design)
plotUMAP(lm.cells, 
         .feature = condition.stats$fit$coefficients[,"ConditionB"],
         .color.label = "log2 FC", 
         .midpoint = 0)

# Interactive with feature hover
lm.cells <- get.features(lm.cells)
plotUMAP(lm.cells, .hover.stats = "marker")
} # }
```
