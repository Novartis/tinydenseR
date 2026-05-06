# Plot PCA

Visualizes landmarks in PCA space with flexible coloring by features,
clusters, or statistical results. Supports both continuous (e.g., gene
expression, fold changes) and categorical (e.g., clusters, cell types)
overlays. Optional interactive hover shows landmark feature signatures.

## Usage

``` r
plotPCA(x, ...)

# S3 method for class 'TDRObj'
plotPCA(
  x,
  .PC.x = 1,
  .PC.y = 2,
  .feature = .tdr.obj@landmark.annot$clustering$ids,
  .cat.feature.color = Color.Palette[1, 1:5],
  .panel.size = if (is.numeric(x = .feature)) 2 else 3,
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
  PCA computed via
  [`get.landmarks()`](https://opensource.nibr.com/tinydenseR/reference/get.landmarks.md).

- ...:

  Additional arguments passed to methods.

- .PC.x:

  Integer specifying x-axis principal component (default 1).

- .PC.y:

  Integer specifying y-axis principal component (default 2).

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

  Numeric panel width/height in inches. Default 2 for numeric, 3 for
  categorical (to accommodate legend).

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

Interactive hover features require the ggiraph package. Install with
`install.packages("ggiraph")`.

## See also

[`plotUMAP`](https://opensource.nibr.com/tinydenseR/reference/plotUMAP.md)
for UMAP visualization,
[`get.features`](https://opensource.nibr.com/tinydenseR/reference/get.features.md)
for hover feature computation

## Examples

``` r
if (FALSE) { # \dontrun{
# From README: Basic PCA visualization
lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks(.nHVG = 500, .nPC = 3) |>
  get.graph()

plotPCA(lm.cells, .point.size = 1, .panel.size = 1.5)

# Color by fold change from get.lm()
condition.stats <- get.lm(lm.cells, .design = design)
plotPCA(lm.cells, 
        .feature = condition.stats$fit$coefficients[,"ConditionB"],
        .color.label = "log2 FC", 
        .midpoint = 0)

# Interactive with feature hover
lm.cells <- get.features(lm.cells)
plotPCA(lm.cells, .hover.stats = "marker")
} # }
```
