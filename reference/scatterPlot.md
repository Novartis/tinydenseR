# Scatter Plot with Feature Coloring

Creates a 2D scatter plot with flexible x/y features and optional
coloring by a third feature. Useful for exploring relationships between
any features in the tinydenseR object (scaled expression, PCA
coordinates, graph embeddings, cluster IDs, metadata, etc.).

## Usage

``` r
scatterPlot(
  .x.feature,
  .y.feature,
  .color.feature = NULL,
  .x.label = "",
  .y.label = "",
  .color.label = "",
  .cat.feature.color = Color.Palette[1, 1:5],
  .panel.size = c(2, 2),
  .midpoint = NULL,
  .plot.title = "",
  .legend.position = "right",
  .point.size = 0.1,
  .seed = 123
)
```

## Arguments

- .x.feature:

  Numeric vector for x-axis values (e.g.,
  `.tdr.obj$scaled.landmarks[,"CD3"]` or `.tdr.obj$pca$embed[,"PC1"]`).

- .y.feature:

  Numeric vector for y-axis values.

- .color.feature:

  Optional vector for point colors. Can be numeric (continuous coloring
  with diverging blue-white-red scale) or categorical (discrete colors).
  If `NULL`, all points receive default ggplot2 coloring.

- .x.label:

  Character x-axis label (default "").

- .y.label:

  Character y-axis label (default "").

- .color.label:

  Character color legend label (default "").

- .cat.feature.color:

  Character vector of colors for categorical features (default
  `Color.Palette[1,1:5]`). Automatically interpolated if more categories
  than colors exist.

- .panel.size:

  Numeric vector `c(width, height)` in inches (default `c(2,2)`).

- .midpoint:

  Numeric midpoint for continuous color gradients. Defaults to median of
  `.color.feature`. Useful for centering diverging scales.

- .plot.title:

  Character plot title (default "").

- .legend.position:

  Character: "right", "left", "top", "bottom", or "none" (default
  "right").

- .point.size:

  Numeric point size (default 0.1). Increase for smaller datasets.

- .seed:

  Integer random seed for point order randomization (default 123).
  Prevents systematic overplotting of one group by another.

## Value

A `ggplot` object.

## Details

This flexible plotting function enables custom visualizations beyond the
standard `plotPCA` and `plotUMAP` interfaces. Use cases include:

- Plotting Laplacian Eigenmap coordinates from `.tdr.obj$graph$LE$embed`

- Exploring relationships between markers (e.g., CD3 vs CD4)

- Overlaying metadata on any 2D embedding

- Creating custom QC plots (e.g., library size vs mitochondrial %)

Points are plotted in randomized order (controlled by `.seed`) to avoid
systematic overplotting bias when groups overlap.

For numeric `.color.feature`, applies a diverging color scale
(blue-white-red) centered at `.midpoint` (defaults to median). For
categorical features, generates discrete colors by interpolating
`.cat.feature.color`.

## See also

[`plotPCA`](https://opensource.nibr.com/tinydenseR/reference/plotPCA.md),
[`plotUMAP`](https://opensource.nibr.com/tinydenseR/reference/plotUMAP.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# After processing
lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks() |> get.graph()

# Plot PC1 vs PC2 colored by cluster
scatterPlot(.x.feature = lm.cells$pca$embed[,"PC1"],
            .y.feature = lm.cells$pca$embed[,"PC2"],
            .color.feature = lm.cells$landmark.annot$clustering$ids,
            .x.label = "PC1", .y.label = "PC2",
            .color.label = "Cluster")

# Marker expression relationship
scatterPlot(.x.feature = lm.cells$scaled.landmarks[,"CD4"],
            .y.feature = lm.cells$scaled.landmarks[,"CD8A"],
            .color.feature = .meta$Condition,
            .x.label = "CD4", .y.label = "CD8A",
            .color.label = "Condition")

# Laplacian Eigenmap coordinates
scatterPlot(.x.feature = lm.cells$graph$LE$embed[,1],
            .y.feature = lm.cells$graph$LE$embed[,2],
            .color.feature = lm.cells$landmark.annot$clustering$ids,
            .x.label = "LE1", .y.label = "LE2")
} # }
```
