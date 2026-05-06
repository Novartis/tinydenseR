# Bidimensional Hexbin Plot for Marker Expression

Creates hexagonal binning plots to visualize the joint distribution of
two markers across landmarks. Useful for exploring marker coexpression
patterns and identifying cell populations. Optionally overlays reference
density from all landmarks when plotting a specific cluster/celltype.

## Usage

``` r
plot2Markers(x, ...)

# S3 method for class 'TDRObj'
plot2Markers(
  x,
  .id = NULL,
  .id.from = "clustering",
  .x.feature = "CD3",
  .y.feature = "CD20",
  .bins = 128,
  .legend.position = "right",
  .plot.title = "",
  .panel.size = 1.5,
  .reference = TRUE,
  .density.bins = 32,
  .sd.range = c(-3, 6),
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object
  processed through
  [`get.landmarks()`](https://opensource.nibr.com/tinydenseR/reference/get.landmarks.md)
  and
  [`get.graph()`](https://opensource.nibr.com/tinydenseR/reference/get.graph.md).

- ...:

  Additional arguments passed to methods.

- .id:

  Optional character: cluster or celltype ID to highlight. If NULL
  (default), plots all landmarks.

- .id.from:

  Character: "clustering" (default) or "celltyping". Source of `.id`.

- .x.feature:

  Character: column name from `.tdr.obj$landmarks` for x-axis (default
  "CD3").

- .y.feature:

  Character: column name from `.tdr.obj$landmarks` for y-axis (default
  "CD20").

- .bins:

  Integer: number of hexagonal bins for main plot (default 128). Higher
  values = finer resolution.

- .legend.position:

  Character: "right" (default), "left", "top", "bottom", or "none".

- .plot.title:

  Character: plot title (default "").

- .panel.size:

  Numeric: panel width/height in inches (default 1.5).

- .reference:

  Logical: if TRUE (default) and `.id` is specified, overlay reference
  density contours from all landmarks for comparison.

- .density.bins:

  Integer: number of bins for reference density contours (default 32).

- .sd.range:

  Numeric vector: range of standard deviations for outlier exclusion
  (default c(-3, 6)). Currently not implemented in the function.

## Value

A `ggplot` object with hexagonal binning showing marker coexpression.

## Details

The function creates a hexagonal heatmap where:

- Each hexagon represents a bin in 2D marker expression space

- Color intensity (red gradient) indicates cell density (log2 scale)

- For cytometry data, expression values are divided by 50 for scaling

- When `.id` is specified and `.reference = TRUE`, gray density contours
  show the distribution of all landmarks for context

This visualization helps identify:

- Marker coexpression patterns (e.g., CD3+CD4+ vs CD3+CD8+ populations)

- Whether a cluster is truly distinct from the background

- Outlier populations in marker expression space

## See also

[`plotPCA`](https://opensource.nibr.com/tinydenseR/reference/plotPCA.md),
[`plotUMAP`](https://opensource.nibr.com/tinydenseR/reference/plotUMAP.md)
for dimensionality reduction visualizations

## Examples

``` r
if (FALSE) { # \dontrun{
# After landmark identification and graph construction
lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks() |>
  get.graph()

# All landmarks
plot2Markers(lm.cells, .x.feature = "CD3", .y.feature = "CD4")

# Specific cluster with reference density
plot2Markers(lm.cells, 
             .id = "1", 
             .x.feature = "CD3", 
             .y.feature = "CD8",
             .reference = TRUE)
} # }
```
