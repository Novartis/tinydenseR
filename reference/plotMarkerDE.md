# Plot Marker DE Results

Visualizes marker identification results from
`get.pbDE(.mode = "marker")` as a colored heatmap showing log fold
changes with significance overlays. Results are stored in
`.tdr.obj$markerDE[[.model.name]][[.comparison.name]]`.

## Usage

``` r
plotMarkerDE(x, ...)

# S3 method for class 'TDRObj'
plotMarkerDE(
  x,
  .de.obj = NULL,
  .model.name = "default",
  .comparison.name = NULL,
  .coefs = NULL,
  .order.by = "none",
  .markers = .tdr.obj@config$markers,
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
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object with
  marker results from `get.pbDE(.mode = "marker")`.

- ...:

  Additional arguments passed to methods.

- .de.obj:

  Optional: Direct DE results object for backward compatibility. If
  NULL, fetches from
  `.tdr.obj$markerDE[[.model.name]][[.comparison.name]]`.

- .model.name:

  Character identifying the model (default "default").

- .comparison.name:

  Character identifying the comparison (e.g., "cluster.1_vs_all").
  Required if `.de.obj` is NULL.

- .coefs:

  Character vector of coefficient names to plot. Default includes
  "(Intercept)" and ".id1" to show overall expression and group-specific
  enrichment.

- .order.by:

  Source for marker ordering: "clustering" or "celltyping".

- .markers:

  Features to include (default: all markers from `.order.by`).

- .q:

  Numeric FDR threshold for significance markers (default 0.1).

- .row.space.scaler:

  Numeric row height scaling (default 0.2).

- .col.space.scaler:

  Numeric column width scaling (default 0.065).

- .label.substr.rm:

  Character pattern to remove from coefficient labels.

## Value

A `ggplot` heatmap showing log fold changes per feature and coefficient.

## Details

The plot shows:

- Fill color: Log fold change (blue = lower in .id1, red = higher in
  .id1)

- Point size: -log10(adjusted p-value)

- Asterisk overlay: Features passing the q-value threshold

For marker analysis, the key coefficient is typically ".id1" which
represents the difference between group 1 and group 2 (or all other
landmarks).

## See also

[`get.pbDE`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md),
[`plotPbDE`](https://opensource.nibr.com/tinydenseR/reference/plotPbDE.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# After running get.pbDE in marker mode
lm.cells <- get.pbDE(lm.cells, .mode = "marker", .id = "cluster.3",
                     .result.name = "cluster3_markers")

# Visualize marker results
plotMarkerDE(lm.cells, .comparison.name = "cluster3_markers",
             .coefs = c("(Intercept)", ".id1"))
} # }
```
