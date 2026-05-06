# Plot plsD Scores (Diagnostic and Component Views)

Visualize plsD results: either a diagnostic overview of all components
or a detailed component-level view showing embedding and Y-vs-score
scatter.

## Usage

``` r
plotPlsD(x, ...)

# S3 method for class 'TDRObj'
plotPlsD(
  x,
  .coef.col,
  .plsD.dim = NULL,
  .embed = "umap",
  .point.size = 1,
  .label.size = 3,
  .panel.size = 2,
  .seed = 123,
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object after
  [`get.plsD()`](https://opensource.nibr.com/tinydenseR/reference/get.plsD.md).

- ...:

  Additional arguments passed to methods.

- .coef.col:

  Character: coefficient name matching a slot in `.tdr.obj$plsD`.

- .plsD.dim:

  Integer or NULL: component to visualize (1-indexed). If NULL, plots Ak
  vs Sk diagnostic scatter for all components.

- .embed:

  Character: either `"umap"` or `"pca"`. Default "umap".

- .point.size:

  Numeric: point size for scatter plots. Default 1.

- .label.size:

  Numeric: label size for diagnostic scatter plots. Default 3. Applies
  only if .plsD.dim is NULL.

- .panel.size:

  Numeric: panel size in inches. Default 2.

- .seed:

  Integer: random seed for point ordering. Default 123.

## Value

A ggplot2 object (if `.plsD.dim = NULL`) or a patchwork composition (if
`.plsD.dim` is integer).

## Details

**Component view** (`.plsD.dim = integer`):

- Left panel: UMAP (or PCA) embedding colored by PLS scores (diverging
  scale)

- Right panel: Scatter of centered Y vs PLS scores, colored by raw
  (uncentered) density contrast coefficient — essential for
  distinguishing genuine contrast signal from structural score balancing

- Right panel scatter note: landmarks are colored by their *raw*
  (uncentered) density contrast coefficient (not the centered Y on the
  x-axis). This is intentional: if a cluster of landmarks with negative
  scores shows warm raw-Y colors (near zero or positive raw
  coefficient), it is likely a structural geometric counterweight rather
  than a genuinely depleted population. The same reasoning applies in
  reverse: large-magnitude positive-score landmarks with near-zero raw Y
  may reflect structural balance from the opposite side, depending on
  contrast direction.

**Diagnostic view** (`.plsD.dim = NULL`):

- X-axis: Smoothness (Sk); higher = large-scale graph-smooth structure

- Y-axis: Y-alignment (Ak); higher = stronger density coupling

- Color: \|q_k\| (Y-loading magnitude; larger = more Y variance
  captured)

- Labels: Component indices (1, 2, ...)

- Warning: high Ak + high Sk (top-left region of the diagnostic plot) is
  the typical pattern for a component dominated by a single extreme
  population. Inspect the corresponding component view and score-vs-Y
  scatter before interpreting gene loadings.

Ak is high by construction for early components: NIPALS PLS1 maximizes
covariance with Y at every deflation step. The diagnostic scatter is
most useful for identifying which components are also graph-smooth (high
Sk).

## See also

[`get.plsD`](https://opensource.nibr.com/tinydenseR/reference/get.plsD.md)
for computing plsD,
[`plotPlsDHeatmap`](https://opensource.nibr.com/tinydenseR/reference/plotPlsDHeatmap.md)
for expression heatmaps

## Examples

``` r
if (FALSE) { # \dontrun{
# After running plsD
lm.obj <- get.plsD(lm.obj, .coef.col = "Infection")

# Diagnostic overview
plotPlsD(lm.obj, .coef.col = "Infection", .plsD.dim = NULL)

# Visualize first component
plotPlsD(lm.obj, .coef.col = "Infection", .plsD.dim = 1)
} # }
```
