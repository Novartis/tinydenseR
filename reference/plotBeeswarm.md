# Bee Swarm Plot of Density Estimate Change

Creates a beeswarm plot showing effect sizes (log fold changes) from
differential density testing, with points colored by significance.
Optionally splits results by cluster or cell type and displays mean cell
percentages alongside for biological context.

## Usage

``` r
plotBeeswarm(x, ...)

# S3 method for class 'TDRObj'
plotBeeswarm(
  x,
  .model.name = "default",
  .coefs,
  .q = 0.1,
  .q.from = "pca.weighted.q",
  .split.by = if (length(x = .coefs) > 1) "none" else "clustering",
  .swarm.title = NULL,
  .label.substr.rm = "",
  .point.size = 0.1,
  .facet.scales = "fixed",
  .row.space.scaler = 0.2,
  .col.space.scaler = 0.1,
  .panel.width = 1.5,
  .legend.position = "right",
  .perc.plot = TRUE,
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

- .coefs:

  Character vector of coefficient names from the design matrix to plot.
  Must match column names in the model fit's coefficients. Can plot
  single or multiple coefficients.

- .q:

  Numeric q-value threshold for significance coloring (default 0.1).
  Points with q \< threshold are colored by direction (red/blue),
  otherwise gray.

- .q.from:

  Character specifying q-value source: "pca.weighted.q" (default,
  PCA-variance weighted) or "density.weighted.bh.fdr" (density
  weighted). See
  [`get.lm()`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md)
  for details.

- .split.by:

  Character controlling plot faceting: "none" (default if multiple
  coefficients), "clustering" (split by clusters), or "celltyping"
  (split by cell types). Default "clustering" for single coefficient.

- .swarm.title:

  Character plot title. If NULL and no facets, uses coefficient name.

- .label.substr.rm:

  Character substring to remove from axis labels (default "").

- .point.size:

  Numeric point size (default 0.1).

- .facet.scales:

  Character facet scales: "fixed" (default), "free", "free_x", "free_y".

- .row.space.scaler:

  Numeric scaling factor for row height when splitting (default 0.2
  inches per row).

- .col.space.scaler:

  Numeric scaling factor for column width when splitting (default 0.1).

- .panel.width:

  Numeric panel width in inches when plotting multiple coefficients
  (default 1.5).

- .legend.position:

  Character: "right" (default), "bottom", "left", "top", "none".

- .perc.plot:

  Logical whether to show cell percentage plot alongside beeswarm
  (default TRUE). Only applies when `.split.by != "none"`.

- .order.ids:

  Logical whether to order IDs based on dendrogram order (default
  FALSE).

## Value

A `patchwork` object combining plots:

- With percentages:

  `.perc.plot = TRUE` and `.split.by != "none"`: percentage plot +
  beeswarm plot side-by-side

- Beeswarm only:

  Otherwise: beeswarm plot alone (still wrapped in patchwork)

## Details

The beeswarm layout prevents overlapping points, making it easier to
assess the distribution of effect sizes across landmarks. Red indicates
significant increases (log FC \> 0), blue indicates significant
decreases, and gray indicates non-significant changes.

When split by clustering/celltyping, the percentage plot shows mean cell
type abundances across samples, helping interpret whether changes occur
in rare or common populations.

## See also

[`get.lm`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md),
[`plotDensity`](https://opensource.nibr.com/tinydenseR/reference/plotDensity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# After statistical testing
lm.obj <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks() |>
  get.graph() |>
  get.map()

design <- model.matrix(~ Condition + Replicate, data = .meta)
lm.obj <- get.lm(lm.obj, .design = design)

# Basic beeswarm split by clusters
plotBeeswarm(lm.obj, .coefs = "ConditionB")

# Multiple coefficients without splitting
plotBeeswarm(lm.obj, .coefs = c("ConditionB", "ConditionC"), .split.by = "none")

# Using a different model fit
plotBeeswarm(lm.obj, .model.name = "full", .coefs = "ConditionB")
} # }
```
