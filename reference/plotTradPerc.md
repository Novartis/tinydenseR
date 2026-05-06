# Plot Traditional Percentages

Creates dot/line plots showing cell percentages per sample for specific
populations. Useful for visually inspecting distribution of cell
abundances across conditions and identifying paired/longitudinal
patterns. Complements statistical test results from
[`plotTradStats()`](https://opensource.nibr.com/tinydenseR/reference/plotTradStats.md).

## Usage

``` r
plotTradPerc(x, ...)

# S3 method for class 'TDRObj'
plotTradPerc(
  x,
  .x.split = colnames(x = .tdr.obj@metadata)[1],
  .x.split.subset = NULL,
  .pop = NULL,
  .pop.from = "clustering",
  .order.pop = FALSE,
  .line.by = NULL,
  .dodge.by = NULL,
  .x.space.scaler = 0.25,
  .height = 1.5,
  .cat.feature.color = Color.Palette[1, 1:5],
  .seed = 123,
  .orientation = "wide",
  .log2.y = FALSE,
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object
  processed through
  [`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md).

- ...:

  Additional arguments passed to methods.

- .x.split:

  Character specifying metadata column for x-axis grouping. Defaults to
  first column.

- .x.split.subset:

  Optional character vector to subset `.x.split` categories. Default
  NULL.

- .pop:

  Character vector of population names to plot. If NULL, plots all
  populations from `.pop.from`.

- .pop.from:

  Character: "clustering" (default) or "celltyping" - which grouping to
  plot.

- .order.pop:

  Logical whether to order populations based on dendrogram order
  (default FALSE).

- .line.by:

  Character metadata column for connecting paired samples with lines
  (e.g., "Subject" for longitudinal data). Default NULL (no lines).

- .dodge.by:

  Character metadata column for coloring/dodging points. Default NULL
  (all black).

- .x.space.scaler:

  Numeric scaling factor for x-axis panel width (default 0.25 inches per
  group).

- .height:

  Numeric plot height in inches (default 1.5).

- .cat.feature.color:

  Character vector of colors for `.dodge.by` categories (default
  `Color.Palette[1,1:5]`).

- .seed:

  Integer random seed for x-axis jitter (default 123).

- .orientation:

  Character: "wide" (default, all populations in one row) or "square"
  (facet grid).

- .log2.y:

  Logical whether to log2-transform y-axis percentages (default FALSE).

## Value

A `ggplot` object showing cell percentages with optional paired
connections.

## Details

This function visualizes the raw cell percentages used in traditional DA
testing. Points show individual samples, and lines (if `.line.by`
specified) connect repeated measures from the same subject/mouse. Useful
for:

- Inspecting data distribution before statistical testing

- Identifying outliers or batch effects

- Visualizing paired/longitudinal designs

- Confirming significant results have biological meaning

## See also

[`plotTradStats`](https://opensource.nibr.com/tinydenseR/reference/plotTradStats.md),
[`get.lm`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# After mapping
lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks() |>
  get.graph() |>
  get.map()

# Plot CD4 T cells across conditions
plotTradPerc(lm.cells, .pop = "CD4.T.cells", .x.split = "Condition")

# Paired design with subject lines
plotTradPerc(lm.cells, 
             .pop = "CD4.T.cells",
             .line.by = "Subject",
             .dodge.by = "Timepoint")
} # }
```
