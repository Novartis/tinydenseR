# Plot Density

Creates dot/line plots showing log2-transformed landmark fuzzy densities
per sample. Unlike
[`plotTradPerc()`](https://opensource.nibr.com/tinydenseR/reference/plotTradPerc.md)
which plots cell percentages at the cluster/celltype level, this
function plots landmark-level densities Useful for inspecting
distributions and paired/longitudinal patterns at the landmark
resolution.

## Usage

``` r
plotDensity(x, ...)

# S3 method for class 'TDRObj'
plotDensity(
  x,
  .x.split = colnames(x = .tdr.obj@metadata)[1],
  .pop = NULL,
  .pop.from = "clustering",
  .subject.id = NULL,
  .color.by = NULL,
  .x.space.scaler = 0.25,
  .height = 1.5,
  .cat.feature.color = Color.Palette[1, 1:5],
  .seed = 123,
  .orientation = "wide",
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

  Character specifying metadata column for x-axis grouping (default
  first column).

- .pop:

  Character vector of population names to plot. If NULL, plots all
  landmarks. If specified, plots only landmarks belonging to that
  population (from `.pop.from`).

- .pop.from:

  Character: "clustering" (default) or "celltyping" - which grouping to
  use for filtering landmarks when `.pop` is specified.

- .subject.id:

  Character metadata column for connecting paired samples with lines
  (e.g., "Subject"). Default NULL.

- .color.by:

  Character metadata column for coloring points. Default NULL (all
  black).

- .x.space.scaler:

  Numeric x-axis width scaling (default 0.25 inches per group).

- .height:

  Numeric plot height in inches (default 1.5).

- .cat.feature.color:

  Character vector of colors (default `Color.Palette[1,1:5]`).

- .seed:

  Integer random seed for jitter (default 123).

- .orientation:

  Character: "wide" (default) or "square" faceting.

## Value

A `ggplot` object showing log2-transformed landmark densities.

## See also

[`plotTradPerc`](https://opensource.nibr.com/tinydenseR/reference/plotTradPerc.md),
[`plotTradStats`](https://opensource.nibr.com/tinydenseR/reference/plotTradStats.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic density plot
plotDensity(lm.cells, .pop = "B.cells")

# With paired subject lines
plotDensity(lm.cells, 
              .pop = "B.cells",
              .subject.id = "MouseID",
              .color.by = "Treatment")
} # }
```
