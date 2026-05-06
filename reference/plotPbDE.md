# Plot Pseudobulk Differential Expression Results

Visualizes pseudobulk differential expression results as a heatmap
showing log fold changes for genes/markers across coefficients. Rows
ordered by cluster/celltype expression patterns, with significance
indicated by asterisks. Helps identify which features drive
population-level changes.

## Usage

``` r
plotPbDE(x, ...)

# S3 method for class 'TDRObj'
plotPbDE(
  x,
  .de.obj = NULL,
  .model.name = "default",
  .population.name = "all",
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
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object
  processed through
  [`get.pbDE()`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md).

- ...:

  Additional arguments passed to methods.

- .de.obj:

  Optional: differential expression results object. If NULL (default),
  retrieves results from
  `.tdr.obj$pbDE[[.model.name]][[.population.name]]`.

- .model.name:

  Character: model name to retrieve from `.tdr.obj$pbDE` (default
  "default"). Ignored if `.de.obj` is provided.

- .population.name:

  Character: population name to retrieve (default "all"). Ignored if
  `.de.obj` is provided.

- .coefs:

  Character vector of coefficient names to plot. Defaults to all
  coefficients.

- .order.by:

  Character: "clustering" (default) or "celltyping" - order rows by mean
  expression in these groups.

- .markers:

  Character vector of feature names (genes/proteins) to plot. Defaults
  to features shown in cluster/celltype heatmap (top PC contributors for
  RNA, all markers for cytometry).

- .q:

  Numeric adjusted p-value threshold for significance marking (default
  0.1).

- .row.space.scaler:

  Numeric row height scaling (default 0.2 inches per feature).

- .col.space.scaler:

  Numeric column width scaling (default 0.065 inches per coefficient).

- .label.substr.rm:

  Character substring to remove from labels (default "").

## Value

A `ggplot` heatmap showing effect sizes (log2 fold changes for RNA,
estimated differences for cytometry) with point size indicating adjusted
p-values and asterisks marking features meeting the significance
threshold.

## Details

This function shows which genes/markers are differentially expressed
between conditions, organized by their expression patterns across
clusters/cell types. The ordering helps identify:

- Marker genes defining specific populations

- Broadly vs. specifically regulated features

- Cell type-specific transcriptional responses

For RNA data, typically shows top PC-loading genes. For cytometry, shows
all markers.

## See also

[`get.pbDE`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md),
[`plotBeeswarm`](https://opensource.nibr.com/tinydenseR/reference/plotBeeswarm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# After pbDE analysis
design <- model.matrix(~ Condition, data = .meta)
lm.cells <- get.pbDE(lm.cells, .design = design)

# Heatmap of DE genes (uses .tdr.obj$pbDE$default$all)
plotPbDE(lm.cells, .coefs = "ConditionB", .order.by = "clustering")

# Plot results from specific population
lm.cells <- get.pbDE(lm.cells, .design = design, .id = "1", 
                     .id.from = "clustering", .population.name = "cluster1")
plotPbDE(lm.cells, .population.name = "cluster1", .coefs = "ConditionB")

# Focus on specific markers
plotPbDE(lm.cells, 
         .coefs = "ConditionB",
         .markers = c("CD4", "CD8A", "CD3D"))
} # }
```
