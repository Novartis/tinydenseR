# Plot Mean Expression Heatmap

Displays a heatmap of median marker/gene expression across clusters or
cell types. Shows the expression patterns that define each population,
helping to validate cluster annotations and identify marker genes. The
heatmap is computed on-the-fly from the current active annotations and
expression data.

## Usage

``` r
plotHeatmap(x, ...)

# S3 method for class 'TDRObj'
plotHeatmap(x, .id.from = "clustering", ...)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object
  processed through
  [`get.graph()`](https://opensource.nibr.com/tinydenseR/reference/get.graph.md).

- ...:

  Additional arguments passed to methods.

- .id.from:

  Character: "clustering" or "celltyping" (default "clustering").
  Determines which population definitions to show. Use "clustering"
  after
  [`get.graph()`](https://opensource.nibr.com/tinydenseR/reference/get.graph.md),
  or "celltyping" after manual annotation with
  [`celltyping()`](https://opensource.nibr.com/tinydenseR/reference/celltyping.md).

## Value

A `gtable`/`grob` object (from `pheatmap`) rendered to the graphics
device. The heatmap includes hierarchical clustering of both features
and populations.

## Details

The heatmap content depends on assay type:

- **Cytometry**: All markers from `.tdr.obj$marker`

- **RNA-seq**: Top 3 positive and 3 negative genes per PC (from highly
  variable genes)

Rows (features) are ordered by hierarchical clustering to group
co-expressed markers/genes. Columns (populations) also clustered to
reveal relationships between cell types.

The heatmap is computed on-the-fly each time this function is called,
using `@landmark.annot[[.id.from]]$ids` and the expression matrix via
[`.compute_annot_pheatmap()`](https://opensource.nibr.com/tinydenseR/reference/dot-compute_annot_pheatmap.md).
This ensures the plot always reflects the current active annotations
(e.g.\\ after
[`set_active_clustering()`](https://opensource.nibr.com/tinydenseR/reference/set_active_clustering.md)
or
[`set_active_celltyping()`](https://opensource.nibr.com/tinydenseR/reference/set_active_celltyping.md)).

## See also

[`get.graph`](https://opensource.nibr.com/tinydenseR/reference/get.graph.md),
[`celltyping`](https://opensource.nibr.com/tinydenseR/reference/celltyping.md),
[`set_active_clustering`](https://opensource.nibr.com/tinydenseR/reference/set_active_clustering.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# After clustering
lm.cells <- setup.tdr.obj(.cells = .cells, .meta = .meta) |>
  get.landmarks() |> get.graph()

# Show cluster marker heatmap
plotHeatmap(lm.cells, .id.from = "clustering")

# After manual annotation
lm.cells <- celltyping(lm.cells, 
                       .celltyping.map = list("CD4_T" = c("cluster.1", "cluster.3"),
                                              "CD8_T" = c("cluster.2")))
plotHeatmap(lm.cells, .id.from = "celltyping")
} # }
```
