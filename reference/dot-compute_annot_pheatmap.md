# Compute annotation pheatmap on-the-fly

Builds the median-expression matrix and pheatmap object for either
clustering or celltyping annotations. Called by `plotHeatmap`,
`plotBeeswarm`, and `plotTradStats` whenever a hierarchically ordered
heatmap is needed.

## Usage

``` r
.compute_annot_pheatmap(.tdr.obj, .id.from = "clustering")
```

## Arguments

- .tdr.obj:

  A TDRObj with `@landmark.annot[[.id.from]]$ids` already set.

- .id.from:

  Character: `"clustering"` or `"celltyping"`.

## Value

A `pheatmap` object (with `$gtable`, `$tree_row`, etc.).
