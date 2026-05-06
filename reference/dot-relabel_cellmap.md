# Re-derive cell-level annotation IDs from fuzzy graph + landmark labels

Reads the cached (on-disk or in-memory) fuzzy graph for each sample and
applies the current landmark-level annotation via the same
weighted-voting logic used in
[`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md).
No expression data is re-read and no UMAP transform is repeated.

## Usage

``` r
.relabel_cellmap(.tdr.obj, .annot.type = "celltyping", .verbose = FALSE)
```

## Arguments

- .tdr.obj:

  A TDRObj where
  [`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
  has already been run.

- .annot.type:

  Character: `"celltyping"` (default) or `"clustering"`.

- .verbose:

  Logical; print progress messages.

## Value

The modified `.tdr.obj`.
