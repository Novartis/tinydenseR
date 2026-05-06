# Access per-cell map data for a single sample

Retrieves per-cell data (cluster IDs, cell type IDs, nearest landmarks,
or fuzzy graph) for a single sample, transparently reading from on-disk
cache when caching is active or from in-memory slots otherwise.

## Usage

``` r
get.cellmap(x, .slot, .sample)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
  processed through
  [`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md).

- .slot:

  Character. One of `"clustering"`, `"celltyping"`, `"nearest.lm"`,
  `"fuzzy.graphs"`. Legacy names (`"clustering.ids"`,
  `"celltyping.ids"`, `"nearest.landmarks"`, `"fuzzy.graph"`) are
  accepted for backward compatibility.

- .sample:

  Character. Sample identifier (must match a name in `names(x@cells)`).

## Value

The per-cell object for that sample and slot: a named character vector
(for `*ids`), an integer matrix (for `nearest.lm`), or a sparse matrix
(for `fuzzy.graphs`).

## See also

[`get.map`](https://opensource.nibr.com/tinydenseR/reference/get.map.md),
[`tdr_cache_validate`](https://opensource.nibr.com/tinydenseR/reference/tdr_cache_validate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Get cell type labels for the first sample
ct <- get.cellmap(x = lm.obj,
                  .slot = "celltyping",
                  .sample = names(lm.obj$cells)[1])
} # }
```
