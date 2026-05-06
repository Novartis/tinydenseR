# Refresh all clustering-dependent downstream slots

Mirrors `.refresh_celltyping` but for clustering annotations. Called
after reclustering or switching the active clustering solution.

## Usage

``` r
.refresh_clustering(.tdr.obj, .verbose = FALSE)
```

## Arguments

- .tdr.obj:

  A TDRObj with `@landmark.annot$clustering$ids` already set to the new
  values.

- .verbose:

  Logical; print progress messages.

## Value

The modified `.tdr.obj`.
